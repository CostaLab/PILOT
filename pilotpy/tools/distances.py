"""
Cost matrices between cell types and optimal transport between samples.

The cost matrix says how expensive it is to move cells from one cell type to
another. The Wasserstein distance between two samples is the cheapest way to
turn the cell type proportions of one sample into those of the other, using
this cost matrix.
"""
import numpy as np
import pandas as pd
import ot
import scipy
import scipy.sparse
from scipy.sparse.csgraph import shortest_path
from scipy.spatial.distance import cdist, pdist, squareform

#: Ways to compute the cost matrix between cell types.
COST_METHODS = ('centroid', 'distribution', 'graph', 'expression', 'jaccard', 'precomputed')

#: Optimal transport methods for the distance between two samples.
OT_METHODS = ('unreg', 'sinkhorn', 'sinkhorn_divergence', 'unbalanced')

#: Metrics with which a k-nearest-neighbor search can use a fast tree index.
_TREE_METRICS = ('euclidean', 'minkowski', 'cityblock', 'manhattan', 'chebyshev')


def _cell_types(annot):
    """Cell types in the order used by PILOT, and the cell type index of every cell."""
    cells = annot['cell_type'].unique()
    codes = pd.Categorical(annot['cell_type'], categories=cells).codes
    return cells, np.asarray(codes)


def _as_cost_df(dis, cells):
    cost = pd.DataFrame(np.asarray(dis, dtype=float), columns=cells)
    cost['cell_types'] = cells
    return cost.set_index('cell_types')


def _check_metric_input(values, metric):
    if metric == 'jensenshannon' and np.nanmin(values) < 0:
        raise ValueError("metric='jensenshannon' needs non-negative values. Use another metric "
                         "(e.g. 'cosine'), or cost_method='expression'.")


def cost_centroid(annot, data, metric='cosine', centroid='median'):
    """
    Distance between the centers of the cell types (the default of PILOT).

    Each cell type is represented by one point: the median (default) or mean
    of its cells, or, with 'sum', the normalized sum of its cells.

    Parameters
    ----------
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    data : DataFrame or numpy.ndarray
        Embedding of the cells, one row per cell (``adata.uns['data']``).
    metric : str, optional
        Metric of ``scipy.spatial.distance.pdist``, by default 'cosine'.
    centroid : {'median', 'mean', 'sum'}, optional
        Center of a cell type, by default 'median'.

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    if centroid not in ('median', 'mean', 'sum'):
        raise ValueError("centroid must be 'median', 'mean' or 'sum'.")
    cells = annot[annot.columns[0]].unique()
    centroids = []
    for cell in cells:
        part = data[annot[annot.columns[0]] == cell]
        if centroid == 'median':
            centroids.append(list(part.median(axis=0)))
        elif centroid == 'mean':
            centroids.append(list(part.mean(axis=0)))
        else:
            centroids.append(list(part.sum(axis=0) / (part.sum(axis=0).sum())))
    _check_metric_input(np.asarray(centroids, dtype=float), metric)
    dis = squareform(pdist(centroids, metric=metric), force='no', checks=True)
    return dis, _as_cost_df(dis, cells)


def cost_distribution(annot, data, metric='euclidean', max_cells=500, random_state=0):
    """
    Distance between all cells of two cell types (energy distance).

    Instead of one center per cell type, all cells are compared: two cell
    types are close when their cells are close to each other, compared with
    the spread inside each cell type. At most ``max_cells`` random cells per
    cell type are used, so it stays fast.

    Parameters
    ----------
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    data : DataFrame or numpy.ndarray
        Embedding of the cells, one row per cell (``adata.uns['data']``).
    metric : str, optional
        Metric of ``scipy.spatial.distance.cdist``, by default 'euclidean'.
    max_cells : int, optional
        Maximum number of random cells per cell type, by default 500.
    random_state : int, optional
        Random seed for choosing the cells, by default 0.

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    cells, codes = _cell_types(annot)
    values = np.asarray(data, dtype=float)
    _check_metric_input(values, metric)
    rng = np.random.default_rng(random_state)
    groups = []
    for k in range(len(cells)):
        idx = np.where(codes == k)[0]
        if max_cells is not None and len(idx) > max_cells:
            idx = rng.choice(idx, size=max_cells, replace=False)
        groups.append(values[idx])
    within = [cdist(g, g, metric=metric).mean() for g in groups]
    n = len(cells)
    dis = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            between = cdist(groups[a], groups[b], metric=metric).mean()
            energy = max(2 * between - within[a] - within[b], 0.0)
            dis[a, b] = dis[b, a] = np.sqrt(energy)
    return dis, _as_cost_df(dis, cells)


def cost_graph(annot, data, metric='euclidean', n_neighbors=15):
    """
    Distance between cell types from the neighbor graph of the cells.

    A k-nearest-neighbor graph of all cells is built. Two cell types are
    connected when many of their cells are neighbors (more than expected by
    chance; the same idea as PAGA). The cost is 1 - connectivity, and the
    cost between two cell types is the cheapest path through the graph of
    cell types, so cell types that are linked through an intermediate cell
    type (e.g. healthy -> damaged -> dead cells) are close.

    Parameters
    ----------
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    data : DataFrame or numpy.ndarray
        Embedding of the cells, one row per cell (``adata.uns['data']``).
    metric : str, optional
        Metric of the neighbor search, by default 'euclidean'.
    n_neighbors : int, optional
        Number of neighbors of every cell, by default 15.

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    from sklearn.neighbors import NearestNeighbors
    cells, codes = _cell_types(annot)
    values = np.asarray(data, dtype=float)
    if metric in ('cosine', 'correlation'):
        if metric == 'correlation':
            values = values - values.mean(axis=1, keepdims=True)
        norm = np.linalg.norm(values, axis=1, keepdims=True)
        values = values / np.where(norm == 0, 1, norm)
        nn_metric = 'euclidean'
    else:
        nn_metric = metric
    algorithm = 'auto' if nn_metric in _TREE_METRICS else 'brute'
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, metric=nn_metric, algorithm=algorithm)
    nn.fit(values)
    neighbors = nn.kneighbors(values, return_distance=False)[:, 1:]
    n = len(cells)
    src = np.repeat(codes, n_neighbors)
    dst = codes[neighbors.ravel()]
    edges = np.zeros((n, n))
    np.add.at(edges, (src, dst), 1)
    edges = edges + edges.T
    sizes = np.bincount(codes, minlength=n).astype(float)
    total = sizes.sum()
    expected = 2 * n_neighbors * np.outer(sizes, sizes) / max(total - 1, 1)
    connectivity = np.minimum(1.0, edges / np.where(expected == 0, 1, expected))
    direct = 1.0 - connectivity
    np.fill_diagonal(direct, 0.0)
    dis = shortest_path(direct, method='D', directed=False)
    return dis, _as_cost_df(dis, cells)


def _pseudobulk(adata, annot, layer=None):
    """log(CPM + 1) of the summed counts of every cell type (cell types x genes)."""
    cells, codes = _cell_types(annot)
    X = adata.layers[layer] if layer is not None else adata.X
    if X.shape[0] != len(codes):
        raise ValueError("The expression matrix and the cell annotation have different numbers of cells.")
    indicator = scipy.sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                                        shape=(len(cells), len(codes)))
    sums = indicator @ X
    sums = np.asarray(sums.todense()) if scipy.sparse.issparse(sums) else np.asarray(sums)
    if np.nanmin(sums) < 0:
        raise ValueError("The expression values are negative (e.g. scaled data). Use raw or "
                         "normalized counts, e.g. layer='counts'.")
    totals = sums.sum(axis=1, keepdims=True)
    return cells, np.log1p(sums / np.where(totals == 0, 1, totals) * 1e6)


def cost_expression(adata, annot, metric='correlation', layer=None, n_genes=2000):
    """
    Distance between the average gene expression (pseudobulk) of the cell types.

    The counts of all cells of a cell type are summed, normalized (log CPM),
    and the ``n_genes`` genes that vary most between the cell types are kept.
    Works best with raw counts in ``adata.X`` or in ``layer``.

    Parameters
    ----------
    adata : AnnData
        AnnData object with the counts in ``adata.X`` or in ``layer``.
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    metric : str, optional
        Metric of ``scipy.spatial.distance.pdist``, by default 'correlation'.
    layer : str, optional
        Layer with the counts. If None (default), ``adata.X``.
    n_genes : int, optional
        Number of most variable genes, by default 2000.

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    cells, expr = _pseudobulk(adata, annot, layer)
    if n_genes is not None and n_genes < expr.shape[1]:
        keep = np.argsort(expr.var(axis=0))[::-1][:n_genes]
        expr = expr[:, keep]
    _check_metric_input(expr, metric)
    dis = squareform(pdist(expr, metric=metric), force='no', checks=True)
    return dis, _as_cost_df(dis, cells)


def cost_jaccard(adata, annot, layer=None, n_markers=100):
    """
    Jaccard distance between the marker genes of the cell types.

    For each cell type, the ``n_markers`` genes with the highest expression
    compared with the other cell types are its markers. The distance is
    1 - (shared markers / all markers of the two cell types): 0 when two cell
    types have the same markers, 1 when they share none.

    Parameters
    ----------
    adata : AnnData
        AnnData object with the counts in ``adata.X`` or in ``layer``.
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    layer : str, optional
        Layer with the counts. If None (default), ``adata.X``.
    n_markers : int, optional
        Number of marker genes per cell type, by default 100.

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    cells, expr = _pseudobulk(adata, annot, layer)
    n = len(cells)
    markers = []
    for k in range(n):
        others = np.delete(expr, k, axis=0).mean(axis=0) if n > 1 else 0
        score = np.where(expr[k] > 0, expr[k] - others, -np.inf)
        markers.append(set(np.argsort(score)[::-1][:n_markers]))
    dis = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            union = len(markers[a] | markers[b])
            dis[a, b] = dis[b, a] = 1 - len(markers[a] & markers[b]) / union if union else 0.0
    return dis, _as_cost_df(dis, cells)


def cost_precomputed(annot, cost):
    """Use a cost matrix given by the user (a DataFrame with cell types as rows and columns).

    Parameters
    ----------
    annot : DataFrame
        Cell annotation with the columns 'cell_type', 'sampleID' and 'status'
        (``adata.uns['annot']`` after ``pl.tl.wasserstein_distance``).
    cost : DataFrame
        Cost matrix with the cell types as index and columns (symmetric,
        non-negative, 0 on the diagonal).

    Returns
    -------
    numpy.ndarray, DataFrame
        The cost matrix as an array and as a DataFrame with the cell types as
        index and columns.
    """
    cells = annot['cell_type'].unique()
    if not isinstance(cost, pd.DataFrame):
        raise ValueError("cost must be a pandas DataFrame with the cell types as index and columns.")
    missing = [c for c in cells if c not in cost.index or c not in cost.columns]
    if missing:
        raise ValueError(f"These cell types are missing in cost: {missing}")
    dis = cost.loc[cells, cells].to_numpy(dtype=float)
    if (dis < 0).any() or not np.allclose(dis, dis.T) or not np.allclose(np.diag(dis), 0):
        raise ValueError("cost must be symmetric, non-negative and 0 on the diagonal.")
    return dis, _as_cost_df(dis, cells)


def compute_cost(cost_method, annot, data, adata=None, metric='cosine', centroid='median',
                 cost=None, layer=None, n_genes=2000, n_markers=100, n_neighbors_cost=15,
                 max_cells=500, random_state=0):
    """Compute the cost matrix with the chosen method; see ``pl.tl.wasserstein_distance``."""
    if cost_method == 'centroid':
        return cost_centroid(annot, data, metric=metric, centroid=centroid)
    if cost_method == 'distribution':
        return cost_distribution(annot, data, metric=metric, max_cells=max_cells,
                                 random_state=random_state)
    if cost_method == 'graph':
        return cost_graph(annot, data, metric=metric, n_neighbors=n_neighbors_cost)
    if cost_method == 'expression':
        return cost_expression(adata, annot, metric=metric, layer=layer, n_genes=n_genes)
    if cost_method == 'jaccard':
        return cost_jaccard(adata, annot, layer=layer, n_markers=n_markers)
    if cost_method == 'precomputed':
        if cost is None:
            raise ValueError("cost_method='precomputed' needs the cost matrix in cost=...")
        return cost_precomputed(annot, cost)
    raise ValueError(f"Unknown cost_method '{cost_method}'. Use one of {list(COST_METHODS)}.")


def transport_distances(proportions, cost, method='sinkhorn_divergence', reg=0.1, reg_m=1.0):
    """
    Distances between all samples with one of the newer optimal transport methods.

    'sinkhorn_divergence': entropic OT, corrected so that the distance of a
    sample to itself is 0. 'unbalanced': the proportions do not have to be
    matched exactly; mass that is not moved is penalised with strength
    ``reg_m`` (smaller values allow more mass to stay, so rare cell types
    matter less).
    """
    samples = list(proportions.keys())
    n = len(samples)
    P = [np.asarray(proportions[s], dtype=float) for s in samples]
    M = np.asarray(cost, dtype=float)
    EMD = np.zeros((n, n))
    if method == 'sinkhorn_divergence':
        def entropic(a, b):
            return float(ot.sinkhorn2(a, b, M, reg, method='sinkhorn_stabilized'))
        self_cost = [entropic(p, p) for p in P]
        for i in range(n):
            for j in range(i + 1, n):
                value = entropic(P[i], P[j]) - 0.5 * (self_cost[i] + self_cost[j])
                EMD[i, j] = EMD[j, i] = max(value, 0.0)
    elif method == 'unbalanced':
        # Total cost: transport cost plus the penalty for the mass that is not moved.
        def unbalanced(a, b):
            try:
                return float(ot.unbalanced.mm_unbalanced2(a, b, M, reg_m, div='kl',
                                                          returnCost='total'))
            except TypeError:  # older POT versions
                G = ot.unbalanced.mm_unbalanced(a, b, M, reg_m, div='kl')
                g1, g2 = G.sum(axis=1), G.sum(axis=0)

                def kl(x, y):
                    x = np.maximum(x, 1e-300)
                    return float(np.sum(x * np.log(x / y)) - x.sum() + y.sum())
                return float(np.sum(G * M) + reg_m * (kl(g1, a) + kl(g2, b)))
        for i in range(n):
            for j in range(i + 1, n):
                EMD[i, j] = EMD[j, i] = max(unbalanced(P[i], P[j]), 0.0)
    else:
        raise ValueError(f"Unknown method '{method}'.")
    emd = pd.DataFrame(EMD, columns=samples)
    emd['sampleID'] = samples
    return EMD, emd.set_index('sampleID')
