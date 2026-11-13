import os
import warnings
from pathlib import Path

# Data Science
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from sklearn import metrics
from statsmodels.stats.multitest import multipletests

# Single-Cell/Bio
import anndata
import elpigraph 
import scanpy as sc
from gprofiler import GProfiler
from pydiffmap import diffusion_map

# Plotting/UI
import matplotlib
import matplotlib.ticker
import seaborn as sns
import textwrap as tw
from adjustText import adjust_text
from cycler import cycler
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.image import imread
from matplotlib.lines import Line2D

from ..tools.Gene_cluster_specific_functions import *

# Configurations
warnings.filterwarnings('ignore')


# ---------------------------------------------------------------------------
# Shared defaults and helpers for the sample-level plots
# ---------------------------------------------------------------------------

#: Default size (inches) of the square sample-level plots (trajectory,
#: principal graph and UMAP), so that they all have the same size.
DEFAULT_FIGSIZE = (8, 7)

#: Default base font size of the sample-level plots.
DEFAULT_FONTSIZE = 14

#: Default colors for the categories of a sample-level variable. The first
#: three are the colors used by earlier PILOT versions.
DEFAULT_COLORS = ['#377eb8', '#ff7f00', '#e41a1c', '#4daf4a', '#984ea3',
                  '#a65628', '#f781bf', '#999999', '#dede00']

#: Status names that are recognised as the healthy / control group.
CONTROL_NAMES = ('control', 'normal', 'healthy', 'ctrl', 'healthy_control',
                 'healthy control', 'non-diseased')


def _array_hash(array):
    """Return a short fingerprint of a numeric array, used to detect changes."""
    import hashlib
    array = np.ascontiguousarray(np.asarray(array, dtype=float))
    return hashlib.sha1(array.tobytes() + str(array.shape).encode()).hexdigest()


def _sample_order(adata):
    """Sample IDs in the row order of ``adata.uns['EMD']``."""
    if 'proportions' not in adata.uns:
        raise ValueError("No PILOT results found. Run pl.tl.wasserstein_distance(adata, ...) first.")
    return list(adata.uns['proportions'].keys())


def _sample_table(adata, sample_col=1):
    """
    One row per sample, in the row order of ``adata.uns['EMD']``.

    The table contains the columns of ``adata.uns['annot']`` and, if
    available, the sample-level metadata stored by ``wasserstein_distance``.
    """
    order = _sample_order(adata)
    annot = adata.uns['annot']
    sample_name = annot.columns[sample_col] if isinstance(sample_col, int) else sample_col
    df = annot.drop_duplicates(subset=[sample_name]).set_index(sample_name)
    missing = [s for s in order if s not in df.index]
    if missing:
        raise ValueError(f"Samples missing from adata.uns['annot']: {missing}")
    df = df.loc[order]
    meta = adata.uns.get('sample_metadata', None)
    if isinstance(meta, pd.DataFrame):
        extra = [c for c in meta.columns if c not in df.columns]
        if extra:
            df = df.join(meta[extra].reindex(order))
    df.index.name = 'sampleID'
    return df.reset_index()


def _find_control(categories, root_status=None):
    """Return the category used as control/root, or None if not found."""
    categories = [c for c in categories]
    if root_status is not None:
        for c in categories:
            if str(c) == str(root_status):
                return c
        raise ValueError(f"root_status='{root_status}' was not found. "
                         f"Available values: {categories}")
    lower = {str(c).strip().lower(): c for c in categories}
    for name in CONTROL_NAMES:
        if name in lower:
            return lower[name]
    return None


def _category_order(values, legend_order=None):
    """Order of the categories: user order first, then order of appearance."""
    present = list(pd.unique(pd.Series(values).astype(object)))
    if legend_order is None:
        return present
    legend_order = list(legend_order)
    unknown = [c for c in legend_order if c not in present]
    if unknown:
        raise ValueError(f"legend_order contains values that are not in the data: {unknown}. "
                         f"Available values: {present}")
    return legend_order + [c for c in present if c not in legend_order]


def _default_palette(n):
    """Return n distinct colors."""
    if n <= len(DEFAULT_COLORS):
        return DEFAULT_COLORS[:n]
    return [matplotlib.colors.to_hex(c) for c in sns.color_palette('husl', n)]


def _get_colors(adata, clusters, categories, colors=None):
    """
    Return a {category: color} dict for ``clusters``.

    The colors are taken, in this order, from the ``colors`` argument (list or
    dict), from the colors saved by an earlier plot, or generated
    automatically. The result is saved in ``adata.uns['sample_colors']`` so
    that the trajectory, the principal graph and the UMAP use the same colors.
    """
    saved = adata.uns.get('sample_colors', {})
    if not isinstance(saved, dict):
        saved = {}
    saved = dict(saved)
    previous = dict(saved.get(clusters, {}))

    if colors is None:
        mapping = {str(k): v for k, v in previous.items()}
        new = [c for c in categories if str(c) not in mapping]
        if new:
            used = set(mapping.values())
            pool = [c for c in _default_palette(len(categories) + len(used)) if c not in used]
            for c, col in zip(new, pool):
                mapping[str(c)] = col
    elif isinstance(colors, dict):
        mapping = {str(k): v for k, v in previous.items()}
        mapping.update({str(k): v for k, v in colors.items()})
        missing = [c for c in categories if str(c) not in mapping]
        if missing:
            raise ValueError(f"No color given for: {missing}")
    else:
        colors = list(colors)
        if len(colors) == 0:
            raise ValueError("colors must not be empty.")
        mapping = {str(c): colors[i % len(colors)] for i, c in enumerate(categories)}

    saved[clusters] = {str(k): v for k, v in mapping.items()}
    adata.uns['sample_colors'] = saved
    return {c: mapping[str(c)] for c in categories}


def _scatter_samples(adata,
                     coords,
                     clusters='status',
                     sample_col=1,
                     colors=None,
                     legend_order=None,
                     legend_loc='best',
                     legend_bbox=None,
                     legend_title=None,
                     legend_fontsize=None,
                     show_legend=True,
                     label_act=False,
                     label_fontsize=None,
                     point_size=None,
                     alpha=1,
                     figsize=None,
                     font_size=None,
                     title=None,
                     title_fontsize=None,
                     xlabel=None,
                     ylabel=None,
                     facecolor='white',
                     axes_color='black',
                     axes_line_width=1,
                     ax=None):
    """Scatter plot of the samples colored by a sample-level variable."""
    df = _sample_table(adata, sample_col)
    if clusters not in df.columns:
        raise ValueError(f"'{clusters}' is not a sample-level column. Available columns: "
                         f"{[c for c in df.columns if c != 'sampleID']}")
    if len(df) != coords.shape[0]:
        raise ValueError(f"{len(df)} samples vs {coords.shape[0]} coordinates")

    font_size = DEFAULT_FONTSIZE if font_size is None else font_size
    legend_fontsize = font_size if legend_fontsize is None else legend_fontsize
    label_fontsize = max(font_size - 4, 6) if label_fontsize is None else label_fontsize
    title_fontsize = font_size + 2 if title_fontsize is None else title_fontsize
    figsize = DEFAULT_FIGSIZE if figsize is None else figsize
    if point_size is None:
        point_size = _auto_point_size(len(df))

    values = df[clusters].astype(object).values
    saved_orders = adata.uns.get('sample_legend_order', {})
    saved_orders = dict(saved_orders) if isinstance(saved_orders, dict) else {}
    if legend_order is None and clusters in saved_orders:
        legend_order = [c for c in pd.unique(pd.Series(values).astype(object))
                        if str(c) in [str(x) for x in saved_orders[clusters]]]
        legend_order = sorted(legend_order, key=lambda c: [str(x) for x in saved_orders[clusters]].index(str(c)))
    elif legend_order is not None:
        saved_orders[clusters] = [str(c) for c in legend_order]
        adata.uns['sample_legend_order'] = saved_orders
    categories = _category_order(values, legend_order)
    color_map = _get_colors(adata, clusters, categories, colors)

    rc = {'font.size': font_size, 'axes.edgecolor': axes_color,
          'axes.linewidth': axes_line_width}
    with plt.rc_context(rc):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure
        ax.set_facecolor(facecolor)
        for category in categories:
            mask = values == category
            ax.scatter(coords[mask, 0], coords[mask, 1], color=color_map[category],
                       alpha=alpha, s=point_size, label=str(category),
                       edgecolors='none')
            if label_act:
                for name, x, y in zip(df['sampleID'][mask], coords[mask, 0], coords[mask, 1]):
                    ax.annotate(str(name), (x, y), fontsize=label_fontsize,
                                xytext=(3, 3), textcoords='offset points')
        if show_legend:
            kwargs = dict(loc=legend_loc, fontsize=legend_fontsize, frameon=False,
                          title=legend_title if legend_title is not None else clusters,
                          title_fontsize=legend_fontsize, markerscale=1)
            if legend_bbox is not None:
                kwargs['bbox_to_anchor'] = legend_bbox
            ax.legend(**kwargs)
        if title:
            ax.set_title(title, fontsize=title_fontsize)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        ax.tick_params(labelsize=max(font_size - 2, 6))
    return fig, ax


def _auto_point_size(n_samples):
    """Point size that gets smaller when there are many samples."""
    return float(100 if n_samples <= 50 else max(15, 100 * np.sqrt(50 / n_samples)))


def _save_and_show(fig, save_path, show, dpi=300):
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)


def _drop_keys(adata, keys):
    for key in keys:
        if key in adata.uns:
            del adata.uns[key]


def _compute_trajectory(adata, n_evecs=2, epsilon=1, alpha=0.5, knn=64,
                        recompute=False, root_status=None, verbose=True):
    """
    Compute (or reuse) the diffusion map of the Wasserstein distances.

    The result is cached in ``adata.uns['embedding']``. It is only recomputed
    when the Wasserstein distances or the diffusion map parameters change, or
    when ``recompute=True``. The sign of every axis is fixed, so the same data
    always gives the same picture.
    """
    if 'EMD' not in adata.uns:
        raise ValueError("No Wasserstein distances found. Run pl.tl.wasserstein_distance(adata, ...) first.")
    EMD = np.asarray(adata.uns['EMD'], dtype=float)
    params = {'n_evecs': int(n_evecs), 'epsilon': str(epsilon), 'alpha': float(alpha),
              'knn': int(knn), 'emd_hash': _array_hash(EMD),
              'root_status': 'None' if root_status is None else str(root_status)}
    old = adata.uns.get('trajectory_params', None)
    if (not recompute and 'embedding' in adata.uns and isinstance(old, dict)
            and all(str(old.get(k)) == str(v) for k, v in params.items())):
        return adata.uns['embedding']

    mydmap = diffusion_map.DiffusionMap.from_sklearn(n_evecs=n_evecs, epsilon=epsilon,
                                                     alpha=alpha, k=knn)
    embedding = np.real(np.asarray(mydmap.fit_transform(EMD / EMD.max()), dtype=complex))

    # The eigen solver returns each axis with a random sign. Fix the sign so
    # that the entry with the largest absolute value is positive.
    for j in range(embedding.shape[1]):
        idx = np.argmax(np.abs(embedding[:, j]))
        if embedding[idx, j] < 0:
            embedding[:, j] = -embedding[:, j]

    # Optionally put the control samples on the left side of the first axis.
    df = _sample_table(adata)
    control = _find_control(df['status'].unique(), root_status) if root_status is not None else None
    if control is not None:
        mask = (df['status'] == control).values
        if mask.any() and embedding[mask, 0].mean() > embedding[:, 0].mean():
            embedding[:, 0] = -embedding[:, 0]

    had_graph = 'principal_graph' in adata.uns or 'pseudotime' in adata.uns
    adata.uns['embedding'] = embedding
    adata.uns['trajectory_params'] = params
    # Results that depend on the old embedding are no longer valid.
    _drop_keys(adata, ['principal_graph', 'principal_graph_params', 'pseudotime',
                       'root_node'])
    if had_graph and verbose:
        print("The trajectory changed. Run pl.pl.fit_pricipla_graph(adata, ...) again "
              "to update the pseudotime.")
    return embedding


def _compute_umap(adata, n_neighbors=None, min_dist=0.5, random_state=0, recompute=False):
    """Compute (or reuse) a UMAP of the samples from the Wasserstein distances."""
    import umap
    if 'EMD' not in adata.uns:
        raise ValueError("No Wasserstein distances found. Run pl.tl.wasserstein_distance(adata, ...) first.")
    EMD = np.asarray(adata.uns['EMD'], dtype=float)
    n = EMD.shape[0]
    if n < 4:
        raise ValueError("UMAP needs at least 4 samples.")
    if n_neighbors is None:
        n_neighbors = max(2, min(15, n - 1))
    params = {'n_neighbors': int(n_neighbors), 'min_dist': float(min_dist),
              'random_state': str(random_state), 'emd_hash': _array_hash(EMD)}
    old = adata.uns.get('umap_samples_params', None)
    if (not recompute and 'umap_samples' in adata.uns and isinstance(old, dict)
            and all(str(old.get(k)) == str(v) for k, v in params.items())):
        return adata.uns['umap_samples']
    D = EMD / EMD.max()
    D = (D + D.T) / 2
    np.fill_diagonal(D, 0)
    reducer = umap.UMAP(n_neighbors=int(n_neighbors), min_dist=min_dist, metric='precomputed',
                        random_state=random_state)
    coords = reducer.fit_transform(D)
    adata.uns['umap_samples'] = np.asarray(coords)
    adata.uns['umap_samples_params'] = params
    return adata.uns['umap_samples']


def _graph_edges(pg):
    edges = pg['Edges']
    return np.asarray(edges[0] if isinstance(edges, (tuple, list)) else edges, dtype=int)


def _fit_principal_graph(adata, NumNodes=20, random_state=0, recompute=False):
    """Fit (or reuse) the elastic principal tree on the trajectory embedding."""
    if 'embedding' not in adata.uns:
        raise ValueError("No trajectory found. Run pl.pl.trajectory(adata) first.")
    emb = np.asarray(adata.uns['embedding'], dtype=float)
    params = {'NumNodes': int(NumNodes), 'random_state': str(random_state),
              'emb_hash': _array_hash(emb)}
    old = adata.uns.get('principal_graph_params', None)
    if (not recompute and 'principal_graph' in adata.uns and isinstance(old, dict)
            and all(str(old.get(k)) == str(v) for k, v in params.items())):
        return adata.uns['principal_graph']

    state = np.random.get_state()
    try:
        if random_state is not None:
            np.random.seed(random_state)
        pg = elpigraph.computeElasticPrincipalTree(emb, NumNodes=NumNodes, verbose=False)[0]
    finally:
        np.random.set_state(state)
    graph = {'NodePositions': np.asarray(pg['NodePositions']), 'Edges': _graph_edges(pg)}
    adata.uns['principal_graph'] = graph
    adata.uns['principal_graph_params'] = params
    _drop_keys(adata, ['pseudotime', 'root_node'])
    return graph


def _pseudotime_from(emb, graph, source):
    """Pseudotime of the samples when the graph starts at node ``source``."""
    pg = {'NodePositions': np.asarray(graph['NodePositions']),
          'Edges': [np.asarray(graph['Edges'])]}
    elpigraph.utils.getPseudotime(emb, pg, source=int(source), target=None)
    return np.asarray(pg['pseudotime'], dtype=float)


def _auto_root_node(adata, root_status=None):
    """
    Choose the root of the principal graph automatically.

    Every end node (leaf) of the graph is tried as root. The chosen root is
    the one for which the control samples come first along the pseudotime
    (highest probability that a control sample has a smaller pseudotime than
    a non-control sample). The control group is ``root_status`` or, if not
    given, a status called e.g. 'control', 'normal' or 'healthy'.
    """
    graph = adata.uns['principal_graph']
    nodes = np.asarray(graph['NodePositions'])
    edges = np.asarray(graph['Edges'])
    df = _sample_table(adata)
    control = _find_control(df['status'].unique(), root_status)
    if control is None:
        raise ValueError("Could not find the control group automatically. Set root_status "
                         f"to one of {list(df['status'].unique())}, or give source_node "
                         "as a number.")
    emb = np.asarray(adata.uns['embedding'], dtype=float)
    is_control = (df['status'] == control).values
    degree = np.bincount(edges.ravel(), minlength=nodes.shape[0])
    leaves = np.where(degree == 1)[0]
    candidates = leaves if len(leaves) > 0 else np.arange(nodes.shape[0])
    if is_control.all() or not is_control.any():
        center = emb[is_control][:, :nodes.shape[1]].mean(axis=0) if is_control.any() \
            else emb[:, :nodes.shape[1]].mean(axis=0)
        return int(candidates[np.argmin(np.linalg.norm(nodes[candidates] - center, axis=1))])
    other_center = emb[~is_control][:, :nodes.shape[1]].mean(axis=0)
    best, best_key = None, None
    for leaf in candidates:
        t = _pseudotime_from(emb, graph, leaf)
        tc, to = t[is_control], t[~is_control]
        auc = ((tc[:, None] < to[None, :]).sum() + 0.5 * (tc[:, None] == to[None, :]).sum()) \
            / (len(tc) * len(to))
        # Ties: prefer the end node that is farther from the non-control samples.
        key = (round(float(auc), 10), float(np.linalg.norm(nodes[leaf] - other_center)))
        if best_key is None or key > best_key:
            best, best_key = int(leaf), key
    return best


def _compute_pseudotime(adata, source_node='auto', root_status=None, verbose=True):
    """Compute the pseudotime of each sample along the principal graph."""
    if 'principal_graph' not in adata.uns:
        raise ValueError("No principal graph found. Run pl.pl.fit_pricipla_graph(adata) first.")
    graph = adata.uns['principal_graph']
    n_nodes = np.asarray(graph['NodePositions']).shape[0]
    if isinstance(source_node, str):
        if source_node != 'auto':
            raise ValueError("source_node must be a node number or 'auto'.")
        source_node = _auto_root_node(adata, root_status)
        if verbose:
            print(f"Root of the trajectory: node {source_node} (end node where the control samples come first).")
    source_node = int(source_node)
    if not 0 <= source_node < n_nodes:
        raise ValueError(f"source_node must be between 0 and {n_nodes - 1}.")
    emb = np.asarray(adata.uns['embedding'], dtype=float)
    adata.uns['pseudotime'] = _pseudotime_from(emb, graph, source_node)
    adata.uns['root_node'] = source_node
    return adata.uns['pseudotime']


def _auto_heatmap_size(n_rows, n_cols=None):
    """Heatmap size that grows with the number of rows and columns."""
    n_cols = n_rows if n_cols is None else n_cols
    width = min(max(0.35 * n_cols + 3, 7), 24)
    height = min(max(0.35 * n_rows + 3, 7), 24)
    return (width, height)


def _clustermap(df, title, save_path, figsize=None, col_cluster=True, row_cluster=True,
                cmap='Blues_r', font_scale=1, show=True, dpi=300, **kwargs):
    figsize = _auto_heatmap_size(df.shape[0], df.shape[1]) if figsize is None else figsize
    with sns.plotting_context('notebook', font_scale=font_scale):
        g = sns.clustermap(df, cmap=cmap, figsize=figsize, col_cluster=col_cluster,
                           row_cluster=row_cluster, xticklabels=True, yticklabels=True,
                           **kwargs)
        g.fig.suptitle(title, y=1.02)
    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        g.savefig(save_path, dpi=dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(g.fig)
    return g


# ---------------------------------------------------------------------------
# Trajectory
# ---------------------------------------------------------------------------

def trajectory(adata,
               n_evecs=2,
               epsilon=1,
               alpha=0.5,
               knn=64,
               sample_col=1,
               clusters='status',
               label_act=False,
               colors=None,
               location_labels='best',
               figsize=None,
               font_size=None,
               axes_line_width=1,
               axes_color='black',
               facecolor='white',
               point_size=None,
               cmap='viridis',
               fontsize_legend=None,
               alpha_trans=1,
               plot_titel="Trajectory of the disease progression",
               legend_order=None,
               legend_bbox=None,
               legend_title=None,
               show_legend=True,
               root_status=None,
               recompute=False,
               save_format='pdf',
               show=True):
    """
    Find the trajectory with a diffusion map and plot it.

    The trajectory is computed only once and stored in
    ``adata.uns['embedding']``. Calling this function again (e.g. with other
    colors) reuses it, so the picture and the principal graph do not change.
    It is recomputed only when the Wasserstein distances or the diffusion map
    parameters (``n_evecs``, ``epsilon``, ``alpha``, ``knn``) change, or when
    ``recompute=True``.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    n_evecs : int, optional
        Number of diffusion map dimensions, by default 2.
    epsilon : float or str, optional
        Kernel bandwidth of the diffusion map, a number or a pydiffmap method
        such as 'bgh', by default 1.
    alpha : float, optional
        Normalization parameter of the diffusion map, by default 0.5.
    knn : int, optional
        Number of nearest neighbors used to build the kernel, by default 64.
    sample_col : int, optional
        Position of the sample column in ``adata.uns['annot']``, by default 1.
    clusters : str, optional
        Sample-level variable used for the colors, by default 'status'.
        Any column of ``adata.uns['sample_metadata']`` can be used, e.g. 'sex'.
    label_act : bool, optional
        If True, write the sample IDs next to the points, by default False.
    colors : list, dict or None, optional
        Colors of the groups. A list is used in the legend order, a dict maps
        each group to a color. If None (default), colors are chosen
        automatically from the number of groups. The colors are remembered and
        reused by ``pl.pl.umap_samples`` and ``pl.pl.fit_pricipla_graph``.
    location_labels : str, optional
        Legend position, any matplotlib ``loc`` value, by default 'best'.
    figsize : tuple, optional
        Figure size in inches. If None, ``DEFAULT_FIGSIZE`` (8, 7) is used.
    font_size : int, optional
        Base font size. If None, 14.
    axes_line_width : float, optional
        Line width of the axes frame, by default 1.
    axes_color : str, optional
        Color of the axes frame, by default 'black'.
    facecolor : str, optional
        Background color of the plot, by default 'white'.
    point_size : float, optional
        Size of the sample points. If None (default), it is chosen from the
        number of samples (100 for up to 50 samples, smaller for more).
    cmap : str, optional
        Not used. Kept for backwards compatibility.
    fontsize_legend : int, optional
        Font size of the legend. If None, the base font size.
    alpha_trans : float, optional
        Opacity of the points (0 to 1), by default 1.
    plot_titel : str, optional
        Title of the plot, also used as the file name.
    legend_order : list, optional
        Order of the groups in the legend, e.g. ``['control', 'disease']``.
        Groups not listed are added at the end. The order is remembered and
        reused by the UMAP and principal graph plots.
    legend_bbox : tuple, optional
        ``bbox_to_anchor`` of the legend. Use e.g. ``(1.02, 1)`` with
        ``location_labels='upper left'`` to place it outside the plot.
    legend_title : str, optional
        Legend title. If None, the name of ``clusters``.
    show_legend : bool, optional
        Whether to show the legend, by default True.
    root_status : str, optional
        If given, the samples of this group are placed on the left side of the
        first axis.
    recompute : bool, optional
        Force a new computation of the trajectory, by default False.
    save_format : str, optional
        File format of the saved plot ('pdf', 'png', ...), by default 'pdf'.
        Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.

    Returns
    -------
    None
        The trajectory is stored in ``adata.uns['embedding']`` and the plot is
        saved in 'Results_PILOT/plots'.
    """
    embedding = _compute_trajectory(adata, n_evecs=n_evecs, epsilon=epsilon, alpha=alpha,
                                    knn=knn, recompute=recompute, root_status=root_status)
    trajectory_plotting(adata, sample_col=sample_col, clusters=clusters, label_act=label_act,
                        colors=colors, location_labels=location_labels, figsize=figsize,
                        font_size=font_size, axes_line_width=axes_line_width,
                        axes_color=axes_color, facecolor=facecolor, point_size=point_size,
                        fontsize_legend=fontsize_legend, alpha_trans=alpha_trans,
                        plot_title=plot_titel, legend_order=legend_order,
                        legend_bbox=legend_bbox, legend_title=legend_title,
                        show_legend=show_legend, save_format=save_format, show=show)


def trajectory_calculation(
    adata,
    n_evecs=2,
    epsilon=1,
    alpha=0.5,
    knn=64,
    root_status=None,
    recompute=True):
    """
    Compute the trajectory (diffusion map) without plotting it.

    Use ``trajectory_plotting`` to plot the result. The sign of every axis is
    fixed, so the same data always gives the same picture.

    Parameters
    ----------
    adata : AnnData
        AnnData object with the sample-by-sample Wasserstein distances in
        ``adata.uns['EMD']``.
    n_evecs : int, optional
        Number of diffusion map dimensions, by default 2.
    epsilon : float or str, optional
        Kernel bandwidth, a number or a pydiffmap method such as 'bgh',
        by default 1.
    alpha : float, optional
        Normalization parameter of the diffusion map, by default 0.5.
    knn : int, optional
        Number of nearest neighbors used to build the kernel, by default 64.
    root_status : str, optional
        If given, the samples of this group are placed on the left side of the
        first axis.
    recompute : bool, optional
        If False, an existing trajectory with the same parameters is reused,
        by default True.

    Returns
    -------
    None
        The trajectory is stored in ``adata.uns['embedding']``.
    """
    _compute_trajectory(adata, n_evecs=n_evecs, epsilon=epsilon, alpha=alpha, knn=knn,
                        recompute=recompute, root_status=root_status)


def trajectory_plotting(
    adata,
    sample_col=1,
    clusters="status",
    label_act=False,
    colors=None,
    location_labels="best",
    figsize=None,
    font_size=None,
    axes_line_width=1,
    axes_color="black",
    facecolor="white",
    point_size=None,
    fontsize_legend=None,
    alpha_trans=1,
    plot_path="Results_PILOT/plots",
    plot_title="Trajectory of the disease progression",
    legend_order=None,
    legend_bbox=None,
    legend_title=None,
    show_legend=True,
    save_format="png",
    show=True):
    """
    Plot the trajectory stored in ``adata.uns['embedding']``.

    The trajectory is not recomputed, so the picture stays the same.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``trajectory_calculation`` or ``trajectory``.
    sample_col : int, optional
        Position of the sample column in ``adata.uns['annot']``, by default 1.
    clusters : str, optional
        Sample-level variable used for the colors, by default 'status'.
    label_act : bool, optional
        If True, write the sample IDs next to the points, by default False.
    colors : list, dict or None, optional
        Colors of the groups, see ``trajectory``. If None, chosen
        automatically.
    location_labels : str, optional
        Legend position, by default 'best'.
    figsize : tuple, optional
        Figure size in inches. If None, (8, 7).
    font_size : int, optional
        Base font size. If None, 14.
    axes_line_width : float, optional
        Line width of the axes frame, by default 1.
    axes_color : str, optional
        Color of the axes frame, by default 'black'.
    facecolor : str, optional
        Background color of the plot, by default 'white'.
    point_size : float, optional
        Size of the sample points. If None (default), it is chosen from the
        number of samples (100 for up to 50 samples, smaller for more).
    fontsize_legend : int, optional
        Font size of the legend. If None, the base font size.
    alpha_trans : float, optional
        Opacity of the points, by default 1.
    plot_path : str or Path, optional
        Folder where the plot is saved, by default 'Results_PILOT/plots'.
    plot_title : str, optional
        Title of the plot, also used as the file name.
    legend_order : list, optional
        Order of the groups in the legend.
    legend_bbox : tuple, optional
        ``bbox_to_anchor`` of the legend.
    legend_title : str, optional
        Legend title. If None, the name of ``clusters``.
    show_legend : bool, optional
        Whether to show the legend, by default True.
    save_format : str, optional
        File format, by default 'png'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.
    """
    if 'embedding' not in adata.uns:
        raise ValueError("No trajectory found. Run pl.pl.trajectory(adata) first.")
    embedding = np.asarray(adata.uns['embedding'])
    fig, ax = _scatter_samples(adata, embedding, clusters=clusters, sample_col=sample_col,
                               colors=colors, legend_order=legend_order,
                               legend_loc=location_labels, legend_bbox=legend_bbox,
                               legend_title=legend_title, legend_fontsize=fontsize_legend,
                               show_legend=show_legend, label_act=label_act,
                               point_size=point_size, alpha=alpha_trans, figsize=figsize,
                               font_size=font_size, title=plot_title,
                               xlabel='Diffusion component 1', ylabel='Diffusion component 2',
                               facecolor=facecolor, axes_color=axes_color,
                               axes_line_width=axes_line_width)
    save_path = None if save_format is None else Path(plot_path) / f"{plot_title}.{save_format}"
    _save_and_show(fig, save_path, show)


def umap_samples(adata,
                 clusters='status',
                 colors=None,
                 n_neighbors=None,
                 min_dist=0.5,
                 random_state=0,
                 recompute=False,
                 sample_col=1,
                 label_act=False,
                 legend_order=None,
                 legend_loc='best',
                 legend_bbox=None,
                 legend_title=None,
                 show_legend=True,
                 figsize=None,
                 font_size=None,
                 fontsize_legend=None,
                 point_size=None,
                 alpha=1,
                 title='UMAP of samples',
                 facecolor='white',
                 plot_path='Results_PILOT/plots',
                 save_format='pdf',
                 show=True):
    """
    UMAP of the samples (patients) based on the Wasserstein distances.

    Each point is a sample. By default the colors are the same as in the
    trajectory plot, so both plots can be compared directly. The UMAP is
    computed once and stored in ``adata.uns['umap_samples']``; it is only
    recomputed when the distances or the UMAP parameters change.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    clusters : str, optional
        Sample-level variable used for the colors, by default 'status'.
    colors : list, dict or None, optional
        Colors of the groups. If None, the colors of the trajectory plot are
        reused (or chosen automatically).
    n_neighbors : int, optional
        Number of neighbors of the UMAP. If None, ``min(15, n_samples - 1)``.
    min_dist : float, optional
        Minimum distance between points in the UMAP, by default 0.5.
    random_state : int, optional
        Random seed, so the UMAP is reproducible, by default 0.
    recompute : bool, optional
        Force a new computation of the UMAP, by default False.
    sample_col : int, optional
        Position of the sample column in ``adata.uns['annot']``, by default 1.
    label_act : bool, optional
        If True, write the sample IDs next to the points, by default False.
    legend_order : list, optional
        Order of the groups in the legend.
    legend_loc : str, optional
        Legend position, by default 'best'.
    legend_bbox : tuple, optional
        ``bbox_to_anchor`` of the legend.
    legend_title : str, optional
        Legend title. If None, the name of ``clusters``.
    show_legend : bool, optional
        Whether to show the legend, by default True.
    figsize : tuple, optional
        Figure size in inches. If None, (8, 7), the same as the trajectory.
    font_size : int, optional
        Base font size. If None, 14.
    fontsize_legend : int, optional
        Font size of the legend. If None, the base font size.
    point_size : float, optional
        Size of the sample points. If None (default), it is chosen from the
        number of samples (100 for up to 50 samples, smaller for more).
    alpha : float, optional
        Opacity of the points, by default 1.
    title : str, optional
        Title of the plot, also used as the file name, by default
        'UMAP of samples'.
    facecolor : str, optional
        Background color of the plot, by default 'white'.
    plot_path : str or Path, optional
        Folder where the plot is saved, by default 'Results_PILOT/plots'.
    save_format : str, optional
        File format, by default 'pdf'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.

    Returns
    -------
    None
        The UMAP coordinates are stored in ``adata.uns['umap_samples']``.
    """
    coords = _compute_umap(adata, n_neighbors=n_neighbors, min_dist=min_dist,
                           random_state=random_state, recompute=recompute)
    fig, ax = _scatter_samples(adata, np.asarray(coords), clusters=clusters,
                               sample_col=sample_col, colors=colors,
                               legend_order=legend_order, legend_loc=legend_loc,
                               legend_bbox=legend_bbox, legend_title=legend_title,
                               legend_fontsize=fontsize_legend, show_legend=show_legend,
                               label_act=label_act, point_size=point_size, alpha=alpha,
                               figsize=figsize, font_size=font_size, title=title,
                               xlabel='UMAP 1', ylabel='UMAP 2', facecolor=facecolor)
    save_path = None if save_format is None else Path(plot_path) / f"{title}.{save_format}"
    _save_and_show(fig, save_path, show)


def plot_trajectory_continuous(
    adata,
    obs_key,
    out_file,
    sample_col="sampleID",
    title=None,
    label=None,
    cmap="viridis",
    na_color="lightgrey",
    figsize=(8, 8),
    dpi=300,
    point_size=100,
    alpha=1,
    font_size=16,
    vmin=None,
    vmax=None,):
    """
    Scatter the sample-level embedding coloured by a variable from adata.obs.

    The cell-level column is collapsed to one value per sample and reordered
    to match the embedding rows.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with the embedding in adata.uns['embedding'],
        the sample order in adata.uns['proportions'] and cell-level
        metadata in adata.obs.
    obs_key : str
        Column of adata.obs holding the value to color by. Must be
        constant within each sample. Non-numeric values are treated as missing.
    out_file : str or Path
        Output path; the extension sets the format (e.g. '.png', '.pdf').
        Parent directories are created if they do not exist.
    sample_col : str, optional
        Column of adata.obs identifying the sample, by default 'sampleID'.
    title : str, optional
        Plot title. No title is shown if None (default).
    label : str, optional
        Colorbar label. Uses 'obs_key' if None (default).
    cmap : str or Colormap, optional
        Colormap for the values, by default 'viridis'.
    na_color : str, optional
        Color for samples without a numeric value, by default 'lightgrey'.
    figsize : tuple of float, optional
        Figure width and height in inches, by default (8, 8).
    dpi : int, optional
        Resolution of the saved figure, mainly relevant for raster formats
        such as PNG, by default 300.
    point_size : float, optional
        Marker size, by default 100.
    alpha : float, optional
        Marker opacity between 0 (transparent) and 1 (opaque), by default 1.
    font_size : int, optional
        Base font size for title, tick labels, legend and colorbar label,
        by default 16.
    vmin, vmax : float, optional
        Lower and upper limits of the color scale. Taken from the data
        if None (default).

    """
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    embedding = adata.uns["embedding"]
    sample_order = list(adata.uns["proportions"].keys())

    if len(sample_order) != embedding.shape[0]:
        raise ValueError(f"{len(sample_order)} samples vs {embedding.shape[0]} embedding rows")

    grouped = adata.obs.groupby(sample_col, observed=True)[obs_key]

    n_unique = grouped.nunique(dropna=True)
    if (n_unique > 1).any():
        bad = list(n_unique[n_unique > 1].index)
        raise ValueError(f"{obs_key} is not constant within these samples: {bad}")

    values = grouped.first().reindex(sample_order)

    values = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    na = np.isnan(values)
    if na.all():
        raise ValueError(f"No sample has a numeric value for {obs_key}")

    if na.any():
        na_samples = [s for s, m in zip(sample_order, na) if m]
        print(f"{na.sum()} of {len(values)} samples have no numeric {obs_key}: {na_samples}")

    with plt.rc_context({"font.size": font_size}):
        fig, ax = plt.subplots(figsize=figsize)

        if na.any():
            ax.scatter(embedding[na, 0], embedding[na, 1],
                       c=na_color, s=point_size, alpha=alpha, label="NA")
            ax.legend(loc="best", fontsize=font_size)

        sc = ax.scatter(embedding[~na, 0], embedding[~na, 1],
                        c=values[~na], cmap=cmap, s=point_size,
                        alpha=alpha, vmin=vmin, vmax=vmax)

        cbar = fig.colorbar(sc, ax=ax)
        cbar.set_label(label or obs_key, fontsize=font_size)

        if title:
            ax.set_title(title)

        fig.savefig(out_file, dpi=dpi, bbox_inches="tight")
        plt.show()

    plt.close(fig)


# ---------------------------------------------------------------------------
# Cell type importance
# ---------------------------------------------------------------------------

def _plot_one_cell_type(ax, name, result, x, y, x_curve, point_size=100, alpha=1,
                        cmap='viridis', color_back=None, title_fontsize=14,
                        xtick_step=None, show_ylabel=True, ylabel='Cell proportion',
                        xlabel=None, label_fontsize=12, tick_fontsize=10):
    """Draw the proportions of one cell type and its fitted curve."""
    from ..tools.Cell_gene_selection import generate_feature as _features
    func, stats = result[0], result[1]
    if color_back is not None:
        ax.set_facecolor(color_back)
    ax.scatter(x, y, c=y, cmap=cmap, alpha=alpha, s=point_size)
    X = np.append(np.ones((len(x_curve), 1)), _features(func, x_curve), axis=1)
    ax.plot(x_curve, X @ np.asarray(stats['params']), color='tab:blue')
    pval = float(result[-1])
    title = f"{name}\nP-val: {pval:.3e}"
    ax.set_title(title, fontsize=title_fontsize,
                 fontweight='bold' if pval <= 0.05 else 'normal')
    ax.set_xlim(x_curve.min(), x_curve.max())
    if xtick_step:
        ax.set_xticks(np.arange(x.min(), x.max() + 1e-9, xtick_step))
    elif np.allclose(x, np.round(x)):
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.tick_params(labelsize=tick_fontsize)
    if show_ylabel and ylabel:
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=label_fontsize)


def _cell_importance_data(adata):
    if 'cell_importance' not in adata.uns:
        raise ValueError("No cell importance results found. Run "
                         "pl.tl.cell_importance(adata) or pl.tl.cell_importance_calculation(adata) first.")
    res = adata.uns['cell_importance']
    x = np.asarray(res['x'], dtype=float)
    return res, x, np.linspace(x.min(), x.max(), 100)


def cell_importance_grid(adata,
                         ncols=4,
                         panel_size=(4.5, 3.8),
                         width=None,
                         height=None,
                         xlim=None,
                         point_size=60,
                         color_back=None,
                         fontsize=None,
                         title_fontsize=None,
                         label_fontsize=None,
                         tick_fontsize=None,
                         xlabel=None,
                         ylabel='Cell proportion',
                         alpha=1,
                         cmap_proportions='viridis',
                         save_path=None,
                         show=True,
                         dpi=300):
    """
    Plot the proportion of every cell type along the trajectory in one figure.

    The size of the figure follows the number of cell types, so the panels are
    never squeezed. Run ``pl.tl.cell_importance_calculation`` (or
    ``pl.tl.cell_importance``) first.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.cell_importance_calculation``.
    ncols : int, optional
        Number of panels per row, by default 4.
    panel_size : tuple, optional
        Size (width, height) in inches of one panel, by default (4.5, 3.8).
    width, height : float, optional
        Total figure size. If None, computed from ``panel_size``.
    xlim : float, optional
        Step between the x-axis ticks. If None, chosen automatically.
    point_size : float, optional
        Size of the points, by default 60.
    color_back : str, optional
        Background color of the panels.
    fontsize : int, optional
        Font size of the panel titles (kept for backwards compatibility).
    title_fontsize, label_fontsize, tick_fontsize : int, optional
        Font sizes of the panel titles, axis labels and tick labels.
    xlabel : str, optional
        Label of the x-axis. If None, 'Sample order' or 'Pseudotime'.
    ylabel : str, optional
        Label of the y-axis, by default 'Cell proportion'.
    alpha : float, optional
        Opacity of the points, by default 1.
    cmap_proportions : str, optional
        Colormap of the points, by default 'viridis'.
    save_path : str or Path, optional
        File to save the figure to. If None, the figure is not saved.
    show : bool, optional
        Whether to show the plot, by default True.
    dpi : int, optional
        Resolution of the saved figure, by default 300.
    """
    res, x, x_curve = _cell_importance_data(adata)
    best = dict(sorted(res['sorted_best'].items(), key=lambda item: item[1][-1]))
    names = list(best.keys())
    if len(names) == 0:
        print("No cell type passed the p-value threshold.")
        return
    target = res['target']
    ncols = max(1, min(ncols, len(names)))
    nrows = int(np.ceil(len(names) / ncols))
    if width is None:
        width = panel_size[0] * ncols
    if height is None:
        height = panel_size[1] * nrows
    title_fontsize = title_fontsize or fontsize or 13
    label_fontsize = label_fontsize or 12
    tick_fontsize = tick_fontsize or 10
    if xlabel is None:
        xlabel = 'Pseudotime' if res.get('x_axis') == 'pseudotime' else 'Sample order'

    fig, axes = plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False)
    for i, name in enumerate(names):
        ax = axes[i // ncols][i % ncols]
        _plot_one_cell_type(ax, name, best[name], x, target.loc[name].values.astype(float),
                            x_curve, point_size=point_size, alpha=alpha,
                            cmap=cmap_proportions, color_back=color_back,
                            title_fontsize=title_fontsize, xtick_step=xlim,
                            show_ylabel=(i % ncols == 0), ylabel=ylabel,
                            xlabel=xlabel if i // ncols == nrows - 1 or i + ncols >= len(names) else None,
                            label_fontsize=label_fontsize, tick_fontsize=tick_fontsize)
    for j in range(len(names), nrows * ncols):
        axes[j // ncols][j % ncols].axis('off')
    fig.tight_layout()
    _save_and_show(fig, save_path, show, dpi=dpi)


def cell_importance_plot(
    adata,
    xlim=None,
    p_val=1,
    point_size=100,
    color_back="white",
    fontsize=None,
    alpha=1,
    cmap_proportions="viridis",
    save_as_pdf=False,
    save_path=None,
    figsize_single=(6, 5),
    dpi=300,
    xlabel=None,
    ylabel='Cell proportion',
    label_fontsize=None,
    tick_fontsize=None,
    show=False):
    """
    Save one plot per cell type with its proportion along the trajectory.

    Run ``pl.tl.cell_importance_calculation`` (or ``pl.tl.cell_importance``)
    first. If it has not been run, it is run here with ``p_val``.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.cell_importance_calculation``.
    xlim : float, optional
        Step between the x-axis ticks. If None, chosen automatically.
    p_val : float, optional
        P-value threshold, only used if the calculation has not been run yet,
        by default 1.
    point_size : int, optional
        Size of the points, by default 100.
    color_back : str, optional
        Background color of the plot, by default 'white'.
    fontsize : int, optional
        Font size of the titles. If None, 14.
    alpha : float, optional
        Opacity of the points, by default 1.
    cmap_proportions : str, optional
        Colormap of the points, by default 'viridis'.
    save_as_pdf : bool, optional
        Save as PDF instead of PNG, by default False.
    save_path : str or None, optional
        Output folder. If None, the PILOT results folder is used.
    figsize_single : tuple, optional
        Size (width, height) in inches of each plot, by default (6, 5).
    dpi : int, optional
        Resolution of the saved figures, by default 300.
    xlabel : str, optional
        Label of the x-axis. If None, 'Sample order' or 'Pseudotime'.
    ylabel : str, optional
        Label of the y-axis, by default 'Cell proportion'.
    label_fontsize, tick_fontsize : int, optional
        Font sizes of the axis labels and tick labels.
    show : bool, optional
        Whether to also show the plots, by default False.

    Returns
    -------
    None
        Saves one file per cell type in '<save_path>/Cell_types_importance'.
    """
    import re
    if 'cell_importance' not in adata.uns:
        from ..tools.Trajectory import cell_importance_calculation
        cell_importance_calculation(adata, p_val=p_val)
    res, x, x_curve = _cell_importance_data(adata)
    if save_path is None:
        save_path = adata.uns['cell_importance'].get('results_path', 'Results_PILOT')
    out_dir = Path(save_path) / "Cell_types_importance"
    out_dir.mkdir(parents=True, exist_ok=True)
    ext = "pdf" if save_as_pdf else "png"
    if xlabel is None:
        xlabel = 'Pseudotime' if res.get('x_axis') == 'pseudotime' else 'Sample order'
    for name, result in res['sorted_best'].items():
        fig, ax = plt.subplots(figsize=figsize_single)
        _plot_one_cell_type(ax, name, result, x, res['target'].loc[name].values.astype(float),
                            x_curve, point_size=point_size, alpha=alpha,
                            cmap=cmap_proportions, color_back=color_back,
                            title_fontsize=fontsize or 14, xtick_step=xlim, ylabel=ylabel,
                            xlabel=xlabel, label_fontsize=label_fontsize or 12,
                            tick_fontsize=tick_fontsize or 10)
        fig.tight_layout()
        fname = re.sub(r"[^\w\-.]+", "_", str(name)).strip("_")
        _save_and_show(fig, out_dir / f"{fname}.{ext}", show, dpi=dpi)


# ---------------------------------------------------------------------------
# Heatmaps
# ---------------------------------------------------------------------------

def heatmap_cost(adata, figsize=None, col_cluster=True, row_cluster=True, cmap='Blues_r',
                 font_scale=1, title='Cost matrix', save_format='pdf', show=True):
    """
    Heatmap of the cost matrix (distances between cell types).

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    figsize : tuple, optional
        Figure size in inches. If None, it grows with the number of cell types.
    col_cluster, row_cluster : bool, optional
        Whether to cluster the columns / rows, by default True.
    cmap : str, optional
        Colormap, by default 'Blues_r'.
    font_scale : float, optional
        Scale of all fonts, by default 1.
    title : str, optional
        Title of the plot, also used as the file name, by default 'Cost matrix'.
    save_format : str, optional
        File format, by default 'pdf'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.
    """
    annot = adata.uns['annot']
    cost = adata.uns['cost'][annot.cell_type.unique()]
    save_path = None if save_format is None else Path('Results_PILOT/plots') / f"{title}.{save_format}"
    _clustermap(cost, title, save_path, figsize=figsize, col_cluster=col_cluster,
                row_cluster=row_cluster, cmap=cmap, font_scale=font_scale, show=show)


def heatmap_wasserstein(adata, figsize=None, col_cluster=True, row_cluster=True, cmap='Blues_r',
                        font_scale=1, title='Wasserstein distance', save_format='pdf',
                        show=True):
    """
    Heatmap of the Wasserstein distances between samples.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    figsize : tuple, optional
        Figure size in inches. If None, it grows with the number of samples.
    col_cluster, row_cluster : bool, optional
        Whether to cluster the columns / rows, by default True.
    cmap : str, optional
        Colormap, by default 'Blues_r'.
    font_scale : float, optional
        Scale of all fonts, by default 1.
    title : str, optional
        Title of the plot, also used as the file name, by default
        'Wasserstein distance'.
    save_format : str, optional
        File format, by default 'pdf'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.
    """
    emd = adata.uns['EMD_df']
    save_path = None if save_format is None else Path('Results_PILOT/plots') / f"{title}.{save_format}"
    _clustermap(emd, title, save_path, figsize=figsize, col_cluster=col_cluster,
                row_cluster=row_cluster, cmap=cmap, font_scale=font_scale, show=show)


def heatmap_proportions(adata, figsize=None, col_cluster=True, row_cluster=False,
                        cmap='Blues', font_scale=1,
                        title='Proportions of cell types for samples over trajectory',
                        save_format='pdf', show=True):
    """
    Heatmap of the cell type proportions of the samples, ordered by pseudotime.

    Rows are samples (from the start to the end of the trajectory), columns
    are cell types. Each column is scaled to the range 0 to 1.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.pl.fit_pricipla_graph``.
    figsize : tuple, optional
        Figure size in inches. If None, it grows with the data.
    col_cluster : bool, optional
        Whether to cluster the cell types, by default True.
    row_cluster : bool, optional
        Whether to cluster the samples, by default False (keep the trajectory
        order).
    cmap : str, optional
        Colormap, by default 'Blues'.
    font_scale : float, optional
        Scale of all fonts, by default 1.
    title : str, optional
        Title of the plot, also used as the file name.
    save_format : str, optional
        File format, by default 'pdf'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.
    """
    if 'pseudotime' not in adata.uns:
        raise ValueError("No pseudotime found. Run pl.pl.fit_pricipla_graph(adata) first.")
    order = _sample_order(adata)
    cell_types = adata.uns['annot']['cell_type'].unique()
    props = pd.DataFrame.from_dict(adata.uns['proportions']).T
    props.columns = cell_types
    props.index.name = 'sampleID'
    props = props.loc[order]
    sorted_samples = pd.Series(np.asarray(adata.uns['pseudotime'], dtype=float),
                               index=order).sort_values(kind='stable').index
    df = props.loc[sorted_samples]
    df = (df - df.min()) / (df.max() - df.min())
    save_path = None if save_format is None else Path('Results_PILOT/plots') / f"{title}.{save_format}"
    _clustermap(df[cell_types], title, save_path, figsize=figsize, col_cluster=col_cluster,
                row_cluster=row_cluster, cmap=cmap, font_scale=font_scale, show=show)


def heatmaps(adata, figsize=None, col_cluster=True, row_cluster=True, cmap='Blues_r',
             font_scale=1, show=True):
    """
    Plot the heatmaps of the cost matrix and of the Wasserstein distances.

    This calls ``heatmap_cost`` and ``heatmap_wasserstein``. Use them directly
    to plot only one of the two.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    figsize : tuple, optional
        Figure size in inches. If None, it grows with the data.
    col_cluster, row_cluster : bool, optional
        Whether to cluster the columns / rows, by default True.
    cmap : str, optional
        Colormap, by default 'Blues_r'.
    font_scale : float, optional
        Scale of all fonts, by default 1.
    show : bool, optional
        Whether to show the plots, by default True.

    Returns
    -------
    None
        Saves 'Cost matrix.pdf' and 'Wasserstein distance.pdf' in
        'Results_PILOT/plots'.
    """
    heatmap_cost(adata, figsize=figsize, col_cluster=col_cluster, row_cluster=row_cluster,
                 cmap=cmap, font_scale=font_scale, show=show)
    heatmap_wasserstein(adata, figsize=figsize, col_cluster=col_cluster,
                        row_cluster=row_cluster, cmap=cmap, font_scale=font_scale, show=show)


def heatmaps_df(df, figsize=(12, 12), col_cluster=True, row_cluster=True, cmap='Blues_r',
                show=True):
    """
    Heatmap of the cell type proportions of the samples.

    Parameters
    ----------
    df : DataFrame
        Samples (rows) by cell types (columns).
    figsize : tuple, optional
        Figure size in inches. If None, it grows with the data.
    col_cluster : bool, optional
        Whether to cluster the columns, by default True.
    row_cluster : bool, optional
        Whether to cluster the rows, by default True.
    cmap : str, optional
        Colormap, by default 'Blues_r'.
    show : bool, optional
        Whether to show the plot, by default True.

    Returns
    -------
    None
        Saves the heatmap in 'Results_PILOT/plots'.
    """
    title = 'Proportions of cell types for samples over trajectory'
    _clustermap(df, title, Path('Results_PILOT/plots') / f"{title}.pdf", figsize=figsize,
                col_cluster=col_cluster, row_cluster=row_cluster, cmap=cmap, show=show)


# ---------------------------------------------------------------------------
# Principal graph and pseudotime
# ---------------------------------------------------------------------------

def principal_graph(adata,
                    clusters='status',
                    colors=None,
                    show_text=True,
                    show_samples=True,
                    figsize=None,
                    font_size=None,
                    point_size=None,
                    node_size=60,
                    node_color='k',
                    edge_color='k',
                    root_color='#e41a1c',
                    DimToPlot=(0, 1),
                    facecolor='white',
                    title='Principal graph',
                    legend_order=None,
                    legend_loc='best',
                    legend_bbox=None,
                    show_legend=True,
                    plot_path='Results_PILOT/plots',
                    save_format='pdf',
                    show=True):
    """
    Plot the principal graph on top of the trajectory.

    The samples are colored like in the trajectory plot, and the nodes are
    numbered, which helps to choose ``source_node``. The root node, if
    already chosen, is shown in red.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.fit_principal_graph``.
    clusters : str, optional
        Sample-level variable used for the colors, by default 'status'.
    colors : list, dict or None, optional
        Colors of the groups. If None, the colors of the trajectory are reused.
    show_text : bool, optional
        Whether to write the node numbers, by default True.
    show_samples : bool, optional
        Whether to draw the samples, by default True.
    figsize : tuple, optional
        Figure size in inches. If None, (8, 7), the same as the trajectory.
    font_size : int, optional
        Base font size. If None, 14.
    point_size : float, optional
        Size of the sample points. If None (default), it is chosen from the
        number of samples (100 for up to 50 samples, smaller for more).
    node_size : float, optional
        Size of the graph nodes, by default 60.
    node_color, edge_color : str, optional
        Colors of the nodes and edges, by default 'k'.
    root_color : str, optional
        Color of the root node, by default red.
    DimToPlot : tuple, optional
        Dimensions of the trajectory to plot, by default (0, 1).
    facecolor : str, optional
        Background color, by default 'white'.
    title : str, optional
        Title of the plot, also used as the file name.
    legend_order : list, optional
        Order of the groups in the legend.
    legend_loc : str, optional
        Legend position, by default 'best'.
    legend_bbox : tuple, optional
        ``bbox_to_anchor`` of the legend.
    show_legend : bool, optional
        Whether to show the legend, by default True.
    plot_path : str or Path, optional
        Folder where the plot is saved, by default 'Results_PILOT/plots'.
    save_format : str, optional
        File format, by default 'pdf'. Use None to not save.
    show : bool, optional
        Whether to show the plot, by default True.
    """
    if 'principal_graph' not in adata.uns:
        raise ValueError("No principal graph found. Run pl.tl.fit_principal_graph(adata) first.")
    d0, d1 = DimToPlot
    graph = adata.uns['principal_graph']
    nodes = np.asarray(graph['NodePositions'])
    edges = np.asarray(graph['Edges'])
    emb = np.asarray(adata.uns['embedding'])
    font_size = DEFAULT_FONTSIZE if font_size is None else font_size
    if show_samples:
        fig, ax = _scatter_samples(adata, emb[:, [d0, d1]], clusters=clusters, colors=colors,
                                   legend_order=legend_order, legend_loc=legend_loc,
                                   legend_bbox=legend_bbox, show_legend=show_legend,
                                   point_size=point_size, alpha=0.8, figsize=figsize,
                                   font_size=font_size, title=title,
                                   xlabel=f'Diffusion component {d0 + 1}',
                                   ylabel=f'Diffusion component {d1 + 1}', facecolor=facecolor)
    else:
        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE if figsize is None else figsize)
        ax.set_facecolor(facecolor)
        ax.set_title(title, fontsize=font_size + 2)
    for a, b in edges:
        ax.plot(nodes[[a, b], d0], nodes[[a, b], d1], color=edge_color, lw=1.5, zorder=3)
    ax.scatter(nodes[:, d0], nodes[:, d1], color=node_color, s=node_size, zorder=4)
    root = adata.uns.get('root_node', None)
    if root is not None:
        ax.scatter(nodes[root, d0], nodes[root, d1], color=root_color, s=node_size * 2.5,
                   zorder=5, marker='*', edgecolors='k', linewidths=0.5)
    if show_text:
        for i, (x, y) in enumerate(nodes[:, [d0, d1]]):
            ax.annotate(str(i), (x, y), fontsize=max(font_size - 3, 6), zorder=6,
                        xytext=(4, 4), textcoords='offset points')
    save_path = None if save_format is None else Path(plot_path) / f"{title}.{save_format}"
    _save_and_show(fig, save_path, show)


def fit_pricipla_graph(adata,
                       NumNodes=20,
                       source_node=0,
                       show_text=True,
                       Do_PCA=False,
                       figsize=None,
                       X_color='r',
                       Node_color='k',
                       DimToPlot=[0, 1],
                       facecolor='white',
                       title='Principal graph',
                       root_status=None,
                       random_state=0,
                       clusters='status',
                       colors=None,
                       recompute=False,
                       show=True):
    """
    Fit an elastic principal graph to the trajectory and compute the pseudotime.

    This runs three steps, which can also be run one by one:
    ``pl.tl.fit_principal_graph`` (fit the graph), ``pl.pl.principal_graph``
    (plot it) and ``pl.tl.compute_pseudotime`` (choose the root and compute
    the pseudotime). The graph is fitted with a fixed random seed, so the node
    numbers stay the same when the function is run again in the same session.
    In a new Python session elpigraph can give a slightly different graph
    (small numerical differences), so a fixed ``source_node`` number may point
    to another place; ``source_node='auto'`` avoids this. Always check the plot.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.pl.trajectory``.
    NumNodes : int, optional
        Number of nodes of the graph, by default 20.
    source_node : int or 'auto', optional
        Root node of the pseudotime, by default 0. Use 'auto' to choose the end
        node from which the control samples come first along the trajectory
        (see ``root_status``).
    show_text : bool, optional
        Whether to show the node numbers, by default True.
    Do_PCA : bool, optional
        Not used. Kept for backwards compatibility.
    figsize : tuple, optional
        Figure size in inches. If None, (8, 7), the same as the trajectory.
    X_color : str, optional
        Not used, the samples have the colors of the trajectory. Kept for
        backwards compatibility.
    Node_color : str, optional
        Color of the graph nodes, by default 'k'.
    DimToPlot : list, optional
        Dimensions of the trajectory to plot, by default [0, 1].
    facecolor : str, optional
        Background color, by default 'white'.
    title : str, optional
        Title of the plot, by default 'Principal graph'.
    root_status : str, optional
        Control group used when ``source_node='auto'``. If None, a status
        called e.g. 'control', 'normal' or 'healthy' is used.
    random_state : int, optional
        Random seed of the graph fit, by default 0.
    clusters : str, optional
        Sample-level variable used for the sample colors, by default 'status'.
    colors : list, dict or None, optional
        Colors of the groups. If None, the colors of the trajectory are reused.
    recompute : bool, optional
        Force a new fit of the graph, by default False.
    show : bool, optional
        Whether to show the plot, by default True.

    Returns
    -------
    None
        The graph is stored in ``adata.uns['principal_graph']``, the root in
        ``adata.uns['root_node']`` and the pseudotime of each sample in
        ``adata.uns['pseudotime']``.
    """
    _fit_principal_graph(adata, NumNodes=NumNodes, random_state=random_state,
                         recompute=recompute)
    _compute_pseudotime(adata, source_node=source_node, root_status=root_status)
    principal_graph(adata, clusters=clusters, colors=colors, show_text=show_text,
                    figsize=figsize, node_color=Node_color, DimToPlot=DimToPlot,
                    facecolor=facecolor, title=title, show=show)


def clustering_emd(
    adata,
    res=0.3,
    metric="cosine",
    groupby_col="Leiden",
    swap_axes=False,
    cmap="Blues_r",
    dendrogram=True,
    show_gene_labels=True,
    var_group_rotation=45,
    figsize=[12,12],
    save=False,
    sorter_leiden=None,
    use_rep="X",
    order_file=None):
    """
    Perform clustering and visualization of EMD (Earth Mover's Distance) data in AnnData object.
    
    Parameters
    ----------
    adata : AnnData
        Input AnnData object containing EMD data.
    res : float, optional
        Resolution parameter for Leiden clustering. Default is 0.3.
    metric : str, optional
        Distance metric for clustering. Default is "cosine".
    groupby_col : str, optional
        Grouping variable for plotting. "Leiden" groups by predicted clusters, "status" groups by real labels. Default is "Leiden".
    swap_axes : bool, optional
        Swap the axes in the heatmap. Default is False.
    cmap : str, optional
        Colormap for the heatmap. Default is "Blues_r".
    dendrogram : bool, optional
        Display dendrograms in the heatmap. Default is True.
    show_gene_labels : bool, optional
        Show gene labels in the heatmap. Default is True.
    var_group_rotation : int, optional
        Rotation angle for gene labels. Default is 45 degrees.
    figsize : list, optional
        Size of the heatmap figure. Default is [12, 12].
    save : bool, optional
        Save the heatmap figure. Default is False.
    sorter_leiden : list or None, optional
        Custom order for Leiden clusters. If not provided, the default order is used.
    use_rep: str or None, data passed to sc.pp.neighbors and sc.tl.dendrogram.
            "X" (default) clusters on adata.X as is.
            "PCA" computes a PCA first and clusters on the principal components.
             None mirrors scanpy's automatic choice: adata.X if it has at most
             50 columns, calculates the PCA otherwise. The 50 columns (features)
             comes from sc.settings.N_PCS.
    order_file: Path where the row/column order should be saved. Default: None
    
    Returns
    -------
    proportion_df : pd.DataFrame
        DataFrame containing proportions of sub-clusters in each sample.
    """
    
    EMD = adata.uns["EMD"]
    proportions = adata.uns["proportions"]
    annot = adata.uns["annot"]
    real_labels = adata.uns["real_labels"]
    
    EMD_df = pd.DataFrame(EMD,columns=proportions.keys())
    EMD_df["sampleID"] = proportions.keys()
    EMD_df["status"] = list(real_labels)
    
    # Get the dimensions of the EMD matrix.
    n_vars = EMD.shape[1]
    # Checking what we should do with the data.
    if use_rep == "X":
        # We use it as is. No PCA.
        rep = "X"
    elif use_rep == "PCA":
        # We calculate the PCA irrespective of the number of columns.
        rep = "X_pca"
    elif use_rep is None:
        # We calculate the PCA if the number of columns is >50, otherwise we
        # use the EMD data as is. This is the original behaviour.
        rep = "X" if n_vars <= 50 else "X_pca"
    else:
        raise ValueError(f"use_rep must be 'X', 'PCA' or None. Given: {use_rep}")
    
    adata_emd = sc.AnnData(EMD)
    if rep == "X_pca":
        sc.pp.pca(adata_emd)
    
    sc.pp.neighbors(adata_emd, metric=metric, use_rep=rep)
    sc.tl.leiden(adata_emd, resolution=res)
    predicted_labels = np.array(adata_emd.obs.leiden)
    
    proportion_df = pd.DataFrame(proportions)
    proportion_df = proportion_df.T
    proportion_df.columns = annot.cell_type.unique()
    
    proportion_df["Predicted_Labels"] = predicted_labels
    proportion_df["sampleID"] = list(proportions.keys())
    
    EMD_df["Leiden"] = predicted_labels
    if groupby_col == "status":
        sorter=np.unique(real_labels)
        EMD_df["status"] = EMD_df.status.astype("category")
        EMD_df["status"] = EMD_df["status"].cat.set_categories(sorter)
        EMD_df = EMD_df.sort_values(["status"])
    elif groupby_col == "Leiden":
        if sorter_leiden is None:
            sorter = EMD_df.Leiden.unique()
        else:
            sorter=sorter_leiden
        
        EMD_df["Leiden"] = EMD_df.Leiden.astype("category")
        EMD_df["Leiden"] = EMD_df["Leiden"].cat.set_categories(sorter)
        EMD_df = EMD_df.sort_values(["Leiden"])
    
    obs = pd.DataFrame()
    obs["sampleID"] = EMD_df.sampleID.astype(str)
    obs["status"] = EMD_df.status.astype(str)
    obs["Leiden"] = EMD_df.Leiden.astype(str)
    df_genes = pd.DataFrame(index = EMD_df.columns[0:EMD_df.shape[0]])
    adata_emd = anndata.AnnData(X = EMD_df[ EMD_df.columns[0:EMD_df.shape[0]]].values, var=df_genes, obs=obs )
    adata_emd.obs[groupby_col] = adata_emd.obs[groupby_col].astype("category")
    if dendrogram:
        if rep == "X_pca":
            sc.pp.pca(adata_emd)
        
        sc.tl.dendrogram(adata_emd, groupby=groupby_col, use_rep=rep)
    
    if order_file is not None:
        col_order = list(adata_emd.obs.sampleID)
        if dendrogram:
            cats = list(adata_emd.uns[f"dendrogram_{groupby_col}"]["categories_ordered"])
        else:
            cats = list(adata_emd.obs[groupby_col].cat.categories)
        row_order = [s for c in cats
                     for s in adata_emd.obs.sampleID[adata_emd.obs[groupby_col] == c]]
        pd.DataFrame({"row_order": row_order, "col_order": col_order}).to_csv(order_file, index=False)
    
    sc.pl.heatmap(adata_emd, adata_emd.obs.sampleID,groupby=[groupby_col],swap_axes=swap_axes,cmap=cmap,dendrogram=dendrogram,show_gene_labels=show_gene_labels,var_group_rotation=var_group_rotation,figsize=figsize,save=save)
    return proportion_df


def Sil_computing(EMD, real_labels, metric='cosine'):
    """
    Compute the Silhouette score based on Wasserstein distances.

    Parameters
    ----------
    EMD : numpy.ndarray
        Wasserstein distances matrix.
    real_labels : list or numpy.ndarray
        True labels or ground truth.
    metric : str, optional
        Metric for calculating pairwise distances, by default 'cosine'.

    Returns
    -------
    float
        Silhouette score indicating cluster quality.
    """
    Silhouette = metrics.silhouette_score(EMD, real_labels, metric =metric)
    #print("Silhouette score: ", Silhouette) 
    return Silhouette


def select_best_sil(adata,marker="o",figsize=(10,10),facecolor="white",metric="cosine",path=None,start=0.2,step=0.1,end=2,use_rep="X"):
    """
    Parameters
    ----------
    adata : adata,
    marker : str, optional (default="o")
            Marker style for data points in the plot.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (10, 10).
    facecolor : str, optional
        Background color of the saved plot image. Default is "white".
    metric: str, metric for leiden clustering and calculating sil.
    path: str, the location were the images should be saved.
    start: float or int, starting resolution for the leiden algorithm. Default 0.2.
    step: float or int, step resolution for the leiden algorithm. Default 0.1.
    end: float or int, final resolution for the leiden algorithm. Default 2.
    use_rep: str or None, data passed to sc.pp.neighbors.
            "X" (default) clusters on adata.X as is.
            "PCA" computes a PCA first and clusters on the principal components.
             None mirrors scanpy's automatic choice: adata.X if it has at most
             50 columns, calculates the PCA otherwise. The 50 columns (features)
             comes from sc.settings.N_PCS.

    Returns
    -------
    None,
    plot Silhouette Score vs. Resolution to figure out the best sil.
    plot Silhouette Score vs. Number of clusters.
    """

    resolutions = np.linspace(start, end, int(round((end - start) / step)) + 1)
    # Create a list to store the Silhouette Scores
    best_res=start
    best_sil=float("-inf")
    if path is None:
        result_path=Path("Results_PILOT/plots")
        result_path.mkdir(parents=True, exist_ok=True)
    else:
        result_path=Path(path)

    sil_scores = []
    number_of_clusters=[]
    EMD=adata.uns["EMD"]

    # Get the dimensions of the EMD matrix.
    n_vars = EMD.shape[1]
    # Checking what we should do with the data.
    if use_rep == "X":
        # We use it as is. No PCA.
        rep = "X"
    elif use_rep == "PCA":
        # We calculate the PCA irrespective of the number of columns.
        rep = "X_pca"
    elif use_rep is None:
        # We calculate the PCA if the number of columns is >50, otherwise we
        # use the EMD data as is. This is the original behaviour.
        rep = "X" if n_vars <= 50 else "X_pca"
    else:
        raise ValueError(f"use_rep must be 'X', 'PCA' or None. Given: {use_rep}")

    adata_emd = sc.AnnData(EMD)
    if rep == "X_pca":
        sc.pp.pca(adata_emd)

    sc.pp.neighbors(adata_emd, metric=metric, use_rep=rep)

    # Define the range of resolutions/number of clusters to test
    # Loop through the resolutions and calculate the Silhouette Score for each
    for resolution in resolutions:
        # Cluster the data using Leiden clustering at the specified resolution
        sc.tl.leiden(adata_emd, resolution = resolution)
        # Calculate the Silhouette Score
        predicted_labels = np.array(adata_emd.obs.leiden)

        Silhouette = metrics.silhouette_score(EMD, predicted_labels, metric=metric)
        if float(Silhouette) > best_sil:
            best_sil=Silhouette
            best_res=resolution

        # Append the Silhouette Score to the list
        sil_scores.append(Silhouette)
        number_of_clusters.append(len(np.unique(predicted_labels)))

    adata.uns["best_res"]=float(best_res)
    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, sil_scores, marker=marker)
    plt.xlabel("Resolution")
    plt.ylabel("Silhouette Score")
    plt.title("Silhouette Score vs. Resolution")
    plt.savefig(result_path / "Silhouette_score_VS_resolution.pdf")
    plt.show()

    # Plot the Silhouette Scores against the resolutions
    plt.figure(figsize = figsize,facecolor=facecolor)
    plt.plot(resolutions, number_of_clusters, marker=marker)
    plt.xlabel("Resolution")
    plt.ylabel("Number of Clusters")
    plt.title("Number of Clusters vs. Resolution")
    plt.savefig(result_path / "Number_of_clusters_VS_resolution.pdf")
    plt.show()


def cell_type_diff_two_sub_patient_groups(
    proportions: pd.DataFrame = None,
    cell_types: list = None,
    labels: str = "Predicted_Labels",
    group1: str = "Tumor 1",
    group2: str = "Tumor 2",
    pval_thr: float = 0.05,
    figsize: tuple = (15, 4),
    file_path: str = None,fontsize:int = 24):
    if file_path == None:
        plot_path = Path("Results_PILOT/plots")
        plot_path.mkdir(parents=True, exist_ok=True)
    else:
        plot_path=Path(file_path)
    
    group1_proportions = proportions[proportions[labels] == group1]
    group2_proportions = proportions[proportions[labels] == group2]
    
    stats_group12 = []
    scores_group12 = []
    # using Welch’s t-test to detect the statistically significant difference
    # between group 1 and group 2
    for cell_type in cell_types:
        ttest_result = ttest_ind(group1_proportions[cell_type], group2_proportions[cell_type], equal_var=False)
        scores_group12.append(ttest_result[0])
        stats_group12.append(ttest_result[1])
    
    # adjust the p-values
    stats_group12 = multipletests(list(stats_group12), method="fdr_bh")[1]
    
    stats_bar_group12 = pd.DataFrame()
    stats_bar_group12["cell_type"] = cell_types
    stats_bar_group12["adjPval"] = stats_group12
    stats_bar_group12["-logPval"] = -np.log(stats_group12)
    stats_bar_group12["score"] = scores_group12
    
    # sort the data based on p-value and statistic score
    stats_bar_group12 = stats_bar_group12.sort_values(by=["score", "-logPval"], ascending=[False, False])
    
    # save the table of statistical differentiation of two groups
    # as positive score shows cell types statistically significant on the first group
    # and negative scrore shows cell types statistically significant on the second group
    
    diff_result_path = Path("Results_PILOT/Diff_Expressions_Results")
    diff_result_path.mkdir(parents=True, exist_ok=True)
    stats_bar_group12.to_csv(diff_result_path / f"Celltype_diff_{group1}_vs_{group2}.csv", header=True, index=None)
    
    # filter data based on a p-value threshold
    stats_bar_group12 = stats_bar_group12[stats_bar_group12["adjPval"] < pval_thr]
    
    fig, ax = plt.subplots(figsize = figsize)
    pl.pl.plot_hor_vs_vert(stats_bar_group12, 1, x="score", y="cell_type", c="type",
                           xlabel="statistic score", ylabel=None,
                           rotation=None, tick_bottom=True, tick_left=False,
                           title=f"Celltype rank {group1} vs {group2}", fontsize=fontsize)
    fig.tight_layout()
    plt.savefig(plot_path / f"Celltype_diff_{group1}_VS_{group2}.png", facecolor="white")
    return stats_bar_group12


def plot_numeric_with_trajectory(
    adata,
    feature,
    sample_col="sampleID",
    ax=None,
    figsize=(7, 6),
    annotate_stats=True,
    fit_line=True,
    fit_color="firebrick",
    point_color="steelblue",
    point_size=60,
    alpha=1,
    font_size=14,
    title=None,
    xlabel="Time score",
    ylabel=None,
    pos_text_x=0.03,
    pos_text_y=0.97):
    """
    Scatter one numeric variable against the trajectory Time_score.
    
    Uses the same preparation as correlation_numeric_with_trajectory, so the
    points shown are exactly the samples that entered the correlation.
    
    Parameters:
    -- adata (AnnData): Input AnnData object containing trajectory and clinical data.
    -- feature (str): Single numeric variable to plot. Not a list.
    -- sample_col (str): Name of the column in adata.obs that contains sample IDs.
    -- ax (matplotlib Axes or None): Draw into an existing axes. A new figure is
       created when None.
    -- annotate_stats (bool): Print Spearman rho, p-value and n in the corner.
    -- fit_line (bool): Overlay a least-squares line. Purely descriptive; the
       reported statistic is a rank correlation, not a linear fit.
    -- title, xlabel, ylabel (str or None): Axis text. Defaults derive from the
       feature name.
    
    Returns:
    -- fig, ax: The figure and axes, for further modification and saving.
       Nothing is written to disk.
    """
    if not isinstance(feature, str):
        raise TypeError(f"feature must be a single column name, got {type(feature).__name__}.")
    
    missing = {sample_col, feature} - set(adata.obs.columns)
    if missing:
        raise KeyError(f"adata.obs is missing required columns: {missing}")
    
    if "orders" not in adata.uns:
        raise KeyError("adata.uns is missing 'orders'. Run the cell importance step first.")
    
    orders = adata.uns["orders"]
    if "Time_score" not in orders.columns:
        raise KeyError("adata.uns['orders'] is missing the Time_score column.")
    
    clinical_data = adata.obs[[sample_col, feature]].rename(columns={sample_col: "sampleID"})
    # Casting out of Categorical, otherwise pd.to_numeric below fails.
    clinical_data[feature] = clinical_data[feature].astype(object)
    
    merged_data = orders.merge(clinical_data, on="sampleID", how="left")
    
    # Detect and remove non-numeric values.
    values_to_remove = detect_values_to_remove(merged_data[feature])
    merged_data = merged_data[~merged_data[feature].isin(set(values_to_remove))].copy()
    
    # Ensuring that the feature is actually numeric.
    merged_data[feature] = pd.to_numeric(merged_data[feature], errors="coerce")
    
    # Dropping missing values before de-duplicating, so a sample is not lost
    # just because its first cell happened to carry the NaN.
    merged_data = merged_data.dropna(subset=["Time_score", feature])
    unique_merged_data = merged_data.drop_duplicates(subset="sampleID")
    
    time_score = unique_merged_data["Time_score"].to_numpy(dtype=float)
    test_feature = unique_merged_data[feature].to_numpy(dtype=float)
    number_samples = len(unique_merged_data)
    
    if number_samples == 0:
        raise ValueError(
            f"No samples left for '{feature}' after removing non-numeric and "
            f"missing values. Nothing to plot."
        )
    
    spearman_corr = spearman_p_value = np.nan
    enough_samples = number_samples >= 3
    has_variation = np.ptp(time_score) > 0 and np.ptp(test_feature) > 0
    if enough_samples and has_variation:
        spearman_corr, spearman_p_value = spearmanr(time_score, test_feature)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()
    
    ax.scatter(time_score, test_feature, s=point_size, alpha=alpha,
               color=point_color, edgecolors="none")
    
    # A straight line needs two distinct x values to be defined.
    if fit_line and number_samples >= 2 and np.ptp(time_score) > 0:
        slope, intercept = np.polyfit(time_score, test_feature, 1)
        x_line = np.array([time_score.min(), time_score.max()])
        ax.plot(x_line, slope * x_line + intercept, color=fit_color, linewidth=2)
    
    if annotate_stats:
        if np.isnan(spearman_corr):
            stats_text = f"n = {number_samples}\nnot testable"
        else:
            p_text = "p < 0.001" if spearman_p_value < 0.001 else f"p = {spearman_p_value:.3f}"
            stats_text = f"rho = {spearman_corr:.2f}\n{p_text}\nn = {number_samples}"
        
        ax.text(pos_text_x, pos_text_y, stats_text, transform=ax.transAxes,
                va="top", ha="left", fontsize=font_size,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="lightgrey"))
    
    ax.set_xlabel(xlabel, fontsize=font_size)
    ax.set_ylabel(ylabel if ylabel is not None else feature, fontsize=font_size)
    ax.set_title(title if title is not None else f"{feature} along the trajectory",
                 fontsize=font_size)
    ax.tick_params(labelsize=font_size)
    
    fig.tight_layout()
    return fig, ax


def plot_cell_types_distributions(
    proportions: pd.DataFrame = None,
    cell_types: list = None,
    plot_name: str = None,
    labels: str = "Predicted_Labels",
    file_path: str = None,
    figsize: tuple = (15, 7),
    label_order=None,
    label_colors=None,
    fontsize=24,
    rotation=45):
    """
    Boxplot of cell type proportions, split by group.
    
    Parameters:
    proportions (DataFrame): Cell type proportions per sample, plus the label column.
    cell_types (list): Columns of proportions to plot.
    plot_name (str): File name for the saved figure. Required. An extension is
                     added if none is given.
    labels (str): Column of `proportions` holding the group assignment.
    file_path (str or Path): Output directory. Default 'Results_PILOT/plots'.
    label_order (list or None): Groups to plot, in the order given. If None, all
                                groups are plotted in sorted order.
    label_colors (dict, list, str or None): Passed to seaborn as `palette`. A dict
                                keyed by group label is safest, since it stays
                                correct when `label_order` selects a subset.
    """
    if plot_name is None:
        raise ValueError("plot_name is required; there is no default file name.")
    
    # Checking if the labels and cell types are present.
    required = {labels, *cell_types}
    missing = required - set(proportions.columns)
    if missing:
        raise KeyError(f"proportions is missing the columns: {missing}")
    
    present = set(map(str, proportions[labels].dropna().unique()))
    
    if label_order is not None:
        unknown = {g for g in label_order if str(g) not in present}
        if unknown:
            raise ValueError(
                f"label_order has groups not in {labels}: {unknown}. "
                f"We have: {sorted(present)}"
            )
        
        label_order = [str(x) for x in label_order]
        if len(label_order) == 0:
            raise ValueError("label_order is empty; nothing would be plotted.")
        
    else:
        label_order = sorted(present)
    
    plot_path = Path("Results_PILOT/plots") if file_path is None else Path(file_path)
    plot_path.mkdir(parents=True, exist_ok=True)
    
    # Keeping the label column last so pd.melt uses it as the id variable.
    col_names = [c for c in cell_types if c != labels] + [labels]
    df_long = pd.melt(proportions.loc[:, col_names], labels)
    df_long[labels] = df_long[labels].astype(str)
    
    # Restricting to the requested groups. When label_order was not given this
    # keeps everything, since label_order is then all present groups.
    df_long = df_long[df_long[labels].isin(label_order)].copy()
    df_long[labels] = pd.Categorical(df_long[labels], categories=label_order, ordered=True)
    df_long = df_long.sort_values(labels)
    
    fig = plt.figure(figsize=figsize)
    ax = plt.gca()
    sns.boxplot(x="variable", hue=labels, y="value", data=df_long, palette=label_colors, ax=ax)
    ax.set_ylabel("Cell proportions", fontsize=fontsize)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelsize=fontsize)
    plt.setp(ax.get_xticklabels(), rotation=rotation, ha="right", rotation_mode="anchor")
    ax.legend(title=labels, fontsize=fontsize)
    fig.tight_layout()
    
    out_file = Path(plot_name)
    if not out_file.suffix:
        out_file = out_file.with_suffix(".png")
    
    fig.savefig(plot_path / out_file.name, facecolor="white")
    plt.show()
    plt.close(fig)


def gene_annotation_cell_type_subgroup(
    data: pd.DataFrame = None,
    cell_type: str = "Unknown",
    group: str = None,
    symbol: str = 'gene',
    sig_col: str = 'significant_gene',
    sources: list = None,
    num_gos: int = 10,
    figsize: tuple = (6, 4),
    font_size: int = 12,
    max_length: int = 50,
    path_to_results: str = 'Results_PILOT',
    my_pal: dict = None,
    save_plot: bool = True,
    organism: str = 'hsapiens',
    marker_size: int = 300
):
    """
    Performs GO enrichment analysis for a specific cell-type subgroup and 
    visualizes the top enriched terms.

    Parameters
    ----------
    data : pandas.DataFrame, optional
        A dataframe containing gene significance labels. If None, the function 
        attempts to load genes from a CSV file based on the `path_to_results`.
    cell_type : str, optional
        The specific cell type being analyzed. Used for path naming and 
        plot titles. Default is "Unknown".
    group : str, optional
        The specific group label (e.g., "Upregulated") to filter for in 
        `sig_col`.
    symbol : str, optional
        The column name in `data` containing gene symbols. Default is 'gene'.
    sig_col : str, optional
        The column name in `data` that defines the significance group. 
        Default is 'significant_gene'.
    sources : list, optional
        The g:Profiler data sources to query (e.g., ["GO:BP", "KEGG"]). 
        If None, queries all default sources.
    num_gos : int, optional
        The maximum number of top GO terms to display in the plot. Default is 10.
    figsize : tuple, optional
        The dimensions of the output plot. Default is (6, 4).
    font_size : int, optional
        The base font size for plot text elements. Default is 12.
    max_length : int, optional
        The character limit for wrapping GO term names on the y-axis. 
        Default is 50.
    path_to_results : str, optional
        The root directory for reading data and saving results. 
        Default is 'Results_PILOT'.
    my_pal : dict, optional
        A dictionary mapping group names to specific colors. 
        Default is None (defaults to 'tab:blue').
    save_plot : bool, optional
        Whether to save the plot as a PDF and the full results as a CSV. 
        Default is True.
    organism : str, optional
        The organism ID for g:Profiler (e.g., 'hsapiens', 'mmusculus'). 
        Default is 'hsapiens'.
    marker_size : int, optional
        The size of the dots in the scatter plot. Default is 300.

    Returns
    -------
    None or str
        Displays a plot and saves files if `save_plot` is True. 
        Returns an error string if no data is found or enrichment fails.
    """
    
    base_path = Path(path_to_results)
    go_path = base_path / 'Diff_Expressions_Results' / cell_type / 'GO_analysis'
    
    # Plotting - color palette 
    if my_pal and group in my_pal:
        plot_color = my_pal[group]
    else:
        plot_color = 'tab:blue'

    # Get query genes
    if data is not None:
        query_genes = data.loc[data[sig_col] == group, symbol].dropna().unique().tolist()
    else:
        # Load from CSV if no DataFrame provided
        file_path = base_path / 'Diff_Expressions_Results' / cell_type / f"Significant_genes_{cell_type}_{group}.csv"
        if not file_path.exists():
            return f"Error: No data provided and file not found at {file_path}"
        
        df_load = pd.read_csv(file_path, index_col=0)
        query_genes = df_load['0'].dropna().tolist() if '0' in df_load.columns else df_load.index.dropna().tolist()

    if not query_genes:
        return f"No significant genes found for group: {group}"

    # GProfiler Analysis
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(
        organism=organism,
        query=query_genes,
        no_evidences=False,
        sources=sources
    )

    if gprofiler_results is None or gprofiler_results.empty:
        return "Not enough GO information found for these genes."

    # Data preparation for plotting
    num_gos = min(num_gos, gprofiler_results.shape[0])
    selected_gps = gprofiler_results.head(num_gos).copy()
    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    # Wrap labels for better fit
    selected_gps['name'] = selected_gps['name'].apply(lambda x: "\n".join(tw.wrap(x, max_length)))

    # Plotting
    plt.figure(figsize=figsize, dpi=100)
    plt.style.use('default')
    sns.scatterplot(data=selected_gps, x="nlog10", y="name", s=marker_size, color=plot_color)

    plt.title(f'GO enrichment in {cell_type} ({group})\n(n={len(query_genes)} genes)', fontsize=font_size + 2)
    plt.xticks(size=font_size)
    plt.yticks(size=font_size)
    plt.ylabel("GO Terms", size=font_size)
    plt.xlabel("$-log_{10}$ (P-value)", size=font_size)
    
    # Saving
    if save_plot:
        go_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(go_path / f"{group}_GO_plot.pdf", bbox_inches='tight', facecolor='white')
        gprofiler_results.to_csv(go_path / f"{group}_full_GO_results.csv")
    
    plt.show()


def exploring_specific_genes(cluster_name='cell_type',
                             font_size=24,
                             gene_list=[],
                             fig_size=(64, 56),
                             p_value=0.01,
                             create_new_plot_folder=True,
                             fc_ther=0.5):
    """
    Explore specific genes within a cluster to analyze their patterns in comparison to other cell types.

    Parameters:
        cluster_name : str
            The name of the cluster you're interested in exploring.
        sort : list
            List of column names for sorting the results.
        number_genes : int
            Number of top genes to display for each pattern (linear, quadratic, etc.).
        cluster_names : list
            List of cluster names for exploration (if empty, all clusters will be considered).
        font_size : int
            Font size for text and labels in plots.
        gene_list : list
            List of specific genes to explore within the cluster.
        fig_size: tuple,optional
              Size of the plot.

        fc_ther: float, optional
            threshold for FC.

    Returns:
        Show the genes for the interested cell types
    """
    result_path=Path('Results_PILOT')
    file_name = "Whole_expressions.csv"
    cluster_names = [f.stem for f in (result_path / 'cells').iterdir() if f.is_file()]
    
    all_stats_extend = pd.read_csv(result_path / "gene_clusters_stats_extend.csv", sep = ",")
    
    with open(result_path / 'genes_dictionary.pkl', 'rb') as handle:
        gene_dict = pickle.load(handle)
    
    pline = np.linspace(1, 20, 20)
    filtered_all_stats_extend=all_stats_extend[all_stats_extend['gene'].isin(gene_list)]
    filtered_all_stats_extend=filtered_all_stats_extend[filtered_all_stats_extend['cluster'].isin([cluster_name])]
    
 
    plot_stats_by_pattern(cluster_names, filtered_all_stats_extend, gene_dict, pline, result_path, file_name,font_size=font_size,p_value=p_value,create_new_plot_folder=create_new_plot_folder,fc_ther=fc_ther)
    
    # Load the PNG image file
    if create_new_plot_folder:
        image_path = result_path+'/plot_genes_for_'+str(cluster_name)+'/'+str(cluster_name) + ".png"  # Replace with the actual path to your PNG image
        image = imread(image_path)
    else:        
        image_path =result_path+'plots_gene_cluster_differentiation/'+cluster_name+'.png'  # Replace with the actual path to your PNG image
        image = imread(image_path)
    
    # Set the size of the figure
    fig, ax = plt.subplots(figsize=fig_size)  # Adjust the width and height as needed

    # Display the PNG image
    ax.imshow(image)
    ax.axis('off')  # Turn off axis labels and ticks
    plt.show()
    
  
def go_enrichment(df,
                  num_gos=20,
                  source=None,
                  cell_type='cell_type',
                  fontsize=32,
                  s=200, 
                  figsize = (15,12),
                  color = 'tab:blue',
                  dpi=100,
                  bbox_inches = 'tight', 
                  facecolor='white', 
                  transparent=False,
                  organism='hsapiens'):
    
    """
    Perform Gene Ontology (GO) enrichment analysis and create a scatterplot of enriched terms.

    Parameters:
    -----------
    df : DataFrame
        Input DataFrame containing gene information.
    num_gos : int, optional
        Number of top enriched GO terms to visualize. Default is 20.
    source : str, optional
        Specify the source of GO terms. The default is None.
    cell_type : str, optional
        Name of the cell type for which the GO enrichment is performed. Default is 'cell_type'.
    fontsize : int, optional
        Font size for the plot title and labels. Default is 32.
    s : int, optional
        Marker size for the scatterplot. Default is 200.
    figsize : tuple, optional
        Figure size (width, height) in inches. Default is (15, 12).
    color : str, optional
        Color of the scatterplot markers. Default is 'tab:blue'.
    dpi : int, optional
        Dots per inch for the saved plot image. Default is 100.
    bbox_inches : str, optional
        Bounding box options for saving the plot. Default is 'tight'.
    facecolor : str, optional
        Background color of the saved plot image. Default is 'white'.
    transparent : bool, optional
        Whether the saved plot image has a transparent background. Default is False.
    organism : str, optional
        Organism for GO enrichment analysis. Default is 'hsapiens'.

    Returns:
    --------
    None
        Saves the scatterplot of enriched GO terms as a PDF file.
    """

    result_path=Path('Results_PILOT/')
    df_sorted =df
    gp = GProfiler(return_dataframe=True)
    gprofiler_results = gp.profile(organism = organism,
                query = list(df_sorted['gene'].values))
    
    go_path = result_path / 'GO' / cell_type
    go_path.mkdir(parents=True, exist_ok=True)
    gprofiler_results.to_csv(result_path+'GO/'+cell_type+'/'+cell_type+"_all_gprofiler_results.csv")
    if(gprofiler_results.shape[0] < num_gos):
        num_gos = gprofiler_results.shape[0]

    #selected_gps = gprofiler_results.loc[0:num_gos,['name', 'p_value']]
    if source:   
        gprofiler_results = gprofiler_results[gprofiler_results['source']==source]
    selected_gps = gprofiler_results.head(num_gos)[['name', 'p_value']]

    selected_gps['nlog10'] = -np.log10(selected_gps['p_value'].values)

    figsize = figsize

    plt.figure(figsize = figsize, dpi = dpi)
    plt.style.use('default')
    sns.scatterplot(data = selected_gps, x="nlog10", y="name", s = s, color =color)

    plt.title('GO enrichment in ' + cell_type, fontsize = fontsize)

    plt.xticks(size = fontsize)
    plt.yticks(size = fontsize)

    plt.ylabel("GO Terms", size = fontsize)
    plt.xlabel("-$log_{10}$ (P-value)", size = fontsize)
    #plt.savefig(path+'GO/'+cell_type+".pdf", bbox_inches = 'tight', facecolor='white', transparent=False)
    save_path = Path(result_path / "GO" / cell_type)
    save_path.mkdir(parents=True, exist_ok=True)
    plt.savefig(save_path / f"{cell_type}.pdf", bbox_inches = bbox_inches, facecolor=facecolor, transparent=transparent)
    

def plt_gene_cluster_differentiation(cellnames=['healthy_CM','Myofib'],font_size=22,p_value=0.01,fc_ther=0.5):
    """
    Generate and save plots showcasing gene expression patterns for selected cell clusters.

    Parameters:
    -----------
    cellnames : list, optional
        List of cluster names for which gene expression patterns will be visualized. Default is ['healthy_CM', 'Myofib'].
    font_size : int, optional
        Font size for labels and titles in the plots. Default is 16.
    p_value : float, optional
        P-value threshold for significance. Default is 0.05.
    fc_ther: float, optional
            threshold for FC.
    Returns:
    --------
    None
        Generates and saves scatterplots of gene expression patterns for selected clusters.
    """
    
    file_name = "Whole_expressions.csv"
    result_path=Path('Results_PILOT')
    all_stats_extend=pd.read_csv(result_path / 'gene_clusters_stats_extend.csv')
    with open(result_path / 'genes_dictionary.pkl', 'rb') as handle:
        gene_dict = pickle.load(handle)
    
    pline = np.linspace(1, 20, 20)
    plot_stats_by_pattern(cluster_names=cellnames, all_stats_extend=all_stats_extend, gene_dict=gene_dict, 
                          pline=pline, path_to_results=result_path, file_name=file_name,font_size=font_size,
                          p_value=p_value,fc_ther=fc_ther)


def qq_plot_gene(target, data, sorted_best, gene_name):
    """
    Generate a QQ plot for a specific gene's performance.

    Parameters
    ----------
        target : pd.DataFrame
            The target DataFrame containing gene labels.
        data : pd.DataFrame
            The data DataFrame containing gene data.
        sorted_best : dict
            A dictionary containing sorted gene results.
        gene_name : str
            The name of the gene to create the QQ plot for.

    Returns:
        None: Displays the QQ plot for the specified gene's performance.
    """

    x = np.array(list(data['label']))

    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    # pattern = sorted_best[gene_name][2]
    data = pd.DataFrame({'x': x, 'y':best_tf})

    X = generate_feature(best_func, x)
    new_X = np.append(np.ones((len(X),1)), X, axis=1)
    # print(X.shape)
    # print(best_results['params'])
    predictions = np.matmul(new_X, best_results['params'])
    
    plt.figure(figsize = (8,6))
    stats.probplot( (best_tf - predictions), dist="norm", plot=pylab)
    pylab.show()


def plot_best_matches_cell_types(target, 
                                 data,
                                 df,
                                 sorted_best, 
                                 scale_name, 
                                 min_target=0, 
                                 max_target=35,
                                 num=11,
                                 width=25,
                                 height=25,
                                 xlim=4,
                                 point_size=100,
                                 color_back=None,
                                 fontsize=28,
                                 alpha=1,
                                 cmap='viridis'):
    """
    Plot the best-fitted models for cell types.

    Parameters:
        target : pd.DataFrame
            The target data for gene activity.
        data : pd.DataFrame
            The data containing cell type labels.
        df : pd.DataFrame
            The data frame containing sample information.
        sorted_best : dict
            A dictionary containing the best-fitted model results for each factor, sorted by R-squared or modified R-squared.
        scale_name : str
            The name of the scale.
        min_target : int, optional
            The minimum target value.
        max_target : int, optional 
            The maximum target value.
        num : int, optional 
            Number of models to plot.
        width : int, optional
            Width of the figure.
        height : int, optional
            Height of the figure.
        xlim : int, optional
            X-axis limit.
        point_size : int, optional
            Size of data points.
        color_back : str, optional
            Background color of the plot.
        fontsize : int, optional
            Font size for the titles and labels.
        alpha : float, optional
            Alpha value for data points.
        cmap : str, optional
            Colormap for data points.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    
    sorted_best = dict(sorted(sorted_best.items(), key=lambda x: x[1][-1]))

    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    start=df[df['Time_score']==min_x]['sampleID'].unique()
    end=df[df['Time_score']==max_x]['sampleID'].unique()
    xlabel='From  '+ start+'  to  '+end
    
    plt.figure(figsize=((width, height)))
    plt.subplots_adjust(wspace = 0.5, hspace = 1 )
   # plt.suptitle(xlabel, fontsize=20, y=0.95)
    plt.tight_layout()
    
    j = 1
    for i in range(num):
        
      
        
        
        best_tf_name = list(sorted_best)[i]
        best_tf = target.loc[[best_tf_name]].values
        best_results = list(sorted_best.items())[i][1][1]
        best_func = list(sorted_best.items())[i][1][0]
        p_vals_best_model = sorted_best[best_tf_name][1]['pvalues']
        polyline = np.linspace(min_target,max_target)
        polyline = generate_feature(best_func, polyline)
        
        ax = plt.subplot(math.ceil(num/4), 4, j)
        plt.xticks(np.arange(min_x,max_x,xlim))
        if color_back!=None:
            ax.set_facecolor(color_back)
        ax.scatter(x, best_tf, c =best_tf ,cmap=cmap,alpha=alpha,s=point_size)
        new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
        curve = np.matmul(new_polyline, best_results['params'])
        ax.plot(np.linspace(min_target,max_target), curve)
        
        formatted_number= "{:.3e}".format(sorted_best[best_tf_name][-1])
        
        line1=str(best_tf_name)+"\n" +"P-val: "+ "\n"
        bold_word=formatted_number                     

        if float(sorted_best[best_tf_name][-1]) <= 0.05:
            title =line1+"{}".format(bold_word)
            font_props = FontProperties(weight='bold')
            ax.set_title(title, fontproperties=font_props,fontsize=fontsize)
        else:
            combined_line = "{} {}".format(line1,bold_word)
            ax.set_title(combined_line,fontsize=fontsize)
            
            
      #  ax.set_title(str(best_tf_name)+ "\n" +"P-val: "+ "\n" +str("{:.3e}".format(sorted_best[best_tf_name][-1])),fontsize=fontsize)
            
       # best_tf_name=str(best_tf_name)+"\n Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj']))
       # ax.set_title(best_tf_name,fontsize=fontsize)
       
        ax.axis(xmin = min_target,xmax = max_target)
    
        if((j+4) % 4 == 1):
            #ax.set(ylabel = scale_name)
            ax.set_ylabel(scale_name, fontsize=fontsize)
            
       # if show_R:
        #    if x_lab==None:
         #       x_lab=np.mean(x)
          #  if y_lab==None:
           #     y_lab=np.mean(best_tf)
            #ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
             #        (x_lab,y_lab),
              #       color=color_label,
               #      size=fontsize)
     
    
        j += 1


def plot_best_matches(target, data,df, sorted_best, scale_name, plot_color='tab:orange',num=16,width=25,height=25,x_lim=4,fontsize=24,alpha=0.5,cmap='viridis',color_back=None):
    """
    Plot the best-fitted models for different patterns.

    Parameters:
        target : pd.DataFrame
            The target data for gene activity.
        data : pd.DataFrame
            The data containing cell type labels.
        df : pd.DataFrame
            The data frame containing sample information.
        sorted_best : dict
            A dictionary containing the best-fitted model results for each factor, sorted by R-squared or modified R-squared.
        scale_name : str
            The name of the scale.
        plot_color : str, optional
            Color of the plots. Default is 'tab:orange'.
        num : int, optional
            Number of models to plot. Default is 16.
        width : int, optional
            Width of the figure. Default is 25.
        height : int, optional
            Height of the figure. Default is 25.
        x_lim : int, optional
            X-axis limit. Drfault is 4.
        fontsize : int, optional
            Font size for titles and labels. Default is 24.
        alpha : float, optional
            Alpha value for data points. Default is 0.5.
        cmap : str, optional
            Colormap for data points. Default is 'viridis'.
        color_back : str, optional
            Background color of the plot. Default is None.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    
    plt.rcParams.update({'font.size': fontsize})

    
    all_motifs = np.array(list(sorted_best.keys()))
    all_motifs = all_motifs.reshape(all_motifs.shape[0],1)
    all_patterns = np.array([i[2] for i in list(sorted_best.values())])
    all_patterns = all_patterns.reshape(all_patterns.shape[0],1)
    mat = np.append(all_motifs, all_patterns, axis = 1)
    patterns=np.unique(mat[:,1])

    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    start=df[df['Time_score']==min_x]['sampleID'].unique()
    end=df[df['Time_score']==max_x]['sampleID'].unique()

    plt_count=0
    plt.figure(figsize=((width, height)))
    xlabel='From  '+ start+'  to  '+end
   # plt.suptitle(xlabel, fontsize=20, y=0.05)
    #plt.tight_layout()
    plt.subplots_adjust(


                    wspace=0.5,
                    hspace=1)

    genes_expressions=patterns
    pltt=0
    for pattern in genes_expressions:
        counter=0
        flag=False
        j=0

        if (pltt+4)%4==3:
            pltt=pltt+1
        elif (pltt+4)%4==2:
            pltt=pltt+2
        elif (pltt+4)%4==1:
            pltt=pltt+3

        for key,values in sorted_best.items():

            if list(sorted_best.values())[counter][2]==pattern:
                pltt+=1
                ax = plt.subplot(len(patterns), 4, pltt)
                if color_back!=None:
                    ax.set_facecolor(color_back)
                plt.xticks(np.arange(min_x,max_x,x_lim))
                #flag=True
                best_tf_name = list(sorted_best)[counter]
                best_tf = target.loc[[best_tf_name]].values
                best_results = list(sorted_best.items())[counter][1][1]
                best_func = list(sorted_best.items())[counter][1][0]

                data = pd.DataFrame({'x': x.tolist(), 'y':best_tf.T.tolist()})

                pline = np.linspace(min_x, max_x)
                polyline = generate_feature(best_func, pline)
                new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
                curve = np.matmul(new_polyline, best_results['params'])


               
                    
                ax.scatter(x, best_tf, c = best_tf,cmap=cmap, alpha = alpha)
                 
                
                    
                formatted_number= "{:.3e}".format(sorted_best[best_tf_name][-1])

                ax.plot(pline, curve, c = plot_color)
                
                line1=str(best_tf_name)+"\n" +"P-val: "+ "\n"
                bold_word=formatted_number 
                
          
                
                
                
                
                if float(sorted_best[best_tf_name][-1]) <= 0.05:
                    
                    title =line1+"{}".format(bold_word)
                    font_props = FontProperties(weight='bold')
                    ax.set_title(title, fontproperties=font_props,fontsize=fontsize)
                    
                   
                else:
                    combined_line = "{} {}".format(line1,bold_word)
                    ax.set_title(combined_line,fontsize=fontsize)
                    
                
                
               # best_tf_name=best_tf_name+"\n Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj']))
               # ax.set_title(best_tf_name,fontsize=fontsize)

               

                ax.axis(xmin = min(x)-0.01,xmax = max(x)+0.01)

                if((j+4) % 4 == 0):
                    ax.set_ylabel(pattern, fontsize=fontsize)
                    #ax[plt_count,j].set(ylabel = pattern,fontsize=8)
                    #ax[plt_count,j].set_xlabel(str(xlabel),fontsize=12)
            #    if show_R:
             #       if x_lab==None:
              #          x_lab=np.mean(x)
               #     if y_lab==None:
                #        y_lab=np.mean(best_tf)
                       
                 #   ax.annotate("Modified adj R2 = " + str("{:.2f}".format(best_results['mod_rsquared_adj'])),
                              
                  #           (x_lab,y_lab),
                   #          color=color_label,
                    #         size=fontsize)

                j += 1
                if j%4==0:
                    break




            counter=counter+1  

 
def plot_two_genes(adata, sorted_best_WT, sorted_best_KO, gene_name, scale_name, plot_color1 = 'tab:blue', plot_color2 = 'tab:red'):
    
    """
    Plot the gene activity for two different conditions (e.g., WT and KO) along with their best-fitted models.

    Parameters:
        adata : AnnData
            An AnnData object containing the data.
        sorted_best_WT : dict
            A dictionary containing the best-fitted model results for the wild-type (WT) condition.
        sorted_best_KO : dict
            A dictionary containing the best-fitted model results for the knock-out (KO) condition.
        gene_name : str
            The name of the gene for which the activity is plotted.
        scale_name : str
            The name of the scale.
        plot_color1 : str, optional
            Color for the WT condition plot. Default is 'tab:blue'.
        plot_color2 : str, optional
            Color for the KO condition plot. Default is 'tab:red'.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    WT_ATAC_data = pd.DataFrame()
    WT_ATAC_target = adata[adata.obs.group == 'WT'].to_df().transpose()
    WT_ATAC_data['label'] = list(adata[adata.obs.group == 'WT'].obs['label'].values)
    
    WT_x = WT_ATAC_data['label']
    WT_tf = list(WT_ATAC_target.loc[[gene_name]].values.ravel())
    WT_results = sorted_best_WT[gene_name][1]
    WT_func = sorted_best_WT[gene_name][0]
    # WT_pattern = sorted_best_WT[gene_name][2]
    
    KO_ATAC_data = pd.DataFrame()
    KO_ATAC_target = adata[adata.obs.group == 'KO'].to_df().transpose()
    KO_ATAC_data['label'] = list(adata[adata.obs.group == 'KO'].obs['label'].values)
    
    KO_x = KO_ATAC_data['label']
    KO_tf = list(KO_ATAC_target.loc[[gene_name]].values.ravel())
    KO_results = sorted_best_KO[gene_name][1]
    KO_func = sorted_best_KO[gene_name][0]
    # KO_pattern = sorted_best_KO[gene_name][2]
    
    pline = np.linspace(min(KO_x), max(KO_x))
    
    polyline = generate_feature(WT_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    WT_curve = np.matmul(new_polyline, WT_results['params'])
    
    polyline = generate_feature(KO_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    KO_curve = np.matmul(new_polyline, KO_results['params'])
    
    ax = plt.subplot(1, 1, 1)
    ax.scatter(WT_x, WT_tf, c = WT_tf, alpha = 0.5, cmap='Blues', norm=colors.CenteredNorm(-1),)
    ax.scatter(KO_x, KO_tf, c = KO_tf, alpha = 0.5, cmap='Reds', norm=colors.CenteredNorm(-1),)
    
    
    
    ax.axis(xmin=min(WT_x),xmax=max(WT_x))

    ax.plot(np.linspace(min(WT_x),max(WT_x)), WT_curve, c = plot_color1)
    ax.plot(np.linspace(min(WT_x),max(WT_x)), KO_curve, c = plot_color2)
    ax.set_title(re.sub('.*?:', '', gene_name))
    ax.set_ylabel(scale_name)
    

def plot_one_gene(target, data, sorted_best, gene_name, scale_name, plot_color):
    
    """
    Plot the gene activity and its best-fitted model for a single gene.

    Parameters:
        target : pd.DataFrame
            A Pandas DataFrame containing the target gene expression data.
        data : pd.DataFrame
            A Pandas DataFrame containing the data.
        sorted_best : dict
            A dictionary containing the best-fitted model results for multiple genes.
        gene_name : str
            The name of the gene for which the activity is plotted.
        scale_name : str
            The name of the scale.
        plot_color : str
            Color for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    # data = pd.DataFrame({'x': x, 'y':best_tf})


    
    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])

    ax = plt.subplot(1, 1, 1)
    ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='plasma')
    ax.axis(xmin=min(x),xmax=max(x))

    ax.plot(np.linspace(min(x),max(x)), curve, c = plot_color)
    ax.set_title(re.sub('.*?:', '', gene_name))
    ax.set_ylabel(scale_name)
    
    
    ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                  (0.2,np.max(best_tf)-0.1),
                  color='black',
                  size=14)
    ax.set_xlabel(pattern)


def plot_gene(target, data, sorted_best, gene_name, scale_name, plot_color):
    """
    Plot gene expression data along with the best-fitted curve and statistical information.

    Parameters:
        target : pd.DataFrame
            A Pandas DataFrame containing gene expression data.
        data : pd.DataFrame
            A Pandas DataFrame containing additional data such as labels.
        sorted_best : dict
            A dictionary containing the best-fitted models for different genes.
        gene_name : str
            The name of the gene to plot.
        scale_name : str
            The name of the scale for the y-axis.
        plot_color : str
            The color to use for plotting.

    Returns:
        None

    Note:
        This function plots gene expression data for a specific gene along with the best-fitted curve,
        statistical information, and other visual elements using the ggplot library.
        It includes the gene name, adjusted R-squared value, p-value, and the fitted function type in the plot.
        The appearance of the plot may vary based on the fitted function type (linear, quadratic, or linear-quadratic).
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())

    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    data = pd.DataFrame({'x': x, 'y':best_tf})

    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])
    
    p = ggplot(aes(x='x', y='y'), data,)
    p = p + geom_pointdensity(size = 3, show_legend=False)
    
    dt = pd.DataFrame({'x': pline, 'y': curve})
    p = p + geom_line(aes(x='x', y='y'), dt, colour = plot_color, size = 1)
    
    p = p + theme(panel_background=element_rect(fill='white'),              
          axis_line_x=element_line(color='black'),
          axis_line_y=element_line(color='black'),
          panel_grid=element_blank(),
          panel_border=element_blank())
    p = p + labs(title = gene_name + "\n\n" + \
                 "adj $R^{2}$ = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
                            + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4]))
                 , y = scale_name, x = pattern)
    # p = p + annotate(geom = "text", 
    #                   label = "adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
    #                         + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4])),
    #                   x = 0, 
    #                   y = np.max(best_tf))
    if(best_func == 'linear'):
        p = p + annotate(geom = "text", 
                          label = "ge = " +
                              str(round(best_results['params'][0],2)) \
                              + ["", "+"][best_results['params'][1] > 0]+ \
                                  str(round(best_results['params'][1],2)) + "t" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " +
                              str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'linear_quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + 
                              str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "t" \
                                   + ["", "+"][best_results['params'][2] > 0]+\
                                       str(round(best_results['params'][2],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    p.draw()
    

def plot_gene_specific(target, data, sorted_best, gene_name, scale_name, plot_color):
    
    """
    Plot gene expression and its best-fitted model for a single gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing the data.
        sorted_best (dict): A dictionary containing the best-fitted model results for multiple genes.
        gene_name (str): The name of the gene for which the expression is plotted.
        scale_name (str): The name of the scale.
        plot_color (str): Color for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = data['label']
    min_x = min(x)
    max_x = max(x)
    
    best_tf = list(target.loc[[gene_name]].values.ravel())
    best_results = sorted_best[gene_name][1]
    best_func = sorted_best[gene_name][0]
    pattern = sorted_best[gene_name][2]
    
    new_x = generate_feature(best_func, np.array(list(data['label'])))
    all_main_polyline = np.append(np.ones((len(new_x),1)), new_x, axis=1)
    new_curve = np.matmul(all_main_polyline, best_results['params'])
    
    my_data = pd.DataFrame({'x': x, 'y':best_tf})
    pline = np.linspace(min_x, max_x)
    polyline = generate_feature(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    curve = np.matmul(new_polyline, best_results['params'])
    
    
    # stdev = np.sqrt(sum((best_tf - target.loc[gene_name])**2) / (len(target.loc[gene_name]) - 2))
    
    
    print("mean value: " + str(np.mean(target.loc[gene_name])))
    ss_res = np.sum( (target.loc[gene_name] - new_curve )**2)
    print("SSres: " + str(ss_res) )
    
    delata = 1.35
    modified_ss_res = 0
    res_e = target.loc[gene_name] - new_curve
    for e in res_e:
        if abs(e) < delata:
            modified_ss_res += 1/2*np.sum(e**2)
        else:
            modified_ss_res += delata *(abs(e) - (1/2)*delata)
    print("Modified SSres: " + str(modified_ss_res))
    
    ss_tot = np.sum( (target.loc[gene_name] - np.mean(target.loc[gene_name]))**2)
    print("SStot: " + str(ss_tot) )
    r2 = 1 - ( ss_res / ss_tot )
    print("R2: " + str(r2) )
    
    r2 = 1 - ( modified_ss_res / ss_tot )
    r2 = 1 - (1-r2)*(len(best_tf)-1)/(len(best_tf)-new_x.shape[1]-1)
    print("Modified R2: " + str(r2) )
    
    mse = mean_absolute_error(new_curve, target.loc[gene_name])
    print("MAE: " + str(mse))
    
    p = ggplot(aes(x='x', y='y'), my_data,)
    p = p + geom_pointdensity(size = 3, show_legend=False)
    
    dt = pd.DataFrame({'x': pline, 'y': curve})
    p = p + geom_line(aes(x='x', y='y'), dt, colour=plot_color, size = 1)
    
    mean_dt = pd.DataFrame({'x': pline, 'y': np.mean(target.loc[gene_name])*np.ones(len(pline))})
    p = p + geom_line(aes(x='x', y='y'), mean_dt, colour='tab:red', size = 1)
    
    p = p + theme(panel_background=element_rect(fill='white'),              
          axis_line_x=element_line(color='black'),
          axis_line_y=element_line(color='black'),
          panel_grid=element_blank(),
          panel_border=element_blank())
    p = p + labs(title = gene_name, y = scale_name, x = pattern)
    p = p + annotate(geom = "text", 
                      label = "\n\n adj $R^2$ = " + str("{:.2f}".format(best_results['rsquared_adj'])) \
                            + "\n Modified adj $R^2$ = " + str("{:.2f}".format(r2))\
                            + "\n p-value = " + str("{:.2f}".format(sorted_best[gene_name][4])),
                      x = np.mean(pline), 
                      y = np.max(best_tf)+0.15)
    if(best_func == 'linear'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                              + ["", "+"][best_results['params'][1] > 0]+ \
                                  str(round(best_results['params'][1],2)) + "t" ,
                          color = "tab:orange",
                          size = 16,
                          x = np.max(pline)-15, 
                          y = np.max(curve)+0.01)
    elif(best_func == 'quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    elif(best_func == 'linear_quadratic'):
        p = p + annotate(geom = "text", 
                          label = "ge = " + str(round(best_results['params'][0],2)) \
                               + ["", "+"][best_results['params'][1] > 0]+ \
                                   str(round(best_results['params'][1],2)) + "t" \
                                   + ["", "+"][best_results['params'][2] > 0]+\
                                       str(round(best_results['params'][2],2)) + "$t^2$" ,
                          x = pline[10], 
                          y = curve[10])
    p = p + annotate(geom = "text", 
                    label = "Mean Model" ,
                        x = np.max(mean_dt['x'])-8, 
                        y = np.max(mean_dt['y'])-0.03,
                    color = 'tab:red',
                     size = 16
                        )
    p.draw()
   

def plot_gene_distribtion(target, gene_name):
    """
    Plot the distribution of gene expression for a specific gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        gene_name (str): The name of the gene for which the expression distribution is plotted.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    data = pd.DataFrame({'gene': target.loc[gene_name]})
    plt.figure(figsize = (8,6))
    sns.histplot(data=data, x='gene', edgecolor='k')
    plt.title(gene_name + " expression distribtion")
    plt.xlabel('gene expression', size = 16)
    plt.ylabel('Frequency', size = 16)
    plt.yticks(size = 12, weight = 'bold')
    plt.xticks(size = 12, weight = 'bold')
    plt.show()


def plot_gene_density(target, data, sorted_best, gene_name, scale_name, plot_color):
    """
    Plot the density distribution of gene expression for a specific gene.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expression data.
        data (pd.DataFrame): A Pandas DataFrame containing the data with 'label' and 'y' columns.
        sorted_best (dict): A dictionary containing sorted best-fitted models for genes.
        gene_name (str): The name of the gene for which the density distribution is plotted.
        scale_name (str): The name of the scale or measurement associated with the gene expression.
        plot_color (str): The color for the density plot.

    Returns:
        None: This function generates the density plot but does not return any value.
    """

    x = data['label']

    best_tf = list(target.loc[[gene_name]].values.ravel())
    data = pd.DataFrame({'x': x, 'y':best_tf})
    
    labels=np.sort([y for y in list(data.x.unique())])
    fig, axes = joypy.joyplot(data, by="x", column="y", labels=labels, 
                          grid="y", linewidth=1, legend=False, overlap=0.5, figsize=(6,5),kind="kde", bins=80,
                          title=gene_name, ylabels=False,
                          colormap=cm.autumn_r)


def plot_pval_rsq_correlation(table, feature1, feature2, show_fit = True, log_transform = False):
    
    """
    Plot the correlation between two features, optionally showing a linear fit.

    Parameters:
        table (pd.DataFrame): A Pandas DataFrame containing the data.
        feature1 (str): The name of the first feature to be plotted on the x-axis.
        feature2 (str): The name of the second feature to be plotted on the y-axis.
        show_fit (bool, optional): If True, a linear fit line will be plotted. Default is True.
        log_transform (bool, optional): If True, the second feature will be log-transformed (-log10) before plotting. Default is False.

    Returns:
        None: This function generates the scatter plot but does not return any value.
    """

    
    x = [float(x) for x in table[feature1].values]
    
    if(log_transform):
        y = [float(x) for x in table[feature2].values]
        y = -np.log10(y)
        data = {"x": x, "y": y}
    
    plt.figure(figsize=(6,6))
    sns.scatterplot(x, y, color='tab:blue')
    
    if(show_fit):
        results = ols("y ~ x", data=data).fit()
        print('ratio: ' + str(results.params[1]))
        pline = np.linspace(min(x),max(x))
        curve = np.matmul(np.array(results.params), [np.ones(len(pline)),pline])
        # y_pred = np.matmul(np.array(results.params), [np.ones(len(x)),x])
    
        sns.lineplot(pline, curve, color = 'tab:red')
    
        # textstr = '\n'.join((
        #         r'RMSE=%.2f' % (compute_metrics(y, y_pred)[1], ),
        #         r'MAE=%.2f' % (compute_metrics(y, y_pred)[2], ),
        #         r'r2=%.2f' % (compute_metrics(y, y_pred)[0], ),
        #         r'Cor=%.2f' % (compute_metrics(y, y_pred)[3], )))
        # plt.text(min(pline), max(y), textstr, fontsize=16,
        #         verticalalignment='top')
        
    plt.xlabel(feature1, fontsize=16)
    if(log_transform):
        plt.ylabel("$-log_{10}($" + str(feature2) + ")", fontsize=16)
    else:
        plt.ylabel(feature2, fontsize=16)
    plt.xticks(size = 14)
    plt.yticks(size = 14)
    plt.show()


def plot_condition(target, data, sorted_best, condition_type, scale_name, plot_color):
    """
    Plot the best-fitted models based on up or downregulation patterns in gene expression.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expressions.
        data (pd.DataFrame): A Pandas DataFrame containing the feature data.
        sorted_best (dict): A dictionary containing the best-fitted models for each gene.
        condition_type (str): The type of condition to plot, either 'increasing' or 'decreasing'.
        scale_name (str): The name of the scale for the y-axis.
        plot_color (str): The color to use for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    assert condition_type in ['increasing', 'decreasing'], 'condition type parameter must be increasing, decreasing'
    
    x = np.array(list(data['label']))
    j = 1
    i = 0
    
    if(condition_type == 'increasing'):
        plt.figure(figsize=(15, 12))
        plt.subplots_adjust(hspace=0.5, wspace=0.3)
        plt.suptitle("Increasing Expression", fontsize=18, y=0.95)
        plt.tight_layout()
   
        while(j<=9):
            
            best_tf_name = list(sorted_best)[i]
            best_tf = target.loc[[best_tf_name]].values
            best_results = list(sorted_best.items())[i][1][1]
            best_func = list(sorted_best.items())[i][1][0]
            
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(best_func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, best_results['params'])
            mit = list(sorted_best.items())[i][1][3]
            if(mit > 0):
                
                
                ax = plt.subplot(3, 3, j)

                ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
                ax.axis(xmin = min(x),xmax = max(x))
                ax.set(ylabel = scale_name)
                ax.plot(pline, curve, c = plot_color)
                ax.set_title(best_tf_name)
                
                ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                      (0,np.max(best_tf)-0.2),
                      color='black',
                      size=14)
        
                j += 1
            i += 1
    else:
        plt.figure(figsize=(15, 12))
        plt.subplots_adjust(hspace=0.5, wspace=0.3)
        plt.suptitle("Decreasing Expression", fontsize=18, y=0.95)
        plt.tight_layout()
        
        while(j<=9):
            
            best_tf_name = list(sorted_best)[i]
            best_tf = target.loc[[best_tf_name]].values
            best_results = list(sorted_best.items())[i][1][1]
            best_func = list(sorted_best.items())[i][1][0]
            
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(best_func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, best_results['params'])
            mit = list(sorted_best.items())[i][1][3]
            if(mit < 0):
                
                
                ax = plt.subplot(3, 3, j)
                ax.scatter(x, best_tf, c = best_tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
                ax.axis(xmin = min(x),xmax = max(x))
                ax.set(ylabel = scale_name)
                ax.plot(pline, curve, c = plot_color)
                ax.set_title(best_tf_name)
                
                ax.annotate("adj R2 = " + str("{:.2f}".format(best_results['rsquared_adj'])), 
                      (1.8,np.max(best_tf)-0.1),
                      color='black',
                      size=14)
        
                j += 1
            i += 1


def plot_6_best(target, data, sorted_best, scale_name, plot_color):
    """
    Plot the 6 best-fitted models for each pattern type (linear, quadratic, linear_quadratic) and regulation direction (up, down) in gene expression.

    Parameters:
        target (pd.DataFrame): A Pandas DataFrame containing the target gene expressions.
        data (pd.DataFrame): A Pandas DataFrame containing the feature data.
        sorted_best (dict): A dictionary containing the best-fitted models for each gene.
        scale_name (str): The name of the scale for the y-axis.
        plot_color (str): The color to use for the plot.

    Returns:
        None: This function generates the plot but does not return any value.
    """

    x = np.array(list(data['label']))
    
    plt.figure(figsize=(15, 12))
    plt.subplots_adjust(hspace=0.5)
    
    pattern_types = ['linear', 'quadratic', 'linear_quadratic']
    mit_types = ['up', 'down']
    
    k = 1
    for mit in mit_types:
        for pattern in pattern_types:
            motifs = data_interst(sorted_best, pattern, mit)
            tf_name = list(motifs.keys())[0]
            tf = target.loc[[tf_name]].values
            results = list(motifs.items())[0][1][1]
            func = list(motifs.items())[0][1][0]
            patt = list(motifs.items())[0][1][2]
            pline = np.linspace(min(x), max(x))
            polyline = generate_feature(func, pline)
            new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
            curve = np.matmul(new_polyline, results['params'])
            
            ax = plt.subplot(3, 3, k)
            ax.scatter(x, tf, c = tf, alpha = 0.5, cmap='RdBu', norm=colors.CenteredNorm(), )
            ax.axis(xmin=min(x),xmax=max(x))
            
            ax.plot(np.linspace(min(x),max(x)), curve, c = plot_color)
            ax.set_title(re.sub('.*?:', '', tf_name))
            if(k % 3 == 1):
                ax.set(ylabel = scale_name)
            
            ax.annotate("adj R2 = " + str("{:.2f}".format(results['rsquared_adj'])), 
                          (0.2,np.max(tf)-0.1),
                          color='black',
                          size=14)
            ax.set_xlabel(patt)
            k += 1
            

def plot_hor_vs_vert(data, subplot, x, y, c, xlabel, ylabel, rotation,
                     tick_bottom, tick_left, title,fontsize=24):
    
    '''
    Plot horizontal and vertical bar charts using Seaborn.
    
    Parameters:
        data : pandas.DataFrame
            The input DataFrame containing the data to be plotted.
        subplot : int
            The subplot number (1 or 2) for the placement of the chart.
        x : str
            The name of the column to be plotted on the x-axis.
        y : str
            The name of the column to be plotted on the y-axis.
        c : str
            The color column used to determine the color of bars.
        xlabel : str
            The label for the x-axis.
        ylabel : str
            The label for the y-axis.
        rotation : float
            The rotation angle for x-axis labels in degrees.
        tick_bottom : bool
            If True, display ticks on the bottom of the plot.
        tick_left : bool
            If True, display ticks on the left side of the plot.
        title : str
            The title of the plot.
        fontsize : int, optional
            The font size for title, labels, and ticks (default is 24).
    
    Returns:
        None
    '''
    ax=plt.subplot(1,2,subplot)
    cols = ['tab:red' if x >= 0 else 'tab:blue' for x in data[x]]
    
    sns.barplot(x=x, y=y, data=data, ci=None, palette=cols)
    plt.title(title, fontsize=fontsize, fontweight = 'bold')
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.xticks(fontsize=fontsize, rotation=rotation)
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    sns.despine(bottom=False, left=True)
    ax.grid(False)
    ax.tick_params(bottom=tick_bottom, left=tick_left)
    return None


def map_color(
    df,
    fc_thr,
    pv_thr,
    selected_labels_all,
    my_pal=None,
    feature1=None,
    feature2=None,
    hue_order=None,
    palette=None
):
    '''
    Assigns color categories, point shapes, and sizes for volcano plot visualization.
    It also identifies "important" genes for specialized shape mapping and calculates point sizes based 
    on significance.
    
    Parameters
    ----------
    df : pd.DataFrame
        The differential expression results. Must contain 'log2FoldChange', 
        'nlog10' (negative log10 of the p-value), and 'symbol' columns.
    fc_thr : float
        The fold change threshold used to define 'higher' and 'lower' expression.
    pv_thr : float
        The significance threshold (in -log10 scale).
    selected_labels_all : list or set
        A collection of gene symbols considered "important" to be assigned 
        a specific shape in the plot.
    my_pal : dict, optional
        A dictionary mapping feature names to specific hex or named colors.
    feature1 : str, optional
        The name of the first experimental group (used to pull colors from my_pal).
    feature2 : str, optional
        The name of the second experimental group (used to pull colors from my_pal).
    hue_order : list, optional
        Explicit ordering of color categories for the plot legend.
    palette : list, optional
        List of colors corresponding to the categories in hue_order.

    Returns
    -------
    df : pd.DataFrame
        The input DataFrame with added 'color', 'shape', and 'baseMean' columns.
    hue_order : list
        The final list of categories used for plot coloring.
    palette : list
        The final list of colors mapped to the hue_order.
    '''
    def map_color(row, fc_thrr, pv_thrr):
        fc = row["log2FoldChange"]
        sig = row["nlog10"] >= pv_thrr
        if not sig:
            if fc >= fc_thrr:
                return "higher"
            
            if fc <= -fc_thrr:
                return "lower"
            
            return "no"
        
        if fc >= fc_thrr:
            if fc > 2 * fc_thrr:
                return "very higher"
            
            return "higher"
        
        if fc <= -fc_thrr:
            if fc < -2 * fc_thrr:
                return "very lower"
            
            return "lower"
        
        return "mix"
    
    def map_shape(symbol):
        if symbol in selected_labels_all:
            return "important"
        
        return "not"
    
    # Map color
    df["color"] = df.apply(lambda r: map_color(r, fc_thrr=fc_thr, pv_thrr=pv_thr), axis=1)
    
    # Use my_pal if provided, else default palette
    if my_pal is not None and feature1 in my_pal and feature2 in my_pal:
        color1 = my_pal[feature1]
        color2 = my_pal[feature2]
        if hue_order is None:
            hue_order = ["no", "very higher", "very lower"]
        if palette is None:
            palette = ["lightgrey", color2, color1]
    else:
        if hue_order is None:
            hue_order = ["no", "very higher", "higher", "mix", "very lower", "lower"]
        if palette is None:
            palette = [
                "lightgrey",          # non‑significant / no
                "#d62a2b",            # very higher
                "#D62A2B7A",          # higher (semi‑transparent)
                "lightgrey",          # mix
                "#1f77b4",            # very lower
                "#1F77B47D"           # lower (semi‑transparent)
            ]
    
    # Map shape and baseMean
    df["shape"] = df["symbol"].map(map_shape)
    df["baseMean"] = df["nlog10"] * 10
    
    return df, hue_order, palette


def volcano_plot(
    # input modes (mutually exclusive)
    data=None,                    # DataFrame mode (provide data, symbol_Col, fc_col, pval_col)
    symbol_col=None,
    fc_col=None,
    pval_col=None,
    scores=None,                  # Series/array mode (provide scores, foldchanges, p_values)
    foldchanges=None,
    p_values=None,
    # context
    cell_type=None,
    feature1=None,
    feature2=None,
    # thresholds
    fc_thr=1.0,
    pv_thr=1.0,                   # on -log10 scale
    # labeling
    label_mode="topN",            # "all" or "topN"
    n_p=5,
    n_n=5,
    my_pal=None,                  # optional dict {feature: color}, otherwise use default palette
    hue_order=None,
    palette=None,
    figsize=(20, 20),
    font_size=18,
    marker='o',
    marker_edge_color='w',
    markersize_legend=8,
    font_weight_legend='normal',
    size_legend=12,
    dpi=100,
    output_path=None,
    save_prefix="Volcano",
    save_csv=True,
):
    """
    Generate a unified volcano plot from either a DataFrame or Series/array input.

    Input modes:
      1) DataFrame: data + symbol_col + fc_col + pval_col
      2) Series/array: scores, foldchanges, p_values

    Parameters
    ----------
    data : pd.DataFrame or None, optional
        Input DataFrame containing gene‑level results. If given, symbol_col, fc_col, and pval_col must also be provided.
        If None, scores, foldchanges, and p_values must be provided instead.
    symbol_col : str or None, optional
        Column name in data containing gene symbols.
    fc_col : str or None, optional
        Column name in data containing log2 fold changes.
    pval_col : str or None, optional
        Column name in data containing p‑values or adjusted p‑values.
    scores : pd.Series or array‑like or None, optional
        Index must correspond to gene symbols; used in Series/array mode instead of data.
    foldchanges : array‑like or None, optional
        Array of log2 fold changes, used in Series/array mode.
    p_values : pd.Series or array‑like or None, optional
        Array or Series of p‑values (or adjusted p‑values), used in Series/array mode.
    cell_type : str or None, optional
        Name of the cell type being analyzed, used for CSV and title context.
    feature1 : str or None, optional
        Name of the first group/condition in the comparison (e.g., reference).
    feature2 : str or None, optional
        Name of the second group/condition in the comparison (e.g., test).
    fc_thr : float, optional
        Threshold on log2 fold change (absolute) for significance. Use ±fc_thr.
    pv_thr : float, optional
        Threshold on -log₁₀(p‑value) for significance.
    label_mode : {'all', 'topN'}, optional
        How to select labels for display.
        'all': label all points that pass the thresholds.
        'topN': label only the top n_p and n_n most significant points per side.
    n_p : int, optional
        Number of highest‑significance points on the positive log2FC side to label (used when label_mode == 'topN').
    n_n : int, optional
        Number of highest‑significance points on the negative log2FC side to label (used when label_mode == 'topN').
    my_pal : dict or None, optional
        Optional color palette as a dict {feature1: color1, feature2: color2}.
        If provided and valid, a two‑color scheme is used for the legend instead of the default palette.
    hue_order : list of str or None, optional
        Order of color categories in the legend; if None, sensible defaults are used.
    palette : list of color‑spec or None, optional
        Color palette to use for hue categories; if None, a default palette with semi‑transparent colors is used.
    figsize : tuple of float, optional
        Figure size (width, height) for the plot.
    font_size : int, optional
        Base font size for titles, axis labels, and tick labels.
    marker : str, optional
        Marker style for all points (e.g., 'o', 's').
    marker_edge_color : str, optional
        Color of the marker edge (often 'w' for white edge).
    markersize_legend : int, optional
        Size of legend markers for the color legend.
    font_weight_legend : str, optional
        Font weight for legend text (e.g., 'normal', 'bold').
    size_legend : int, optional
        Font size for legend text.
    dpi : int, optional
        Dots per inch for the saved figure.
    output_path : str or None, optional
        Directory path where the plot (and optionally CSVs) should be saved. If None, only show the plot.
    save_prefix : str, optional
        Prefix used when saving the PDF file (e.g., "Volcano").
    save_csv : bool, optional
        Whether to save CSV files of significant genes per feature (used when output_path is not None).
    """
    # Build DataFrame (log2FoldChange, nlog10, symbol)
    # ------------------------------------------------------------------
    if data is not None:
        if any(v is None for v in (symbol_col, fc_col, pval_col)):
            raise ValueError("In DataFrame mode, symbol_col, fc_col, and pval_col must be provided.")
        df = pd.DataFrame({
            "log2FoldChange": data[fc_col].values,
            "nlog10": -np.log10(data[pval_col].values),
            "symbol": data[symbol_col].values
        })
    else:
        if scores is None or foldchanges is None or p_values is None:
            raise ValueError("In Series mode, scores, foldchanges, and p_values must all be provided.")
        df = pd.DataFrame({
            "log2FoldChange": np.asarray(foldchanges),
            "nlog10": -np.log10(np.asarray(p_values)),
            "symbol": np.asarray(scores.index)
        })
    
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(subset=["nlog10"], how="all", inplace=True)
    
    # Significance masks, CSV export
    # ------------------------------------------------------------------
    sig_mask = (np.abs(df["log2FoldChange"]) >= fc_thr) & (df["nlog10"] >= pv_thr)
    selected_labels_all = df.loc[sig_mask, "symbol"].values
    
    group1_mask = (df["log2FoldChange"] <= -fc_thr) & (df["nlog10"] >= pv_thr)
    group2_mask = (df["log2FoldChange"] >=  fc_thr) & (df["nlog10"] >= pv_thr)
    
    group1_selected = df.loc[group1_mask, "symbol"].values
    group2_selected = df.loc[group2_mask, "symbol"].values
    
    if output_path is not None and save_csv:
        if cell_type is None:
            cell_type = "NA"
        
        if feature1 is None:
            feature1 = "group1"
        
        if feature2 is None:
            feature2 = "group2"
        
        pd.DataFrame(group1_selected, columns=["symbol"]).to_csv(
            Path(output_path) / f"Significant_genes_{cell_type}_{feature1}.csv", index=False
        )
        pd.DataFrame(group2_selected, columns=["symbol"]).to_csv(
            Path(output_path) / f"Significant_genes_{cell_type}_{feature2}.csv", index=False
        )
    
    # Color mapping
    df, hue_order, palette = map_color(
        df=df,
        fc_thr=fc_thr,
        pv_thr=pv_thr,
        selected_labels_all=selected_labels_all,
        my_pal=my_pal,
        feature1=feature1,
        feature2=feature2,
        hue_order=None,      # or pass a custom one if needed
        palette=None
    )
    
    # Label selection
    # ------------------------------------------------------------------
    if label_mode == "all":
        subset_labels = selected_labels_all
    elif label_mode == "topN":
        filtered_df = df.loc[df["nlog10"] >= pv_thr]
        pos_df = filtered_df.loc[filtered_df["log2FoldChange"] >= fc_thr]
        pos_df = pos_df.sort_values(by="nlog10", ascending=False).head(n_p)
        neg_df = filtered_df.loc[filtered_df["log2FoldChange"] <= -fc_thr]
        neg_df = neg_df.sort_values(by="nlog10", ascending=False).head(n_n)
        subset_labels = np.concatenate([pos_df["symbol"].values, neg_df["symbol"].values])
    else:
        raise ValueError("label_mode must be 'all' or 'topN'.")
    
    # Plotting
    # ------------------------------------------------------------------
    plt.figure(figsize=figsize, frameon=False, dpi=dpi)
    plt.style.use("default")
    
    ax = sns.scatterplot(
        data=df,
        x="log2FoldChange",
        y="nlog10",
        hue="color",
        hue_order=hue_order,
        palette=palette,
        style="shape",
        style_order=["not", "important"],
        markers=["o", "o"],
        size="baseMean",
        sizes=(40, 800)
    )
    
    # Lines
    ax.axhline(pv_thr, zorder=0, c="k", lw=2, ls="--")
    ax.axvline(fc_thr,  zorder=0, c="k", lw=2, ls="--")
    ax.axvline(-fc_thr, zorder=0, c="k", lw=2, ls="--")
    
    # Labels
    texts = []
    for i in range(len(df)):
        if df.iloc[i].symbol in subset_labels:
            if df.iloc[i].nlog10 >= pv_thr and df.iloc[i].log2FoldChange >= fc_thr:
                texts.append(
                    plt.text(
                        x=df.iloc[i].log2FoldChange,
                        y=df.iloc[i].nlog10,
                        s=df.iloc[i].symbol,
                        fontsize=font_size,
                        weight="bold",
                        family="sans-serif",
                    )
                )
            if df.iloc[i].nlog10 >= pv_thr and df.iloc[i].log2FoldChange <= -fc_thr:
                texts.append(
                    plt.text(
                        x=df.iloc[i].log2FoldChange,
                        y=df.iloc[i].nlog10,
                        s=df.iloc[i].symbol,
                        fontsize=font_size,
                        weight="bold",
                        family="sans-serif",
                    )
                )
    if texts:
        #adjust_text(texts, arrowprops=dict(arrowstyle="->", color='black', lw=1))
        adjust_text(texts)
    
    # Legend
    if feature1 is None:
        feature1 = "group1"
    
    if feature2 is None:
        feature2 = "group2"
    
    col_pos = palette[1] if len(palette) > 1 else "#d62a2b"
    col_neg = palette[-1] if len(palette) > 1 else "#1f77b4"
    
    custom_lines = [
        Line2D([0], [0], marker=marker, color=marker_edge_color,
               markerfacecolor=col_pos, markersize=markersize_legend),
        Line2D([0], [0], marker=marker, color=marker_edge_color,
               markerfacecolor=col_neg, markersize=markersize_legend),
    ]
    plt.legend(
        custom_lines,
        [f"Higher expressions in {feature2}", f"Higher expressions in {feature1}"],
        loc=1,
        bbox_to_anchor=(1, 1.1),
        frameon=False,
        prop={"weight": font_weight_legend, "size": size_legend},
    )
    
    # Axis styling
    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=2)
    ax.set_ylim(bottom=0)
    
    plt.title(f"Expression Score\n{feature1} - {feature2}", fontsize=font_size)
    plt.xticks(size=font_size, weight="bold")
    plt.yticks(size=font_size, weight="bold")
    plt.xlabel("$log_{2}$ (Fold Change)", size=font_size)
    plt.ylabel("-$log_{10}$ (P-value)", size=font_size)
    
    # Save
    if output_path is not None:
        plt.savefig(
            Path(output_path) / f"{save_prefix} {feature1}-{feature2} FC.pdf",
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
        )
    
    plt.show()


def plot_stats_by_pattern(cluster_names: list = None,
                          all_stats_extend: pd.DataFrame = None,
                          gene_dict: dict = None,
                          pline: list = None,
                          path_to_results: str = None,
                          file_name: str = "/Whole_expressions.csv",
                          font_size: int = 24,p_value=0.01,create_new_plot_folder=False,fc_ther=0.5):
    """
    Parameters
    ----------
    cluster_names : list, optional
        List of cluster names to be considered.
    all_stats_extend : DataFrame, optional
        DataFrame containing stats for each gene and the cluster its belong.
    gene_dict : dict, optional
        Dictionary map gene to cluster its belong.
    pline : list, optional
        Values for axis x to plot curve. 
    path_to_results : str, optional
        path to save result contained Markers and cells.
    file_name : str, optional
        File name same for all clusters have fitting information for genes.
        The default is "/Whole_expressions.csv".
    font_size : int, optional
        Specify font size for labels and names. The default is 12.
    fc_ther: float, optional
        Threshold for FC.

    Returns
    -------
    return figures for each cluster ploting top 4 genes for each pattern.

    """
    plt.rcParams.update({'font.size': font_size})
    # plot results
    for cluster in cluster_names:
        my_data = all_stats_extend[ (all_stats_extend['cluster'] == str(cluster)) & (all_stats_extend['pvalue'] < p_value)]
        my_data['FC']=my_data['FC'].astype(float)
        my_data=my_data[my_data['FC'] > fc_ther]
        sort_my_data = my_data.sort_values(['pvalue'],
                ascending=[True]).groupby('Expression pattern').head(4)
        expression_patterns = np.unique(sort_my_data['Expression pattern'])
        if(len(expression_patterns)):
            
            n_col = 4
            n_row = len(expression_patterns)
            n_px = 10
            
            plt.rcParams["figure.facecolor"] = 'w'
                            
            #plt.figure(figsize=(80, 80))
            #f, axs = plt.subplots(n_row, n_col, figsize=(16, 8))
            plt.figure(figsize=(80, 80))  # Set the overall figure size
            
            # Adjust the size of individual subplots
            subplot_width = 8  # Choose an appropriate value
            subplot_height = 6  # Choose an appropriate value
            
            f, axs = plt.subplots(n_row, n_col, figsize=(n_col * subplot_width, n_row * subplot_height))
            axs = np.atleast_2d(axs)
           
            p = 0
            for pattern in expression_patterns:
                pattern_stats = sort_my_data[sort_my_data['Expression pattern'] == pattern]
                k = 0
                for i, row in pattern_stats.iterrows():
                    cluster_n = row['cluster']
                    gene_name = row['gene']
                    WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve = get_two_tables(gene_name, cluster_n, gene_dict,
                                                                                  file_name, pline, path_to_results)
                    
                    if(KO_x is not None):
                        axs[p, k].scatter([x+0.2 for x in KO_x], KO_tf, alpha = 0.2, c = "tab:gray")
                    
                    axs[p, k].scatter(WT_x, WT_tf, alpha = 0.5, c = WT_tf, cmap = 'viridis',
                              norm = colors.CenteredNorm(np.mean(WT_tf)))
                    axs[p, k].axis(xmin = min(WT_x), xmax = max(WT_x))
                    if(KO_x is not None):
                        axs[p, k].plot(np.linspace(min(WT_x),max(WT_x)), KO_curve,
                                       c = "dimgray", linewidth = 4.0)
                    
                    axs[p, k].plot(np.linspace(min(WT_x),max(WT_x)), WT_curve,
                                   c = "tab:orange", linewidth = 4.0)
                    
                    axs[p, k].set_title(gene_name, size = font_size, weight = 'bold')
                    if( k == 0):
                        plt.rcParams.update({'font.size': font_size-4})
                        axs[p, k].set_ylabel(pattern, size = font_size-4)
                    #else:
                        
                     #   axs[p, k].set_ylabel('Gene expression', size = 16)
                    plt.rcParams.update({'font.size': font_size})
                    for item in (axs[p, k].get_xticklabels() + axs[p, k].get_yticklabels()):
                        item.set_fontsize(font_size)
                    
                    plt.text(.01, .99, 'p-value = %.2e ' % + Decimal(str(row['pvalue'])),
                             ha = 'left', va = 'top',
                             transform=axs[p, k].transAxes, size = font_size)
                    axs[p, k].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    k += 1
                while(k != 4):
                    axs[p, k].set_axis_off()
                    k += 1
                p += 1
            plt.subplots_adjust(wspace = 0.5, hspace = 0.7)
            
            if create_new_plot_folder:
                if not os.path.exists(path_to_results+'/plot_genes_for_'+str(cluster)+'/'):
                    os.makedirs(path_to_results+'/plot_genes_for_'+str(cluster)+'/')
                
                save_path = path_to_results+'/plot_genes_for_'+str(cluster)+'/'+str(cluster) + ".png"
                plt.savefig(save_path)
                plt.close()
            else:   
                if not os.path.exists(path_to_results+'/plots_gene_cluster_differentiation/'):
                    os.makedirs(path_to_results+'/plots_gene_cluster_differentiation/')
                
                save_path = path_to_results+'/plots_gene_cluster_differentiation/'+str(cluster) + ".png"
                plt.savefig(save_path)
                print('Plot for '+str(cluster))
                plt.show()
                plt.close()


def create_sample_grid(df, pattern, image_dir=' figures/', output_path="Results_PILOT/plots/sample_grid_cores.png",zoom=0.1,num_rows = 10):
    """
    Create a grid of sample images based on ordering in a DataFrame.
    
    Args:
        df: DataFrame containing at least:
            - 'sampleID': unique sample identifiers
            - 'order': numerical ordering for sorting
        image_dir: Directory containing images (named as {sampleID}.png)
        output_path: Where to save the resulting grid
    """
    # Sort samples by the specified order
    sorted_df = df.sort_values('Time_score')
    sample_names = sorted_df['sampleID'].tolist()
    
    # Calculate grid dimensions (10 rows, auto columns)
    num_samples = len(sample_names)
    num_rows = num_rows
    num_cols = math.ceil(num_samples / num_rows)
    
    # Prepare figure
    fig, axes = plt.subplots(num_rows, num_cols, 
                           figsize=(num_cols * 2, num_rows * 2))
    axes = axes.flatten()
    
    # Image loading function
    def get_image(path, zoom=zoom):
        if os.path.exists(path):
            return OffsetImage(plt.imread(path), zoom=zoom)
        return None
    
    # Plot each sample
    for i, sample in enumerate(sample_names):
        img_path = os.path.join(image_dir, f"show_{sample}_{pattern}.png")
        ax = axes[i]
        
        # Add image if exists
        img = get_image(img_path)
        if img:
            ab = AnnotationBbox(img, (0.5, 0.5), 
                              frameon=False, 
                              box_alignment=(0.5, 0.5))
            ax.add_artist(ab)
        
        # Minimal labeling
        ax.set_title(sample, fontsize=6)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
    
    # Hide unused subplots
    for i in range(num_samples, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    fig.canvas.draw()  # Force render
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    print(f"Saved grid to {output_path}")
