"""
Storage of the per cell type expression tables (``Results_PILOT/cells``).

``pl.tl.genes_importance`` needs, for each cell type, the expression of its
cells together with the sample and the pseudotime (``Time_score``) of each
cell. Up to version 2.1.1 this table was written as a dense CSV file
(``cells/<cell type>.csv``), which is slow and large. It is now written as a
sparse ``.h5ad`` file (``cells/<cell type>.h5ad``):

- ``X``: expression (cells x genes, sparse),
- ``obs``: ``sampleID`` and ``Time_score`` of each cell.

All functions that read these tables accept both formats, so results made
with an older version (CSV files) keep working without re-running anything.
When a table is read, it is returned in the same layout as the old CSV file:
columns ``sampleID``, ``Time_score`` and then the genes.
"""

import os
from functools import lru_cache

import numpy as np
import pandas as pd
import scipy.sparse as sp

CELL_FORMATS = ('h5ad', 'csv')
META_COLUMNS = ('sampleID', 'Time_score')


def cells_path(path_to_results='Results_PILOT'):
    """Folder of the per cell type tables."""
    return os.path.join(str(path_to_results), 'cells')


def _file(path, name, fmt):
    return os.path.join(str(path), f"{name}.{fmt}")


def cell_table_file(path, name):
    """
    Path of the stored table of ``name`` in the folder ``path``, or None.

    If both a ``.h5ad`` and a ``.csv`` file exist, the newest one is used.
    """
    found = [_file(path, name, fmt) for fmt in CELL_FORMATS
             if os.path.isfile(_file(path, name, fmt))]
    if not found:
        return None
    return max(found, key=os.path.getmtime)


def list_cell_types(path):
    """Names of the cell types stored in the folder ``path`` (both formats)."""
    if not os.path.isdir(str(path)):
        return []
    names = []
    for f in sorted(os.listdir(str(path))):
        stem, ext = os.path.splitext(f)
        if ext.lstrip('.') in CELL_FORMATS and os.path.isfile(os.path.join(str(path), f)) \
                and stem not in names:
            names.append(stem)
    return names


def save_cell_table(path, name, X, genes, sample_ids, time_scores, cells_format='h5ad'):
    """
    Save the table of one cell type.

    Parameters
    ----------
    path : str
        Output folder (created if needed).
    name : str
        Name of the cell type (file name).
    X : array or sparse matrix
        Expression, cells x genes.
    genes : list
        Gene names.
    sample_ids, time_scores : array-like
        Sample and pseudotime of each cell.
    cells_format : {'h5ad', 'csv'}
        'h5ad' (default): sparse and fast. 'csv': the format of PILOT <= 2.1.1.

    Returns
    -------
    str
        Path of the written file.
    """
    import anndata
    os.makedirs(str(path), exist_ok=True)
    genes = [str(g) for g in genes]
    if cells_format == 'h5ad':
        X = sp.csr_matrix(X) if not sp.issparse(X) else X.tocsr()
        obs = pd.DataFrame({'sampleID': np.asarray(sample_ids).astype(str),
                            'Time_score': np.asarray(time_scores, dtype=float)},
                           index=[str(i) for i in range(X.shape[0])])
        var = pd.DataFrame(index=pd.Index(genes))
        out = _file(path, name, 'h5ad')
        anndata.AnnData(X=X, obs=obs, var=var).write_h5ad(out)
    elif cells_format == 'csv':
        dense = X.toarray() if sp.issparse(X) else np.asarray(X)
        df = pd.DataFrame(dense, columns=genes)
        df.insert(0, 'Time_score', np.asarray(time_scores, dtype=float))
        df.insert(0, 'sampleID', np.asarray(sample_ids))
        out = _file(path, name, 'csv')
        df.to_csv(out)
    else:
        raise ValueError("cells_format must be 'h5ad' or 'csv'.")
    _read_h5ad.cache_clear()
    return out


@lru_cache(maxsize=4)
def _read_h5ad(file, mtime):
    """Read (and keep in memory) a stored table; ``mtime`` invalidates the cache."""
    import anndata
    return anndata.read_h5ad(file)


def load_cell_table(path, name, columns=None, index_col=True):
    """
    Read the table of one cell type (``.h5ad`` or old ``.csv``).

    Parameters
    ----------
    path : str
        Folder of the tables (e.g. 'Results_PILOT/cells').
    name : str
        Name of the cell type.
    columns : list, optional
        Only these columns (genes, 'sampleID', 'Time_score'). They are
        returned in the order of the file (like ``pd.read_csv(usecols=...)``).
        A missing column raises a ValueError.
    index_col : bool, optional
        Keep the first column of an old CSV file as index (True, like
        ``loadTarget``). Ignored for ``.h5ad`` files.

    Returns
    -------
    DataFrame or None
        None if the cell type has no stored table.
    """
    file = cell_table_file(path, name)
    if file is None:
        return None
    if file.endswith('.csv'):
        if columns is not None:
            return pd.read_csv(file, usecols=list(columns))
        return pd.read_csv(file, index_col=0 if index_col else None)
    adata = _read_h5ad(file, os.path.getmtime(file))
    genes = list(adata.var_names)
    if columns is None:
        meta = [c for c in adata.obs.columns]
        wanted_genes = genes
    else:
        columns = [str(c) for c in columns]
        known = set(genes) | set(adata.obs.columns)
        missing = [c for c in columns if c not in known]
        if missing:
            raise ValueError(f"Usecols do not match columns, columns expected but not found: {missing}")
        meta = [c for c in adata.obs.columns if c in columns]
        chosen = set(columns)
        wanted_genes = [g for g in genes if g in chosen]
    X = adata[:, wanted_genes].X if len(wanted_genes) else np.zeros((adata.n_obs, 0))
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    df = pd.DataFrame(X.astype(float), columns=wanted_genes)
    for i, c in enumerate(meta):
        values = adata.obs[c].values
        df.insert(i, c, values.astype(float) if c == 'Time_score' else values.astype(object))
    return df
