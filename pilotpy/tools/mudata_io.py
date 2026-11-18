"""
MuData support: save the PILOT results at the patient (sample) level, and use
MuData objects as input.

A MuData object holds several AnnData objects ("modalities") in one file
(``.h5mu``). ``to_mudata`` creates one with

- 'cells': the single-cell data (one row per cell, as before), and
- 'samples': one row per patient/sample with the PILOT results: the cell type
  proportions (X), the clinical variables, status, pseudotime, trajectory
  position and clusters (obs), the trajectory and UMAP coordinates (obsm), the
  Wasserstein distances between the samples (obsp) and the cost matrix
  between the cell types (varp).

The 'samples' modality is a normal AnnData, so other tools (scanpy, sklearn,
...) can use the patient-level results directly. ``from_mudata`` reads such an
object back, so the PILOT analysis can be continued without computing the
distances again.
"""

import numpy as np
import pandas as pd

__all__ = ['sample_adata', 'to_mudata', 'from_mudata', 'get_modality']

#: PILOT keys of ``adata.uns`` that are moved to the 'samples' modality.
PILOT_UNS_KEYS = ('data', 'annot', 'proportions', 'cost', 'EMD', 'EMD_df', 'real_labels',
                  'sample_metadata', 'wasserstein_params', 'embedding', 'trajectory_params',
                  'umap_samples', 'umap_samples_params', 'principal_graph',
                  'principal_graph_params', 'pseudotime', 'root_node', 'branch',
                  'cell_importance', 'cellnames', 'orders', 'sample_colors', 'Sil', 'ARI',
                  'trajectory_associations', 'trajectory_associations_skipped',
                  'trajectory_associations_values', 'cluster_associations',
                  'cluster_associations_skipped', 'cluster_associations_values',
                  'sample_clusters', 'cluster_names', 'best_res')


def _require_mudata():
    try:
        import mudata
    except ImportError as err:
        raise ImportError("This function needs the 'mudata' package: pip install mudata") from err
    return mudata


def _plain(value):
    """Values that can be written to an .h5ad/.h5mu file (dicts, arrays, strings, numbers)."""
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items() if v is not None}
    if isinstance(value, (list, tuple)):
        try:
            return np.asarray(value)
        except Exception:
            return np.asarray([str(v) for v in value])
    if isinstance(value, (np.ndarray, str, int, float, bool, np.integer, np.floating, np.bool_)):
        return value
    return str(value)


def _clean_frame(df):
    """DataFrame with columns that can be written to .h5ad (objects as strings)."""
    df = df.copy()
    df.columns = [str(c) for c in df.columns]
    for col in df.columns:
        if df[col].dtype == object:
            df[col] = df[col].map(lambda v: v if pd.isna(v) else str(v))
    return df


def get_modality(data, modality=None):
    """
    Return the AnnData to work with: ``data`` itself, or one modality of a MuData.

    Parameters
    ----------
    data : AnnData or MuData
        Input object.
    modality : str, optional
        Name of the modality (e.g. 'rna'). If None: the only modality, else
        'cells', 'rna' or 'gex' if present.

    Returns
    -------
    AnnData
        The AnnData object (not a copy: PILOT results written to it are kept
        in the MuData). Columns of the global ``mdata.obs`` (e.g. clinical
        variables) that are missing in the modality are copied to it.
    """
    if not hasattr(data, 'mod'):
        return data
    names = list(data.mod.keys())
    if modality is None:
        candidates = [m for m in ('cells', 'rna', 'gex', 'RNA', 'GEX') if m in names]
        if len(names) == 1:
            modality = names[0]
        elif candidates:
            modality = candidates[0]
        else:
            raise ValueError(f"The MuData has several modalities {names}; choose one with "
                             "modality='...'.")
    if modality not in names:
        raise KeyError(f"'{modality}' is not a modality of the MuData. Available: {names}")
    adata = data.mod[modality]
    glob = data.obs
    prefix = f"{modality}:"
    for col in glob.columns:
        name = col[len(prefix):] if col.startswith(prefix) else col
        if ':' in name or name in adata.obs.columns:
            continue
        values = glob[col].reindex(adata.obs_names)
        if values.notna().any():
            adata.obs[name] = values.values
    return adata


def sample_adata(adata):
    """
    The PILOT results at the patient (sample) level as an AnnData.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance`` (and optionally
        the trajectory, principal graph, cell importance and clustering).

    Returns
    -------
    AnnData
        One row per sample (same order as the Wasserstein distances):

        - ``X``: cell type proportions (samples x cell types);
        - ``obs``: status and sample-level variables, and if available
          'pseudotime', 'trajectory_order', 'graph_node', 'in_branch' and
          'cluster';
        - ``obsm``: 'X_trajectory' (diffusion map), 'X_umap' (sample UMAP);
        - ``obsp``: 'wasserstein' (distances between the samples);
        - ``varp``: 'cost' (cost matrix between the cell types);
        - ``uns['pilot']``: the settings and the other PILOT results.
    """
    import anndata
    from ..plot.ploting import _sample_order, _sample_table, _sample_nodes, _branch_mask
    if 'EMD' not in adata.uns:
        raise ValueError("No PILOT results found. Run pl.tl.wasserstein_distance(adata, ...) first.")
    order = [str(s) for s in _sample_order(adata)]
    cell_types = list(pd.unique(adata.uns['annot']['cell_type']))
    props = np.vstack([np.asarray(v, dtype=float) for v in adata.uns['proportions'].values()])
    obs = _sample_table(adata).drop(columns=['cell_type'], errors='ignore')
    obs.index = order
    obs = obs.drop(columns=['sampleID'], errors='ignore')
    uns = adata.uns
    if 'pseudotime' in uns:
        obs['pseudotime'] = np.asarray(uns['pseudotime'], dtype=float)
    if 'principal_graph' in uns and 'embedding' in uns:
        obs['graph_node'] = _sample_nodes(adata)
    if 'branch' in uns:
        obs['in_branch'] = _branch_mask(adata)
    if 'orders' in uns:
        orders = uns['orders']
        obs['trajectory_order'] = pd.Series(orders['Time_score'].values,
                                            index=orders['sampleID'].astype(str)).reindex(order).values
    if 'sample_clusters' in uns:
        cl = uns['sample_clusters']
        obs['cluster'] = pd.Series(cl['Predicted_Labels'].astype(str).values,
                                   index=cl['sampleID'].astype(str)).reindex(order).values
    samples = anndata.AnnData(X=props, obs=_clean_frame(obs),
                              var=pd.DataFrame(index=[str(c) for c in cell_types]))
    samples.obsp['wasserstein'] = np.asarray(uns['EMD'], dtype=float)
    cost = uns.get('cost', None)
    if isinstance(cost, pd.DataFrame):
        samples.varp['cost'] = cost.reindex(index=cell_types, columns=cell_types).values.astype(float)
    if 'embedding' in uns:
        samples.obsm['X_trajectory'] = np.asarray(uns['embedding'], dtype=float)
    if 'umap_samples' in uns:
        samples.obsm['X_umap'] = np.asarray(uns['umap_samples'], dtype=float)
    pilot = {}
    for key in ('wasserstein_params', 'trajectory_params', 'umap_samples_params',
                'principal_graph', 'principal_graph_params', 'root_node', 'branch',
                'sample_colors', 'cellnames', 'cluster_names'):
        if key in uns and uns[key] is not None:
            pilot[key] = _plain(uns[key])
    for key in ('trajectory_associations', 'cluster_associations'):
        if isinstance(uns.get(key, None), pd.DataFrame):
            pilot[key] = _clean_frame(uns[key])
    res = uns.get('cell_importance', None)
    if isinstance(res, dict) and 'sorted_best' in res:
        rows = []
        for name, r in res['sorted_best'].items():
            # r = [fitted function, model results, pattern, slope, p-value]
            fit = r[1] if isinstance(r[1], dict) else {}
            rows.append({'cell_type': str(name), 'fitted_function': str(r[0]),
                         'pattern': str(r[2]), 'slope': float(np.ravel(r[3])[0]),
                         'p_value': float(np.min(np.ravel(r[4]))),
                         'r_squared_adj': float(fit.get('rsquared_adj', np.nan))})
        pilot['cell_importance'] = pd.DataFrame(rows)
    samples.uns['pilot'] = pilot
    return samples


def to_mudata(adata, include_cells=True, cells_key='cells', samples_key='samples', path=None):
    """
    Save the PILOT results as a MuData object with a patient-level modality.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance`` (and optionally
        the next PILOT steps).
    include_cells : bool, optional
        Also put the single-cell data in the MuData (modality ``cells_key``),
        by default True. With False only the patient level is kept (small).
    cells_key, samples_key : str, optional
        Names of the two modalities, by default 'cells' and 'samples'.
    path : str, optional
        If given, the MuData is also written to this ``.h5mu`` file.

    Returns
    -------
    MuData
        With the modalities 'cells' (a copy of ``adata`` without the PILOT
        results in ``uns``, they are in 'samples') and 'samples' (see
        ``pl.tl.sample_adata``).

    Example
    -------
    >>> mdata = pl.tl.to_mudata(adata, path='PILOT_results.h5mu')
    >>> mdata['samples']            # one row per patient
    >>> adata = pl.tl.from_mudata('PILOT_results.h5mu')
    """
    mudata = _require_mudata()
    samples = sample_adata(adata)
    mods = {samples_key: samples}
    if include_cells:
        cells = adata.copy()
        for key in PILOT_UNS_KEYS:
            cells.uns.pop(key, None)
        mods = {cells_key: cells, samples_key: samples}
    mdata = mudata.MuData(mods)
    mdata.uns['pilot_modalities'] = {'cells': cells_key if include_cells else '',
                                     'samples': samples_key}
    if path is not None:
        mdata.write(str(path))
    return mdata


def from_mudata(mdata, modality=None, samples_key=None):
    """
    Read a MuData object (or ``.h5mu`` file) for PILOT.

    - If it was made by ``pl.tl.to_mudata``, the cells are returned with the
      PILOT results restored (Wasserstein distances, trajectory, principal
      graph, pseudotime, ...). The analysis can be continued directly, e.g.
      with ``pl.pl.trajectory``, ``pl.tl.cell_importance`` or
      ``pl.tl.correlation_with_trajectory``, without computing the distances
      again.
    - Otherwise the chosen modality (e.g. 'rna') is returned, with the
      columns of the global ``mdata.obs`` added, ready for
      ``pl.tl.wasserstein_distance``.

    Parameters
    ----------
    mdata : MuData or str
        MuData object or path of a ``.h5mu`` file.
    modality : str, optional
        Modality with the cells. If None, chosen automatically.
    samples_key : str, optional
        Modality with the PILOT sample-level results. If None, the one written
        by ``to_mudata`` (if any).

    Returns
    -------
    AnnData
        The cells, with the PILOT results in ``uns`` when available.
    """
    mudata = _require_mudata()
    if isinstance(mdata, str):
        mdata = mudata.read_h5mu(mdata)
    names = dict(mdata.uns.get('pilot_modalities', {}))
    if samples_key is None:
        samples_key = names.get('samples', None)
    if modality is None and names.get('cells'):
        modality = names['cells']
    if samples_key and modality is None:
        others = [m for m in mdata.mod.keys() if m != samples_key]
        if len(others) == 1:
            modality = others[0]
    if samples_key and samples_key in mdata.mod and (modality is None or modality not in mdata.mod):
        raise ValueError("This MuData has only the patient-level results (made with "
                         "include_cells=False); use mdata['" + str(samples_key) + "'] directly.")
    adata = get_modality(mdata, modality)
    if samples_key and samples_key in mdata.mod:
        _restore(adata, mdata.mod[samples_key])
    return adata


def _restore(adata, samples):
    """Put the PILOT results of the 'samples' modality back into ``adata.uns``."""
    from .Trajectory import (extract_data_anno_scRNA_from_h5ad, extract_data_anno_pathomics_from_h5ad,
                             Cluster_Representations, return_real_labels, _sample_metadata)
    pilot = dict(samples.uns.get('pilot', {}))
    params = dict(pilot.get('wasserstein_params', {}))
    if not params:
        raise ValueError("The 'samples' modality has no PILOT settings (uns['pilot']).")
    for key in ('regulizer', 'reg', 'reg_m', 'seconds_cost', 'seconds_transport'):
        if key in params:
            params[key] = float(params[key])
    params['normalization'] = str(params.get('normalization', 'True')) in ('True', '1', 'true')
    if params.get('data_type', 'scRNA') == 'scRNA':
        data, annot = extract_data_anno_scRNA_from_h5ad(
            adata, emb_matrix=params['emb_matrix'], clusters_col=params['clusters_col'],
            sample_col=params['sample_col'], status=params['status'])
    else:
        data, annot = extract_data_anno_pathomics_from_h5ad(
            adata, var_names=list(adata.var_names), clusters_col=params['clusters_col'],
            sample_col=params['sample_col'], status=params['status'])
    n_pcs = params.get('n_pcs', 'all')
    if str(n_pcs) != 'all':
        data = data.iloc[:, :int(n_pcs)]
    proportions = Cluster_Representations(annot, regulizer=params.get('regulizer', 0.2),
                                          normalization=params['normalization'])
    order = [str(s) for s in samples.obs_names]
    if [str(s) for s in proportions.keys()] != order:
        raise ValueError("The samples of the cells and of the 'samples' modality do not match.")
    EMD = np.asarray(samples.obsp['wasserstein'], dtype=float)
    adata.uns['data'] = data
    adata.uns['annot'] = annot
    adata.uns['proportions'] = proportions
    adata.uns['EMD'] = EMD
    adata.uns['EMD_df'] = pd.DataFrame(EMD, index=list(proportions.keys()),
                                       columns=list(proportions.keys()))
    if 'cost' in samples.varp:
        names = list(samples.var_names)
        adata.uns['cost'] = pd.DataFrame(np.asarray(samples.varp['cost']), index=names,
                                         columns=names)
    adata.uns['real_labels'] = return_real_labels(annot)
    adata.uns['sample_metadata'] = _sample_metadata(adata, params['sample_col'],
                                                    list(proportions.keys()), 'auto')
    adata.uns['wasserstein_params'] = params
    if 'X_trajectory' in samples.obsm:
        adata.uns['embedding'] = np.asarray(samples.obsm['X_trajectory'], dtype=float)
    if 'X_umap' in samples.obsm:
        adata.uns['umap_samples'] = np.asarray(samples.obsm['X_umap'], dtype=float)
    for key in ('trajectory_params', 'umap_samples_params', 'principal_graph_params',
                'sample_colors', 'cluster_names'):
        if key in pilot:
            adata.uns[key] = dict(pilot[key])
    if 'principal_graph' in pilot:
        g = pilot['principal_graph']
        adata.uns['principal_graph'] = {'NodePositions': np.asarray(g['NodePositions'], dtype=float),
                                        'Edges': np.asarray(g['Edges'], dtype=int)}
    if 'pseudotime' in samples.obs:
        adata.uns['pseudotime'] = samples.obs['pseudotime'].values.astype(float)
    if 'root_node' in pilot:
        adata.uns['root_node'] = int(pilot['root_node'])
    if 'branch' in pilot:
        b = pilot['branch']
        adata.uns['branch'] = {'name': str(b['name']),
                               'nodes': [int(x) for x in np.atleast_1d(b['nodes'])],
                               'samples': [str(x) for x in np.atleast_1d(b['samples'])]}
    if 'cluster' in samples.obs:
        adata.uns['sample_clusters'] = pd.DataFrame({'sampleID': order,
                                                     'Predicted_Labels': samples.obs['cluster'].astype(str).values})
    return adata
