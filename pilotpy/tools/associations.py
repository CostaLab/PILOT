"""
Automatic tests of the sample variables (clinical, technical, batch) against
the PILOT trajectory or the PILOT clusters of samples.

The user does not need to say which variables are categorical or numeric:
all columns of ``adata.obs`` are checked, the usable ones are kept, their type
is detected, and every variable is tested with the right test. Columns that
cannot be used (one value only, a different value for every sample such as
IDs, too many missing values, ...) are skipped and the reason is reported.
"""

import warnings

import numpy as np
import pandas as pd
from scipy.stats import f_oneway, kruskal, spearmanr, chi2_contingency
from statsmodels.stats.multitest import multipletests

from .Trajectory import _sample_metadata
from ..plot.ploting import _sample_order


__all__ = ['sample_variables', 'correlation_with_trajectory', 'correlation_with_clustering',
           'name_clusters', 'MISSING_VALUES']

#: Values that are treated as missing (case is ignored).
MISSING_VALUES = {'', 'nan', 'na', 'n/a', 'none', 'null', 'unknown', 'not reported',
                  'not available', 'not applicable', '-', '?'}


def _params(adata):
    params = adata.uns.get('wasserstein_params', None)
    return params if isinstance(params, dict) else {}


def _default_sample_col(adata, sample_col):
    if sample_col is not None:
        return sample_col
    sample_col = _params(adata).get('sample_col', None)
    if sample_col is None:
        raise ValueError("Please give sample_col (the column of adata.obs with the samples).")
    return sample_col


def _clean_missing(values):
    """Values as objects, with the usual missing codes ('unknown', 'NA', ...) as NaN."""
    values = pd.Series(values, dtype=object).copy()
    text = values.astype(str).str.strip().str.lower()
    values[values.isna() | text.isin(MISSING_VALUES)] = np.nan
    return values


def _as_number(values):
    """Numbers of a variable, NaN where a value is not a number."""
    if pd.api.types.is_bool_dtype(values):
        return pd.Series(np.nan, index=values.index)
    return pd.to_numeric(values, errors='coerce')


def _same_grouping(a, b, kind_a='categorical', kind_b='categorical'):
    """
    True if two variables carry the same information: the same groups
    (categorical) or the same order of the samples (numeric).
    """
    if kind_a != kind_b:
        return False
    if kind_a == 'numeric':
        mask = a.notna() & b.notna()
        if mask.sum() < 3 or mask.sum() != a.notna().sum() or mask.sum() != b.notna().sum():
            return False
        ra = a[mask].astype(float).rank().values
        rb = b[mask].astype(float).rank().values
        return bool(np.array_equal(ra, rb) or np.array_equal(ra, len(rb) + 1 - rb))
    mask = a.notna() & b.notna()
    if mask.sum() == 0 or mask.sum() != a.notna().sum() or mask.sum() != b.notna().sum():
        return False
    pairs = pd.DataFrame({'a': a[mask].astype(str), 'b': b[mask].astype(str)}).drop_duplicates()
    return pairs['a'].is_unique and pairs['b'].is_unique


def sample_variables(adata,
                     sample_col=None,
                     features='auto',
                     include_cell_level=True,
                     agg='mean',
                     max_missing=0.5,
                     min_group_size=2,
                     max_levels=None,
                     numeric_min_unique=6,
                     exclude=None,
                     samples=None,
                     max_mixed=0.5):
    """
    Find the sample variables that can be tested and detect their type.

    Every column of ``adata.obs`` is checked (or only ``features``):

    - a column with one value per sample is a sample variable (e.g. sex, age,
      disease, batch, hospital); if only some samples have more than one
      value (``max_mixed``), the most frequent value of each sample is used;
    - a numeric column with different values inside a sample (e.g. number of
      genes or % mitochondrial reads per cell) is summarised per sample with
      ``agg`` if ``include_cell_level`` is True;
    - other columns (e.g. the cell type) are skipped.

    Then the type is detected: numbers with at least ``numeric_min_unique``
    different values are 'numeric', everything else is 'categorical'. Values
    such as 'unknown', 'NA' or '' are treated as missing. A variable is
    skipped if it has only one value, a different value for every sample
    (IDs), too many missing values, too many groups, or if it splits the
    samples exactly like a variable already kept (e.g. 'disease' and
    'disease_ontology_term_id').

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    sample_col : str, optional
        Column of ``adata.obs`` with the samples. If None, the one used in
        ``pl.tl.wasserstein_distance``.
    features : 'auto' or list, optional
        'auto' (default): all columns. A list: only these columns (they are
        still checked and their type is still detected).
    include_cell_level : bool, optional
        Also use numeric cell-level columns, summarised per sample with
        ``agg``, by default True.
    agg : str, optional
        How a numeric cell-level column is summarised, e.g. 'mean' (default)
        or 'median'.
    max_missing : float, optional
        Skip a variable if more than this fraction of the samples has no
        value, by default 0.5.
    min_group_size : int, optional
        Groups of a categorical variable with fewer samples are not tested
        (their samples are treated as missing), by default 2.
    max_levels : int, optional
        Skip a categorical variable with more groups. If None, a variable is
        skipped when it has more groups than half of its samples.
    numeric_min_unique : int, optional
        A numeric variable needs at least this many different values to be
        treated as numeric; otherwise (e.g. 0/1, stage 1-4) it is treated as
        categorical, by default 6.
    exclude : list, optional
        Columns that should not be used.
    samples : list, optional
        Samples to use (default: all samples of the PILOT results).
    max_mixed : float, optional
        A categorical column with more than one value in some samples (e.g.
        the tissue, when a few patients have biopsies from two regions) is
        used with the most frequent value of each sample if at most this
        fraction of the samples has more than one value, by default 0.5.
        Columns that differ in most samples (e.g. the cell type) are skipped.

    Returns
    -------
    values : DataFrame
        One row per sample, one column per usable variable.
    info : DataFrame
        One row per checked column: 'Feature', 'Type' ('numeric',
        'categorical' or 'skipped'), 'Level' ('sample' or 'cell (<agg>)'),
        'N_samples', 'N_groups' and 'Reason' (why it was skipped).
    """
    sample_col = _default_sample_col(adata, sample_col)
    if sample_col not in adata.obs.columns:
        raise KeyError(f"'{sample_col}' is not a column of adata.obs.")
    if samples is None:
        samples = [str(s) for s in _sample_order(adata)]
    samples = [str(s) for s in samples]
    obs = adata.obs.copy()
    obs[sample_col] = obs[sample_col].astype(str)
    params = _params(adata)
    skip = {sample_col}
    if exclude is not None:
        skip |= set([exclude] if isinstance(exclude, str) else exclude)
    if isinstance(features, str) and features == 'auto':
        candidates = [c for c in obs.columns if c not in skip]
        # Readable columns first, so that e.g. 'sex' is kept and
        # 'sex_ontology_term_id' is reported as its duplicate.
        status = params.get('status', None)

        def _rank(c):
            name = str(c).lower()
            id_like = ('ontology' in name or name == 'id' or name.endswith('_id')
                       or name.endswith('.id'))
            return (c != status, id_like)

        # The status first, then readable columns, then IDs, so that e.g.
        # 'sex' is kept and 'sex_ontology_term_id' is reported as its copy.
        candidates = sorted(candidates, key=_rank)
        given = False
    else:
        candidates = [features] if isinstance(features, str) else list(features)
        missing = [c for c in candidates if c not in obs.columns]
        if missing:
            raise KeyError(f"These columns are not in adata.obs: {missing}")
        given = True
    # Columns with one value per sample.
    groups = obs.groupby(sample_col, observed=True)
    columns, rows, kept_values, kinds = {}, [], [], {}
    for col in candidates:
        level = 'sample'
        try:
            n_unique = groups[col].nunique(dropna=True)
        except TypeError:
            rows.append(dict(Feature=col, Type='skipped', Level='', N_samples=0, N_groups=0,
                             Reason='values cannot be compared'))
            continue
        column = obs[col]
        numeric_cells = (pd.api.types.is_numeric_dtype(column)
                         and not pd.api.types.is_bool_dtype(column))
        mixed = float((n_unique > 1).mean())
        if (n_unique <= 1).all():
            values = groups[col].first().astype(object)
        elif not numeric_cells and mixed <= max_mixed:
            # Mostly one value per sample (e.g. a few patients with biopsies
            # from two regions): use the most frequent value of each sample.
            values = groups[col].agg(lambda x: x.value_counts().index[0]
                                     if x.notna().any() else np.nan).astype(object)
            level = 'sample (most frequent)'
        else:
            if numeric_cells and include_cell_level:
                values = groups[col].agg(agg).astype(object)
                level = f'cell ({agg})'
            else:
                reason = ('differs between the cells of a sample (cell-level)'
                          if not numeric_cells else 'cell-level (include_cell_level=False)')
                rows.append(dict(Feature=col, Type='skipped', Level='cell', N_samples=0,
                                 N_groups=0, Reason=reason))
                continue
        values.index = values.index.astype(str)
        values = _clean_missing(values.reindex(samples))
        n_valid = int(values.notna().sum())
        row = dict(Feature=col, Type='skipped', Level=level, N_samples=n_valid, N_groups=0,
                   Reason='')
        if n_valid == 0 or 1 - n_valid / len(samples) > max_missing:
            row['Reason'] = f'too many missing values ({len(samples) - n_valid} of {len(samples)})'
            rows.append(row)
            continue
        numbers = _as_number(values)
        n_numbers = int(numbers.notna().sum())
        is_number = n_numbers >= max(3, 0.9 * n_valid)
        if is_number and numbers.nunique() >= numeric_min_unique:
            if n_numbers < n_valid:
                row['Reason'] = f'{n_valid - n_numbers} non-numeric values ignored'
            values = numbers.astype(float)
            row.update(Type='numeric', N_samples=n_numbers, N_groups=int(numbers.nunique()))
        else:
            if is_number and n_numbers < n_valid:
                values = numbers.astype(object)
            n_distinct = int(values.nunique())
            limit = max_levels if max_levels is not None else max(2, n_valid // 2)
            counts = values.value_counts()
            small = counts[counts < min_group_size].index
            if len(small):
                values = values.where(~values.isin(small), np.nan)
                counts = values.value_counts()
            n_levels = len(counts)
            row.update(N_samples=int(counts.sum()), N_groups=n_levels)
            if n_distinct <= 1:
                row['Reason'] = 'only one value'
            elif n_distinct > limit:
                row['Reason'] = (f'too many different values ({n_distinct} for {n_valid} '
                                 'samples, e.g. an ID)')
            elif n_levels < 2:
                row['Reason'] = f'fewer than 2 groups with >= {min_group_size} samples'
            else:
                row['Type'] = 'categorical'
                if len(small):
                    row['Reason'] = f'{len(small)} group(s) with < {min_group_size} samples not tested'
            if row['Type'] == 'skipped':
                rows.append(row)
                continue
        duplicate = None
        if not given:
            for other in kept_values:
                if _same_grouping(values, columns[other], row['Type'], kinds[other]):
                    duplicate = other
                    break
        if duplicate is not None:
            row.update(Type='skipped', Reason=f"same groups as '{duplicate}'")
            rows.append(row)
            continue
        columns[col] = values
        kinds[col] = row['Type']
        kept_values.append(col)
        rows.append(row)
    table = pd.DataFrame(columns, index=pd.Index(samples, name='sampleID'))
    info = pd.DataFrame(rows, columns=['Feature', 'Type', 'Level', 'N_samples', 'N_groups',
                                       'Reason'])
    return table, info


def _adjust(p_values, method='fdr_bh'):
    p = np.asarray(p_values, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method=method)[1]
    return out


def _trajectory_position(adata):
    """Position of each sample on the trajectory (as used by the old functions)."""
    if 'orders' in adata.uns:
        orders = adata.uns['orders']
        return pd.Series(orders['Time_score'].values.astype(float),
                         index=orders['sampleID'].astype(str).values)
    if 'pseudotime' in adata.uns:
        return pd.Series(np.asarray(adata.uns['pseudotime'], dtype=float),
                         index=[str(s) for s in _sample_order(adata)])
    raise ValueError("No trajectory found. Run pl.pl.fit_pricipla_graph(adata) (pseudotime) "
                     "or pl.tl.cell_importance(adata) first.")


def correlation_with_trajectory(adata,
                                sample_col=None,
                                features='auto',
                                include_cell_level=True,
                                agg='mean',
                                max_missing=0.5,
                                min_group_size=2,
                                max_levels=None,
                                numeric_min_unique=6,
                                exclude=None,
                                max_mixed=0.5,
                                p_adjust='fdr_bh',
                                show_skipped=True):
    """
    Test all usable sample variables against the trajectory, in one table.

    You do not need to list the variables or say their type: the columns of
    ``adata.obs`` are checked automatically (see ``pl.tl.sample_variables``)
    and each usable variable is tested:

    - categorical (e.g. sex, hospital, batch, disease): ANOVA (and
      Kruskal-Wallis) of the trajectory position between the groups;
    - numeric (e.g. age, BMI, eGFR, % mitochondrial reads): Spearman
      correlation with the trajectory position.

    The trajectory position is the one of ``pl.tl.cell_importance`` (sample
    order) if available, otherwise the pseudotime. If a branch was selected
    (``pl.tl.select_branch``) only its samples are used.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.pl.fit_pricipla_graph`` or
        ``pl.tl.cell_importance``.
    sample_col : str, optional
        Column with the samples; default: the one of ``wasserstein_distance``.
    features : 'auto' or list, optional
        'auto' (default) checks all columns; or a list of columns.
    include_cell_level, agg, max_missing, min_group_size, max_levels, numeric_min_unique, exclude, max_mixed
        How the variables are chosen; see ``pl.tl.sample_variables``.
    p_adjust : str, optional
        Correction for multiple testing (``statsmodels`` method), by default
        'fdr_bh' (Benjamini-Hochberg).
    show_skipped : bool, optional
        Print the skipped columns and why, by default True.

    Returns
    -------
    DataFrame
        One row per tested variable, sorted by p-value: 'Feature', 'Type',
        'Level', 'Test', 'Statistic', 'P_value', 'Adjusted_P_value',
        'N_samples', 'N_groups', 'Note'. Also stored in
        ``adata.uns['trajectory_associations']``; the skipped columns are in
        ``adata.uns['trajectory_associations_skipped']``.

    Example
    -------
    >>> res = pl.tl.correlation_with_trajectory(adata)
    >>> res[res['P_value'] < 0.05]
    >>> pl.pl.trajectory_associations(adata)
    """
    position = _trajectory_position(adata)
    table, info = sample_variables(adata, sample_col=sample_col, features=features,
                                   include_cell_level=include_cell_level, agg=agg,
                                   max_missing=max_missing, min_group_size=min_group_size,
                                   max_levels=max_levels, numeric_min_unique=numeric_min_unique,
                                   exclude=exclude, samples=list(position.index),
                                   max_mixed=max_mixed)
    rows = []
    for _, meta in info[info['Type'] != 'skipped'].iterrows():
        feature = meta['Feature']
        values = table[feature]
        ok = values.notna() & position.notna()
        x, t = values[ok], position[ok]
        row = dict(Feature=feature, Type=meta['Type'], Level=meta['Level'], Test='',
                   Statistic=np.nan, P_value=np.nan, N_samples=int(ok.sum()), N_groups=0,
                   Note=meta['Reason'])
        if meta['Type'] == 'numeric':
            row.update(Test='Spearman', N_groups=int(x.nunique()))
            if len(x) >= 3 and np.ptp(t.values) > 0 and np.ptp(x.values.astype(float)) > 0:
                stat, p = spearmanr(t.values, x.values.astype(float))
                row.update(Statistic=float(stat), P_value=float(p))
        else:
            groups = [t[x == g].values for g in pd.unique(x)]
            row.update(Test='ANOVA', N_groups=len(groups))
            if len(groups) >= 2 and sum(len(g) for g in groups) > len(groups):
                stat, p = f_oneway(*groups)
                row.update(Statistic=float(stat), P_value=float(p))
                try:
                    row['Kruskal_P_value'] = float(kruskal(*groups).pvalue)
                except ValueError:
                    row['Kruskal_P_value'] = np.nan
        rows.append(row)
    columns = ['Feature', 'Type', 'Level', 'Test', 'Statistic', 'P_value', 'Adjusted_P_value',
               'Kruskal_P_value', 'N_samples', 'N_groups', 'Note']
    res = pd.DataFrame(rows)
    if res.empty:
        res = pd.DataFrame(columns=columns)
    else:
        res['Adjusted_P_value'] = _adjust(res['P_value'], p_adjust)
        res = res.reindex(columns=columns)
        res = res.sort_values('P_value', na_position='last', kind='stable').reset_index(drop=True)
    skipped = info[info['Type'] == 'skipped'].reset_index(drop=True)
    adata.uns['trajectory_associations'] = res
    adata.uns['trajectory_associations_skipped'] = skipped
    adata.uns['trajectory_associations_values'] = table
    if show_skipped:
        _print_summary(res, skipped)
    return res


def _print_summary(res, skipped):
    n_cat = int((res['Type'] == 'categorical').sum()) if len(res) else 0
    n_num = int((res['Type'] == 'numeric').sum()) if len(res) else 0
    print(f"Tested {len(res)} variables ({n_cat} categorical, {n_num} numeric); "
          f"skipped {len(skipped)} columns.")
    if len(skipped):
        by_reason = skipped.groupby('Reason')['Feature'].apply(list)
        for reason, feats in by_reason.items():
            shown = ', '.join(map(str, feats[:8])) + (', ...' if len(feats) > 8 else '')
            print(f"  skipped ({reason}): {shown}")


def _cluster_labels(adata, proportion_df=None):
    """Cluster of each sample, from ``proportion_df`` or the last clustering."""
    if proportion_df is None:
        stored = adata.uns.get('sample_clusters', None)
        if stored is None:
            raise ValueError("Give proportion_df (the output of pl.pl.clustering_emd).")
        proportion_df = stored
    if not {'sampleID', 'Predicted_Labels'} <= set(proportion_df.columns):
        raise KeyError("proportion_df must contain 'sampleID' and 'Predicted_Labels'.")
    return pd.Series(proportion_df['Predicted_Labels'].astype(str).values,
                     index=proportion_df['sampleID'].astype(str).values)


def correlation_with_clustering(adata,
                                proportion_df=None,
                                sample_col=None,
                                features='auto',
                                include_cell_level=True,
                                agg='mean',
                                max_missing=0.5,
                                min_group_size=2,
                                max_levels=None,
                                numeric_min_unique=6,
                                exclude=None,
                                max_mixed=0.5,
                                p_adjust='fdr_bh',
                                show_skipped=True):
    """
    Test all usable sample variables against the clusters of samples, in one table.

    Like ``pl.tl.correlation_with_trajectory``, but for the sample clusters of
    ``pl.pl.clustering_emd``:

    - categorical variables: chi-squared test of the table clusters x groups;
    - numeric variables: ANOVA (and Kruskal-Wallis) between the clusters.

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    proportion_df : DataFrame, optional
        Output of ``pl.pl.clustering_emd`` (columns 'sampleID' and
        'Predicted_Labels'). If None, the last clustering is used.
    sample_col, features, include_cell_level, agg, max_missing, min_group_size, max_levels, numeric_min_unique, exclude, max_mixed, p_adjust, show_skipped
        See ``pl.tl.correlation_with_trajectory``.

    Returns
    -------
    DataFrame
        One row per tested variable, sorted by p-value: 'Feature', 'Type',
        'Level', 'Test', 'Statistic', 'P_value', 'Adjusted_P_value',
        'N_samples', 'N_groups' (groups of the variable), 'N_clusters',
        'Note'. Also stored in ``adata.uns['cluster_associations']``.
    """
    clusters = _cluster_labels(adata, proportion_df)
    table, info = sample_variables(adata, sample_col=sample_col, features=features,
                                   include_cell_level=include_cell_level, agg=agg,
                                   max_missing=max_missing, min_group_size=min_group_size,
                                   max_levels=max_levels, numeric_min_unique=numeric_min_unique,
                                   exclude=exclude, samples=list(clusters.index),
                                   max_mixed=max_mixed)
    rows = []
    for _, meta in info[info['Type'] != 'skipped'].iterrows():
        feature = meta['Feature']
        values = table[feature]
        ok = values.notna()
        x, c = values[ok], clusters[ok]
        row = dict(Feature=feature, Type=meta['Type'], Level=meta['Level'], Test='',
                   Statistic=np.nan, P_value=np.nan, N_samples=int(ok.sum()),
                   N_groups=int(x.nunique()), N_clusters=int(c.nunique()), Note=meta['Reason'])
        if meta['Type'] == 'numeric':
            groups = [x[c == g].values.astype(float) for g in pd.unique(c)]
            groups = [g for g in groups if len(g) >= 2]
            row['Test'] = 'ANOVA'
            if len(groups) >= 2:
                stat, p = f_oneway(*groups)
                row.update(Statistic=float(stat), P_value=float(p))
                try:
                    row['Kruskal_P_value'] = float(kruskal(*groups).pvalue)
                except ValueError:
                    row['Kruskal_P_value'] = np.nan
        else:
            row['Test'] = 'Chi-squared'
            ct = pd.crosstab(c, x.astype(str))
            ct = ct.loc[ct.sum(axis=1) > 0, ct.sum(axis=0) > 0]
            if ct.shape[0] >= 2 and ct.shape[1] >= 2:
                stat, p, _, _ = chi2_contingency(ct)
                row.update(Statistic=float(stat), P_value=float(p))
        rows.append(row)
    columns = ['Feature', 'Type', 'Level', 'Test', 'Statistic', 'P_value', 'Adjusted_P_value',
               'Kruskal_P_value', 'N_samples', 'N_groups', 'N_clusters', 'Note']
    res = pd.DataFrame(rows)
    if res.empty:
        res = pd.DataFrame(columns=columns)
    else:
        res['Adjusted_P_value'] = _adjust(res['P_value'], p_adjust)
        res = res.reindex(columns=columns)
        res = res.sort_values('P_value', na_position='last', kind='stable').reset_index(drop=True)
    skipped = info[info['Type'] == 'skipped'].reset_index(drop=True)
    adata.uns['cluster_associations'] = res
    adata.uns['cluster_associations_skipped'] = skipped
    adata.uns['cluster_associations_values'] = table
    adata.uns['sample_clusters'] = pd.DataFrame({'sampleID': clusters.index,
                                                 'Predicted_Labels': clusters.values})
    if show_skipped:
        _print_summary(res, skipped)
    return res


def name_clusters(adata, proportion_df, by=None, sample_col=None, min_fraction=0.0,
                  mixed_name='mixed'):
    """
    Give the sample clusters names from the majority group of their samples.

    Leiden numbers the clusters 0, 1, 2, ... This function names each cluster
    after the most frequent group (by default the status) of its samples,
    e.g. 'control' or 'disease'. If several clusters have the same majority
    group they are numbered by size: 'control_1' (largest), 'control_2', ...

    Parameters
    ----------
    adata : AnnData
        AnnData object after ``pl.tl.wasserstein_distance``.
    proportion_df : DataFrame
        Output of ``pl.pl.clustering_emd`` (columns 'sampleID' and
        'Predicted_Labels').
    by : str, optional
        Column of ``adata.obs`` used for the names. If None (default), the
        status of ``pl.tl.wasserstein_distance``.
    sample_col : str, optional
        Column with the samples; default: the one of ``wasserstein_distance``.
    min_fraction : float, optional
        If the majority group has a smaller fraction of the samples of a
        cluster, the cluster is called ``mixed_name``, by default 0 (always
        use the majority).
    mixed_name : str, optional
        Name of clusters without a clear majority, by default 'mixed'.

    Returns
    -------
    DataFrame
        Copy of ``proportion_df`` where 'Predicted_Labels' holds the new names
        and the old numbers are kept in 'Leiden'. The mapping is stored in
        ``adata.uns['cluster_names']`` and the table in
        ``adata.uns['sample_clusters']``.
    """
    labels = _cluster_labels(adata, proportion_df)
    sample_col = _default_sample_col(adata, sample_col)
    if by is None:
        by = _params(adata).get('status', None)
    if by is None or by not in adata.obs.columns:
        raise KeyError("Please give 'by', a column of adata.obs (e.g. the status).")
    groups = (adata.obs[[sample_col, by]].astype({sample_col: str})
              .drop_duplicates(subset=sample_col).set_index(sample_col)[by].astype(str))
    group = groups.reindex(labels.index)
    summary = []
    for cluster in pd.unique(labels):
        members = group[labels == cluster]
        counts = members.value_counts()
        top = counts.index[0] if len(counts) else mixed_name
        frac = counts.iloc[0] / counts.sum() if len(counts) else 0
        if len(counts) > 1 and counts.iloc[0] == counts.iloc[1]:
            top = '/'.join(sorted(counts.index[counts == counts.iloc[0]].astype(str)))
        if frac < min_fraction:
            top = mixed_name
        summary.append((cluster, top, int(len(members)), frac))
    summary = pd.DataFrame(summary, columns=['Leiden', 'base', 'size', 'fraction'])
    mapping = {}
    for base, sub in summary.groupby('base', sort=False):
        sub = sub.sort_values(['size', 'fraction'], ascending=False, kind='stable')
        if len(sub) == 1:
            mapping[sub['Leiden'].iloc[0]] = base
        else:
            for i, cluster in enumerate(sub['Leiden'], start=1):
                mapping[cluster] = f"{base}_{i}"
    out = proportion_df.copy()
    out['Leiden'] = out['Predicted_Labels'].astype(str).values
    out['Predicted_Labels'] = out['Leiden'].map(mapping).values
    summary['name'] = summary['Leiden'].map(mapping)
    adata.uns['cluster_names'] = {str(k): v for k, v in mapping.items()}
    adata.uns['sample_clusters'] = pd.DataFrame({'sampleID': out['sampleID'].astype(str).values,
                                                 'Predicted_Labels': out['Predicted_Labels'].values})
    print("Cluster names (majority of '" + str(by) + "'):")
    for _, r in summary.sort_values('name').iterrows():
        print(f"  {r['Leiden']} -> {r['name']}  ({r['size']} samples, "
              f"{100 * r['fraction']:.0f}% {r['base']})")
    return out
