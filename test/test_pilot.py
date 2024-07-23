import pilotpy as pl
import scanpy as sc

def test_pilot():

    adata_T=sc.read_h5ad('Tutorial/Datasets/Kidney_IgAN_T.h5ad') 
    adata_G=sc.read_h5ad('Tutorial/Datasets/Kidney_IgAN_G.h5ad')
    
    pl.tl.wasserstein_distance(
    adata_G,
    clusters_col = 'Cell_type',
    sample_col = 'sampleID',
    status = 'status',
    data_type = 'Pathomics'
    )
    
    pl.tl.wasserstein_distance(
        adata_T,
        clusters_col = 'Cell_type',
        sample_col = 'sampleID',
        status = 'status', 
        data_type = 'Pathomics'
        )


    assert adata_T.uns['EMD'].shape[0]==adata_T.uns['EMD'].shape[1]
    assert adata_G.uns['EMD'].shape[0]==adata_G.uns['EMD'].shape[1]
    assert len(adata_T.uns['real_labels'])==adata_T.uns['EMD'].shape[1]
    adata_Com = adata_G
    adata_Com.uns['EMD'] = adata_G.uns['EMD'] + adata_T.uns['EMD']
    pl.pl.fit_pricipla_graph(adata_Com, source_node = 2)
    pl.pl.heatmaps(adata_T)
    pl.pl.trajectory(adata_T, colors = ['red','blue','orange'])
    pl.tl.cell_importance(adata_Com, xlim = 125)
    pl.tl.extract_cells_from_pathomics(adata_Com)
    data=pl.tl.norm_morphological_features(
     column_names = ['glom_sizes', 'glom_distance_to_closest_glom','glom_diameters','glom_tuft_sizes',
                     'glom_bowman_sizes'],
     name_cell = 'All'
        )
    pl.tl.morphological_features_importance(data, height = 10, x_lim = 160, width = 20) 
    assert len(adata_T.uns['orders']) == adata_Com.uns['orders'].shape[1]
    
