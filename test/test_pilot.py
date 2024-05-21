import pilotpy as pl
import scanpy as sc

def test_pilot():

    adata_T=sc.read_h5ad('../Tutorial/Datasets/Kidney_IgAN_T.h5ad') 

    pl.tl.wasserstein_distance(
        adata_T,
        clusters_col = 'Cell_type',
        sample_col = 'sampleID',
        status = 'status', 
        data_type = 'Pathomics'
        )


    assert adata_T.uns['EMD'].shape[0]==adata_T.uns['EMD'].shape[1]

    pl.pl.heatmaps(adata_T)
    pl.pl.trajectory(adata_T, colors = ['red','blue','orange'])

    assert len(adata_T.uns['real_labels'])==adata_T.uns['EMD'].shape[1]
