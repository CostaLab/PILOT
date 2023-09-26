
from .Gene_cluster_specific_functions import *
from .Gene_cluster_specific_functions import _fit_best_model_,_fit_model_
from ..plot.ploting import *



def infer_gene_cluster_differentiation(gene_list: list = None,
                                       cluster_names: list = None,
                                       col: str = 'Time_score',
                                       fc_thr: float = 1.5,
                                       eigenThresh: float = 1e-8,
                                       n_points: int = 20,
                                       start: int = 1,
                                       end: int = 20,
                                       path_to_results: str = None,
                                       file_name: str = "/Whole_expressions.csv",font_size=24,cell_show_plot=None):
    """
    

    Parameters
    ----------
    gene_list : list, optional
        list of gene you want to check their differentiation.
        The default is None would consider all.
    cluster_names : list, optional
        list of clusters to be considered for defferentiation.
        The default is None would condsider all.
    col : str, optional
        name of column containing time points. The default is 'Time_score'.
    fc_thr : float, optional
        define threshold for fold change as how much defferentitation you 
        expect to have between the cluster and all others. The default is 1.5.
    eigenThresh : float, optional
        defien eigen threshold to how many eigen vector show be consider based 
        on eigenThresh x highest eigen value consider . The default is 1e-8.
    n_points : int, optional
        number of points to be considered on the curve (better to be more than 
        number of time points). The default is 20.
    start : int, optional
        define start time point. The default is 1.
    end : int, optional
        define end time point. The default is 20.
    path_to_results : str, optional
        define path storing Markers, cells, and results.
    file_name : str, optional
        name of file storing fitted model for each gene same for all clusters.
        The default is "/Whole_expressions.csv".
    font_size:float, font size
    Returns
    -------
    Plot significant genes for each pattern for each cell types.

    """
    l2fc = np.log2(fc_thr)
    logFCCutoff = np.log(np.power(2,l2fc))
    pline = np.linspace(start, end, n_points)
    n_bootstraps = 50
    
    if(cluster_names is None):
        cluster_names = [os.path.splitext(f)[0] for f in listdir(path_to_results + '/cells/') \
                         if isfile(join(path_to_results + '/cells/', f))]
    else:
        if(len(cluster_names) < 2):
            warn("Number of clusters should be at least two to be compared.")
    
    if(gene_list is None):
        gene_list = []
        for c in cluster_names:
            data = pd.read_csv(path_to_results + '/Markers/' + str(c) + file_name, index_col = 0)
            gene_list.extend(data['Gene ID'].tolist())
        gene_list = np.unique(gene_list)
        
    gene_dict = get_gene_dict(gene_list, cluster_names, path_to_results, file_name)

    all_stats = pd.DataFrame(columns = ['gene', 'cluster', 'waldStat', 'df', 'pvalue', 'FC'])
    for gene_name in gene_list:
       # print(gene_name)
        gene_clusters = gene_dict[gene_name]
        if(len(gene_clusters) > 1):
            for i in range(len(gene_clusters)):
            
                # get points of fitted model
                #   for intrested gene with its specific cluster
                table1 = pd.read_csv(path_to_results + '/Markers/' + str(gene_clusters[i]) + file_name, index_col = 0)
                polyline1 = make_linearCombination(table1, gene_name, n_points = n_points, start = start, end = end)
    
                # get curve to compute fold change
                polyline = generate_feature(table1[table1['Gene ID'] == gene_name].iloc[0, 3], pline)
                new_polyline = np.append(np.ones((len(polyline), 1)), polyline, axis = 1)
                curve1 = np.matmul(new_polyline, get_params_table(table1, gene_name))
    
                # get points of fitted model
                #   for rest of clusters the intrested gene belongs
                #   and get their average distribution
                RNA_target_clusters = pd.DataFrame(columns=list(range(len(pline))))
                for j in range(len(gene_clusters)):
                    if( j != i ):
                        table2 = pd.read_csv(path_to_results + '/Markers/' + str(gene_clusters[j]) + file_name, index_col = 0)
                        func_type = table2[table2['Gene ID'] == gene_name]['Fitted function'].values[0]
                        polyline = generate_feature(func_type, pline)
                        new_polyline = np.append(np.ones((len(polyline), 1)), 
                                                 polyline, axis=1)
                        params = [x for x in table2[table2['Gene ID'] == gene_name][['Intercept', 'Treat', 'Treat2']].values[0] if (math.isnan(x) == False)]
                        curve = np.matmul(new_polyline, params)
                        RNA_target_clusters.loc[len(RNA_target_clusters)] = list(curve)     
                RNA_data = pd.DataFrame()
                RNA_data['label'] = pline
                table2 = _fit_best_model_(pd.DataFrame({gene_name: RNA_target_clusters.mean()}).transpose(),
                                        RNA_data, pval_thr = 1)
                polyline2 = make_linearCombination(table2, gene_name, 'Feature ID', n_points, start, end)
                
                # get curve to compute fold change
                polyline = generate_feature(table2[table2['Feature ID'] == gene_name].iloc[0, 3], pline)
                new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
                curve2 = np.matmul(new_polyline, get_params_table(table2, gene_name, 'Feature ID'))
    
                # create betas with putting the coefficients of two models together
                gene_betas = pd.DataFrame(columns=['Intercept', 'Feature', 'Feature2'])
                gene_beta = get_betas_by_table(table1, gene_name)
           
                gene_beta = pd.DataFrame([gene_beta])
               # gene_betas = gene_betas.append(gene_beta, ignore_index=True)
                gene_betas = pd.concat([gene_betas, gene_beta],ignore_index=True)
                gene_beta = get_betas_by_table(table2, gene_name, 'Feature ID')
                gene_beta = pd.DataFrame([gene_beta])
               # gene_betas = gene_betas.append(gene_beta, ignore_index=True)
                gene_betas = pd.concat([gene_betas, gene_beta],ignore_index=True)
                
                # estimate variance-covariance matrix using bootstraping method
                data = load_data(path_to_results + '/cells/', str(gene_clusters[i]), [gene_name, 'Time_score'])
                data = data.sort_values('Time_score', ignore_index = True)
                RNA_target = np.transpose(data).loc[gene_name]
                RNA_data = np.array(data[col])
                func_type = table1[table1['Gene ID'] == gene_name]['Fitted function'].values[0]
                RNA_data_boots = make_bootstraps(RNA_data.reshape(-1,1), n_bootstraps)
                boot_betas1 = pd.DataFrame(columns=['Intercept', 'Feature', 'Feature2'])
                for key in RNA_data_boots.keys():
                    RNA_data_boot = RNA_data_boots[key]['boot'].ravel()
                    md = _fit_model_(func_type, RNA_data_boot, RNA_target, 'HuberRegressor')
                    betas = get_betas_by_model(md, func_type)
                    betas = pd.DataFrame([betas])
                   # boot_betas1 = boot_betas1.append(betas, ignore_index = True)
                    boot_betas1 = pd.concat([boot_betas1, betas],ignore_index=True)
                    
                RNA_target = RNA_target_clusters.mean()
                RNA_data = pline
                func_type = table2[table2['Feature ID'] == gene_name]['Fitted function'].values[0]
                RNA_data_boots = make_bootstraps(RNA_data.reshape(-1,1), n_bootstraps)
                boot_betas2 = pd.DataFrame(columns = ['Intercept', 'Feature', 'Feature2'])
                for key in RNA_data_boots.keys():
                    RNA_data_boot = RNA_data_boots[key]['boot'].ravel()
                    md = _fit_model_(func_type, RNA_data_boot, RNA_target, 'HuberRegressor')
                    betas = get_betas_by_model(md, func_type)
                    betas = pd.DataFrame([betas])
                   # boot_betas2 = boot_betas2.append(betas, ignore_index = True)
                    boot_betas2 = pd.concat([boot_betas2, betas],ignore_index=True)
                cov_features = pd.concat([boot_betas1, boot_betas2], axis = 1)
                gene_covs = cov_features.cov()
                
                # make linearCombination matrix to have pairwise comparison 
                #   between two time points of two different clusters
                CC = np.zeros((n_points, 2 * 3))
                CC = np.concatenate((polyline1, polyline2 * -1), axis = 1)
            
                estFC = np.matmul(CC, np.array(gene_betas.values.flatten()).reshape(2 * 3, 1))
                est = np.sign(estFC) * np.maximum(0, np.abs(estFC) - logFCCutoff)
            
                part_sigma = np.matmul(CC, gene_covs)
                sigma = np.matmul(part_sigma, CC.transpose())
            
                eSigma_values, eSigma_vectors = eig(sigma)
                r = np.sum(eSigma_values / np.max(eSigma_values) > eigenThresh)
                halfCovInv = np.matmul(eSigma_vectors[eSigma_values / np.max(eSigma_values) > eigenThresh].transpose(), 
                                       np.diag(1/ np.sqrt(eSigma_values[eSigma_values / np.max(eSigma_values) > eigenThresh])))
                halfStat = np.matmul(est.transpose(), halfCovInv)
            
                stat = np.matmul(halfStat, halfStat.transpose())
            
                waldStat_score = stat.ravel()[0].real
                waldStat_pval = stats.chi2.sf(waldStat_score, r)
                log_fold_change = np.log2(curve1.mean()) - np.log2(curve2.mean())
                all_stats.loc[len(all_stats)] = [str(gene_name),
                                                   str(gene_clusters[i]),
                                                   waldStat_score, r,
                                                   waldStat_pval,
                                                   log_fold_change]
              #  print([str(gene_name), str(gene_clusters[i]), waldStat_score, r,
               #        waldStat_pval, log_fold_change])
        else:
            if(len(gene_clusters) == 1):
                all_stats.loc[len(all_stats)] = [str(gene_name), str(gene_clusters[0]),
                                                   1, 1, 0.0, 0.0]
    # all_stats['pvalue'].fillna(1.0, inplace = True)
    # all_pvals = smf.multitest.multipletests(list(all_stats['pvalue'].values), method='fdr_bh')[1]
    # all_stats['pvalue'] = all_pvals
    all_stats_extend = extend_stats(all_stats, path_to_results)
   # plot_stats_by_pattern(cluster_names, all_stats_extend, gene_dict, pline, path_to_results, file_name,font_size=font_size)



