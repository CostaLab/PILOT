import os
from genericpath import isfile
from logging import info, warn
import pandas as pd
import numpy as np
from typing import Tuple, Dict
from collections import defaultdict
from sklearn.covariance import empirical_covariance
from numpy.linalg import eig
from scipy import stats
from sklearn.linear_model import HuberRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
import statsmodels
from statsmodels.stats.multitest import multipletests
import statsmodels.stats as smf
from os import listdir
from os.path import isfile, join
import math
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.lines import Line2D
from decimal import Decimal
from matplotlib.ticker import FormatStrFormatter
import pickle


def get_gene_dict(gene_list: list = None,
              cluster_names:list = None,
              path_to_results: str = None,
              file_name: str = "/Whole_expressions.csv"):
    """
    Map genes to clusters it belongs
    

    Parameters
    ----------
    gene_list : list
        List of genes to make the dictionary
    cluster_names: list
        List of cell types to consider for checking their differentiation
    path_to_results : str
        Path where 'Markers' folder had been created
    file_name: str
        File name same in all cell types tables created

    Returns
    -------
    Dictionary of genes map to clusters

    """
    gene_dict = defaultdict(list)
    for gene in gene_list:
        for c in cluster_names:
            data = pd.read_csv(path_to_results + '/Markers/' + str(c) + \
                               file_name, index_col = 0)
            if gene in data['Gene ID'].tolist():
                gene_dict[gene].append(c)
                
    with open(path_to_results + '/genes_dictionary.pkl', 'wb') as f:
        pickle.dump(gene_dict, f)
    return gene_dict


def load_data(path: str = None,
              name: str = None,
              col_names: list = None):
    """
    load matrix of samples(cells) x genes
    

    Parameters
    ----------
    path : str
        Path to get the data
    name: str
        Name of the file
    col_names: list
        get specific gene along the time_scores 

    Returns
    -------
    Desired data read as Dataframe

    """
    for format in [".csv"]:
        p = os.path.join(path, name + format)
        if isfile(p):
            info("loading " + p)
            try: 
                if format == ".csv":
                    if(col_names is not None):
                        target = pd.read_csv(p, usecols = col_names)
                    else:
                        target = pd.read_csv(p)
                return target
            except (Exception):
                warn("loading " + format + " failed")

def get_betas_by_table(table: pd.DataFrame = None,
                       gene_name: str = None,
                       feature_name: str = 'Gene ID'):
    """
    get coefficients of specific gene from table
    

    Parameters
    ----------
    table: DataFrame
        containing fitted model of each gene
    gene_name: str
        gene name to get coefficients
    feature_name: str
        column name to containing gene names

    Returns
    -------
    Dictionary of coefficients which can be at most three

    """
    table = table[table[feature_name] == gene_name]
    if(table['Fitted function'].values[0] == 'linear'):
        gene_beta = {'Intercept': float(table['Intercept'].values[0]),
                     'Feature': float(table['Treat'].values[0]),
                     'Feature2': float(0) }
    elif(table['Fitted function'].values[0] == 'quadratic'):
        gene_beta = {'Intercept': float(table['Intercept'].values[0]),
                     'Feature': float(0),
                     'Feature2': float(table['Treat'].values[0]) }
    else:
        gene_beta = {'Intercept': float(table['Intercept'].values[0]),
                     'Feature': float(table['Treat'].values[0]),
                     'Feature2': float(table['Treat2'].values[0]) }
    return gene_beta

def get_betas_by_model(md: object = None,
                       func_type: str = None):
    """
    get coefficients of specific gene from fitted model with specific function
    

    Parameters
    ----------
    md: object
        Fitted LinearRegression / HuberRegressor estimator
    func_type: str
        Used function to fit can be linear, quadratic, linear_quadratic

    Returns
    -------
    Dictionary of coefficients which can be at most three

    """
    if(func_type == 'linear'):
        gene_beta = {'Intercept': float(md['params'][0]),
                     'Feature': float(md['params'][1]),
                     'Feature2': float(0) }
    elif(func_type == 'quadratic'):
        gene_beta = {'Intercept': float(md['params'][0]),
                     'Feature': float(0),
                     'Feature2': float(md['params'][1]) }
    else:
        gene_beta = {'Intercept': float(md['params'][0]),
                     'Feature': float(md['params'][1]),
                     'Feature2': float(md['params'][2]) }
    return gene_beta

def make_bootstraps(data: np.array = None,
                    n_bootstraps: int = 100):
    """
    Function to generate bootstrapped samples
    
    Parameters
    ----------
    data: np.array
        Array of input data
    n_bootstraps: int
        Integer number of bootstraps to produce

    Returns
    -------
    {'boot_n': {'boot': np.array, 'test': np.array}}:
    dictionary of dictionaries containing 
        the bootstrap samples & out-of-bag test sets
    """

    # initialize output dictionary, unique value count, sample size, & list of indices
    dc       = {}
    n_unival = 0
    sample_size = data.shape[0]
    idx = [i for i in range(sample_size)]
    # loop through the required number of bootstraps
    for b in range(n_bootstraps):
        # obtain boostrap samples with replacement
        sidx   = np.random.choice(idx, replace = True, size = sample_size)
        b_samp = data[sidx,:]
        # compute number of unique values contained in the bootstrap sample
        n_unival += len(set(sidx))
        # obtain out-of-bag samples for the current b
        oob_idx = list(set(idx) - set(sidx))
        t_samp = np.array([])
        if oob_idx:
            t_samp = data[oob_idx,:]
        # store results
        dc['boot_' + str(b)] = {'boot': b_samp, 'test': t_samp}
    # state the mean number of unique values in the bootstraps
    # print('Mean number of unique values in each bootstrap: {:.2f}'.format(n_unival / n_bootstraps))
    # return the bootstrap results
    return(dc)

def generate_feature_by_fill_empty(func_type: str = None,
                                   X: list = None):
    """
    create features list of the model with size two and fill non exist feature
    with zero in the same position
    
    Parameters
    ----------
    func_type: str 
        Function type of fitted model
    X: list
        List of time scores
    
    Returns
    -------
    numpy array of size num_points x 2 for any function
    """
    assert func_type in ['linear', 'quadratic', 'linear_quadratic'], \
        'func type parameter must be linear, quadratic, linear_quadratic'
    
    if func_type == 'linear':
        return np.column_stack((X, np.zeros(X.shape[0])))
        
    if func_type == 'linear_quadratic':
        return np.column_stack((X, np.power(X, 2)))
    
    if func_type == 'quadratic':
        return np.column_stack((np.zeros(X.shape[0]), np.power(X, 2)))
    
def make_linearCombination(table: pd.DataFrame = None,
                           feature_name: str = None,
                           feature_type: str = 'Gene ID',
                           n_points: int = 100,
                           start: int = 1,
                           end: int = 20):
    """
    Create linear combination of features of the model to be used for Wald test
    
    Parameters
    ----------
    table: DataFrame
        containing fitted model of each gene
    feature_name: str
        name of gene
    feature_type: str
        name of column containing gene names
    n_points: int
        number of points to consider in fitted curve
    start: int
        start point of curve
    end: int
        end point of curve

    Returns
    -------
    numpy array of size num_points x 3 for any function
    """
    i = 0
    filter_table = table.copy()
    filter_table = filter_table[filter_table[feature_type] == feature_name]
    
    best_func = filter_table.iloc[i, 3]
    pline = np.linspace(start, end, n_points)
    polyline = generate_feature_by_fill_empty(best_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    return new_polyline

def generate_feature(func_type: str = None,
                     X: list = None):
    """
    Create features list of the model with size two and fill non exist feature
    with zero in the same position
    
    Parameters
    ----------
    func_type: str 
        Function type of fitted model
    X: list
        List of time scores
    
    Returns
    -------
    numpy array of size num_points x number function parameters
    """
    assert func_type in ['linear', 'quadratic', 'linear_quadratic'], \
        'func type parameter must be linear, quadratic, linear_quadratic'
    
    if func_type == 'linear':
        return X.reshape(X.shape[0],1)
        
    if func_type == 'linear_quadratic':
        return np.column_stack((X, np.power(X,2) ))
    
    if func_type == 'quadratic':
        double_X = np.power(X,2)
        return double_X.reshape(double_X.shape[0],1)

def _fit_model_(func_type: str = None,
                       X: np.array = None,
                       y: np.array = None,
                       model_type: str = None):
    """
    Fit linearRegression / HuberRegressor model on the given data
    
    Parameters
    ----------
    func_type: str 
        Function type of fitted model
    X: np.array
        Training vector, n_samples x one gene
    y: np.array
        Target vector relative to X
    model_type: str
        Whether you want to user LienarRegression or HuberRegressor
    
    Returns
    -------
    Dictionary of fitted model parameters and attributes
    """
    X = generate_feature(func_type, X)
    
    results = dict.fromkeys(['pvalues', 'rsquared_adj', 'params'], [])
    
    
    if(model_type == 'LinearRegression'):
        model = LinearRegression().fit(X, y)
    elif(model_type == 'HuberRegressor'):
        epsilon = 1.35
        for step in list(np.arange(0, 1, 0.05)):
            try:
                model = HuberRegressor(epsilon = epsilon + step).fit(X, y)
            except ValueError:
                continue
            break
    params = np.append(model.intercept_,model.coef_)
    predictions = model.predict(X)
    results['rsquared_adj'] = 1 - (1 - model.score(X, y)) * (len(y) - 1) / \
                                    (len(y) - X.shape[1] - 1)
    results['rmse'] = mean_squared_error(y, predictions, squared=False)

    new_X = np.append(np.ones((len(X),1)), X, axis=1)
    M_S_E = (sum((y-predictions)**2)) / (len(new_X)-len(new_X[0]))
    v_b = M_S_E * (np.linalg.inv(np.dot(new_X.T,new_X)).diagonal())
    s_b = np.sqrt(v_b)
    t_b = params/ s_b
    p_val =[2 * (1 - stats.t.cdf(np.abs(i),
                                 (len(new_X)-len(new_X[0])))) for i in t_b]

    p_val_adj = smf.multitest.multipletests(p_val, method='fdr_bh')[1]

    results['pvalues'] = p_val_adj
    results['params'] = params
    
    return results

def _fit_best_model_(target: pd.DataFrame = None,
                   data: pd.DataFrame = None,
                   pval_thr: float = 0.05,
                   model_type: str = 'HuberRegressor',
                   fun_types: list = ['linear', 'linear_quadratic', 'quadratic']):
    """
    Find best fitted model for each gene
    
    Parameters
    ----------
    target: DataFrame
        matrix of samples(cells) x genes
    data: DataFrame
        matrix with column name label to infer       
    pval_thr: float
        p-value threshold for checking significancy of p-value
    model_type: str
        model type to fit data (LinearRegression, or HuberRegressor)        
    fun_types: list
            list of function to use (['linear', 'linear_quadratic', 'quadratic'])
    
    Returns
    -------
    Dataframe of genes with their information of fitted model
    """
    
    assert model_type in ['LinearRegression', 'HuberRegressor'], \
        'model type must be LinearRegression, or HuberRegressor'
    
    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    
    all_results = []
    best_r2 = {}
    for index, tf in target.iterrows():
        max_r2 = -1000
        best_results = 0
        best_func = 0
        
        try:
            for func_type in fun_types:
                
                results = _fit_model_(func_type, x, tf, model_type)
            
                if(all(i <= pval_thr for i in list(results['pvalues']))):
                    
                    r2 = results['rsquared_adj']
                    if(r2 > max_r2):
                        best_func = func_type
                        best_results = results
                        max_r2 = r2                  
          
            if(best_func):
            
                pline = np.linspace(min_x, max_x)
                polyline = generate_feature(best_func, pline)
                new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
                curve = np.matmul(new_polyline, best_results['params'])
                slope = (curve[-1] - curve[0]) 
                
                rho, pval = stats.pearsonr(x, tf)
            
                if best_func == 'linear':
                    if best_results['params'][1] >= 0:
                        exp_pattern = 'linear up'
                    else:
                        exp_pattern = 'linear down'
                
                elif best_func == 'quadratic':
                    if best_results['params'][1] >= 0:
                            exp_pattern = 'quadratic up'
                    else:
                        exp_pattern = 'quadratic down'
        
                elif best_func == 'linear_quadratic':
                    if best_results['params'][1] >= 0 and best_results['params'][2] >= 0:
                        exp_pattern = 'linear up quadratic up'
                    elif best_results['params'][1] >= 0 and best_results['params'][2] < 0:
                        exp_pattern = 'linear up quadratic down'
                    elif best_results['params'][1] < 0 and best_results['params'][2] >= 0:
                        exp_pattern = 'linear down quadratic up'
                    else:
                        exp_pattern = 'linear down quadratic down'
            
            
                all_results.append(best_results)
                best_r2.update({index: [best_func, best_results, exp_pattern, slope, pval]})
                
        except BaseException:
            
            print(index)
        
    sorted_best = {k: v for k, v in sorted(best_r2.items(),
                                           key = lambda item: item[1][1]['rsquared_adj'], reverse=True)}
    try: 
        if(bool(sorted_best)):
            table = make_table(sorted_best)
            return table
    except (Exception):
        warn("No data available!")

def make_table(dictionary: dict = None):
    """
    Convert dictionary of the fitted model to a table
    
    Parameters
    ----------
    dictionary: dict
        Dictionary of fitted model parameters and attributes for each gene
    
    Returns
    -------
    Dataframe of genes with their information of fitted model
    """
    try: 
        if(bool(dictionary)):
            all_motifs = np.array(list(dictionary.keys()))
            all_motifs = all_motifs.reshape(all_motifs.shape[0],1)
            
            all_patterns = np.array([i[2] for i in list(dictionary.values())])
            all_patterns = all_patterns.reshape(all_patterns.shape[0],1)
            
            
            all_slopes = np.array([i[3] for i in list(dictionary.values())])
            all_slopes = all_slopes.reshape(all_slopes.shape[0],1)
            
            mat = np.append(all_motifs, all_patterns, axis = 1)
            mat = np.append(mat, all_slopes, axis = 1)
        
            tmp = [i[1] for i in list(dictionary.values())]
            
            all_func_types = np.array([i[0] for i in list(dictionary.values())])
            all_func_types = all_func_types.reshape(len(dictionary.keys()),1)
            mat = np.append(mat, all_func_types, axis = 1)
        
            my_stats = pd.Series([i['params'] for i in [j[1] for j in list(dictionary.values())] ])
            df = pd.DataFrame(my_stats.values.tolist(), index = my_stats.index)
        
            # intercept1
            my_stats = np.array(df.iloc[:,0])
            my_stats = my_stats.reshape(my_stats.shape[0],1)
            mat = np.append(mat, my_stats, axis = 1)
        
            # feature2
            my_stats = np.array(df.iloc[:,1])
            my_stats = my_stats.reshape(my_stats.shape[0],1)
            mat = np.append(mat, my_stats, axis = 1)
        
            # feature3
            if(df.shape[1] >= 3):
                my_stats = np.array(df.iloc[:,2])
                my_stats = my_stats.reshape(my_stats.shape[0],1)
                mat = np.append(mat, my_stats, axis = 1)
                
                column_names = ['Feature ID', 'Pattern', 'Slope', 'Fitted function',
                                'Intercept', 'Treat', 'Treat2', 'adjusted P-value',
                                'adjusted R-squared', 'RMSE']
            else:
                column_names = ['Feature ID', 'Pattern', 'Slope', 'Fitted function',
                                'Intercept', 'Treat', 'adjusted P-value',
                                'adjusted R-squared', 'RMSE']
        
            # compute adjusted p-values
            all_pvals = np.array([i[4] for i in list(dictionary.values())])
            all_pvals = multipletests(list(all_pvals), method='fdr_bh')[1]
            all_pvals = all_pvals.reshape(all_pvals.shape[0],1)
            mat = np.append(mat, all_pvals, axis = 1)
        
            my_stats = np.array([i['rsquared_adj'] for i in tmp])
            my_stats = my_stats.reshape(my_stats.shape[0],1)
            mat = np.append(mat, my_stats, axis = 1)
            
            my_stats = np.array([i['rmse'] for i in tmp])
            my_stats = my_stats.reshape(my_stats.shape[0],1)
            mat = np.append(mat, my_stats, axis = 1)
            
            mat_df = pd.DataFrame(mat , columns = column_names)
            return mat_df
    except (Exception):
        warn("No data available!")

def get_params_table(table: pd.DataFrame = None,
                     feature_name: str = None,
                     feature_type: str = 'Gene ID'):
    """
    get coefficients of specific gene from table
    

    Parameters
    ----------
    table: DataFrame
        containing fitted model of each gene
    feature_name: str
        gene name to get coefficients
    feature_type: str
        column name to containing gene names

    Returns
    -------
    List of coefficients which can be at most three

    """
    table = table[table[feature_type] == feature_name]
    func = table.iloc[0,3]
    if(table.shape[1] > 8):
        if( (func == 'linear') | (func == 'quadratic') ):
            results_params = [float(table.iloc[0,4]), 
                              float(table.iloc[0,5])]
        else:
            results_params = [float(table.iloc[0,4]), 
                              float(table.iloc[0,5]),
                              float(table.iloc[0,6])]    
    else:
        results_params = [float(table.iloc[0,4]), 
                          float(table.iloc[0,5])]
    return results_params

def extend_stats(all_stats: pd.DataFrame = None,
                 path_to_results: str = None,
                 col_names: list = ['adjusted P-value', 'R-squared', 'mod_rsquared_adj']):
    """
    Extend Wald stats table by including the fitted model attributes

    Parameters
    ----------
    all_stat : pd.DataFrame
        matrix with statistical information generated by Wald stats
        for each gene and for each cluster. 
    path_to_results : str, optional
        define path storing Markers, cells, and results.. 
    col_names : list, optional
        set the column names to add from fitted model.

    Returns
    -------
    Dataframe of statistical test of differentitiation and fitted model
    for each gene with specific clsuter.

    """
    all_stats_extend = pd.DataFrame(columns = ['gene', 'cluster', 'waldStat', 'df', 'pvalue', 'FC', 
                                           'Expression pattern', 'fit-pvalue', 'fit-rsquared', 'fit-mod-rsquared'])
    
    all_stats_extend = pd.concat([all_stats_extend,all_stats
       [['gene', 'cluster', 'waldStat', 'df', 'pvalue', 'FC']]])
   # all_stats_extend = all_stats_extend.append(all_stats[
    #    ['gene', 'cluster', 'waldStat', 'df', 'pvalue', 'FC']])
    for i, row in all_stats_extend.iterrows():
        data = pd.read_csv(path_to_results + "/Markers/"+ str(row['cluster']) + "/Whole_expressions.csv", index_col=0)
        all_stats_extend.at[i,'Expression pattern'] = data[data['Gene ID'] == row['gene']]['Expression pattern'].values[0]
        all_stats_extend.at[i,'fit-pvalue'] = data[data['Gene ID'] == row['gene']]['adjusted P-value'].values[0]
        all_stats_extend.at[i,'fit-rsquared'] = data[data['Gene ID'] == row['gene']]['R-squared'].values[0]
        all_stats_extend.at[i,'fit-mod-rsquared'] = data[data['Gene ID'] == row['gene']]['mod_rsquared_adj'].values[0]
        # for col in col_names:
        #     all_stats_extend.at[i, 'fit-'+str(col)] = data[data['Gene ID'] == row['gene']][str(col)].values[0]
       
    if  os.path.exists(path_to_results + "/gene_clusters_stats_extend.csv"):
        all_stats_extend.to_csv(path_to_results + "/gene_clusters_stats_extend_exploring_specific_genes.csv", index = None)
    else:
        all_stats_extend.to_csv(path_to_results + "/gene_clusters_stats_extend.csv", index = None)
    return all_stats_extend

def get_coordinates(WT_adata: pd.DataFrame = None,
                    KO_tf: pd.Series = None,
                    KO_x: pd.DataFrame = None,
                    WT_table: pd.DataFrame = None,
                    KO_table: pd.DataFrame = None,
                    feature_name: str = None, 
                    target_name: str = 'Time_score'):
    """
    Make the curve for two curve and get the data for plotting

    Parameters
    ----------
    WT_adata : pd.DataFrame, optional
        Store time points for the intrested gene with specific cluster
    KO_tf : pd.Series, optional
        Average distribution of expression for other clusters.
    KO_x : pd.DataFrame, optional
        Time ponts of other clusters.
    WT_table : pd.DataFrame, optional
        Table stored fitted model for the gene intrested.
    KO_table : pd.DataFrame, optional
        Table stored fitted model for other clusters.
    feature_name : str, optional
        Gene name of instrest.
    target_name : str, optional
        name of column store time points. The default is 'Time_score'.

    Returns
    -------
    WT_x : TYPE
        Time ponts of other clusters.
    WT_tf : TYPE
        Average distribution of expression for specific clusters.
    WT_curve : TYPE
        Fitted model for gene instrest with its cluster 
    KO_x : TYPE
        Time ponts of other clusters.
    KO_tf : TYPE
        Average distribution of expression for other clusters.
    KO_curve : TYPE
        Fitted model for other clusters.

    """
    WT_adata_copy = WT_adata.copy()
    WT_data = pd.DataFrame()
    WT_target = WT_adata_copy.transpose()
    
    if target_name in WT_adata_copy.columns:
        WT_data['label'] = list(WT_adata_copy[target_name].values)
    else:
        raise Exception("target_name are not in the adata obs columns")
    
    
    WT_x = WT_data['label']
    WT_tf = list(WT_target.loc[[feature_name]].values.ravel())
    
    WT_table = WT_table[WT_table['Gene ID'] == feature_name]
    WT_func = WT_table.iloc[0,3]
    if(WT_table.shape[1]>8):
        if( (WT_func == 'linear') | (WT_func == 'quadratic') ):
            WT_results_params = [float(WT_table.iloc[0,4]), 
                                   float(WT_table.iloc[0,5])]
        else:
            WT_results_params = [float(WT_table.iloc[0,4]), 
                                   float(WT_table.iloc[0,5]),
                                   float(WT_table.iloc[0,6])]
            
    else:
        WT_results_params = [float(WT_table.iloc[0,4]), 
                             float(WT_table.iloc[0,5])]
    
    if(KO_tf is not None):
        KO_table = KO_table[KO_table['Feature ID'] == feature_name]
        KO_func = KO_table.iloc[0,3]
        if(KO_table.shape[1]>8):
            if( (KO_func == 'linear') | (KO_func == 'quadratic') ):
                KO_results_params = [float(KO_table.iloc[0,4]), 
                                       float(KO_table.iloc[0,5])]
            else:
                KO_results_params = [float(KO_table.iloc[0,4]), 
                                       float(KO_table.iloc[0,5]),
                                       float(KO_table.iloc[0,6])]

        else:
            KO_results_params = [float(KO_table.iloc[0,4]), 
                                 float(KO_table.iloc[0,5])]

    pline = np.linspace(min(WT_x), max(WT_x))
    
    polyline = generate_feature(WT_func, pline)
    new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
    WT_curve = np.matmul(new_polyline, WT_results_params)
    
    if(KO_tf is not None):
        polyline = generate_feature(KO_func, pline)
        new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
        KO_curve = np.matmul(new_polyline, KO_results_params)
    else:
        KO_curve = None
    return (WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve)

def get_two_tables(gene_name: str = None,
                   cluster_n: str = None,
                   gene_dict: dict = None,
                   file_name: str = "/Whole_expressions.csv",
                   pline: list = None,
                   path_to_results: str = None):
    """
    

    Parameters
    ----------
    gene_name : str, optional
        name of the gene that is intersted.
    cluster_n : str, optional
        name of the cluster the gene belongs.
    gene_dict : dict, optional
        dictionary mapping the genes to the clusters its belongs.
    file_name : str, optional
        Name of the file the fitted data store.
        The default is "/Whole_expressions.csv".
    pline : list, optional
        default line to be used to plot curve. The default is None.
    path_to_results : str, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    WT_x : TYPE
        Time ponts of other clusters.
    WT_tf : TYPE
        Average distribution of expression for specific clusters.
    WT_curve : TYPE
        Fitted model for gene instrest with its cluster 
    KO_x : TYPE
        Time ponts of other clusters.
    KO_tf : TYPE
        Average distribution of expression for other clusters.
    KO_curve : TYPE
        Fitted model for other clusters.

    """
    
    WT_adata = load_data(path_to_results + '/cells/', str(cluster_n), [gene_name, 'Time_score'])
    WT_adata = WT_adata.sort_values('Time_score', ignore_index=True)
    gene_clusters = gene_dict[gene_name]
    
    RNA_target_clusters_mean = []
    RNA_data_list = []
    for gc in gene_clusters:
        if( gc != cluster_n):
            data = load_data(path_to_results + '/cells/', str(gc), [gene_name, 'Time_score', 'sampleID'])
            data = data.sort_values('Time_score', ignore_index=True)
            data_group = data.groupby('Time_score').sample(n = min(data.groupby('Time_score').count()['sampleID']))
            RNA_target_clusters_mean.extend(list(np.transpose(data_group).loc[gene_name]))
            RNA_data_list.extend(data_group['Time_score'].values)
    
    table1 = pd.read_csv(path_to_results + '/Markers/' + str(cluster_n) + file_name, index_col = 0)
    
    if(len(gene_clusters) == 1):
        WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve = get_coordinates(WT_adata, WT_table = table1, feature_name = gene_name)
    else:
        
        RNA_target_clusters = pd.DataFrame(columns=list(range(len(pline))))
        for gc in gene_clusters:
            if( gc != cluster_n):             
                table2 = pd.read_csv(path_to_results + '/Markers/' + str(gc) + file_name, index_col = 0)
                func_type = table2[table2['Gene ID'] == gene_name]['Fitted function'].values[0]
                polyline = generate_feature(func_type, pline)
                new_polyline = np.append(np.ones((len(polyline),1)), polyline, axis=1)
                params = [x for x in table2[table2['Gene ID'] == gene_name][['Intercept', 'Treat', 'Treat2']].values[0]
                          if (math.isnan(x) == False)]
                curve = np.matmul(new_polyline, params)
                RNA_target_clusters.loc[len(RNA_target_clusters)] = list(curve)
        RNA_data = pd.DataFrame()
        RNA_data['label'] = pline
        table2 = _fit_best_model_(pd.DataFrame({gene_name: RNA_target_clusters.mean()}).transpose(), RNA_data,
                                         pval_thr = 1)
        WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve = get_coordinates(WT_adata, RNA_target_clusters_mean,
                                                                       RNA_data_list, table1, table2, gene_name)
    return (WT_x, WT_tf, WT_curve, KO_x, KO_tf, KO_curve)


