#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mina shaigan
"""

#%%% load libraries
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

#%%% load functions
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

def infer_gene_cluster_differentiation(gene_list: list = None,
                                       cluster_names: list = None,
                                       col: str = 'Time_score',
                                       fc_thr: float = 1.5,
                                       eigenThresh: float = 1e-8,
                                       n_points: int = 20,
                                       start: int = 1,
                                       end: int = 20,
                                       path_to_results: str = None,
                                       file_name: str = "/Whole_expressions.csv",font_size=12):
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
    plot_stats_by_pattern(cluster_names, all_stats_extend, gene_dict, pline, path_to_results, file_name,font_size=font_size)
    
    
def plot_stats_by_pattern(cluster_names: list = None,
                          all_stats_extend: pd.DataFrame = None,
                          gene_dict: dict = None,
                          pline: list = None,
                          path_to_results: str = None,
                          file_name: str = "/Whole_expressions.csv",
                          font_size: int = 16):
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

    Returns
    -------
    return figures for each cluster ploting top 4 genes for each pattern.

    """
    plt.rcParams.update({'font.size': font_size})
    # plot results
    for cluster in cluster_names:
        my_data = all_stats_extend[ (all_stats_extend['cluster'] == str(cluster)) & (all_stats_extend['pvalue'] < 0.01)]
        sort_my_data = my_data.sort_values(['Expression pattern', 'FC', 'fit-pvalue'],
                ascending=[True, False, True]).groupby('Expression pattern').head(4)
        expression_patterns = np.unique(sort_my_data['Expression pattern'])
        if(len(expression_patterns)):

            n_col = 4
            n_row = len(expression_patterns)
            n_px = 10

            plt.rcParams["figure.facecolor"] = 'w'
            plt.figure()
            f, axs = plt.subplots(n_row, n_col, figsize=(n_col * n_px * 2, n_row * n_px))
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
                        axs[p, k].set_ylabel(pattern, size = font_size)
                    else:
                        axs[p, k].set_ylabel('Gene expression', size = font_size)
                    for item in (axs[p, k].get_xticklabels() + axs[p, k].get_yticklabels()):
                        item.set_fontsize(font_size)

                    plt.text(.01, .99, 'p-value = %.2e ' % + Decimal(str(row['fit-pvalue'])),
                             ha = 'left', va = 'top',
                             transform=axs[p, k].transAxes, size = font_size)
                    axs[p, k].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    k += 1
                while(k != 4):
                    axs[p, k].set_axis_off()
                    k += 1
                p += 1
            plt.subplots_adjust(wspace = 0.5, hspace = 0.7)
            if not os.path.exists(path_to_results+'/plots_gene_cluster_differentiation/'):  
                    os.makedirs(path_to_results+'/plots_gene_cluster_differentiation/')
            save_path = path_to_results+'/plots_gene_cluster_differentiation/'+str(cluster) + ".png"
            plt.savefig(save_path)
            plt.show()
