
import os
from genericpath import isfile
from logging import info, warn
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import re
from scipy import stats
from sklearn.impute import SimpleImputer
from plotnine import *
import joypy
import statsmodels.stats as smf
import statsmodels.api as sm
import plotly.express as px
from matplotlib import cm
from sklearn.linear_model import LinearRegression, HuberRegressor
from sklearn.metrics import r2_score
from gprofiler import GProfiler
import plotly.io as pio
import math
import matplotlib as mpl
import elpigraph
from matplotlib.font_manager import FontProperties
from ..plot.ploting import *

   
def compute_metrics(y_test, y_pred):
    """
    Calculate regression evaluation metrics.

    Parameters:
        y_test (array-like): Ground truth target values.
        y_pred (array-like): Predicted target values.

    Returns:
        tuple: A tuple containing the following regression metrics:
            - R-squared (R2) Score
            - Root Mean Squared Error (RMSE)
            - Mean Absolute Error (MAE)
            - Pearson Correlation Coefficient
    """
    r2 = r2_score(y_test, y_pred)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    mae = mean_absolute_error(y_test, y_pred)
    cor = stats.pearsonr(y_test, y_pred)[0]
    return (r2, rmse, mae, cor)


def loadTarget(path, name):
    """
    Load a target file from a specified path and name.

    Parameters:
        path (str): The directory path where the target file is located.
        name (str): The name of the target file (without the file extension).

    Returns:
        pd.DataFrame: A Pandas DataFrame containing the loaded target data.

    Note:
        This function attempts to load the target data from different formats (currently supports ".csv"). 
        It returns the loaded data as a Pandas DataFrame or None if loading fails.
    """

    
    for format in [".csv"]:
        p = os.path.join(path, name + format)
        if isfile(p):
            info("loading " + p)
            try: 
                if format == ".csv":
                    target = pd.read_csv(p, index_col=0)
                 
                    # target = target.fillna(0)
                    # target = target.loc[(target!=0).any(axis=1)]
                 
                return target
            except (Exception):
                warn("loading " + format + " failed")
                
def loadData(path, name):
    """
    Load data from a file in the specified format.

    Parameters:
        path (str): The directory path where the data file is located.
        name (str): The base name of the data file (without the file extension).

    Returns:
        pd.DataFrame or None: A Pandas DataFrame containing the loaded data if successful, or None if loading failed.
    """
    for format in [".csv"]:
        p = os.path.join(path, name + format)
        if isfile(p):
            info("loading " + p)
            try: 
                if format == ".csv":
                    data = pd.read_csv(p)
                    # data = pd.Series(np.array(inf['label']), index=inf['sample_id'])
                return data
            except (Exception):
                warn("loading " + format + " failed")
                
def convert_nan(target, convert_type):
    """
    Convert NaN (missing) values in a DataFrame to a specified format.

    Parameters:
        target (pd.DataFrame): The DataFrame containing NaN values to be converted.
        convert_type (str): The type of conversion to apply, either 'impute' or 'zero'.

    Returns:
        pd.DataFrame: A new DataFrame with NaN values converted based on the specified type.
    """
    assert convert_type in ['impute', 'zero'], 'func type parameter must be impute, zero'
    
    if convert_type == 'impute':
        imp = SimpleImputer(missing_values=np.nan, strategy='mean')
        idf = pd.DataFrame(imp.fit_transform(target))
        idf.columns = target.columns
        idf.index = target.index
    elif convert_type == 'zero':
        idf = target.fillna(0)
        idf = idf.loc[(idf!=0).any(axis=1)]
    return idf
      
def calculate_sparsity(data):
    """
    Calculate the sparsity of a given dataset.

    Parameters:
        data (numpy.ndarray or scipy.sparse matrix): The input data for which to calculate sparsity.

    Returns:
        float: The sparsity of the input data as a decimal value between 0 and 1, where 0 indicates no sparsity (dense) and 1 indicates full sparsity (completely empty).
    """
    non_zero = np.count_nonzero(data)
    total_val = np.product(data.shape)
    sparsity = (total_val - non_zero) / total_val
    return(sparsity)

def generate_feature_cell_types(func_type, X):
    """
    Generate a feature matrix based on the specified functional type.

    Parameters:
        func_type (str): The type of feature construction to perform. 
            - 'linear': Linear feature construction (constant, X).
            - 'quadratic': Quadratic feature construction (constant, X^2).
            - 'linear_quadratic': Linear and quadratic feature construction (constant, X, X^2).
        X (numpy.ndarray): The input data or feature vector.

    Returns:
        numpy.ndarray: The feature matrix constructed based on the specified functional type.
    """
    assert func_type in ['linear', 'quadratic', 'linear_quadratic'], 'func type parameter must be linear, quadratic, cubic, sqrt, reciprocal'
    
    if func_type == 'linear':
        X = np.column_stack((np.ones(len(X)), X))
        
    if func_type == 'linear_quadratic':
        X = np.column_stack((np.ones(len(X)), X, np.power(X,2) ))
    
    if func_type == 'quadratic':
        X = np.column_stack((np.ones(len(X)), np.power(X,2) ))
        
    return X

                
def generate_feature(func_type, X):
    """
    Generate a feature matrix based on the specified functional type.

    Parameters:
        func_type (str): The type of feature construction to perform. 
            - 'linear': Linear feature construction (X).
            - 'quadratic': Quadratic feature construction (X^2).
            - 'linear_quadratic': Linear and quadratic feature construction (X, X^2).
        X (numpy.ndarray): The input data or feature vector.

    Returns:
        numpy.ndarray: The feature matrix constructed based on the specified functional type.
    """
    assert func_type in ['linear', 'quadratic', 'linear_quadratic'], 'func type parameter must be linear, quadratic, cubic, sqrt, reciprocal'
    
    if func_type == 'linear':
        return X.reshape(X.shape[0],1)
        
    if func_type == 'linear_quadratic':
        return np.column_stack((X, np.power(X,2) ))
    
    if func_type == 'quadratic':
        double_X = np.power(X,2)
        return double_X.reshape(double_X.shape[0],1)
            
def fit_model_activity(func_type, X, y,max_iter_huber,epsilon_huber,model_type):
    """
    Fit a linear regression or Huber regression model and calculate various model statistics.

    Parameters:
        func_type (str): The type of feature construction to perform. 
            - 'linear': Linear feature construction (X).
            - 'quadratic': Quadratic feature construction (X^2).
            - 'linear_quadratic': Linear and quadratic feature construction (X, X^2).
        X (numpy.ndarray): The input data or feature vector.
        y (numpy.ndarray): The target variable.
        max_iter_huber (int): The maximum number of iterations for HuberRegressor.
        epsilon_huber (float): The value of epsilon for HuberRegressor.
        model_type (str): The type of regression model to fit.
            - 'LinearRegression': Ordinary Least Squares (OLS) linear regression.
            - 'HuberRegressor': Huber regression with specified epsilon.

    Returns:
        dict: A dictionary containing various model statistics.
            - 'pvalues' (list): The p-values for each feature.
            - 'rsquared_adj' (float): The adjusted R-squared value.
            - 'params' (numpy.ndarray): The model coefficients (intercept and coefficients).
            - 'mod_rsquared_adj' (float): The modified adjusted R-squared value (Huber).
    """
    X = generate_feature(func_type, X)
    
    results = dict.fromkeys(['pvalues', 'rsquared_adj', 'params'], [])
    
    if(model_type == 'LinearRegression'):
          model = LinearRegression().fit(X, y)
    elif(model_type == 'HuberRegressor'):
        epsilon = epsilon_huber
        for step in list(np.arange(0, 1, 0.05)):
            try:
                model = HuberRegressor(epsilon = epsilon + step).fit(X, y)
            except ValueError:
                continue
            break
    
    params = np.append(model.intercept_,model.coef_)
    predictions = model.predict(X)
    results['rsquared_adj'] = 1 - (1-model.score(X, y))*(len(y)-1)/(len(y)-X.shape[1]-1)

    ss_res = np.sum( (y - predictions )**2)
    delata = 1.35
    modified_ss_res = 0
    res_e = y - predictions
    for e in res_e:
        if abs(e) < delata:
            modified_ss_res += 1/2*np.sum(e**2)
        else:
            modified_ss_res += delata *(abs(e) - (1/2)*delata)
            
    ss_tot = np.sum( (y - np.mean(y))**2)

    r2 = 1 - ( ss_res / ss_tot )

    r2 = 1 - ( modified_ss_res / ss_tot )
    results['mod_rsquared_adj'] = 1 - (1-r2)*(len(y)-1)/(len(y)-X.shape[1]-1)
    
    # [elem for elem in my_list if elem == 'two']
    # sum(1/2*np.sum(e**2) for e in res_e if abs(e) < delata else delata *(abs(e) - (1/2)*delata))

    new_X = np.append(np.ones((len(X),1)), X, axis=1)
    M_S_E = (sum((y-predictions)**2))/(len(new_X)-len(new_X[0]))
    v_b = M_S_E*(np.linalg.inv(np.dot(new_X.T,new_X)).diagonal())
    s_b = np.sqrt(v_b)
    t_b = params/ s_b
    p_val =[2*(1-stats.t.cdf(np.abs(i),(len(new_X)-len(new_X[0])))) for i in t_b]
    
    results['pvalues'] = p_val
    results['params'] = params
    
    # if(func_type == 'linear'):
    #     plt.figure(figsize = (8,6))
    #     # sns.histplot(y - predictions)
    #     # plt.show()
    #     stats.probplot( (y - predictions), dist="norm", plot=pylab)
    #     pylab.show()
    
    return results


def fit_model_activity_cell_types(func_type, X, y ):
    """
    Fit a linear model for gene activity with cell types as features.

    Parameters:
        func_type (str): The type of function to construct features ('linear', 'quadratic', or 'linear_quadratic').
        X (np.ndarray): The input feature matrix with cell types.
        y (np.ndarray): The target gene activity values.

    Returns:
        statsmodels.regression.linear_model.RegressionResults: The results of the linear regression model.
    """

    X = generate_feature_cell_types(func_type, X)
    
    model = sm.OLS(y, X)
    results = model.fit()
 
    return results

    
def fit_best_model_cell_types(target, data, sort_type = 'rsquare', min_target=0, max_target=35):
    """
    Find the best-fitted model for each factor based on cell types as features.

    Parameters:
        target (pd.DataFrame): The target data for gene activity.
        data (pd.DataFrame): The data containing cell type labels.
        sort_type (str, optional): The type of sorting ('pvalue' or 'rsquare'). Default is 'rsquare'.
        min_target (float, optional): The minimum target value. Default is 0.
        max_target (float, optional): The maximum target value. Default is 35.

    Returns:
        dict: A dictionary containing the best-fitted model results for each factor, sorted by the specified sort type.
    """

    assert sort_type in ['pvalue', 'rsquare'], 'sort type parameter must be pvalue, rsquare'
    
    
    
    
    x = np.array(list(data['label']))
    
    all_results = []
    best_r2 = {}
    for index, tf in target.iterrows():
        max_r2 = -1000
        best_results = 0
        best_func = 0
        
        for func_type in ['linear', 'linear_quadratic', 'quadratic']:
        
            results = fit_model_activity_cell_types(func_type, x, tf)
 
            if sort_type == 'rsquare':
                r2 = results.rsquared_adj
                if(r2 > max_r2):
                    best_func = func_type
                    best_results = results
                    max_r2 = r2
            else:
                if(all(i <= 0.05 for i in list(results.pvalues))):
                    r2 = results.rsquared_adj
                    if(r2 > max_r2):
                        best_func = func_type
                        best_results = results
                        max_r2 = r2
                     
        if(best_func):
        
            polyline = np.linspace(min_target,max_target)
    
            polyline = generate_feature_cell_types(best_func, polyline)
            curve = np.matmul(polyline, best_results.params)
            slope = (curve[-1] - curve[0]) / (polyline[-1,1] - polyline[0,1])
            
        
            if best_func == 'linear':
                if best_results.params[1] >= 0:
                    exp_pattern = 'linear up'
                else:
                    exp_pattern = 'linear down'
            
            if best_func == 'quadratic':
                if best_results.params[1] >= 0:
                        exp_pattern = 'quadratic up'
                else:
                    exp_pattern = 'quadratic down'
    
            if best_func == 'linear_quadratic':
                if best_results.params[1] >= 0 and best_results.params[2] >= 0:
                    exp_pattern = 'linear up quadratic up'
                elif best_results.params[1] >= 0 and best_results.params[2] < 0:
                    exp_pattern = 'linear up quadratic down'
                elif best_results.params[1] < 0 and best_results.params[2] >= 0:
                    exp_pattern = 'linear down quadratic up'
                else:
                    exp_pattern = 'linear down quadratic down'
        
        
            all_results.append(best_results)
            best_r2.update({index: [best_func, best_results, exp_pattern, slope]})
            
    if sort_type == 'pvalue':
        # sorted_best = {k: v for k, v in sorted(best_r2.items(), key=lambda item: item[1][1].pvalues[-1])}
        sorted_best = {k: v for k, v in sorted(best_r2.items(), key=lambda item: item[1][1].rsquared_adj, reverse=True)}
    else:
        sorted_best = {k: v for k, v in sorted(best_r2.items(), key=lambda item: item[1][1].rsquared_adj, reverse=True)}
    return sorted_best  
    
  


    
def fit_best_model(target, data, model_type,max_iter_huber,epsilon_huber,pval_thr,modify_r2):
    """
    Find the best-fitted model for each factor.

    Parameters:
        target (pd.DataFrame): The target data for gene activity.
        data (pd.DataFrame): The data containing cell type labels.
        model_type (str): The type of model to fit ('LinearRegression' or 'HuberRegressor').
        max_iter_huber (int): Maximum number of iterations for HuberRegressor.
        epsilon_huber (float): Epsilon parameter for HuberRegressor.
        pval_thr (float): Threshold for p-values of model coefficients.
        modify_r2 (bool): Whether to use modified R-squared. If True, uses modified R-squared; if False, uses regular R-squared.

    Returns:
        dict: A dictionary containing the best-fitted model results for each factor, sorted by R-squared or modified R-squared, depending on the `modify_r2` parameter.
    """

    
    assert model_type in ['LinearRegression', 'HuberRegressor'], 'model type must be LinearRegression, or HuberRegressor'
    
    x = np.array(list(data['label']))
    min_x = min(x)
    max_x = max(x)
    
    all_results = []
    best_r2 = {}
    for index, tf in target.iterrows():
        # if(index == 'ALDH2'):
        max_r2 = -1000
        best_results = 0
        best_func = 0
        
        for func_type in ['linear', 'linear_quadratic', 'quadratic']:
        
            results = fit_model_activity(func_type, x, tf,max_iter_huber,epsilon_huber, model_type)
            # print(results['pvalues'])
            if(all(i <= pval_thr for i in list(results['pvalues']))):
                
                if(modify_r2):
                    r2 = results['mod_rsquared_adj']
                else:
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
            slope = (curve[-1] - curve[0]) / (pline[-1] - pline[0])
            
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
        
    if(modify_r2):
        sorted_best = {k: v for k, v in sorted(best_r2.items(), key=lambda item: item[1][1]['mod_rsquared_adj'], reverse=True)}
    else:
        sorted_best = {k: v for k, v in sorted(best_r2.items(), key=lambda item: item[1][1]['rsquared_adj'], reverse=True)}
    return sorted_best    


   
def data_interst(data, func_type, slope):
    """
    Filter a dictionary of data based on pattern type and slope direction.

    Parameters:
        data (dict): A dictionary containing data entries with associated patterns and slopes.
        func_type (str): The desired pattern type ('linear', 'quadratic', 'linear_quadratic', or 'none').
        slope (str): The desired slope direction ('up', 'down', or 'none').

    Returns:
        dict: A filtered dictionary containing data entries that match the specified pattern type and slope direction.
    """

    dictionary = {}
    if func_type == 'linear' and slope == 'up':
        for key, value in list(data.items()):
            if value[2] == 'linear up':
                dictionary[key] = value
    elif func_type == 'linear' and slope == 'down':
        for key, value in list(data.items()):
            if value[2] == 'linear down':
                dictionary[key] = value
    elif func_type == 'linear_quadratic' and slope == 'up':
        for key, value in list(data.items()):
            if (value[2] == 'linear up quadratic up' or value[2] == 'linear down quadratic up'):
                dictionary[key] = value
    elif func_type == 'linear_quadratic' and slope == 'down':
        for key, value in list(data.items()):
            if (value[2] == 'linear down quadratic down' or value[2] == 'linear up quadratic down'):
                dictionary[key] = value
    elif func_type == 'quadratic' and slope == 'up':
        for key, value in list(data.items()):
            if value[2] == 'quadratic up':
                dictionary[key] = value
    elif func_type == 'quadratic' and slope == 'down':
        for key, value in list(data.items()):
            if value[2] == 'quadratic down':
                dictionary[key] = value     
    elif func_type == 'linear' and slope == 'none':
        for key, value in list(data.items()):
            if value[0] == 'linear':
                dictionary[key] = value
    elif func_type == 'quadratic' and slope == 'none':
        for key, value in list(data.items()):
            if value[0] == 'quadratic':
                dictionary[key] = value
    elif func_type == 'linear_quadratic' and slope == 'none':
        for key, value in list(data.items()):
            if value[0] == 'linear_quadratic':
                dictionary[key] = value
                
    return(dictionary)

def display_count_patterns(data, name):
    """
    Display the count of statistically significant genes and their expression patterns for a given cell type.

    Parameters:
        data (dict): A dictionary containing gene data for a specific cell type.
        name (str): The name or identifier of the cell type.

    Returns:
        None: This function prints the information but does not return a value.
    """

    print("For cell_type " + str(name) + ", " + str(len(data.keys())) + " genes are statistically significant.")
    df = pd.DataFrame({'Expression pattern': [item[2] for item in data.values()]})
    df2 = df.groupby(['Expression pattern'])['Expression pattern'].count().reset_index(name='counts')
    print(df2)

def data_interst_pattern(data, pattern):
    """
    Filter data based on a more specific expression pattern.

    Parameters:
        data (dict): A dictionary containing gene data.
        pattern (str): The specific expression pattern to filter for.

    Returns:
        dict: A dictionary containing genes that match the specified expression pattern.
    """

    dictionary = {}
    
    for key, value in list(data.items()):
        if value[2] == pattern:
            dictionary[key] = value
       
    return(dictionary)

def make_table(dictionary, column_names):
    """
    Create a table from a dictionary containing gene data.

    Parameters:
        dictionary (dict): A dictionary containing gene data.
        column_names (list): A list of column names for the table.

    Returns:
        pd.DataFrame: A Pandas DataFrame representing the table.
        """

    
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

    stats = pd.Series([i['params'] for i in [i[1] for i in list(dictionary.values())] ])
    df = pd.DataFrame(stats.values.tolist(), index=stats.index)

    # intercept1
    stats = np.array(df.iloc[:,0])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    # feature2
    stats = np.array(df.iloc[:,1])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    # feature3
    stats = np.array(df.iloc[:,2])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    # compute adjusted p-values
    all_pvals = np.array([i[4] for i in list(dictionary.values())])
    all_pvals = smf.multitest.multipletests(list(all_pvals), method='fdr_bh')[1]
    all_pvals = all_pvals.reshape(all_pvals.shape[0],1)
    mat = np.append(mat, all_pvals, axis = 1)

    stats = np.array([i['rsquared_adj'] for i in tmp])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    # stats = np.array([i.f_pvalue for i in tmp])
    # stats = stats.reshape(stats.shape[0],1)
    # mat = np.append(mat, stats, axis = 1)

    # stats = np.array([i.fvalue for i in tmp])
    # stats = stats.reshape(stats.shape[0],1)
    # mat = np.append(mat, stats, axis = 1)
    
    mat_df = pd.DataFrame(mat ,columns=column_names)
    return mat_df

def filter_table(mat_df, feature, thr = 0.05, get_type = 'low'):
    """
    filter table based on multiple testing corrections
    

    Parameters
    ----------
    mat_df : Dataframe
        Data table of statistic results of each feature.
    pval_thr : float
        A threshold for adjusted p-vallue to remove features.
    filter_type : string
        Specify filter type.

    Returns
    -------
    filtered table as Dataframe.

    """
    assert get_type in ['high', 'low'], 'filter type parameter must be high or low'
    
    if get_type == 'high':
        return mat_df[pd.to_numeric(mat_df[feature]) > thr]
    elif get_type == 'low':
        return mat_df[pd.to_numeric(mat_df[feature]) < thr]
    
    
def statistics(mat_df):
    df=mat_df.groupby(["Expression pattern"])['Expression pattern'].count().reset_index(name='count').sort_values(['count'], ascending=False)
                             
    return df[["Expression pattern",'count']]


def save_data(dictionary, column_names,save_path,name,p_val,pro,gprofil=False):
    """
    Generate a table from a dictionary of gene data and save it along with associated files.

    Parameters:
        dictionary (dict): A dictionary containing gene data.
        column_names (list): A list of column names for the table.
        save_path (str): The path where the data and associated files will be saved.
        name (str): The name for the saved files and folder.
        p_val (float): The threshold for the adjusted p-value for significance.
        pro (pd.DataFrame): A Pandas DataFrame containing additional protein-related data.
        gprofil (bool): If True, perform GProfiler analysis and save the results.

    Returns:
        None
    """
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

    stats = pd.Series([i['params'] for i in [i[1] for i in list(dictionary.values())] ])
    df = pd.DataFrame(stats.values.tolist(), index=stats.index)

    # intercept1
    stats = np.array(df.iloc[:,0])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    # feature2
    stats = np.array(df.iloc[:,1])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)
   

    # feature3
    if(df.shape[1]>=3):   
        stats = np.array(df.iloc[:,2])
        stats = stats.reshape(stats.shape[0],1)
        mat = np.append(mat, stats, axis = 1)
    else:
        r=column_names.pop(-4)   
    # compute adjusted p-values
    all_pvals = np.array([i[4] for i in list(dictionary.values())])
    all_pvals = smf.multitest.multipletests(list(all_pvals), method='fdr_bh')[1]
    all_pvals = all_pvals.reshape(all_pvals.shape[0],1)
    mat = np.append(mat, all_pvals, axis = 1)
    
    stats = np.array([i['rsquared_adj'] for i in tmp])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)
    
    
    stats = np.array([i['mod_rsquared_adj'] for i in tmp])
    stats = stats.reshape(stats.shape[0],1)
    mat = np.append(mat, stats, axis = 1)

    
    mat_df = pd.DataFrame(mat ,columns=column_names)
    
   
    
    if not os.path.exists(save_path+'/'+name):
        os.makedirs(save_path+'/'+name)
        
    mat_df['adjusted P-value'] = mat_df['adjusted P-value'].astype(float)
    mat_df=mat_df[mat_df['adjusted P-value'] <= p_val]
    
    
    if column_names[0]=='Cell name':
        mat_df.to_csv(save_path+'/'+name+'_Report.csv')
        
    else:
        mat_df=pd.merge(mat_df, pro, on='Gene ID')

        mat_df.to_csv(save_path+'/'+name+'/Whole_expressions.csv')  


        df=statistics(mat_df)
        print ('For this cell_type, p-value of  {} genes are statistically significant.'.format(len(mat_df)))
        print(df)


        lists_names= list(mat_df['Expression pattern'].unique())
        for exp in lists_names:
            df_exp=mat_df[mat_df['Expression pattern']==exp]
            if not os.path.exists(save_path+name+'/'+'PILOT'):
                os.makedirs(save_path+name+'/'+'PILOT')
            df_exp.to_csv(save_path+name+'/'+'PILOT/'+exp+'.csv')

            if gprofil:
                gp = GProfiler(return_dataframe=True,)
                genes_lists=list(df_exp['Gene ID'][0:])
                df=gp.profile(organism='hsapiens',
                        query=genes_lists,no_evidences=False)
                if not os.path.exists(save_path+name+'/'+'Gprofiler'):
                    os.makedirs(save_path+name+'/'+'Gprofiler')
                df.to_csv(save_path+name+'/'+'Gprofiler'+'/'+exp+'.csv') 




        print("data saved successfully")


    
   
