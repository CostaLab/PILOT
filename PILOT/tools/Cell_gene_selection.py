

#%%% load libraries
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


#%%% load functions    
def compute_metrics(y_test, y_pred):
    r2 = r2_score(y_test, y_pred)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    mae = mean_absolute_error(y_test, y_pred)
    cor = stats.pearsonr(y_test, y_pred)[0]
    return (r2, rmse, mae, cor)


def loadTarget(path, name):
    """
    loads the file count matrix
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
    loads the info
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
    calculate sparsity of the data
    """
    non_zero = np.count_nonzero(data)
    total_val = np.product(data.shape)
    sparsity = (total_val - non_zero) / total_val
    return(sparsity)

def generate_feature_cell_types(func_type, X):
    """
    construct the features matrix
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
    construct the features matrix
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
    fit the linear model
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

def qq_plot_gene(target, data, sorted_best, gene_name):
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

def fit_model_activity_cell_types(func_type, X, y ):
    """
    fit the linear model
    """
    X = generate_feature_cell_types(func_type, X)
    
    model = sm.OLS(y, X)
    results = model.fit()
 
    return results

    
def fit_best_model_cell_types(target, data, sort_type = 'rsquare', min_target=0, max_target=35):
    """
    find best fitted model for each factor 
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
    find best fitted model for each factor 
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


def plot_best_matches_cell_types(target, data,df,sorted_best, scale_name, min_target=0, max_target=35,num=11,width=25,height=25,xlim=4,point_size=100,color_back=None,fontsize=28,alpha=1,cmap='viridis'):
    """
    plot 16 best fitted model
    """
    

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
    plot 4 of each best fitted pattern
    
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
    x = data['label']

    best_tf = list(target.loc[[gene_name]].values.ravel())
    data = pd.DataFrame({'x': x, 'y':best_tf})
    
    labels=np.sort([y for y in list(data.x.unique())])
    fig, axes = joypy.joyplot(data, by="x", column="y", labels=labels, 
                          grid="y", linewidth=1, legend=False, overlap=0.5, figsize=(6,5),kind="kde", bins=80,
                          title=gene_name, ylabels=False,
                          colormap=cm.autumn_r)

def plot_pval_rsq_correlation(table, feature1, feature2, show_fit = True, log_transform = False):
    
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
    plot 9 best fitted model based on up or down regulation pattern
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
    plot 6 best fitted model for each pattern
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
            
   
def data_interst(data, func_type, slope):
    """
    filter data based on the pattern
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
    print("For cell_type " + str(name) + ", " + str(len(data.keys())) + " genes are statistically significant.")
    df = pd.DataFrame({'Expression pattern': [item[2] for item in data.values()]})
    df2 = df.groupby(['Expression pattern'])['Expression pattern'].count().reset_index(name='counts')
    print(df2)

def data_interst_pattern(data, pattern):
    """
    filter data based on the more specific pattern
    """
    dictionary = {}
    
    for key, value in list(data.items()):
        if value[2] == pattern:
            dictionary[key] = value
       
    return(dictionary)

def make_table(dictionary, column_names):
    """
    make table of the data
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
    make table and save the data
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


    
   
