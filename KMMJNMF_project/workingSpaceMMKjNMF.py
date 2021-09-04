# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:40:08 2019

@author: da.salazarb
"""
# %% In[*** Loading packages, data and parameters + MMjNMF algorithm ***]
# %reset_selective -f "^(?!merged)"
import os
import sys
import numpy as np
import pandas as pd
from csv import writer
import fusion_CCLE_TCGA
from datetime import datetime
from sklearn.model_selection import ParameterGrid
# %%
arguments = sys.argv[1:];
if len(arguments) == 0:
    print("------------------------------------------------------------------------------------")
    print("At least include the path/to/kmmjnmf, e.g., python workingSpaceOmic.py path/to/kmmjnmf")
    print("")
    print("")
    # break
elif len(arguments) == 1:
    print("---------------------------------------------------------------------------------------------------")
    print("It will run M&M-KjNMF with the default parameters, please check how to modify them in the README.md.")
    print("")
    print("")
    path = arguments[0]
    param_grid = {'K': [20],
                  'r1': [5e-2], 'r2': [1e-3],
                  'L1': [10], 'L2': [10],
                  'd1': [1e-3], 'o1': [1e-7], 'o2': [1e-7],
                  'sigma_ccle': [1], 'sigma_tcga': [1],
                  'sigma_ccle_diff': [2], 'sigma_tcga_diff': [2]
                  };
elif len(arguments) > 1 and len(arguments) < 11:
    print("---------------------------------------------------------------------------------------------")
    print("There are seven parameters required to execute M&M-KjNMF: path/to/mmjnmf K r1 r2 L1 L2 d1")
    print("For example include them as follows: path/to/mmjnmf [60] [3.5e-6] [3.5e-6] [10] [10] [3.5e-3] [1e-3] [1e-3] [1e-3] [1e-3]")
    print("")
    print("")
elif len(arguments) == 11:
    print("--------------------------------------------------")
    print("Ok ... M&M-KjNMF will be executing in a moment ... ")
    print("")
    print("")
    path = arguments[0]
    param_grid = {'K': [int(item) for item in arguments[1].strip('][').split(',')],
                  'r1': [float(item) for item in arguments[2].strip('][').split(',')], 'r2': [float(item) for item in arguments[3].strip('][').split(',')],
                  'L1': [float(item) for item in arguments[4].strip('][').split(',')], 'L2': [float(item) for item in arguments[5].strip('][').split(',')], 
                  'd1': [float(item) for item in arguments[6].strip('][').split(',')],
                  'sigma_ccle': [float(item) for item in arguments[7].strip('][').split(',')], 'sigma_tcga': [float(item) for item in arguments[8].strip('][').split(',')],
                  'sigma_ccle_diff': [float(item) for item in arguments[9].strip('][').split(',')], 'sigma_tcga_diff': [float(item) for item in arguments[10].strip('][').split(',')]
                  };
#%%
## profiles!!
# path = os.path.dirname(os.getcwd()).replace("\\", "/")+"/"; #print(path)## e.g. path = "D:/"
merged = fusion_CCLE_TCGA.fusion_CCLE_TCGA(path);
merged.loadData();

## Constraints!!
merged.constraints_theta_method(path);
merged.constraints_r_method(path);
# %%
from KMMJNMF import kmmjnmf
tabla = pd.DataFrame(); numero = 0;
np.random.seed(13)

grid = ParameterGrid(param_grid); #len(grid)

for params in grid:
    print("\n                ************************************************* ")
    print("                *                                               * ")
    print('                                 '+' i = '+str(numero+1)+' of '+str(len(grid))+'                            ')
    print("                *                                               * ")
    print("                ************************************************* ")
    param = dict();
    
    param['K']=params['K'];
    param['r1']=params['r1']; param['r2']=params['r2'];
    
    param['L1']=params['L1']; param['L2']=params['L2'];
    param['d1']=params['d1']; 
    param['o1']=params['o1']; param['o2']=params['o2'];
    
    param['sigma_ccle']=params['sigma_ccle']; param['sigma_tcga']=params['sigma_tcga'];
    param['sigma_ccle_diff']=params['sigma_ccle_diff']; param['sigma_tcga_diff']=params['sigma_tcga_diff'];
    
    param['index_train_tcga']=merged.index_train_tcga; param['index_test_tcga']=merged.index_test_tcga;
    param['index_train_ccle']=merged.index_train_ccle; param['index_test_ccle']=merged.index_test_ccle;
    
    param['tol']=1e-7; param['maxiter']=500; param['repeat']=3; param['enrichment']=True;
    
    paramResultKernels = kmmjnmf(merged, param);
    
    if tabla.empty:
        for valor in ['DateAndHour','K', 'r1', 'L1', 'L2', 'r2', 'd1','o1','o2','sigma_ccle','sigma_tcga','sigma_ccle_diff','sigma_tcga_diff', 'repetition', 'last_stop_control' ,'last_delta_control', 'objectiveFunction', 'WeightedAverage(0.2*rho+0.4AUC+0.2*BioScore)']:
            tabla.insert(tabla.shape[1],valor, None)
        for perfil in paramResultKernels['rho_C'].keys():
            tabla.insert(tabla.shape[1],'rho_'+perfil, None)
        for perfil in paramResultKernels['AUC_H_average'].keys():
            tabla.insert(tabla.shape[1],'AUC_H_'+perfil,None)
  
        tempLista = [j+"_"+i for j in paramResultKernels['rawCod_stats'].columns for i in paramResultKernels['rawCod_stats'].index];
        for item in tempLista:
            tabla.insert(tabla.shape[1],item, None)
        for valor in ['BioRelevanceScore', 'codEnrich_over_codTotal', 'densityPenalization', 'TotalDots/DotsRef', 'geneRatio_avg', 'geneRatio_sd']:
            tabla.insert(tabla.shape[1],valor, None)
            
    #--------------------M&M-KjNMF RESULTS---------------------------#
    parcial_kernels = []
    for v in paramResultKernels['rho_C'].values():
        parcial_kernels.append(v)

    raw_cod_metrics_kernel = []
    for c in paramResultKernels['rawCod_stats'].columns:
        raw_cod_metrics_kernel = raw_cod_metrics_kernel  + [v for v in paramResultKernels['rawCod_stats'][c]]
                
    cod_metrics_kernel = []
    for c in range(len(paramResultKernels['dotScore_stats'])):
        cod_metrics_kernel.append(paramResultKernels['dotScore_stats'][c])

    now = datetime.now()
    tabla.loc[numero] = [now, param['K'], param['r1'], param['L1'], param['L2'], param['r2'], param['d1'],param['o1'],param['o2'],param['sigma_ccle'],
                         param['sigma_tcga'],param['sigma_ccle_diff'],param['sigma_tcga_diff'],param['repeat'],
                         paramResultKernels['best_stop_control'][0][len(np.where(paramResultKernels['best_stop_control'] > 0)[1])-1],
                         paramResultKernels['best_delta_control'][0][len(np.where(paramResultKernels['best_delta_control'] > 0)[1])-1], paramResultKernels['FO'],
                         paramResultKernels['weightedAverage']] + parcial_kernels +[v for v in paramResultKernels['AUC_H_average'].values()]+raw_cod_metrics_kernel+cod_metrics_kernel;
    
    if os.path.isfile(path+"/pathFeatureLabel/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv") == False:
        tabla.to_csv(path+"/pathFeatureLabel/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv")
    else:
        with open(path+"/pathFeatureLabel/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv", 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj)
            csv_writer.writerow(np.hstack([numero,tabla.iloc[numero,:]]))
     
    numero += 1;
    
print("Ok, Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv is saved in /pathFeatureLabel/co-mod_tabulated_results!")
print("Ok, Variable clusters are saved in individual folders in pathFeatureLabel folder, e.g., 0_drug folder!")
