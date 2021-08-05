# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:40:08 2019

@author: da.salazarb
"""
# In[packages]
import os
import dataTCGA
import numpy as np
import pandas as pd
from csv import writer
from sklearn.model_selection import ParameterGrid
from func_JNMF_oneStep import hp_metrics_plots, new_syntheticData, LoadSyntheticData

# In[Datasets]: Load or create a synthetic datasets
dim_features = {'mirna': 750, 'mrna': 1400, 'cnv': 1250}; K = 5;
dim_samples = [150,150]; numero_grupos_W = 5; numero_grupos_H = 5;
path = os.path.dirname(os.getcwd()).replace("\\", "/")+"/"; print(path)## e.g. path = "D:/"
dataTCGA.dataTCGA.pathFeatureLabel = path+"pathFeatureLabel";
ruta = dataTCGA.dataTCGA.pathFeatureLabel; 

## For the USER:
## Please, modify this path to execute R scripts
dataTCGA.dataTCGA.command = 'C:/Program Files/R/R-3.6.2/bin/Rscript.exe';
## End for the USER

## First, generate new synthetic datasets, then run LoadSyntheticData.
# new_syntheticData(dim_features, dim_samples, numero_grupos_W, numero_grupos_H, ruta, K)

## To load synthetic datasets
SynData = LoadSyntheticData(ruta, dim_features)

# In[algorithm]: M&M-jNMF + Synthetic data
from MMJNMF import mjnmf

tabla = pd.DataFrame()
numero = 0
useMalaCard=False
save = True

# param_grid = {'r1': [1e-6, 0.1, 1, 10], 'L1': [1e-6, 0.1, 1, 10], 'L2': [1e-6, 0.1, 1, 10], 
#               'r2': [1e-6, 0.1, 1, 10], 'd1': [1e-6, 0.1, 1, 10], 
#               'noise_proportion': [0], 'tol': [1e-6,1e-8]};

param_grid = {'r1': [1e-6], 'L1': [1e-6], 'L2': [1e-6], 
              'r2': [1e-6], 'd1': [1e-6], 
              'noise_proportion': [0], 'tol': [1e-6]};

grid = ParameterGrid(param_grid); i=0;

for params in grid:
# for i in range(100):
    print("                ************************************************* ")
    print("                *                                               * ")
    print('                                     '+' i = '+str(i)+'                            ')
    print("                *                                               * ")
    print("                ************************************************* ")
    
    np.random.seed(i+1)
    param = dict();
    
    param['K']=5;
    # param['r1']=0;
    # param['r2']=0;
    # param['L1']=0;
    # param['L2']=0;
    # param['d1']=0;
    # param['noise_proportion'] = 0;
    
    # param['r1']=np.random.choice([1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000]);
    # param['r2']=np.random.choice([1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000]);
    # param['L1']=np.random.choice([1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000]);
    # param['L2']=np.random.choice([1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000]);
    # param['d1']=np.random.choice([1e-6, 1e-5, 1e-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000]);
    
    # param['K']=params['K'];
    param['r1']=params['r1']; param['r2']=params['r2'];
    param['L1']=params['L1']; param['L2']=params['L2'];
    param['d1']=params['d1'];
    param['noise_proportion'] = params['noise_proportion']
    
    # param['K']=np.random.randint(3, 20);
    # param['r1']=np.round(np.random.uniform(1e-8,2.1), 6);
    # param['r2']=np.round(np.random.uniform(1e-8,2.1), 6);
    # param['L1']=np.round(np.random.uniform(1e-8,2.1), 6);
    # param['L2']=np.round(np.random.uniform(1e-8,2.1), 6);
    # param['d1']=np.round(np.random.uniform(1e-8,2.1), 6);
    # param['noise_proportion'] = np.random.choice([0,1]);
    
    param['tol']=params['tol']; param['maxiter']=600; param['sample']="sample"; param['command']="";
    param['repeat']=2;param['path']=ruta; param['drug_to_heatmapComplex']=""; param["quick_version"]="";
    param["W_simulated"]=SynData['W_simulated']; param["H_record_simulated"]=SynData['H_record_simulated'];
    param["theta_record"]=SynData['theta_record'];param["R_record"]=SynData['R_record']; param['complete_enrichment']="";
    param["numero_grupos_W"]=numero_grupos_W; param["numero_grupos_H"]=numero_grupos_H;
    param["pac_por_comod"]=""; param["pac_por_comod_ccle"]="";
    param['max_columna_H_real']=SynData['max_columna_H_'];param['max_columna_W_real']=SynData['max_columna_W_'];
    param['saveReplicas']=False; param['loadReplicas']=False; param['saveBest_H_W']=False; param['useMalaCard']=False;
    param['noise_data_0']=SynData['noise_data_0'];param['noise_data_1']=SynData['noise_data_1'];
    
    # ccle, tcga, best_W, best_H, comodule, module_patient_cell, rho_C, rec_surv, delta_control, best_stop_control, silhouette_avg;
    paramResult = mjnmf([], param);
    
    ## Residual sum of squares
    rss_ccle = {'rss_ccle_'+k: sum(sum( (v-(np.dot(paramResult['best_W']["W0_ccle"],paramResult['best_H'][k])))**2 )) for k, v in paramResult['ccle'].items()};
    rss_tcga = {'rss_tcga_'+k: sum(sum( (v-(np.dot(paramResult['best_W']["W0_tcga"],paramResult['best_H'][k])))**2 )) for k, v in paramResult['tcga'].items()};
    
    r2_ccle = {'r2_ccle_'+k: 1-((rss_ccle['rss_ccle_'+k] * v.shape[0] * v.shape[1])/(sum(sum( v**2 )) * (v.shape[0]*v.shape[1]) - param['K'] * (v.shape[0] + v.shape[1])) ) for k, v in paramResult['ccle'].items()};
    r2_tcga = {'r2_tcga_'+k: 1-((rss_tcga['rss_tcga_'+k] * v.shape[0] * v.shape[1])/(sum(sum( v**2 )) * (v.shape[0]*v.shape[1]) - param['K'] * (v.shape[0] + v.shape[1])) ) for k, v in paramResult['tcga'].items()};       
    
    parcial = []
    for v in paramResult['rho_C'].values():
        parcial.append(v)
        
    if tabla.empty:
        for valor in ['K', 'r1', 'L1', 'L2', 'r2', 'd1', 'No. repetitions', 'tolerance', 'last_stop_control', 'noise_proportion']:
            tabla.insert(tabla.shape[1],valor, None)
        for perfil in paramResult['rho_C'].keys():
            tabla.insert(tabla.shape[1],'CoephCoeff '+perfil, None)
            
        tempLista = [rss_ccle, rss_tcga ,r2_ccle, r2_tcga]
        for l in tempLista:
            temp_l = [k for k in l.keys()]; ## rss_ccle
            for valor in temp_l:
                tabla.insert(tabla.shape[1], valor, None)
        # for perfil in paramResult['silhouette_avg'].keys():
        #     tabla.insert(tabla.shape[1],'silhouette '+perfil, None)
        tempLista = [i+"_"+j for j in paramResult['rec_surv'].columns for i in paramResult['rec_surv'].index];
        for item in tempLista:
            tabla.insert(tabla.shape[1],item, None)
        tempLista = [j+"_"+i for j in paramResult['rec_auc_XY'].columns for i in paramResult['rec_auc_XY'].index];
        for item in tempLista:
            tabla.insert(tabla.shape[1],item, None)
    
    class_metrics = []
    for c in paramResult['rec_surv'].columns:
        class_metrics = class_metrics  + [v for v in paramResult['rec_surv'][c]]
        
    class_metrics_XY = []
    for c in paramResult['rec_auc_XY'].columns:
        class_metrics_XY = class_metrics_XY + [v for v in paramResult['rec_auc_XY'][c]]
    
    tabla.loc[numero] = [param['K'], param['r1'], param['L1'], param['L2'], param['r2'], param['d1'], param['repeat'], param['tol'], paramResult['best_stop_control'][0,-1],  
                         param['noise_proportion']] + parcial + [v for v in rss_ccle.values()] + [v for v in rss_tcga.values()] + [v for v in r2_ccle.values()] + [v for v in r2_tcga.values()] + class_metrics + class_metrics_XY; #  [v for v in paramResult['silhouette_avg'].values()] + 
    
    if os.path.isfile(dataTCGA.dataTCGA.pathFeatureLabel+"/co-mod_tabulated_results/Supplementary_File_S1.csv") == False:
        tabla.to_csv(dataTCGA.dataTCGA.pathFeatureLabel+"/co-mod_tabulated_results/Supplementary_File_S1.csv")
    else:
        with open(dataTCGA.dataTCGA.pathFeatureLabel+"/co-mod_tabulated_results/Supplementary_File_S1.csv", 'a+', newline='') as write_obj:
            csv_writer = writer(write_obj)
            csv_writer.writerow(np.hstack([numero,tabla.iloc[numero,:]]))
    
    numero+=1; i+=1;
    
    print("    ---------------------------------------------------------- ")
    print("    * Classification metrics")
    print(paramResult['rec_surv'])
    print("    ---------------------------------------------------------- ")
    
    print("    ---------------------------------------------------------- ")
    print("    * Classification metrics for X_tcga and Y_ccle")
    print(paramResult['rec_auc_XY'])
    print("    ---------------------------------------------------------- \n")

print("Ok, Supp_File_S1 is saved in co-mod_tabla_metricas!")

## plots for grid results
# hp_metrics_plots(param['sample'], dataTCGA.dataTCGA.pathFeatureLabel, path, dataTCGA.dataTCGA.command, "")
