# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 21:07:43 2019

@author: da.salazarb
"""

#%% Funciones JNMF.py
# import re
import os
import subprocess
import numpy as np
import pandas as pd
from sklearn import metrics
from munkres import Munkres
from sklearn.metrics import pairwise
from collections import defaultdict
from itertools import combinations, permutations, product

def roc_curve_MjNMF(y_true, y_pred):
    ## https://stackoverflow.com/questions/61321778/how-to-calculate-tpr-and-fpr-in-python-without-using-sklearn
    # from sklearn.metrics import confusion_matrix
    # cnf_matrix = confusion_matrix(y_true, y_pred)
    FP = np.sum((y_pred == 1) & (y_true == 0))
    FN = np.sum((y_pred == 0) & (y_true == 1))
    TP = np.sum((y_pred == 1) & (y_true == 1))
    TN = np.sum((y_pred == 0) & (y_true == 0))
    
    # Sensitivity, hit rate, recall, or true positive rate
    TPR = TP/(TP+FN)
    FPR = FP/(FP+TN)
    
    return TPR, FPR
    
def norm_min_max(profile):
    from sklearn.preprocessing import MinMaxScaler
    scaler = MinMaxScaler(copy=False, feature_range=(0,1))
    profile = scaler.fit_transform(profile)
    return profile

# %% KMMJNMF Functions
def normalize_AH_rowH(A0,H0,value=1):
    """
    Normaliza A y H. Normaliza filas de H (gran H) y columnas de A 

    Args
    -------
    A (dict) 
    H (dict)
    value (int)

    Returns
    -------
    Matrices A y H normalizadas
    """
    vec_sumN = {k: v.shape[1] for k, v in H0.items()}; H = np.concatenate([v for v in H0.values()],axis=1);#gran H kx(p1+..+pn), una HI al lado de la otra
    sh = np.sum(H, axis=1); sh[sh == 0] = value; psh = np.divide(value, sh);
    A = defaultdict(lambda: np.ndarray(0)); A['A0_ccle']={k: np.dot(v,np.diag(sh)) for k, v in A0['A0_ccle'].items()};  A['A0_tcga']={k: np.dot(v,np.diag(sh)) for k, v in A0['A0_tcga'].items()}; H=np.dot(np.diag(psh),H);    
    H_init = defaultdict(lambda: dict());
    for i, t in enumerate(H0.keys()):
        H_init[t] = H[:,0:vec_sumN[t]]
        H = np.delete(H, obj=slice(vec_sumN[t]), axis=1)
    return A, H_init,sh

def compute_diff(H_record,Hsumtheta_record,sumHRt_record,AtKA,KIJ_ccle,KIJ_tcga,A,L1,L2,r1,r2,d1,o1,o2,K,miss_profile):
    numN = len(H_record);
    D2 = 0;
    #sum_{i}
    for i in range(0,numN):
        Hsumtheta = Hsumtheta_record[list(H_record.keys())[i]];
        HtheHt = np.dot(Hsumtheta, H_record[list(H_record.keys())[i]].T);
        D2 = D2 + np.trace(HtheHt);
    diff2 = -(L1/2)*D2
    D3 = 0;
    for i in range(0,numN):
        sumHRt = sumHRt_record[list(H_record.keys())[i]];
        HRHt = np.dot(H_record[list(H_record.keys())[i]],sumHRt.T);
        D3 = D3 + np.trace(HRHt);
    diff3 = -L2*D3;
    D4_ccle=0
    D4_tcga=0
    for i in H_record.keys():
        #D4_ccle=D4_ccle+np.trace(AtKA['AtKA_ccle'][i])
        #Alt. FO
        D4_ccle=D4_ccle+np.trace(np.dot(A['A0_ccle'][i].T,A['A0_ccle'][i]))
        if not miss_profile[0] in i:
            #D4_tcga=D4_tcga+np.trace(AtKA['AtKA_tcga'][i])  
            #Alt. FO
            D4_tcga=D4_tcga+np.trace(np.dot(A['A0_tcga'][i].T,A['A0_tcga'][i]))
    diff4_ccle = r1*D4_ccle
    diff4_tcga = r2*D4_tcga
    D5 = 0;
    for i in H_record.keys():
        #Transpuesta (HHt)t=HHt
        HHt= np.dot(H_record[i],H_record[i].T).T;
        D5 = D5 + np.dot(np.dot(np.ones([1,K]),HHt),np.ones([K,1]));
    diff5 = d1*D5;
    D6_ccle=0
    D6_tcga=0
    for i in H_record.keys():
        ci=int(float(i[0:1]))
        for j in H_record.keys():
            cj=int(float(j[0:1]))
            if i!=j:
                D6_ccle=D6_ccle+np.trace(AtKA['AtKA_ccle'][i]-(2*np.dot(np.dot(A['A0_ccle'][i].T,KIJ_ccle[ci][cj]),A['A0_ccle'][j]))+AtKA['AtKA_ccle'][j])
                if not miss_profile[0] in i and not miss_profile[0] in j:
                    #print('i',i)
                    #print('j',j)
                    D6_tcga=D6_tcga+np.trace(AtKA['AtKA_tcga'][i]-(2*np.dot(np.dot(A['A0_tcga'][i].T,KIJ_tcga[ci][cj]),A['A0_tcga'][j]))+AtKA['AtKA_tcga'][j])
                    # print('1Sum',np.trace(AtKA['AtKA_tcga'][i]))
                    # print('2sum',np.trace((2*np.dot(np.dot(A['A0_tcga'][i].T,KIJ_tcga[ci][cj]),A['A0_tcga'][j]))))
                    # print('3sum',np.trace(AtKA['AtKA_tcga'][j]))
                    #print('sumaActual',np.trace(AtKA['AtKA_tcga'][i]-(2*np.dot(np.dot(A['A0_tcga'][i].T,KIJ_tcga[ci][cj]),A['A0_tcga'][j]))+AtKA['AtKA_tcga'][j]))
                    #print('D6_tcga',D6_tcga)
    diff6_ccle=o1*D6_ccle
    diff6_tcga=o2*D6_tcga                    
    return diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga

def StopCriterion_rule1(Kn,A,H,AtKA,diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga,miss_profile):
    diff1_ccle=0
    diff1_tcga=0
    for i in H.keys():
        diff1_ccle=diff1_ccle+np.trace(Kn['K_ccle'][i]-(2*np.dot(np.dot(Kn['K_ccle'][i],A['A0_ccle'][i]),H[i]))+(np.dot(np.dot(H[i].T,AtKA['AtKA_ccle'][i]),H[i])));
        if not miss_profile[0] in i:
            diff1_tcga=diff1_tcga+np.trace(Kn['K_tcga'][i]-(2*np.dot(np.dot(Kn['K_tcga'][i],A['A0_tcga'][i]),H[i]))+(np.dot(np.dot(H[i].T,AtKA['AtKA_tcga'][i]),H[i])));
    retVal = diff1_ccle+diff1_tcga+diff2+diff3+diff4_ccle+diff4_tcga+diff5+diff6_ccle+diff6_tcga;
    return retVal,diff1_ccle, diff1_tcga

# %%
def compute_XHt_HHt(X_record_ccle, X_record_tcga, H_record,K):
    """
    Computa XHt y HHt.

    Args
    -------
    X_record (dict): data de omicas orignales.

    H_record (dict): Hi donde i corresponde a cada una de las omicas.

    Returns
    -------
    sumXHt (array) calcula

    .. math::
        (\sum_{I}X_{I}H_{I}^T)_{ij}

    sumHHt (dict de array) calcula

    .. math::
        (\sum_{I}H_{I}H_{I}^T)_{ij}

    XHt_record (dict) calcula para cada omica

    .. math:: 
        X_{I}H_{I}^T

    HHt_record (dict) calcula para cada omica

    .. math::
        H_{I}H_{I}^T

    """
    numN = len(X_record_ccle); M=X_record_ccle[list(X_record_ccle.keys())[0]].shape[0]; N=X_record_tcga[list(X_record_tcga.keys())[0]].shape[0];
    XHt_record_ccle=defaultdict(lambda: np.ndarray(0)); XHt_record_tcga=defaultdict(lambda: np.ndarray(0));
    HHt_record = defaultdict(lambda: np.ndarray(0)); sumHHt = defaultdict(lambda: np.zeros([K,K]));
    sumXHt_ccle = np.zeros([M,K]); sumXHt_tcga = np.zeros([N,K]);
    for i in range(0,numN):
        ### XHt_record para ccle y tcga (try: except:)
        XHt_record_ccle[list(X_record_ccle.keys())[i]] = np.dot(X_record_ccle[list(X_record_ccle.keys())[i]],H_record[list(X_record_ccle.keys())[i]].T);
        try: ## cuando haya perfil de drugs, habra un perfil más y saldra error, por lo que solo suma los que hay para tcga
            XHt_record_tcga[list(X_record_tcga.keys())[i]] = np.dot(X_record_tcga[list(X_record_tcga.keys())[i]],H_record[list(X_record_tcga.keys())[i]].T);
        except:
            pass
        ### HHt_record usando ccle index porque tienen todos los perfiles, contando drug
        HHt_record[list(X_record_ccle.keys())[i]] = np.dot(H_record[list(X_record_ccle.keys())[i]],H_record[list(X_record_ccle.keys())[i]].T);
        ### sumXHt para ccle y tcga (try: except:)
        sumXHt_ccle = sumXHt_ccle + XHt_record_ccle[list(X_record_ccle.keys())[i]]; 
        try:
            sumXHt_tcga = sumXHt_tcga + XHt_record_tcga[list(X_record_tcga.keys())[i]];
        except:
            pass
        ### sumHHt para ccle y tcga (try: except:)
        sumHHt["HHt_ccle"] = sumHHt["HHt_ccle"] + HHt_record[list(X_record_ccle.keys())[i]];
        # try: 
        #     sumHHt["HHt_tcga"] = sumHHt["HHt_tcga"] + HHt_record[list(X_record_ccle.keys())[i]];
        # except:
        #     pass
        if i!=0:
            sumHHt["HHt_tcga"] = sumHHt["HHt_tcga"] + HHt_record[list(X_record_ccle.keys())[i]];
        
    return sumXHt_ccle,sumXHt_tcga,sumHHt,XHt_record_ccle,XHt_record_tcga,HHt_record

def compute_HsumTheta_sumHRt_record(H_record, theta_record, R_record):
    """
    Computa Hsumtheta_record y sumHRt_record para cada omica

    Args
    -------
    H_record (dict)

    theta_record (dict)

    sumHRt_record (dict)

    Returns
    -------
    Hsumtheta_record (dict) calcula

    .. math::
        \sum_{t}H_{I}[\\theta_{I}+(\\theta_{I}^{(t)})^{T}]
    
    sumHRt_record (dict) calculo

    .. math::
        \sum_{I\\neq J}H_{J}R_{IJ}^{T}
    """
    numN = len(H_record); Hsumtheta_record = defaultdict(lambda: np.ndarray(0));
    sumHRt_record = defaultdict(lambda: np.ndarray(0));
    sumtheta_record = compute_sumtheta_record(H_record,theta_record);
    for i in range(0,numN):
        Hsumtheta_record[list(H_record.keys())[i]] = np.dot(H_record[list(H_record.keys())[i]],sumtheta_record[list(H_record.keys())[i]]);
        sumHRt_record[list(H_record.keys())[i]] = compute_sumHRt(H_record,R_record,i);
    return Hsumtheta_record, sumHRt_record

def compute_sumtheta_record(H_record,theta_record):
    """
    Computa sumtheta para cada omica.

    Args
    -------
    H_record (dict)

    theta_record (dict)

    Returns
    -------
    sumtheta_record (dict) calcula para cada omica

    .. math:: 
        \sum_{t}\\theta_{I}+(\\theta_{I}^{(t)})^{T}
    """
    numN = len(H_record); sumtheta_record = defaultdict(lambda: np.ndarray(0));
    for i in range(0,numN):
        N = H_record[list(H_record.keys())[i]].shape[1]; sumtheta = np.zeros([N,N]);
        # print("    ** matriz inicial de Theta")
        # print("    -> shape sumtheta: "+str(sumtheta.shape))
        # print("    -> sum sumtheta: "+str(sum(sum(sumtheta))))
        for value in theta_record[list(H_record.keys())[i]]:
            sumtheta = sumtheta+value.T;
            # print("        ** sumando Theta")
            # print("        -> shape sumtheta: "+str(sumtheta.shape))
            # print("        -> sum sumtheta: "+str(sum(sum(sumtheta))))
        sumtheta_record[list(H_record.keys())[i]] = sumtheta;
    return sumtheta_record

def compute_sumHRt(H_record,R_record,index):
    """
    Computa sumHRt.

    Args
    -------
    H_record (dict)

    R_record (dict)

    index (int): indice de la omica a calcular

    Returns
    -------
    sumHRt (array) calcula

    .. math::
        \sum_{I\\neq J}H_{J}R_{IJ}^{T}
    """
    numN = len(H_record); H = H_record[list(H_record.keys())[index]]; K,N_index = H.shape;
    sumHRt = np.zeros([K,N_index]);
    for i in range(0,numN):
        if (str(i) != str(index)) & (R_record[index,i].shape[0] != 0):
            sumHRt = sumHRt + np.dot(H_record[list(H_record.keys())[i]],R_record[index,i].T);
    return sumHRt

def compute_actualize_H_record(H_record, Hsumtheta_record, sumHRt_record,AtKAH,Kn,A0, L1, L2, d1, K, miss_profile):
    """
    Actualiza cada uno de los Hi

    Args
    -------
    H_record (dict)

    Hsumtheta_record (dict)

    sumHRt_record (dict)
    
    AtKAH (dict)
    
    Kn (dict)
    
    A0 (dict)

    L1, L2, r2 (int)

    WtW (array)

    K (int) 

    Returns
    -------
    H_record (dict) 

    """
    factorH={}
    for i in H_record.keys():
        if miss_profile[0] in i:
            numerator=np.dot(A0['A0_ccle'][i].T,Kn['K_ccle'][i].T)+((L1/2)*Hsumtheta_record[i])+((L2/2)*sumHRt_record[i])
            denominator=AtKAH['AtKAH_ccle'][i]+(d1*np.dot(np.ones([K,K]),H_record[i]))+np.finfo(float).eps
        else:
            numerator=np.dot(A0['A0_ccle'][i].T,Kn['K_ccle'][i].T)+np.dot(A0['A0_tcga'][i].T,Kn['K_tcga'][i].T)+((L1/2)*Hsumtheta_record[i])+((L2/2)*sumHRt_record[i])
            denominator=AtKAH['AtKAH_ccle'][i]+AtKAH['AtKAH_tcga'][i]+(d1*np.dot(np.ones([K,K]),H_record[i]))+np.finfo(float).eps
        H_record[i] =  H_record[i]*(numerator / denominator);
        factorH[i]=(numerator / denominator)
        
    return H_record,factorH

def jNMF_module(H,t,path,featureLabel, print_savePath,is_best_H,method_clustering,enrichment,nameDotPlot,merged):
    from itertools import product
    from collections import defaultdict 
    comodule = defaultdict(lambda: defaultdict(dict));
    comodule_count = defaultdict(lambda: defaultdict(int));
    connectivity_matrix = defaultdict(lambda: defaultdict(0));
    
    for perfil, value in H.items():
        if value.shape[1] > 1000:
            pass
        else:
            coord_matrix = list(range(0,H[perfil].shape[1]))
            coord_matrix = product(coord_matrix, coord_matrix)
            connectivity_matrix[perfil] = {k: 0 for k in coord_matrix}
    
    if method_clustering == "weight": 
            ##quantiles para cada H
            H_q3 = {k: np.quantile(v, .75, axis=1) for k, v in H.items()};
            
            for perfil, value in H.items():
                if print_savePath == True:
                    try:
                        # Create target Directory
                        os.mkdir(path+'/'+perfil+'_cluster_weight')
                    except FileExistsError:
                        pass
                
                # limpiar carpeta para evitar incluir resultados de experimentos anteriores
                folder=path+"/"+perfil+'_cluster_weight'
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    os.unlink(file_path)
                    
                for i in range(0, value.shape[0]):
                    comodule[perfil].update({'co-md_'+str(i): list(np.where(value[i,:] > H_q3[perfil][i] ))});
                    comodule_count[perfil]['co-md_'+str(i)] = len(comodule[perfil]['co-md_'+str(i)][0])
                    
                    ## Connectivity matrix
                    if value.shape[1] > 1000:
                        pass
                    else:
                        matriz = product(list(comodule[perfil]['co-md_'+str(i)][0]), list(comodule[perfil]['co-md_'+str(i)][0]));
                        matriz = list(matriz);
                        for j in range(0,len(matriz)):
                            connectivity_matrix[perfil][matriz[j]] = 1
                    if print_savePath == True:
                        listaFeatures = [featureLabel[perfil][v] for v in comodule[perfil]['co-md_'+str(i)][0]]
                        with open(path+'/'+perfil+'/'+perfil+'_cluster_weight'+'/'+perfil+'_co-md_'+str(i)+'.txt', 'w') as f:
                            f.writelines("%s\n" % l for l in listaFeatures)
                            
    elif method_clustering == "H_first_second_max":
        ## maximos_first
        #para cada feature, cluster al que pertenece (con max)
        H_maximos = {perfil: np.argmax(v, axis=0) for perfil, v in H.items()};
        
        #valor de la pertenencia al cluster
        H_maximos_valor = {perfil: np.max(v, axis=0) for perfil, v in H.items()}
        
        # Outliers
        #clusters unicos para los features de la omica
        unicos = {perfil: np.unique(H_maximos[perfil]) for perfil in H_maximos.keys()}
        
        def unicos_H_max_Q2(unicos, H, H_maximos, H_maximos_valor, perfil):
            #para cada cluster i
            for i in unicos[perfil]:
                #valor de la pertenencia al cluster i, dist. de valor de la pertenancia al cluster i (solo con features asignados al cluster i )
                partial_H = H[perfil][i,np.where(H_maximos[perfil] == i)[0]]
                #cuartil 1 y 2 el cuartil es un tipo de cuantil. Hay diferentes tipos de cuantiles (dividen la distribucion en partes iguales): cuartiles, quintiles, deciles y percentiles
                Q1 = np.quantile(partial_H, .25)
                Q2 = np.quantile(partial_H, .50)
                
                # los valores maximos proximos a cero (debajo de Q1, debajo de Q2 en mrna)se omiten por lo que se fija su valor en 9999
                # al ser llamados mas abajo en la formacion de co-modulos, no se encontraran y no 
                # entraran al comodulo.Recuerda H_maximos->para cada feature, cluster al que pertenece (con max)
                if "mrna" in perfil:
                    H_maximos[perfil][np.where((H_maximos_valor[perfil] < Q2) & (H_maximos[perfil] == i))] = 9999
                else:
                    H_maximos[perfil][np.where((H_maximos_valor[perfil] < Q1) & (H_maximos[perfil] == i))] = 9999
                
            return H_maximos[perfil]
        
        for perfil in H_maximos.keys():
            H_maximos[perfil] = unicos_H_max_Q2(unicos, H, H_maximos, H_maximos_valor, perfil)
    
        ## maximos_second
        #para cada feature, cluster al que pertenece (con 2do max)
        H_maximos_second = {perfil: np.argsort(-v, axis=0)[1] for perfil, v in H.items()} #np.argsort: Returns the indices that would sort(ascendent)an array .-v descendent, toma segunda pos, SEGUNDO MAXIMO
        
        H_maximos_valor = {perfil: np.abs(np.sort(-v, axis=0)[1]) for perfil, v in H.items()}
        
        # Outliers
        unicos = {perfil: np.unique(H_maximos_second[perfil]) for perfil in H_maximos_second.keys()}
        
        def unicos_H_max_second_Q2(unicos, H, H_maximos_second, H_maximos_valor, perfil):
            for i in unicos[perfil]:
                #Con 2do maximo.valor de la pertenencia al cluster i, dist. de valor de la pertenancia al cluster i            
                partial_H = H[perfil][i,np.where(H_maximos_second[perfil] == i)[0]]
                Q1 = np.quantile(partial_H, .25)
                Q3 = np.quantile(partial_H, .75)
                IQR = Q3 - Q1
               
                # los valores maximos proximos a cero se omiten por lo que se fija su valor en 9999
                # al ser llamados mas abajo en la formacion de co-modulos, no se encontraran y no 
                # entraran al comodulo.Recuerda H_maximos_second->para cada feature, cluster al que pertenece (con 2do max)
                H_maximos_second[perfil][np.where((H_maximos_valor[perfil] < Q3 + 1.5 * IQR) & (H_maximos_second[perfil] == i))] = 9999
            
            return H_maximos_second[perfil]
        
        for perfil in H_maximos_second.keys():
            H_maximos_second[perfil] = unicos_H_max_second_Q2(unicos, H, H_maximos_second, H_maximos_valor, perfil)
        
        for perfil, value in H.items():
             if print_savePath == True:
                 try:
                     # Create target Directory
                     os.mkdir(path+'/'+perfil)
                 except FileExistsError:
                     pass
             
             # limpiar carpeta para evitar incluir resultados de experimentos anteriores
             folder=path+'/'+perfil
             for filename in os.listdir(folder):
                 file_path = os.path.join(folder, filename)
                 os.unlink(file_path)
             
             #recorre clusters.cluster i 
             for i in range(0, value.shape[0]):
                 #diccionario, para cada omica muestra los K comodulos con los features que pertenecen a dichos comodulos 
                 comodule[perfil].update({'co-md_'+str(i): list(np.where(i == H_maximos[perfil]))});
                 adicional_segundo = np.where(i == H_maximos_second[perfil]);
                 comodule[perfil]['co-md_'+str(i)][0] = np.append(comodule[perfil]['co-md_'+str(i)][0], adicional_segundo[0]);
                 #diccionario, para cada omica muestra los K comodulos con el # de features que pertenecen a dichos comodulos 
                 comodule_count[perfil]['co-md_'+str(i)] = len(comodule[perfil]['co-md_'+str(i)][0]);
                 
                 ## Connectivity matrix.
                 #product->cartesian product from the given iterator. Parejas de features en el mismo comodulo (incluye i=j) (38,38),(38,50),.. etc
                 if value.shape[1] > 1000:
                     pass
                 else:
                     matriz = product(list(comodule[perfil]['co-md_'+str(i)][0]), list(comodule[perfil]['co-md_'+str(i)][0]));
                     matriz = list(matriz);
                     for j in range(0,len(matriz)):
                         connectivity_matrix[perfil][matriz[j]] = 1
                 if print_savePath == True:
                     listaFeatures = [featureLabel[perfil][v] for v in comodule[perfil]['co-md_'+str(i)][0]]
                     with open(path+'/'+perfil+'/'+perfil+'_co-md_'+str(i)+'.txt', 'w') as f:
                         f.writelines("%s\n" % l for l in listaFeatures)
                         
        if is_best_H:
            rango = H[list(H.keys())[0]].shape[0];
            if enrichment:
                # command = 'C:/Program Files/R/R-3.6.2/bin/Rscript.exe'
                # command = 'C:/Users/da.salazarb/Documents/R/R-3.6.3/bin/Rscript.exe'
                #path2script = '--vanilla 00_gene-gene_BioGRID.R'
                
                if os.path.isfile(path+"/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv") == False:
                    iteration_number=str(1)
                else:
                    iteration_number = pd.read_csv(path+"/co-mod_tabulated_results/Tabulated_results_MM-KjNMF_for_"+str(merged.listaPerfiles)+".csv")
                    iteration_number = str(iteration_number.shape[0])
                    
                command='C:/Program Files/R/R-4.1.1/bin/Rscript.exe'
                                
                scriptName = "13_co-module_interpretation.R"
                path2script = path + "/co-mod_R_scripts/" + scriptName
                
                cmd = [command, path2script] + [path+"__"+nameDotPlot+"__"+iteration_number] #cmd = [command, path2script] + args. Ver pagina rPubs.
                x = subprocess.Popen(cmd).wait()
                if x == 0:
                    print("    * All fine with " + scriptName + '\n')
                else:
                    print("    * Something is wrong with " + scriptName + '\n')
                    
                record_mrna = pd.read_csv(path+"/co-mod_records/record_mrna.csv");
                if record_mrna.shape[0] == 0:
                    BioScore=0
                    codEnrich_over_codTotal=0; T4=0; saraPenalization=0; geneRatio_avg=0; geneRatio_sd=0;
                    cod_stats2 = [BioScore, codEnrich_over_codTotal, T4, saraPenalization, geneRatio_avg, geneRatio_sd];
                else:
                    codEnrich_over_codTotal = len(np.unique(record_mrna.Cluster)) / rango; # e/K
                    # num_total_ptos_over_totalEnrich_totalTerms
                    P = record_mrna.shape[0]; e = len(np.unique(record_mrna.Cluster)); terms = len(np.unique(record_mrna.ID));
                    T4 = (P - e) / ((e * (terms-1))+np.finfo(float).eps)
                    T4 = [1 if T4 > 0.5 else 0];
                    # R2 = (record_mrna.shape[0] - len(np.unique(record_mrna.Cluster))) / (len(np.unique(record_mrna.ID)) * len(np.unique(record_mrna.Cluster)));
                    saraPenalization = P / 1000; #Pref=1000
                    geneRatio_avg = np.mean([eval(t) for t in record_mrna.GeneRatio]);
                    geneRatio_sd = np.std([eval(t) for t in record_mrna.GeneRatio]);
                    BioScore = 0.4*codEnrich_over_codTotal + 0.2*geneRatio_avg + 0.4*terms / P - T4[0]; #metrica GDS
                    # num_total_ptos_over_totalEnrich = record_mrna.shape[0] / len(np.unique(record_mrna.ID)); #sum(record_mrna.groupby("Cluster")['ID'].count()); #.to_dict();
                    # record_mrna_ = codEnrich_over_codTotal - num_total_ptos_over_totalEnrich;
                    geneRatio_avg = np.mean([eval(t) for t in record_mrna.GeneRatio]);
                    geneRatio_sd = np.std([eval(t) for t in record_mrna.GeneRatio]);
                    cod_stats2 = [BioScore, codEnrich_over_codTotal, T4, saraPenalization, geneRatio_avg, geneRatio_sd];
            else:
                BioScore=0
                codEnrich_over_codTotal=0; T4=0; saraPenalization=0; geneRatio_avg=0; geneRatio_sd=0;
                cod_stats2 = [BioScore, codEnrich_over_codTotal, T4, saraPenalization, geneRatio_avg, geneRatio_sd];
            
            # numero de comodulos no vacios / K(comodulos totales)     -> % comodulos no vacios
            cod_no_vacios = {k: sum(1 for j in v.values() if j > 0) / rango for k,v in comodule_count.items()};
            # suma variables incluidas en comodulos no vacios / variables totales ->%variables incluidas en los comodulos
            cod_no_vacios_suma = {k: sum(j for j in v.values() if j > 0) / H[k].shape[1] for k,v in comodule_count.items()};
            # min, max y promedio de las variables incluidas de todos los comodulos (comodule_count)  de la omica
            # min, max y promedio de comodule_count de TODOS los comodulos de la omica
            cod_stats = {k: [j for j in v.values() if j > 0] for k,v in comodule_count.items()};
            cod_stats = {k: [k, min(v), max(v), np.mean(v), cod_no_vacios[k], cod_no_vacios_suma[k]] for k,v in cod_stats.items()};
        else: # si es otra matriz H que no sea best_H
            cod_stats = {};
            cod_stats2 = [];
    return comodule, comodule_count, connectivity_matrix, cod_stats, cod_stats2
## ------------------------------------------------------------------------------------------------------------------------------------
    
def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)

    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / (np.sqrt(np.dot(ssA[:, None],ssB[None])) + np.finfo(float).eps)

def jNMF_module_patients_cellLines(best_W, infoFeatures_tcga, infoFeatures_ccle, param, path, sample):
    module_patient_cell = defaultdict(list)
    ccle_dict = defaultdict(int)
    ccle_dict = {k: 0 for k in infoFeatures_ccle}
    tcga_dict = defaultdict(int)
    tcga_dict = {k: 0 for k in infoFeatures_tcga}
    for i in range(0,best_W["W0_ccle"].shape[1]):
        Q3 = np.quantile(best_W["W0_ccle"][:,i], .75)
        ccle_ = list(np.array(infoFeatures_ccle)[np.where(best_W["W0_ccle"][:,i] > Q3),][0])
        for v in ccle_:
            ccle_dict[v] += 1
        
        Q3 = np.quantile(best_W["W0_tcga"][:,i], .75)
        tcga_ = list(np.array(infoFeatures_tcga)[np.where(best_W["W0_tcga"][:,i] > Q3),][0])
        for v in tcga_:
            tcga_dict[v] += 1
        
        module_patient_cell['co-md_'+str(i)] = ccle_ + tcga_
        
    ccle_ = pd.DataFrame.from_dict(ccle_dict, orient='index').sort_values(by=[0], ascending=False).T
    tcga_ = pd.DataFrame.from_dict(tcga_dict, orient='index').sort_values(by=[0], ascending=False).iloc[0:best_W["W0_ccle"].shape[0]].T
    tcga_.columns = list(map(lambda x: x.replace("TCGA-", ""), tcga_.columns))
    
    if sample != "sample":        
        porc_modulos_no_vacios_ccle_=np.sum(ccle_.values != 0) / best_W["W0_ccle"].shape[0];
    else:
        porc_modulos_no_vacios_ccle_=0
    return module_patient_cell, porc_modulos_no_vacios_ccle_

def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)


def silueta(K, sample_silhouette_values, cluster_labels, ax):
    # import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    y_lower = 10
    for i in range(K):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = \
            sample_silhouette_values[cluster_labels == i]

        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i) / K)
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples
        
    return ax

def func_ccle_groups(best_W, mrna_rows_ccle, pac_por_comod_ccle, path, projects):
    rec_surv = pd.DataFrame(data=np.zeros([best_W["W0_ccle"].shape[0], best_W["W0_ccle"].shape[1]]), index=list(map(lambda x: x, mrna_rows_ccle)))
    for i in range(0,best_W["W0_ccle"].shape[0]):
        rec_surv.iloc[i,np.argsort(-best_W["W0_ccle"][i,:], axis=0)[0]] = 1 ## Selecciona la fila de paciente y busca donde es el maximo
        
    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx
    
    cuantos_co_md = pac_por_comod_ccle ## numero de pacientes permitidos en cada cluster, orignalmente 30
    W_max = np.max(best_W["W0_ccle"], axis=0) # los maximos por K clusters
    W_max[np.where(np.sum(rec_surv)<cuantos_co_md)[0]] = -9999 ## cuales de ellos son menores a cuantos_co_md
    b = mrna_rows_ccle
    menores_Xpac = np.where(np.sum(rec_surv)<cuantos_co_md)[0]
    for i in range(0,menores_Xpac.shape[0]): ## recorrer todos los co-md que tienen menos de X pacientes
        a = rec_surv[menores_Xpac[i]][np.where(rec_surv[menores_Xpac[i]]>0)[0]].index
        if len(a) > 0:
            for j in range(0,len(a)):
                d = np.where(rec_surv.loc[a[j]]>0)[0][0] ## buscar posicion actual de a[j]
                rec_surv.loc[a[j], d] = 0 ## quitar esa posicion
                c = np.max(best_W["W0_ccle"][b.index(a[j]),]) ## maximo en best_W["W0_tcga"] para esa fila
                # print(W_max[find_nearest(W_max, c)])
                rec_surv.loc[a[j],find_nearest(W_max, c)] = 1 ## nueva posicion en co-md
    np.where(np.sum(rec_surv,axis=0)>0)
    
    rec_surv = pd.melt(rec_surv.reset_index(), id_vars="index")
    rec_surv  = rec_surv.loc[rec_surv.value != 0,:]
    rec_surv.rename(columns={'variable': 'cluster'}, inplace=True)
    rec_surv.drop(["value"], axis=1, inplace=True)
    rec_surv.to_csv(path+'/co-mod_observations_clusters/observations_clusters_for_project_'+projects[0]+'_using_best_W.csv')
    
def hp_metrics_plots(sample, path, central_path, command, merged):
    # ---------------------------------------------------------------------------------------- #
    # ### ........... Biological data selection of hyperparameters: metric plots ........... ####
    # ---------------------------------------------------------------------------------------- #
    if sample != "sample":
        print("")
        print('----------------------------------------------------------')
        print("* Generating Bio data metric plots: enrich_over_K, geneRatio_avg, puntosTotales, No.CurvasSurivial_p, porc_Cluster_CCLE, FuncObj")
        path2script = '--vanilla ' + path + "/25_result_selectionHyperParams.R"
        cmd = [command, path2script] + ['empty'+"__"+merged.path] #cmd = [command, path2script] + args. Ver pagina rPubs.K-15-r1_10-r2_10-L1_10-L2_10-d1_10
        x = subprocess.Popen(cmd).wait()
        if x == 0:
            print("Ok, metric plots is saved in co-mod_suplementarios!")
        else:
            print("An error ocurred in 25_result_selectionHyperParams.R script")
        print('---------------------------------------------------------- \n')
        
    # ---------------------------------------------------------------------------------------- #
    # ### ........... Synthetic data selection of hyperparameters: metric plots ........... ####
    # ---------------------------------------------------------------------------------------- #
    if sample == "sample":
        print("")
        print('----------------------------------------------------------')
        print("* Generating Synthetic data metric plots: rss, r2, No. asoc. detected / No. asoc. real, perc_true_identified")
        path2script = '--vanilla ' + path + "/26_result_SimulData_r2_adjusted.R"
        cmd = [command, path2script] + ['empty'+"__"+path] #cmd = [command, path2script] + args. Ver pagina rPubs.K-15-r1_10-r2_10-L1_10-L2_10-d1_10
        x = subprocess.Popen(cmd).wait()
        if x == 0:
            print("Ok, metric plots is saved in co-mod_suplementarios!")
        else:
            print("An error ocurred in 26_result_SimulData_r2_adjusted.R script")
        print('---------------------------------------------------------- \n')
        
    # ---------------------------------------------------------------------------------------- #
    # ### ........... Synthetic data selection of hyperparameters: metric plots ........... ####
    # ---------------------------------------------------------------------------------------- #
    if sample == "sample":
        print("")
        print('----------------------------------------------------------')
        print("* Generating Synthetic data metric plots: average metrics")
        path2script = '--vanilla ' + path + "/28_result_SimulData_r2_adjusted_secondTry.R"
        cmd = [command, path2script] + ['empty'+"__"+path] #cmd = [command, path2script] + args. Ver pagina rPubs.K-15-r1_10-r2_10-L1_10-L2_10-d1_10
        x = subprocess.Popen(cmd).wait()
        if x == 0:
            print("Ok, metric plots is saved in co-mod_suplementarios!")
        else:
            print("An error ocurred in 28_result_SimulData_r2_adjusted_secondTry.R script")
        print('---------------------------------------------------------- \n')
        
    # ---------------------------------------------------------------------------------------- #
    # ### ........... Synthetic data selection of hyperparameters: impact hyperparams plots ........... ####
    # ---------------------------------------------------------------------------------------- #
    if sample == "sample":
        print("")
        print('----------------------------------------------------------')
        print("* Generating Synthetic data metric plots: average metrics")
        path2script = '--vanilla ' + path + "/29_result_SimulData_ImpactHyperParams.R"
        cmd = [command, path2script] + ['empty'+"__"+path] #cmd = [command, path2script] + args. Ver pagina rPubs.K-15-r1_10-r2_10-L1_10-L2_10-d1_10
        x = subprocess.Popen(cmd).wait()
        if x == 0:
            print("Ok, hyperparams impact plot is saved in co-mod_suplementarios!")
        else:
            print("An error ocurred in 29_result_SimulData_ImpactHyperParams.R script")
        print('---------------------------------------------------------- \n')  
    # np.where(np.sum(rec_surv,axis=0)>0)

def testSection(index,A_learned,H_learned,ccle,ccle_test,tcga,tcga_test,R_record,theta_record,param, miss_profile,newNamesProfiles,projects):

    K=param['K']; r1=param['r1']; L1=param['L1']; L2=param['L2']; r2=param['r2']; d1=param['d1']; o1=param['o1']; o2=param['o2']
    tol=param['tol']; maxiter=param['maxiter']; sigma_ccle=param['sigma_ccle']; sigma_tcga=param['sigma_tcga'];sigma_ccle_diff=param['sigma_ccle_diff']; sigma_tcga_diff=param['sigma_tcga_diff']; stop_control = np.zeros([1,maxiter+1]); stop_control[0,0] =1; #stop control de las iteraciones
    delta_control = np.zeros([1,maxiter+1]); delta_control[0,0] =1;
        
    # print('        --------------------------------------------------------')
    print('        ** Fix learned A_I and find: H_I and A_I on testing dataset')
    print('')
    # KERNEL ARRAYS
    Kn=defaultdict(lambda: np.ndarray(0));
    Kn['K_ccle']={k: pairwise.rbf_kernel(v.T, Y=None, gamma=1/(2*(sigma_ccle**2))) for k,v in ccle_test.items()};
    Kn['K_tcga']={k: pairwise.rbf_kernel(v.T, Y=None, gamma=1/(2*(sigma_tcga**2))) for k,v in tcga_test.items()};

    
    KIJ_ccle=[]
    KIJ_tcga=[]
    for i in ccle.keys():
        otraMatriz_ccle=[]
        otraMatriz_tcga=[]
        for j in ccle.keys():
            matriz_ccle=pairwise.rbf_kernel(ccle_test[i].T, Y=ccle_test[j].T, gamma=1/(2*(sigma_ccle_diff**2)))
            otraMatriz_ccle.append(matriz_ccle)
            if i!=newNamesProfiles[miss_profile[0]] and j!=newNamesProfiles[miss_profile[0]]:
                matriz_tcga=pairwise.rbf_kernel(tcga_test[i].T, Y=tcga_test[j].T, gamma=1/(2*(sigma_tcga_diff**2)))
                otraMatriz_tcga.append(matriz_tcga)
            else: #Para TCGA NO hay DRUG, lleno el arreglo con basura (ceros) para conservar orden en codigos 0:0_drug, 1:1_mirna .. que seran usados posteriormente para leer esta matriz!
                matriz_tcga=np.zeros([ccle_test[i].shape[1],ccle_test[j].shape[1]])
                otraMatriz_tcga.append(matriz_tcga)               
        KIJ_ccle.append(otraMatriz_ccle)
        KIJ_tcga.append(otraMatriz_tcga)
    
    #print('Kernel Arrays loaded')       
    #print('* Initializing A and H matrices')
    
    #Fix learned A
    Ainit =A_learned.copy()
    H_record_init = defaultdict(lambda: np.ndarray(0));
    H_record_init = {k: np.random.uniform(0,1,[K, v.shape[1]]) for k, v in ccle.items()};
    
    ## normalizacion de A y H
    A, H,sh = normalize_AH_rowH(Ainit, H_record_init,1);
        
    #print('Listo OK ...')
    #print('---------------------------------------------------------- \n')
    
    
    #%% Computo del error con los parametros iniciales ('MUR')
    #print('----------------------------------------------------------')
    #print('* Computing the error with the initial parameters')
    
    AtKA= defaultdict(lambda: np.ndarray(0));
    AtKA['AtKA_ccle']={k: np.dot(np.dot(A['A0_ccle'][k].T,Kn['K_ccle'][k]),A['A0_ccle'][k]) for k in ccle.keys()}
    AtKA['AtKA_tcga']={k: np.dot(np.dot(A['A0_tcga'][k].T,Kn['K_tcga'][k]),A['A0_tcga'][k]) for k in tcga.keys()}        
    Hsumtheta_record, sumHRt_record = compute_HsumTheta_sumHRt_record(H, theta_record, R_record);
    #matrices omicas, una al lado de la otra en un solo numpy (tcga NO tiene drug)
    #terminos de la FO 
    diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga = compute_diff(H,Hsumtheta_record,sumHRt_record,AtKA,KIJ_ccle,KIJ_tcga,A,L1,L2,r1,r2,d1,o1,o2,K,miss_profile);
    delta1_init,diff1_ccle, diff1_tcga = StopCriterion_rule1(Kn,A,H,AtKA,diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga,miss_profile);
    #print('Listo OK ...')
    #print('---------------------------------------------------------- \n')
    #%% KJNMF algorithm
    # print('                    KJNMF algorithm')
    # print('----------------------------------------------------------')
    iteracion = 0;
    A0=defaultdict(lambda: np.ndarray(0));
    sumKA_ccle = defaultdict(lambda: np.ndarray(0));
    sumKA_tcga = defaultdict(lambda: np.ndarray(0));
    KA=defaultdict(lambda: np.ndarray(0));
    AtKAH=defaultdict(lambda: np.ndarray(0));
    #N
    numN = len(H)
    # print('* Running the Multiplicative Update Rule')
    #Mientras que no cumpla tolerancia minima o haya alcanzado el numero maximo de iteraciones 
    while stop_control[0, iteracion] > tol and iteracion < maxiter:
    
        A0=A.copy()
        
        #% update Hi with A fixed
        AtKAH['AtKAH_ccle']={k: np.dot(np.dot(A0['A0_ccle'][k].T,Kn['K_ccle'][k]),np.dot(A0['A0_ccle'][k],H[k])) for k in ccle.keys()}
        AtKAH['AtKAH_tcga']={k: np.dot(np.dot(A0['A0_tcga'][k].T,Kn['K_tcga'][k]),np.dot(A0['A0_tcga'][k],H[k])) for k in tcga.keys()}
        Hsumtheta_record, sumHRt_record = compute_HsumTheta_sumHRt_record(H, theta_record, R_record);
        #H y H0 son variables ligadas!!!
        H0,factorH=compute_actualize_H_record(H, Hsumtheta_record, sumHRt_record,AtKAH,Kn,A0, L1, L2, d1, K, miss_profile)
        
        #% normalization. H y H0 ya NO estan ligadas!!!            
        A,H,sh =normalize_AH_rowH(A0,H0,1);
    
        #stop criterion
        AtKA['AtKA_ccle']={k: np.dot(np.dot(A['A0_ccle'][k].T,Kn['K_ccle'][k]),A['A0_ccle'][k]) for k in ccle.keys()}
        AtKA['AtKA_tcga']={k: np.dot(np.dot(A['A0_tcga'][k].T,Kn['K_tcga'][k]),A['A0_tcga'][k]) for k in tcga.keys()}        
        Hsumtheta_record, sumHRt_record = compute_HsumTheta_sumHRt_record(H, theta_record, R_record);
        #matrices omicas, una al lado de la otra en un solo numpy (tcga NO tiene drug)
        #terminos de la FO 
        diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga = compute_diff(H,Hsumtheta_record,sumHRt_record,AtKA,KIJ_ccle,KIJ_tcga,A,L1,L2,r1,r2,d1,o1,o2,K,miss_profile);
        delta,diff1_ccle, diff1_tcga = StopCriterion_rule1(Kn,A,H,AtKA,diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga,miss_profile);
    
        if iteracion == 0:
            delta_old = delta1_init;
        iteracion += 1;
    			#FO
        delta_control[0, iteracion] = np.around(delta, decimals=5);
        #stop control de las iteraciones
        #stop control es DECRECIENTE.MUR garantiza FO decrecreciente
        stop_control[0, iteracion] = np.around(np.abs((delta_old - delta)/(delta1_init - delta)), decimals=12);
        delta_old  = delta;
        # if sample != "sample":
        #     if iteracion % 10 == 0:
        #         print('We are in iteration {: d} - stop_control: {}'.format(iteracion, stop_control[0, iteracion]))
        #         #print('We are in iteration {: d} - delta_control: {}'.format(iteracion, delta_control[0, iteracion]))
        # else:
        #     if iteracion % 100 == 0:
        #         print('We are in iteration {: d} - stop_control: {}'.format(iteracion, stop_control[0, iteracion]))
                
    print('        + Repeat {:d} on testing data. Optimum value found in iteration: {:d} -> stop_control: {}'.format(index, iteracion, stop_control[0, iteracion]))        
    print('            -> diff1_{}: {} - diff1_{}: {} - diff2: {} - diff3: {}'.format(projects[0], np.round(diff1_ccle,4), projects[1], np.round(diff1_tcga,4), np.round(diff2,4), np.round(diff3,4)))
    print('            -> diff4_{}: {} - diff4_{}: {} - diff5: {}'.format(projects[0], np.round(diff4_ccle,4), projects[1], np.round(diff4_tcga,4), np.round(diff5, 4)))
    print('            -> diff6_{}: {} - diff6_{}: {} '.format(projects[0], np.round(diff6_ccle,4), projects[1], np.round(diff6_tcga,4)))
    print('            -> Objective Function value: {} \n'.format(delta))
    print('        --------------------------------------------------------')
    print('')
    
    # print('----------------------------------------------------------')    
    # print('Compare H_I and predH_I')
    # print('----------------------------------------------------------') 
    #best_H-> pred_HI. H_record_simulated->H_learned
    pred_HI=H.copy()
    
    # Hungarian algorithm: http://software.clapper.org/munkres/#introduction
    B = {k: 1-corr2_coeff(v, pred_HI[k]) for k,v in H_learned.items()}; #bipartite graph G dim k x k (Gij 1-PCC between the i-th row of H0 and the j-th row of H)
    m = Munkres()
    B = {k: m.compute(v) for k,v in B.items()}  #tupla optima (filaH0,filaH) segun algoritmo hungaro (por defecto busca min)
    
    B = {k: [x[1] for x in v] for k,v in B.items()}; #segundo elemento de la tupla (filaH)
    
    pred_HI = {k: v[B[k],:] for k,v in pred_HI.items()}; #H ordered according H0 
    
    H_learned_binary={k: ((v > np.quantile(v,.75)) * 1).flatten() for k,v in H_learned.items()}
    #Para cada omica 50 vectores binarios. compara con 50 cuantiles diferentes desde el 0/50 hasta el 50/50 
    pred_HI_binary={k: [(v.flatten()>np.quantile(v,c/50))*1 for c in range(50)] for k,v in pred_HI.items()}
    #se calculan 50 puntos tpr_fpr. se compara el vector real con cada uno de los 50 estimados 
    tpr_fpr={k: [roc_curve_MjNMF(H_learned_binary[k],i) for i in v] for k,v in pred_HI_binary.items()}
    #ordena de menor a mayor (ambos valores en la tupla siguen el mismo patron de ascendencia)
    tpr_fpr={k: sorted(v) for k,v in tpr_fpr.items()}
    tpr={k: np.array([i[0] for i in v]) for k,v in tpr_fpr.items()}
    fpr={k: np.array([i[1] for i in v]) for k,v in tpr_fpr.items()}
    AUC_H={k: [metrics.auc(fpr[k],tpr[k])] for k in ccle.keys()}
    
    return AUC_H, stop_control, delta_control

## ---------------------------------------------------------------------------------------------------------------------------------------

def LoadSyntheticData(ruta, dim_features):
    ### Cargar datos
    try:
        W_simulated = dict()
        W_simulated["W0_ccle"] = np.loadtxt(ruta+'/co-mod_simulated_H_W/simulated_W0_ccle.csv', delimiter=',')#+(noise_proportion*np.random.normal(0,1,[dim_samples[0],K]))
        W_simulated["W0_tcga"] = np.loadtxt(ruta+'/co-mod_simulated_H_W/simulated_W0_tcga.csv', delimiter=',')#+(noise_proportion*np.random.normal(0,1,[dim_samples[1],K]))
        
        H_record_simulated = dict()
        H_record_simulated = {k: np.loadtxt(ruta+'/co-mod_simulated_H_W/simulated_H_'+k+'.csv', delimiter=',') for k in dim_features.keys()}

        noise_data_0 = dict()
        noise_data_0 = {k: np.loadtxt(ruta+'/co-mod_simulated_H_W/noise_data_0_'+k+'.csv', delimiter=',') for k in dim_features.keys()}
        noise_data_1 = dict()
        noise_data_1 = {k: np.loadtxt(ruta+'/co-mod_simulated_H_W/noise_data_1_'+k+'.csv', delimiter=',') for k in dim_features.keys()}
        
        # W - asociaciones pacientes
        max_W_simulated = {k: np.argmax(v, axis=1) for k,v in W_simulated.items()}
        max_columna_W = {k: pd.DataFrame([range(len(v)),v], index=['indice','grupo']).T for k,v in max_W_simulated.items()};
        max_columna_W_ = dict()
        for k, v in max_columna_W.items():
            W_asociaciones_por_paciente=[]
            for j in np.unique(v.grupo):
                W_asociaciones_por_paciente = W_asociaciones_por_paciente  + list(combinations(v.indice[v.grupo == j],2))
            max_columna_W_[k] = W_asociaciones_por_paciente
            
        # H - asociaciones moleculas
        max_H_record_simulated = {k: np.argmax(v, axis=0) for k,v in H_record_simulated.items()}
        max_columna_H = {k: pd.DataFrame([range(len(v)),v], index=['indice','grupo']).T for k,v in max_H_record_simulated.items()};
        max_columna_H_ = dict()
        for k, v in max_columna_H.items():
            H_asociaciones_por_moleculas=[]
            for j in np.unique(v.grupo):
                H_asociaciones_por_moleculas = H_asociaciones_por_moleculas  + list(combinations(v.indice[v.grupo == j],2))
            max_columna_H_[k] = H_asociaciones_por_moleculas
            
        ## Theta constraints using H_I
        mayor_de_10 = {k: np.argmax(v>=10, axis=0) for k,v in H_record_simulated.items()}
        theta_permutaciones = dict()
        for k, v in mayor_de_10.items():
            lista_permutaciones =[]
            for j in range(v.shape[0]):
                lista_permutaciones = lista_permutaciones  + list(permutations(np.where(v==j)[0],2))
            theta_permutaciones[k] = lista_permutaciones
            
        # matriz theta
        theta_record = {k: np.zeros([v.shape[1], v.shape[1]]) for k,v in H_record_simulated.items()}
        for k,v in theta_record.items():
            for i in theta_permutaciones[k]:
                v[i] = 1;
        
        ## R constraints using H_I
        combinacion_perfiles = list(combinations(list(H_record_simulated.keys()), 2));
        comb_perfiles_dict = dict()
        for i in combinacion_perfiles:
            comb_perfiles_dict[i] = np.zeros([H_record_simulated[i[0]].shape[1], H_record_simulated[i[1]].shape[1]])
        compilado_dict=dict()
        for k,v in comb_perfiles_dict.items():
            a = mayor_de_10[k[0]]
            b = mayor_de_10[k[1]]
            c_range = np.intersect1d(a,b)
            d_list = []
            for i in c_range:
                d_list = d_list + list(product(list(np.where(mayor_de_10[k[0]]==i)[0]), list(np.where(mayor_de_10[k[1]]==i)[0])))
            compilado_dict[k] = d_list
            
        matrizFinal = []
        for k in H_record_simulated.keys():
            otraMatriz = []
            for k_ in H_record_simulated.keys():
                if (k,k_) in compilado_dict.keys():
                    for j in compilado_dict[(k,k_)]:
                        comb_perfiles_dict[(k,k_)][j] = 1
                    matriz = comb_perfiles_dict[(k,k_)].T
                elif k==k_:
                    matriz = np.zeros([H_record_simulated[k].shape[1], H_record_simulated[k_].shape[1]])
                else:
                    for j in compilado_dict[(k_,k)]:
                        comb_perfiles_dict[(k_,k)][j] = 1
                    matriz = comb_perfiles_dict[(k_,k)]
                otraMatriz.append(matriz);
            matrizFinal.append(np.array(otraMatriz, dtype=object)) ## ir guardando los arrays en un lista
            R_record = np.stack(matrizFinal, axis=-1)
            
        SynData = {'W_simulated': W_simulated, 'H_record_simulated': H_record_simulated, 'theta_record': theta_record, 
                   'R_record': R_record, 'max_columna_H_': max_columna_H_, 'max_columna_W_': max_columna_W_, 
                   'noise_data_0': noise_data_0, 'noise_data_1': noise_data_1};
        return SynData
        
    except IOError:
        print("You must first create simulated_W0_ccle, simulated_W0_tcga and the simulated_H_I")

def new_syntheticData(dim_features, dim_samples, numero_grupos_W, numero_grupos_H, ruta, K):
    ## descomentar para generar nuevos datos!
    ## W
    add_W_ccle = np.zeros([dim_samples[0], K])
    add_W_tcga = np.zeros([dim_samples[1], K])
    groups_W  = np.random.randint(numero_grupos_W, size=dim_samples[0]); ## se asignan las muestras a numero_grupos_W (ej: 10 grupos)
    where_col_max_W_ccle = np.random.choice(K, size=numero_grupos_W, replace=False); ## se escojen 10 columnas que tendran los maximos
    where_col_max_W_tcga = np.random.choice(K, size=numero_grupos_W, replace=False); ## se escojen 10 columnas que tendran los maximos
    
    for j in np.unique(groups_W):
        # for i in range(add_W_ccle.shape[1]):
        add_W_ccle[np.where(groups_W == j)[0], where_col_max_W_ccle[j]] = 10 ## se asignan los maximos
        add_W_tcga[np.where(groups_W == j)[0], where_col_max_W_tcga[j]] = 10 ## se asignan los maximos
    
    W0_ccle = np.random.uniform(0,1,[dim_samples[0],K])+add_W_ccle;
    W0_tcga = np.random.uniform(0,1,[dim_samples[1],K])+add_W_tcga;
    
    np.savetxt(ruta+'/co-mod_simulated_H_W/simulated_W0_ccle.csv', W0_ccle, delimiter=',')
    np.savetxt(ruta+'/co-mod_simulated_H_W/simulated_W0_tcga.csv', W0_tcga, delimiter=',')
    
    ## H_I
    add_H_I = {k: np.zeros([K,v]) for k, v in dim_features.items()}
    groups_H_I = {k: np.random.randint(numero_grupos_H, size=v) for k,v in dim_features.items()}
    where_row_max_H = {k: np.random.choice(K, size=numero_grupos_H, replace=False) for k in dim_features.keys()}; ## se escojen 35 columnas que tendran los maximos
    for k, v in add_H_I.items():
        for j in np.unique(groups_H_I[k]):
            # for i in range(v.shape[0]):
            v[where_row_max_H[k][j], np.where(groups_H_I[k] == j)[0]] = 10
    
    H_record_simulated = defaultdict(lambda: np.ndarray(0));
    # H_record_simulated = {k: np.random.uniform(0,1,[K, v])+add_H_I[k] for k, v in dim_features.items()};
    H_record_simulated = {k: np.random.uniform(0,1,[K, v])+add_H_I[k] for k, v in dim_features.items()};
    
    for k, v in H_record_simulated.items():
        np.savetxt(ruta+"/co-mod_simulated_H_W/simulated_H_"+k+'.csv', v, delimiter=',')
    
    ## noise_simulation
    noise_data_0 = {k: np.random.normal(0,5e-3,[dim_samples[0],v]) for k, v in dim_features.items()}
    noise_data_1 = {k: np.random.normal(0,5e-1,[dim_samples[0],v]) for k, v in dim_features.items()}
    
    for k, v in noise_data_0.items():
        np.savetxt(ruta+"/co-mod_simulated_H_W/noise_data_0_"+k+'.csv', v, delimiter=',')
    
    for k, v in noise_data_1.items():
        np.savetxt(ruta+"/co-mod_simulated_H_W/noise_data_1_"+k+'.csv', v, delimiter=',')
        
### Cargar datos
def explore(X):
    for k,v in X.items():
        print(k)
        print(v.shape)
        print(v[0:3,0:3])