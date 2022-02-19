# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 21:07:16 2019

@author: da.salazarb
"""
import os
import os.path
import numpy as np
import pandas as pd
from func_KMMJNMF import *
from collections import defaultdict
from sklearn.metrics import pairwise

#%% Funcion MMJNMF
def kmmjnmf(merged,param):
    """
    Returns W and H.

    Returns
    -------
    A (dict(array))
    H (dict(array))

    """
#%% Load data and parameters
    np.random.seed(0)
    
    # --------------------------------------------- #
    #### ........... Load parameters ........... ####
    # --------------------------------------------- #
    ## SARA
    K=param['K']; r1=param['r1']; L1=param['L1']; L2=param['L2']; r2=param['r2']; d1=param['d1']; o1=param['o1']; o2=param['o2']
    tol=param['tol']; maxiter=param['maxiter'];
    repeat=param['repeat']; path = merged.path+"/pathFeatureLabel/"; enrichment=param['enrichment'];
    sigma_ccle=param['sigma_ccle']; sigma_tcga=param['sigma_tcga'];
    sigma_ccle_diff=param['sigma_ccle_diff']; sigma_tcga_diff=param['sigma_tcga_diff'];
    index_train_tcga=param['index_train_tcga']; index_test_tcga=param['index_test_tcga'];
    index_train_ccle=param['index_train_ccle']; index_test_ccle=param['index_test_ccle'];
    nameDotPlot="K_{}_r1_{}_r2_{}_L1_{}_L2_{}_d1_{}_o1_{}_o2_{}_sigma_ccle_{}_sigma_tcga_{}_sigma_ccle_diff_{}_sigma_tcga_diff_{}".format(param['K'], param['r1'], param['r2'], param['L1'], param['L2'], param['d1'], param['o1'], param['o2'], param['sigma_ccle'], param['sigma_tcga'], param['sigma_ccle_diff'],param['sigma_tcga_diff'])
    ## END SARA
    
    ## *** Real-world data ***
    ## Verify if there are 2 projects
    project_profile = os.listdir(merged.path+"INPUT DATA");
    
    project1 = sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0]; project2 = sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1];
    profiles_project1 = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if project1 in p]));
    profiles_project2 = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if project2 in p]));
    
    if len(profiles_project1) > len(profiles_project2):
        ccle_tcga = "high"
    elif len(profiles_project1) < len(profiles_project2):
        ccle_tcga = "less"
    else:
        ccle_tcga = "equal"
    
    if ccle_tcga == "high":
        projects = [sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0], sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1]];
    elif ccle_tcga == "less":
        projects = [sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1], sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0]];
    elif ccle_tcga == "equal":
        projects = sorted(list(set([p[0:p.find("_")] for p in project_profile])));
        
    print("    * The projects are {} and {} ".format(projects[1], projects[0]))
    profiles_project_TCGA = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if projects[1] in p]));
    profiles_project_CCLE = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if projects[0] in p]));
    print("    ** For {}, the profiles are {}".format(projects[1], profiles_project_TCGA))
    print("    ** For {}, the profiles are {}".format(projects[0], profiles_project_CCLE))
    
    if len(profiles_project_CCLE) > len(profiles_project_TCGA):
        miss_profile = list(set(profiles_project_CCLE) - set(profiles_project_TCGA))
        newNamesProfiles = {item[0]: str(i)+"_"+item[0] for i, item in enumerate(merged.perfiles_CCLE.items())};
        featureLabel = {v: merged.perfiles_CCLE[k].infoFeatures["columns"] for k, v in newNamesProfiles.items()};
        # predictProject = projects[1]
    elif len(profiles_project_TCGA) > len(profiles_project_CCLE):
        miss_profile = list(set(profiles_project_TCGA) - set(profiles_project_CCLE))
        newNamesProfiles = {item[0]: str(i)+"_"+item[0] for i, item in enumerate(merged.perfiles_TCGA.items())};
        featureLabel = {v: merged.perfiles_TCGA[k].infoFeatures["columns"] for k, v in newNamesProfiles.items()};
        # predictProject = projects[0]
    else:
        miss_profile = ["None"]
        newNamesProfiles = {item[0]: str(i)+"_"+item[0] for i, item in enumerate(merged.perfiles_TCGA.items())};
        featureLabel = {v: merged.perfiles_TCGA[k].infoFeatures["columns"] for k, v in newNamesProfiles.items()};
        # predictProject = ""
        
    ccle = {newNamesProfiles[k]: v.profile for k,v in merged.perfiles_CCLE.items() if any(k in s for s in profiles_project_CCLE)};
    tcga = {newNamesProfiles[k]: v.profile for k,v in merged.perfiles_TCGA.items() if any(k in s for s in profiles_project_TCGA)};
    
    if len(profiles_project_CCLE) > len(profiles_project_TCGA):
        profile2use = ccle
    elif len(profiles_project_TCGA) > len(profiles_project_CCLE):
        profile2use = tcga
    else:
        profile2use = tcga
    
    theta_record = {newNamesProfiles[key]: value for key, value in merged.constraints_theta.items()};
    R_record = merged.constraints_r;

    for perfil in newNamesProfiles.values():
        try:
            # Create target Directory
            os.mkdir(path+'/'+perfil)
            print("    * Directory " , path+'/'+perfil ,  " Created ")
        except FileExistsError:
            print("    * Directory " , path+'/'+perfil ,  " already exists")
    
    print("\n    ________________________________________________________________________")
    print("    Running M&M-KjNMF : \n")
    print("     1. Hyperparameters: ")
    print("           * K: {} - r1: {} - r2: {} - d1: {} - L1: {} - L2: {}".format(K, r1, r2, d1, L1, L2))
    print("           * o1: {} - o2: {} \n".format(o1, o2))
    print("     2. Kernel parameters: ")
    print("           * sigma_{}: {} - sigma_{}: {}".format(project1, sigma_tcga, project2, sigma_ccle))
    print("           * sigma_{}_diff: {} - sigma_{}_diff: {} \n".format(project1, sigma_tcga_diff, project2, sigma_ccle_diff))
    print("     3. Some parameters: ")
    print("           * tol: {} - maxiter: {} - repeat: {} \n".format(tol, maxiter, repeat))
    
    ## SARA
    AUC_H_acum={k: 0 for k in ccle.keys()}
    ## END SARA
    
    ## Consensus Matrix
    consensus_matrix = defaultdict(lambda: defaultdict(np.ndarray(0)))
    
    for perfil, value in profile2use.items():
        if value.shape[1] > 1000:
            print("")
            print("    * For the profile {}, the coephinetic coefficient will not be calculated since the dimensions are very large. It would take too much time to calculate it.".format(perfil))
            print("")
        else:
            consensus_matrix[perfil] = np.zeros([value.shape[1], value.shape[1]])

# %% The M&M-jNMF method
    # --------------------------------------------- #
    #### ........... M&M-jNMF method ........... ####
    # --------------------------------------------- #
    for index in range(0,repeat):
        print('    ----------------------------------------------------------------------------------------------------')
        ## SARA
        ## Initialization stop_control and delta_control of each repetition
        stop_control = np.zeros([1,maxiter+1]); stop_control[0,0] =1;
        delta_control = np.zeros([1,maxiter+1]); delta_control[0,0] =1;
        
        ## Split into train and test dataset
        ccle_train={}; ccle_test={}; tcga_train={}; tcga_test={}
        ccle_train={k: v[index_train_ccle['repeat_'+str(index)],:] for k,v in ccle.items()};
        ccle_test={k: v[index_test_ccle['repeat_'+str(index)],:] for k,v in ccle.items()}; 
        tcga_train={k: v[index_train_tcga['repeat_'+str(index)],:] for k,v in tcga.items()};
        tcga_test={k: v[index_test_tcga['repeat_'+str(index)],:] for k,v in tcga.items()};
        
        # KERNEL ARRAYS, other option for sigmas. sigma_equalData, sigma_difData
        Kn=defaultdict(lambda: np.ndarray(0));
        Kn['K_ccle']={k: pairwise.rbf_kernel(v.T, Y=None, gamma=1/(2*(sigma_ccle**2))) for k,v in ccle_train.items()};
        Kn['K_tcga']={k: pairwise.rbf_kernel(v.T, Y=None, gamma=1/(2*(sigma_tcga**2))) for k,v in tcga_train.items()};
        
        KIJ_ccle=[]
        KIJ_tcga=[]
        for i in ccle.keys():
            otraMatriz_ccle=[]
            otraMatriz_tcga=[]
            for j in ccle.keys():
                matriz_ccle=pairwise.rbf_kernel(ccle_train[i].T, Y=ccle_train[j].T, gamma=1/(2*(sigma_ccle_diff**2)))
                otraMatriz_ccle.append(matriz_ccle)
                if i!=newNamesProfiles[miss_profile[0]] and j!=newNamesProfiles[miss_profile[0]]:
                    matriz_tcga=pairwise.rbf_kernel(tcga_train[i].T, Y=tcga_train[j].T, gamma=1/(2*(sigma_tcga_diff**2)))
                    otraMatriz_tcga.append(matriz_tcga)
                else: #Para TCGA NO hay DRUG, lleno el arreglo con basura (ceros) para conservar orden en codigos 0:0_drug, 1:1_mirna .. que seran usados posteriormente para leer esta matriz!
                    matriz_tcga=np.zeros([ccle_train[i].shape[1],ccle_train[j].shape[1]])
                    otraMatriz_tcga.append(matriz_tcga)               
            KIJ_ccle.append(otraMatriz_ccle)
            KIJ_tcga.append(otraMatriz_tcga)
    ## END SARA
        
        # ---------------------------------------------- #
        #### ........... Init H_I and A_I ........... ####
        # ---------------------------------------------- #
        ## SARA
        Ainit = defaultdict(lambda: np.ndarray(0));
        #INICIALIZA MATRICES A y H con valores estandar (literatura), (se podria revisar cosas mas robustas,articulo: DataFusion, autor: Marinka Zitnik)
        #AccleI pI x K . AtcgaI pI x K
        Ainit['A0_ccle']={k: np.random.uniform(0,1,[v.shape[1], K]) for k,v in ccle.items()}
        Ainit['A0_tcga']={k: np.random.uniform(0,1,[v.shape[1], K]) for k,v in tcga.items()}
        
        ## H_record_init (profile x rango x biomolecula)
        #HI.KxpI
        H_record_init = defaultdict(lambda: np.ndarray(0));
        H_record_init = {k: np.random.uniform(0,1,[K, v.shape[1]]) for k, v in ccle.items()};
        
        ## normalizacion de A y H
        A, H,sh = normalize_AH_rowH(Ainit, H_record_init,1);
        ## END SARA
        
# %% Error computation with initial parameters ('MUR')
        # ------------------------------------------------------- #
        #### ........... Initial error calculation ........... ####
        # ------------------------------------------------------- #
        AtKA= {};
        AtKA['AtKA_ccle']={k: np.dot(np.dot(A['A0_ccle'][k].T,Kn['K_ccle'][k]),A['A0_ccle'][k]) for k in ccle.keys()}
        AtKA['AtKA_tcga']={k: np.dot(np.dot(A['A0_tcga'][k].T,Kn['K_tcga'][k]),A['A0_tcga'][k]) for k in tcga.keys()}        
        Hsumtheta_record, sumHRt_record = compute_HsumTheta_sumHRt_record(H, theta_record, R_record);
        #matrices omicas, una al lado de la otra en un solo numpy (tcga NO tiene drug)
        #terminos de la FO 
        diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga = compute_diff(H,Hsumtheta_record,sumHRt_record,AtKA,KIJ_ccle,KIJ_tcga,A,L1,L2,r1,r2,d1,o1,o2,K,miss_profile);
        delta1_init,diff1_ccle, diff1_tcga = StopCriterion_rule1(Kn,A,H,AtKA,diff2,diff3,diff4_ccle,diff4_tcga,diff5,diff6_ccle,diff6_tcga,miss_profile);
        
#%% MM-KjNMF algorithm
        # -------------------------------------------------------------- #
        #### ........... Multiplicative Update Rule (MUR) ........... ####
        # -------------------------------------------------------------- #
        iteracion = 0;
        ## SARA
        A0=defaultdict(lambda: np.ndarray(0));
        sumKA_ccle = defaultdict(lambda: np.ndarray(0));
        sumKA_tcga = defaultdict(lambda: np.ndarray(0));
        KA=defaultdict(lambda: np.ndarray(0));
        AtKAH=defaultdict(lambda: np.ndarray(0));
        #N
        numN = len(H)
        ## END SARA
        
        while stop_control[0, iteracion] > tol and iteracion < maxiter:
            ## SARA
            ## *** Actualizacion de los A_ccleI y A_tcgaI fijando Hi. MUR            
            #% update A_ccleI with any Hi fixed            
            for i in ccle.keys():
                sumKA_ccle[i]=np.zeros([ccle[i].shape[1],K])
                ci=int(float(i[0:1]))
                for j in ccle.keys():
                    cj=int(float(j[0:1]))
                    if i!=j:
                        sumKA_ccle[i]=sumKA_ccle[i]+(np.dot(KIJ_ccle[ci][cj],A['A0_ccle'][j]))
            
            KA['KA_ccle']={k: np.dot(Kn['K_ccle'][k],A['A0_ccle'][k]) for k in ccle.keys()};
            A0['A0_ccle']={k: A['A0_ccle'][k]*np.divide((np.dot(Kn['K_ccle'][k].T,H[k].T)+(o1*sumKA_ccle[k])),(np.dot(np.dot(KA['KA_ccle'][k],H[k]),H[k].T)+(r1*A['A0_ccle'][k])+(o1*(numN-1)*KA['KA_ccle'][k])+np.finfo(float).eps)) for k in ccle.keys()};
            
            #% update A_tcgaI with any Hi fixed
            for i in tcga.keys():
                sumKA_tcga[i]=np.zeros([tcga[i].shape[1],K])
                ci=int(float(i[0:1]))
                for j in tcga.keys():
                    cj=int(float(j[0:1]))
                    if i!=j:
                        sumKA_tcga[i]=sumKA_tcga[i]+(np.dot(KIJ_tcga[ci][cj],A['A0_tcga'][j]))
            
            KA['KA_tcga']={k: np.dot(Kn['K_tcga'][k],A['A0_tcga'][k]) for k in tcga.keys()};
            A0['A0_tcga']={k: A['A0_tcga'][k]*np.divide((np.dot(Kn['K_tcga'][k].T,H[k].T)+(o2*sumKA_tcga[k])),(np.dot(np.dot(KA['KA_tcga'][k],H[k]),H[k].T)+(r2*A['A0_tcga'][k])+(o2*(numN-2)*KA['KA_tcga'][k])+np.finfo(float).eps)) for k in tcga.keys()};
            
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
            ## END SARA
            
            if iteracion == 0:
                delta_old = delta1_init;
            iteracion += 1;
            delta_control[0, iteracion] = np.around(delta, decimals=5);
            stop_control[0, iteracion] = np.around(np.abs((delta_old - delta)/(delta1_init - delta)), decimals=12);
            delta_old  = delta;

            if iteracion % 10 == 0:
                print('    * We are in iteration {: d} - stop_control: {}'.format(iteracion, stop_control[0, iteracion]))
                    
        print('')
        print('    * Repeat {:d}. Optimum value found in iteration: {:d} -> stop_control: {}'.format(index, iteracion, stop_control[0, iteracion]))        
        print('        -> diff1_{}: {} - diff1_{}: {} - diff2: {} - diff3: {}'.format(projects[0], np.round(diff1_ccle,4), projects[1], np.round(diff1_tcga,4), np.round(diff2,4), np.round(diff3,4)))
        print('        -> diff4_{}: {} - diff4_{}: {} - diff5: {}'.format(projects[0], np.round(diff4_ccle,4), projects[1], np.round(diff4_tcga,4), np.round(diff5, 4)))
        print('        -> diff6_{}: {} - diff6_{}: {} '.format(projects[0], np.round(diff6_ccle,4), projects[1], np.round(diff6_tcga,4)))
        print('        -> Objective Function value: {} \n'.format(delta))
        
        try:
            if stop_control[0][len(np.where(stop_control > 0)[1])-1] < best_stop_control[0][len(np.where(best_stop_control > 0)[1])-1]:
        
                best_A = A.copy();
                best_H = H.copy();
                best_delta_control = np.copy(delta_control);
                best_stop_control = np.copy(stop_control);
                best_repeat = index;
        
            else:
                pass
        except:
            #la primera vez, aca crea best_stop_control y demas
            best_A = A.copy();
            best_H = H.copy();
            best_delta_control = np.copy(delta_control);
            best_stop_control = np.copy(stop_control);
            best_repeat = index;
                
        ## Get co-clusters?
        _, _, connectivity_matrix,_,_ = jNMF_module(H,1.5, path, featureLabel, False, False, "H_first_second_max", False,nameDotPlot,merged);
        
        for perfil, value in profile2use.items():
            if value.shape[1] > 1000:
                pass
            else:
                for key, val in connectivity_matrix[perfil].items():
                    consensus_matrix[perfil][key] += val
        
        ## SARA
        # ------------------------------------------- #
        #### ........... TEST SECTION ........... ####
        # ------------------------------------------- #    
        A_learned=A.copy()
        H_learned=H.copy()    
        AUC_H,stop_control_test, delta_control_test=testSection(index,A_learned,H_learned,ccle,ccle_test,tcga,tcga_test,R_record,theta_record,param,miss_profile,newNamesProfiles,projects)
        # AUC_acum para calculo de AUC_average
        for k in profile2use.keys():
            AUC_H_acum[k]=AUC_H_acum[k]+AUC_H[k][0]
        ## ENDA SARA

    #### 
        
# %% Cophenetic coefficient calculation
    # ------------------------------------------------------------------------- #
    #### ........... Consensus matrix and cophenetic coefficient ........... ####
    # ------------------------------------------------------------------------- #
    for perfil in consensus_matrix.keys():
        consensus_matrix[perfil] = (consensus_matrix[perfil] / (repeat))
    
    rho_C = defaultdict(int)
    for perfil in consensus_matrix.keys():
        rho_C[perfil] = sum(sum(4*(consensus_matrix[perfil]-(1/2))**2)) / (consensus_matrix[perfil].shape[1] * consensus_matrix[perfil].shape[1])
    
    rho_C_avg = np.array([v for v in rho_C.values()]).mean()
    
    for perfil in profile2use.keys():
        if not perfil in list(rho_C.keys()):
            rho_C[perfil] = "NA"

# %%
    ## SARA
    # ---------------------------------------------------------------- #
    #### ........... Co-modulos con el mejor H (best_H) ........... ####
    # ---------------------------------------------------------------- #
    # ## Correr (comodulos) con el mejor H guardado, 
    # #comodule->diccionario, para cada omica (dict) muestra los K comodulos (key) con los features que pertenecen a dichos comodulos 
    # #comodule_count->diccionario, para cada omica (dict) muestra los K comodulos (key) con el # de features que pertenecen a dichos comodulos 
    # #connectivity_matrix->diccionario.Para cada omica tiene dict (tipo producto cartesiano), key(i,j), val(1 en el mismo comodulo,0 dlc)
    # #remove current co-module info. from folders 0_drug, 1_mirna, 2_mrna, 3_cnv
    # #print_savePath=True. Save current co-module info. from folders 0_drug, 1_mirna, 2_mrna, 3_cnv
    # """
    # Para todas las omicas, guarda info. de comodulos (nombre de moleculas en cada comod)
    # 0_drug/0_drug_comodX.txt
    # """
    #best_H=True & Enrichment=True. Read comodules from  2_mrna and returns number of enriched comodules. Save this info in co-md_records as record_mrna.csv
    comodule,comodule_count,_,cod_stats,cod_stats2 = jNMF_module(best_H, 1.5, path, featureLabel, True, True, "H_first_second_max", enrichment,nameDotPlot,merged);
    
    # -------------------------------------------------------------------------- #
    #### ........... Average AUC ........... ####
    # -------------------------------------------------------------------------- #  
    AUC_H_average={}    
    AUC_H_average={k: AUC_H_acum[k]/repeat  for k in ccle.keys()} 
    
    nOmics=len(AUC_H_average)
    AUC_acumm=0
    for i in profile2use.keys():
        AUC_acumm=AUC_acumm+AUC_H_average[i]
    AUC_H_average['average']=AUC_acumm/nOmics

    # -------------------------------------------------------------------------- #
    #### ........... Co-md_Stats........... ####
    # -------------------------------------------------------------------------- #  
    
    rawCod_stats = list()
    for j in cod_stats.values():
        rawCod_stats.append(j)
    
    # [min(v), max(v), np.mean(v), cod_no_vacios[k], cod_no_vacios_suma[k]]
    rawCod_stats = pd.DataFrame(rawCod_stats,columns=['profile', 'min_var_in_cod', 'max_var_in_cod', 'mean_var_in_cod', 'perc_cod_no_vacios', 'perc_var_incluidas']);
    rawCod_stats.index = rawCod_stats.profile
    rawCod_stats.drop(['profile'], axis=1, inplace=True)
    ## ENDA SARA
    
# %% Save the best H_I and W's low-rank matrices
    # ----------------------------------------------------- #
    #### ........... Save best_H and best_W ........... ####
    # ----------------------------------------------------- #
    folder=path+"/co-mod_bestW_and_bestH/"
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        os.unlink(file_path)
        
    for k, v in best_A.items():
        if k == "A0_ccle":
            p = projects[0]
        else:
            p = projects[1]
        for k1, v1 in v.items():
            np.savetxt(path+"/co-mod_bestW_and_bestH/best_A_"+p+"_"+k1+'_'+str(K)+'-r1_'+str(r1)+'-r2_'+str(r2)+'-L1_'+str(L1)+'-L2_'+str(L2)+'-d1_'+str(d1)+'.csv', v1, delimiter=',')
        
    for k, v in best_H.items():
        np.savetxt(path+"/co-mod_bestW_and_bestH/best_H_"+k+'_'+str(K)+'-r1_'+str(r1)+'-r2_'+str(r2)+'-L1_'+str(L1)+'-L2_'+str(L2)+'-d1_'+str(d1)+'.csv', v, delimiter=',')
    print('----------------------------------------------------------') 
    print("")
    print("")
    
    weightedAverage = 0.2 * rho_C_avg + 0.4 * AUC_H_average['average'] + 0.4 * cod_stats2[0]
   
# %% Pre-image
    eta=0.000001
    
    F = {k: defaultdict(lambda: np.array(0)) for k in best_A.keys()}
    for key, value in best_A.items(): ## Ej: dict_keys(['A0_ccle', 'A0_tcga'])
        for key2, value2 in value.items(): # Ej: dict_keys(['0_cnv', '1_drug', '2_mirna', '3_mrna'])
            T=[]
            if not "cnv" in key2:
                # F2 = []
                if "ccle" in key:
                    kernel = pairwise.rbf_kernel(ccle[key2].T, Y=None, gamma=1/(2*(sigma_ccle**2)))
                    X = ccle
                else:
                    kernel = pairwise.rbf_kernel(tcga[key2].T, Y=None, gamma=1/(2*(sigma_tcga**2)))
                    X = tcga
                    
                    
                for i in range(value2.shape[1]):
                    #Alpha
                    alpha=value2[:,i]
                    #alpha=np.matmul(np.matmul(W,H),w_t)
                    
                    p1=np.linalg.inv(np.matmul(X[key2],np.transpose(X[key2])))
                    
                    ### calcular de una vez -> KIJ_tcga!!!!
                    p2=np.matmul(X[key2],np.matmul(np.transpose(X[key2]),X[key2]) - eta*np.linalg.inv(kernel))
                    
                    preimage=np.matmul(np.matmul(p1,p2),alpha)
                    
                    # F2.append(preimage)
                    T.append(preimage)
                    
                # print(preimage)
                # T.append(preimage)
                
                ### obtiene matriz n x k. Escala cada columna entre [0, +1]
                # np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T)
                def q3(x, q):
                    q = np.quantile(x, q)
                    x = (x > q) * 1
                    return x
                F[key][key2] = np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T)
                # F[key][key2] = np.apply_along_axis(lambda x: q3(x, 0.75), 0, np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T))
                # F[key][key2] = np.apply_along_axis(lambda x: q3(x, 0.75), 0, np.array(T).T)
                # F[key][key2] = np.array(T).T
            else:
                pass
            
    ####
    module_patient_cell = {k: defaultdict(lambda: []) for k in best_A.keys()}
    for key, value in F.items():
        for key2, value2 in value.items():
            for i in range(0,value2.shape[1]):
                Q3 = np.quantile(value2[:,i], .75)
                ccle_ = list(np.where(value2[:,i] > Q3))
                
                module_patient_cell[key]['cluster_'+str(i)] += ccle_
    
    module_patient_cell_v2 = {k: defaultdict(lambda: float) for k in best_A.keys()}
    for key, value in module_patient_cell.items():
        
        for key2, value2 in value.items():
            if "ccle" in key:
                a = len(set(value2[0]) & set(value2[1]) & set(value2[2]))
                b = len(set(value2[0]) | set(value2[1]) | set(value2[2]))
                module_patient_cell_v2[key][key2] = a / b
                # print(a / b)
            else:
                a = len(set(value2[0]) & set(value2[1]))
                b = len(set(value2[0]) | set(value2[1]))
                module_patient_cell_v2[key][key2] = a / b
                # print(a / b)
    

# %% Pre-image
    eta=0.000001
    
    F = {k: defaultdict(lambda: np.array(0)) for k in best_A.keys()}
    for key, value in best_A.items(): ## Ej: dict_keys(['A0_ccle', 'A0_tcga'])
        for key2, value2 in value.items(): # Ej: dict_keys(['0_cnv', '1_drug', '2_mirna', '3_mrna'])
            T=[]
            if not "cnv" in key2:
                # F2 = []
                if "ccle" in key:
                    kernel = pairwise.rbf_kernel(ccle[key2].T, Y=None, gamma=1/(2*(sigma_ccle**2)))
                    X = ccle
                else:
                    kernel = pairwise.rbf_kernel(tcga[key2].T, Y=None, gamma=1/(2*(sigma_tcga**2)))
                    X = tcga
                    
                    
                for i in range(value2.shape[1]):
                    #Alpha
                    alpha=value2[:,i]
                    #alpha=np.matmul(np.matmul(W,H),w_t)
                    
                    p1=np.linalg.inv(np.matmul(X[key2],np.transpose(X[key2])))
                    
                    ### calcular de una vez -> KIJ_tcga!!!!
                    p2=np.matmul(X[key2],np.matmul(np.transpose(X[key2]),X[key2]) - eta*np.linalg.inv(kernel))
                    
                    preimage=np.matmul(np.matmul(p1,p2),alpha)
                    
                    # F2.append(preimage)
                    T.append(preimage)
                    
                # print(preimage)
                # T.append(preimage)
                
                ### obtiene matriz n x k. Escala cada columna entre [0, +1]
                # np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T)
                def q3(x, q):
                    q = np.quantile(x, q)
                    x = (x > q) * 1
                    return x
                F[key][key2] = np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T)
                # F[key][key2] = np.apply_along_axis(lambda x: q3(x, 0.75), 0, np.apply_along_axis(lambda a: np.interp(a, (a.min(), a.max()), (0, +1)), 0, np.array(T).T))
                # F[key][key2] = np.apply_along_axis(lambda x: q3(x, 0.75), 0, np.array(T).T)
                # F[key][key2] = np.array(T).T
            else:
                pass
            
    ####
    module_patient_cell = {k: defaultdict(lambda: []) for k in best_A.keys()}
    for key, value in F.items():
        for key2, value2 in value.items():
            for i in range(0,value2.shape[1]):
                Q3 = np.quantile(value2[:,i], .90)
                ccle_ = list(np.where(value2[:,i] > Q3))
                
                module_patient_cell[key]['cluster_'+str(i)] += ccle_
    
    module_patient_cell_v2 = {k: defaultdict(lambda: float) for k in best_A.keys()}
    for key, value in module_patient_cell.items():
        
        for key2, value2 in value.items():
            if "ccle" in key:
                a = len(set(value2[0]) & set(value2[1]) & set(value2[2]))
                b = len(set(value2[0]) | set(value2[1]) | set(value2[2]))
                module_patient_cell_v2[key][key2] = a / b
                # print(a / b)
            else:
                a = len(set(value2[0]) & set(value2[1]))
                b = len(set(value2[0]) | set(value2[1]))
                module_patient_cell_v2[key][key2] = a / b
                # print(a / b)
                

# %% Outputs
    # ------------------------------------- #
    #### ........... Outputs ........... ####
    # ------------------------------------- #
    paramResults = {'kernel': Kn, 'best_A': best_A, 'best_H': best_H, 
                  'rho_C': rho_C, 'AUC_H_average': AUC_H_average, 'weightedAverage':weightedAverage,
                  'best_stop_control': best_stop_control,'best_delta_control': best_delta_control, 
                  'rawCod_stats':rawCod_stats,'dotScore_stats':cod_stats2, 'FO': delta, 'F': F, 
                  "module_patient_cell_v2":module_patient_cell_v2}; 
    
    return paramResults
