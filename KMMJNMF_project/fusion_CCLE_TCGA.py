# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 20:13:21 2019

@author: da.salazarb
"""
import json
import random
import os.path
import dataTCGA
import dataCCLE
import numpy as np
import pandas as pd
from collections import defaultdict

class fusion_CCLE_TCGA(dataCCLE.dataCCLE, dataTCGA.dataTCGA):
        
    def __init__(self, _path):
        """
        Crea objeto fusion_CCLE_TCGA, donde se almacena informacion de perfiles existentes, objetos perfiles_TCGA y perfiles_CCLE
        restricciones de tipo Theta y R y la ruta "path" donde estan los archivos y donde guarda archivos.
        """
        self.listaPerfiles = "" ## ["mrna", "protein", "mut", "drug"]
        self.project1 = ""
        self.project2 = ""
        self.featureLabel = "" # solo se debe hacer para uno, ya que quedan iguales para CCLE y TCGA
        self.perfiles_CCLE = defaultdict() ## diccionario con objetos dataTCGA
        self.perfiles_TCGA = defaultdict() ## diccionario con objetos dataCCLE
        self.constraints_theta = defaultdict() ## matriz de restricciones tipo theta
        self.constraints_r = np.array(0) ## matriz de matrices de restricciones tipo R
        self.path = _path
        self.index_train_tcga = dict()
        self.index_test_tcga = dict()
        self.index_train_ccle = dict()
        self.index_test_ccle = dict()
        self.command = ""
            
    def loadData(self):
        """
        Carga los perfiles y los guarda en el atributo self.perfiles_TCGA y self.perfiles.CCLE
        Crea perfil drug para TCGA en caso de ser llamado
        Identifica barcodes unicos en TCGA y modifica self.perfiles_TCGA
        Realiza estabilizacion de varianza
        """
        ## Verify if there are 2 projects
        project_profile = os.listdir(self.path+"INPUT DATA");
        
        self.project1 = sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0]; self.project2 = sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1];
        profiles_project1 = sorted(list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if self.project1 in p])));
        profiles_project2 = sorted(list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if self.project2 in p])));
        
        if len(profiles_project1) > len(profiles_project2):
            ccle_tcga = "high"
        elif len(profiles_project1) < len(profiles_project2):
            ccle_tcga = "less"
        else:
            ccle_tcga = "equal"
        
        if ccle_tcga == "high":
            projects = [sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0], sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1]];
            self.project1 = projects[0]; self.project2 = projects[1];
            print("    * The projects are {} and {} ".format(projects[0], projects[1]))
        elif ccle_tcga == "less":
            projects = [sorted(list(set([p[0:p.find("_")] for p in project_profile])))[1], sorted(list(set([p[0:p.find("_")] for p in project_profile])))[0]];
            self.project1 = projects[0]; self.project2 = projects[1];
            print("    * The projects are {} and {} ".format(projects[0], projects[1]))
        elif ccle_tcga == "equal":
            projects = sorted(list(set([p[0:p.find("_")] for p in project_profile])))
            self.project1 = projects[0]; self.project2 = projects[1];
            print("    * The projects are {} and {} ".format(projects[0], projects[1]))
            
        profiles_project_TCGA = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if projects[1] in p]));
        profiles_project_CCLE = list(set([p.replace(".csv", "")[p.replace(".csv", "").find("_")+1:] for p in project_profile if projects[0] in p]));
        print("        ** For {}, the profiles are {}".format(projects[1], profiles_project_TCGA))
        print("        ** For {}, the profiles are {}".format(projects[0], profiles_project_CCLE))        
        
        if len(profiles_project_CCLE) > len(profiles_project_TCGA):
            self.listaPerfiles = sorted(profiles_project_CCLE)
            miss_profile = list(set(profiles_project_CCLE) - set(profiles_project_TCGA))
            print("        ** The project {} does not have the profile {}".format(projects[1], miss_profile))
            print("           We will create the profile as a empty array")
        elif len(profiles_project_TCGA) > len(profiles_project_CCLE):
            self.listaPerfiles = sorted(profiles_project_TCGA)
            miss_profile = list(set(profiles_project_TCGA) - set(profiles_project_CCLE))
            print("        ** The project {} does not have the profile(s) {}".format(projects[0], miss_profile))
            print("           We will create the profile as a empty array")
        else:
            self.listaPerfiles = sorted(profiles_project_TCGA)
            miss_profile = []
        
        print("")
        ## profiles_project_TCGA
        self.perfiles_TCGA = {perfil.replace(".csv", "").replace(projects[1]+"_", ""): dataTCGA.dataTCGA(self.path+"INPUT DATA/"+perfil,"", "", "", True) for perfil in project_profile if projects[1] in perfil};
        ## profiles_project_CCLE
        self.perfiles_CCLE = {perfil.replace(".csv", "").replace(projects[0]+"_", ""): dataTCGA.dataTCGA(self.path+"INPUT DATA/"+perfil,"", "", "", True) for perfil in project_profile if projects[0] in perfil};
        
        ## If there are a miss profile then create a empty
        if len(miss_profile) > 0:
            if len(profiles_project_CCLE) > len(profiles_project_TCGA):
                for item in miss_profile:
                    self.perfiles_TCGA[item] = dataTCGA.dataTCGA(self.path+"INPUT DATA/"+projects[0]+"_"+item+".csv","", "", "", True);
                    self.perfiles_TCGA[item].profileSobrante = [] ## remover informacion del perfil
                    self.perfiles_TCGA[item].info = list(["", "", ""]) ## remover informacion del perfil
                    self.perfiles_TCGA[item].infoFeatures["rows"] = self.perfiles_TCGA[profiles_project_TCGA[0]].infoFeatures["rows"];
                    self.perfiles_TCGA[item].infoFeatures["columns"] = self.perfiles_CCLE[item].infoFeatures["columns"];
                    self.perfiles_TCGA[item].profile = self.perfiles_TCGA[item].profile.copy()
                    self.perfiles_TCGA[item].profile = np.empty([len(self.perfiles_TCGA[profiles_project_TCGA[0]].infoFeatures["rows"]),len(self.perfiles_TCGA[item].infoFeatures["columns"])]) ## crear matriz vacia con dimensiones nuevas
            elif len(profiles_project_TCGA) > len(profiles_project_CCLE): 
                for item in miss_profile:
                    self.perfiles_CCLE[item] = dataTCGA.dataTCGA(self.path+"INPUT DATA/"+projects[1]+"_"+item+".csv","", "", "", True);
                    self.perfiles_CCLE[item].profileSobrante = [] ## remover informacion del perfil
                    self.perfiles_CCLE[item].info = list(["", "", ""]) ## remover informacion del perfil
                    self.perfiles_CCLE[item].infoFeatures["rows"] = self.perfiles_CCLE[profiles_project_TCGA[0]].infoFeatures["rows"]
                    self.perfiles_CCLE[item].infoFeatures["columns"] = self.perfiles_TCGA[item].infoFeatures["columns"]
                    self.perfiles_CCLE[item].profile = self.perfiles_CCLE[item].profile.copy()
                    self.perfiles_CCLE[item].profile = np.empty([len(self.perfiles_CCLE[profiles_project_TCGA[0]].infoFeatures["rows"]),len(self.perfiles_CCLE[item].infoFeatures["columns"])]) ## crear matriz vacia con dimensiones nuevas
            else:
                pass
        else:
            pass
        
        ## Guardar index_data e index_test
        if os.path.isfile(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project1+'.json') == True and os.path.isfile(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project1+'.json') == True: ## ¿ya existe la restrccion?
            
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project1+'.json', 'r') as f:
                self.index_train_tcga = json.load(f)
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project1+'.json', 'r') as f:
                self.index_test_tcga = json.load(f)
        else:
            nrows=len(self.perfiles_TCGA[list(self.perfiles_TCGA.keys())[0]].infoFeatures["rows"])
            if nrows % 2 != 0:
                self.index_train_tcga = {'repeat_'+str(i): random.sample(list(range(nrows)),int(np.round((nrows+1)*.5, 0))) for i in range(100)}
                self.index_test_tcga = {'repeat_'+str(i): list(set(list(range(nrows))) - set(self.index_train_tcga['repeat_'+str(i)])) + random.sample(list(range(nrows)),1) for i in range(100)}
            else:
                self.index_train_tcga = {'repeat_'+str(i): random.sample(list(range(nrows)),int(np.round(nrows*.5, 0))) for i in range(100)}
                self.index_test_tcga = {'repeat_'+str(i):  list(set(list(range(nrows))) - set(self.index_train_tcga['repeat_'+str(i)])) for i in range(100)}
                
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project1+'.json', 'w') as fp:
                json.dump(self.index_train_tcga, fp)
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project1+'.json', 'w') as fp:
                json.dump(self.index_test_tcga, fp)
          
        if os.path.isfile(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project2+'.json') == True and os.path.isfile(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project2+'.json') == True: ## ¿ya existe la restrccion?
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project2+'.json', 'r') as f:
                self.index_train_ccle = json.load(f)
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project2+'.json', 'r') as f:
                self.index_test_ccle = json.load(f)
      
        else:
            nrows=len(self.perfiles_CCLE[list(self.perfiles_CCLE.keys())[0]].infoFeatures["rows"])
            if nrows % 2 != 0:
                self.index_train_ccle = {'repeat_'+str(i): random.sample(list(range(nrows)),int(np.round((nrows+1)*.5, 0))) for i in range(100)}
                self.index_test_ccle = {'repeat_'+str(i):  list(set(list(range(nrows))) - set(self.index_train_ccle['repeat_'+str(i)])) + random.sample(list(range(nrows)),1) for i in range(100)}
            else:
                self.index_train_ccle = {'repeat_'+str(i): random.sample(list(range(nrows)),int(np.round(nrows*.5, 0))) for i in range(100)}
                self.index_test_ccle = {'repeat_'+str(i):  list(set(list(range(nrows))) - set(self.index_train_ccle['repeat_'+str(i)])) for i in range(100)}
                
            ## guardar            
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_train_'+self.project2+'.json', 'w') as fp:
                json.dump(self.index_train_ccle, fp)
                
            with open(self.path+'/pathFeatureLabel/co-mod_train_test_split/index_test_'+self.project2+'.json', 'w') as fp:
                json.dump(self.index_test_ccle, fp)
        
        self.featureLabel = dataTCGA.dataTCGA.concatFeatureLabel(self.perfiles_CCLE, self.path, "") ## Modifica y guarda el featureLabel de los perfiles actualizados. Este es util para crear las restricciones (usando scripts en R).
        dataTCGA.dataTCGA.saveBarcode_projects(self.perfiles_TCGA[list(self.perfiles_TCGA.keys())[0]].infoFeatures["rows"], self.path,projects[1])
        dataTCGA.dataTCGA.saveBarcode_projects(self.perfiles_CCLE[list(self.perfiles_CCLE.keys())[0]].infoFeatures["rows"], self.path,projects[0])
        
        print("    * The dimensions are: ")
        print(" .................................. ")
        for perfil, v in self.perfiles_CCLE.items():
            print("{} dimensions for {}: {} x {}".format(perfil, projects[0], v.profile.shape[0], v.profile.shape[1]))
            print("{} dimensions for {}: {} x {}".format(perfil, projects[1], self.perfiles_TCGA[perfil].profile.shape[0], self.perfiles_TCGA[perfil].profile.shape[1]))
            print("--''--")
        
        print('---------------------------------------------------------- \n')
        
    def constraints_theta_method(self, path):
        """
        Identifica si existe la restriccion en self.path+pathFeatureLabel
        Si existe solo lo carga
        Si no existe, o crea una matriz de ceros o ejecuta Script en R
        Guarda las restricciones como un array en self.constraints_theta
        """
        print("    * We are creating or loading Theta-type constraints")
        import os.path
        
        profile1_profile2 = [];
        for file in os.listdir(self.path+"CONSTRAINTS DATA"):
            if file.endswith("_theta.csv"):
                profile1_profile2.append(file)
        
        theta = defaultdict(list)
        
        for perfil in self.listaPerfiles: ## recorre perfiles para identificar si existe o debe correr script R
            theta[perfil] = []
            if not any(perfil in s for s in profile1_profile2):
                print("It doesn't exist ... I'm creating "+perfil+" theta constraints")
                theta[perfil].append(np.zeros([len(self.perfiles_TCGA[perfil].infoFeatures["columns"]), len(self.perfiles_TCGA[perfil].infoFeatures["columns"])])) ## generar matriz de ceros
            else:
                print("oK ... I'm reading from CONSTRAINTS DATA the file "+perfil+"_"+perfil+"_theta.csv") 
                matriz = pd.read_csv(self.path+"CONSTRAINTS DATA/"+perfil+"_"+perfil+"_theta.csv") ## cargarla en caso de q exista
                # matriz = matriz.values ## solo interesan los valores no los nombres
                matriz = matriz.iloc[:,1:] ## solo interesan los valores no los nombres
                theta[perfil].append(matriz.values) ## guardar valores (numpy.array)
        self.constraints_theta = theta
        print("Ok, ready to use Theta constraints!")
        print("")
    
    def constraints_r_method(self, path):
        """
        Idnetifica las relaciones de perfiles que deben ser creadas, ejecutadas (script R) o buscadas en self.constraints_r
        Genera una matriz de matrices de las restricciones. Las diagonales son matrices de un perfil contra el mismo (i=j), 
        y los demas espacios son relaciones perfil vs otro perfil (i != j).
        """
        print("    * We are creating or loading R-type constraints")
        numero = 0 ## para asignar un numero a cada perfil
        codigo={} ## para guardar relacion perfil - codigo
        registro = dict() ## registro de si ya se ejecuto un script en R
        
        profile1_profile2 = [];
        for file in os.listdir(self.path+"CONSTRAINTS DATA"):
            if file.endswith("_r.csv"):
                profile1_profile2.append(file)
        
        for i in self.listaPerfiles: # da un numero a cada perfil
            codigo[i] = numero
            numero+=1
        
        matrizFinal = [] ## matriz que contiene el array final de restricciones
        for col in self.listaPerfiles: ## doble loop para comparar cada perfil entre si
            otraMatriz = [] ## matriz parcial para guardar una porcion de restriccion. Ej: la columna de (mrna vs mrna), (mut vs mrna) y (protein vs mrna)
            for row in self.listaPerfiles:
                if os.path.isfile(self.path+"CONSTRAINTS DATA/"+row+"_"+col+"_r.csv") == True: ## ¿ya existe la restrccion?
                    print("oK ... I'm reading from CONSTRAINTS DATA the file "+row+"_"+col+"_r.csv") 
                    matriz = pd.read_csv(self.path+"CONSTRAINTS DATA/"+row+"_"+col+"_r.csv") ## cargarla en caso de q exista
                    matriz = matriz.iloc[:,1:] ## solo interesan los valores no los nombres
                    matriz = matriz.values ## guardar valores (numpy.array)
                    registro[row+"_"+col+"_"] = [codigo[row], codigo[col]] ## Guarde el registro de que ya se corrio el script
                elif row == col: ## ¿es igual el perfil? crear matriz de ceros
                    matriz = np.zeros([len(self.perfiles_TCGA[row].infoFeatures["columns"]), len(self.perfiles_TCGA[col].infoFeatures["columns"])])
                elif any(col+"_"+row+"_" in s for s in profile1_profile2): ## esta la relacion de perfiles en el diccionario constraints_r_dict?
                    matriz = np.zeros([len(self.perfiles_TCGA[row].infoFeatures["columns"]), len(self.perfiles_TCGA[col].infoFeatures["columns"])])
                    registro[row+"_"+col+"_"] = [codigo[row], codigo[col]] ## Guarde el registro de que ya se corrio el script
                else: ## no hay nada, cree una matriz de ceros para la relacion
                    matriz = np.zeros([len(self.perfiles_TCGA[row].infoFeatures["columns"]), len(self.perfiles_TCGA[col].infoFeatures["columns"])])

                otraMatriz.append(matriz) ## se va guardando
                
            matrizFinal.append(np.array(otraMatriz, dtype=object)) ## ir guardando los arrays en un lista
            
            self.constraints_r = np.stack(matrizFinal, axis=-1) ## ir concatenando las restricciones
        
        for reg in registro.values():
            self.constraints_r[reg[1]][reg[0]] = self.constraints_r[reg[0]][reg[1]].T ## transponer matriz ubicada en self.constraints_r
        print("Ok, ready to use R constraints")
        print("")