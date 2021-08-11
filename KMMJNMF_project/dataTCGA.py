# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 16:01:39 2019

@author: da.salazarb
"""
#%% tcga
import numpy as np
import pandas as pd
# np.random.seed(11)
class dataTCGA:
    """
    
    """
    pathFeatureLabel="D:/pathFeatureLabel"
    # pathFeatureLabel="E:/pathFeatureLabel"
    command = 'C:/Program Files/R/R-4.0.5/bin/Rscript.exe'
    # command = 'C:/Users/da.salazarb/Documents/R/R-3.6.3/bin/Rscript.exe'
    
    def __init__(self, _ruta, _tipoPerfil, _estudio, _techSeq, _preTreatedData):
        self.tipoPerfil = _tipoPerfil ## "mrna", "protein", "cnv", "mut", "drug", "mirna", "meth"
        self.info = list([_ruta, _estudio, _techSeq]) ## informacion adicional
        if _preTreatedData != True:
            self.profile = pd.read_csv(_ruta, sep="\t")    
            
            if _tipoPerfil == "cnv": ## para cnv solo seleccionar la banda positiva
                self.profile = self.profile[self.profile["Strand"] == "+"]
            elif _tipoPerfil == "meth": ## para meth cambiar el nombre de genes a GeneSymbol
                self.profile = self.profile.drop(self.profile.index[0])
                self.profile.rename(columns={"Hybridization REF":"GeneSymbol"}, inplace=True)
            elif _tipoPerfil == "protein": ## para protein dejar anticuerpo como GeneSymbol (para emparejamiento facil con CCLE)
                self.profile.rename(columns={"GeneSymbol":"GeneSymbol_1", "ProteinAntibody":"GeneSymbol"}, inplace=True)
                self.profile["GeneSymbol"] = list(map(lambda x: x[:-4], self.profile["GeneSymbol"]))
            elif _tipoPerfil == "drug":
                self.profile = pd.read_csv(_ruta, index_col=0) ## cargar archivos
                self.profile.drop(columns="drug_CCLE_name", inplace=True)
                #self.profile.rename(columns={"drug_CCLE_name":"GeneSymbol"}, inplace=True)
                self.profile = self.profile.T
                self.profile["GeneSymbol"] = list(map(lambda x: x[5:], self.profile.index))
            else:
                pass
            
            self.profileSobrante = self.profile.filter(regex="^(?!TCGA)") ## algunos perfiles traen columnas de mas
            self.profile = self.profile.loc[self.profile["GeneSymbol"] != "?",:] ## eliminar variables cuyo nombre sea: ?
            self.profile.columns = list(map(lambda x: x[0:16], self.profile.columns)) ## tomar los primeros 15 caracteres del barcode de TCGA
            self.infoFeatures = {"columns": list(map(lambda x: self.tipoPerfil + "_" + x.upper(), self.profile["GeneSymbol"])), "rows": self.profile.filter(regex="^TCGA").columns} ## guardar barcodes y nombres moleculas
            self.profile = self.profile.filter(regex="^TCGA").values.T ## transponer para que quede de la forma (paciente vs molecula)
            self.profile = self.profile.astype('float32') ## agilizar computo al fijar el tipo del array
            
            if _tipoPerfil == "mirna": ## modificar mirna porque se cuantificaron los stem-loops y salen varias mediciones para un solo mirna
                self.process_mirna_TCGA()
            else:
                pass
        else:
            self.profile = pd.read_csv(_ruta, index_col=0)
            self.profileSobrante = ""
            self.infoFeatures = {"columns": list(self.profile.columns), "rows": list(self.profile.index)} ## guardar barcodes y nombres moleculas
            self.profile = self.profile.values
            self.profile = self.profile.astype('float32')
        
    def norm_pseudo_log2(self): ## log2 de la pseudo suma de cada dato
        self.profile = np.log2(self.profile + 1)
        #return self.profile
        
    def norm_min_max(self): ## min-max escala
        from sklearn.preprocessing import MinMaxScaler
        scaler = MinMaxScaler(copy=False, feature_range=(0, 1))
        scaler.fit_transform(self.profile)
        #return self.profile
        
    def preProcess_sofImpute(self):
        from softImpute_func import SoftImpute
        impute = SoftImpute()
        impute.fit(self.profile)
        self.profile = impute.predict(self.profile)
        
    def preProcess_IterativeImputer(self):
        from sklearn.experimental import enable_iterative_imputer
        from sklearn.impute import IterativeImputer
        from sklearn.ensemble import ExtraTreesRegressor
        impute_estimator=ExtraTreesRegressor(n_estimators=10, random_state=0)
        imput=IterativeImputer(random_state=0, estimator=impute_estimator)
        imput.fit_transform(self.profile)
        
    def preProcess_KNNImputer(self):
        from sklearn.impute import KNNImputer
        imputer = KNNImputer(n_neighbors=2, weights="distance")
        imputer.fit_transform(self.profile)
        
    def preProcess_Bayesianridge(self):
        from sklearn.experimental import enable_iterative_imputer
        from sklearn.impute import IterativeImputer
        from sklearn.linear_model import BayesianRidge
        impute_estimator=BayesianRidge()
        imput=IterativeImputer(random_state=0, estimator=impute_estimator)
        imput.fit_transform(self.profile)
        
    def preProcess_sofImpute_R(self, pathFeatureLabel, scriptName, command):
        import subprocess
        
        np.savetxt(pathFeatureLabel+"/co-mod_best_results_last_run/profile2R_softimpute.txt", self.profile)
        # command = 'C:/Program Files/R/R-3.6.2/bin/Rscript.exe'
        # command = 'C:/Users/da.salazarb/Documents/R/R-3.6.3/bin/Rscript.exe'
        path2script = '--vanilla ' + pathFeatureLabel + "/" + scriptName
        cmd = [command, path2script] + [pathFeatureLabel]
        x = subprocess.Popen(cmd).wait()
        self.profile = np.loadtxt(pathFeatureLabel+"/co-mod_best_results_last_run/profile2R_softimpute.txt")
    
    def to_data_frame(self): ## tomar la info del objeto y convertir en data frame, util para procesar la informacion
        return pd.DataFrame(self.profile, index=self.infoFeatures["rows"], columns=self.infoFeatures["columns"])
    
    @staticmethod
    def concatFeatureLabel(perfiles, pathFeatureLabel, prefijo): ## concatena todos los nombres de todas las moleculas en una lista
        ## util para los scripts que se ejecuten en R
        featureLabel = []
        for i in range(0,len(perfiles.keys())):
            featureLabel.extend(perfiles[list(perfiles.keys())[i]].infoFeatures["columns"])
            
        with open(pathFeatureLabel+'/pathFeatureLabel/co-mod_records/featureLabel'+prefijo+'.txt', 'w') as filehandle:
            for listitem in featureLabel:
                filehandle.write('%s\n' % listitem)
        return featureLabel
    
    @staticmethod
    def saveBarcode_projects(barcodes,pathFeatureLabel,project): ## guarda los identificadores de los pacientes TCGA
        with open(pathFeatureLabel+'/pathFeatureLabel/co-mod_records/barcodes_'+project+'.txt', 'w') as filehandle:
            for listitem in barcodes:
                filehandle.write('%s\n' % listitem)
    
    @staticmethod
    def constraints(pathFeatureLabel, scriptName, perfil_1, perfil_2, numero, command): ## ejecuta un script de R, retorna la restriccion
        ## ver mas INFO: 
        ## https://www.r-bloggers.com/integrating-python-and-r-part-ii-executing-r-from-python-and-vice-versa/
        import subprocess
        
        # command = 'C:/Program Files/R/R-3.6.2/bin/Rscript.exe'
        # command = 'C:/Users/da.salazarb/Documents/R/R-3.6.3/bin/Rscript.exe'
        #path2script = '--vanilla 00_gene-gene_BioGRID.R'
        path2script = '--vanilla ' + pathFeatureLabel + "/" + scriptName
        
        cmd = [command, path2script] + [pathFeatureLabel] #cmd = [command, path2script] + args. Ver pagina rPubs.
        x = subprocess.Popen(cmd).wait()
        if x == 0:
            print("Todo Ok ... ")
        else:
            print("Something is wrong with " + scriptName)
        
        if perfil_1 == perfil_2:
            df_const = pd.read_csv(pathFeatureLabel+"/co-mod_interact_R_Theta/interact_"+perfil_1+"_"+perfil_2+"_"+numero+"_theta.csv", index_col=0)
        else:
            df_const = pd.read_csv(pathFeatureLabel+"/co-mod_interact_R_Theta/interact_"+perfil_1+"_"+perfil_2+"_r.csv", index_col=0)
        df_const = df_const.values
        return df_const
    
    def process_mirna_TCGA(self): ## procesa mirna porque existen varias repeticiones de un solo mirna (stem-loops)
        df = pd.DataFrame(self.profile, index=self.infoFeatures["rows"], columns=self.infoFeatures["columns"])
        df_1 = df.drop(list(df.filter(regex="-[0-9]$").columns), axis=1) #seleccionar mirna sin repaticiones
        df_names = list(np.unique(list(map(lambda  x: x[:-2],df.filter(regex="-[0-9]$").columns)))) ## nombres sin -#
        df_temp = df.filter(regex="-[0-9]$")
        for names in df_names:
            df_1[names] = df_temp.filter(regex="^"+names).mean(axis=1)
        self.profile = df_1.values
        self.infoFeatures = {"columns": list(df_1.columns), "rows": list(df_1.index)}
        
    def metric_rss(self, W, H):
        rss = sum((self.profile-W*H)**2)
        return rss
    
    @staticmethod
    def variance_threshold(data, perfil):
        # tomado de https://gist.github.com/calippo/9f726ef197b36d6e9f39
        from sklearn.feature_selection import VarianceThreshold
        if perfil == "cnv" or perfil == "mirna":
            selector  = VarianceThreshold(threshold=(.99 * (1 - .99)))
        else:
            selector  = VarianceThreshold(threshold=(.80 * (1 - .80)))
        selector.fit_transform(data)
        #return data.loc[:,data.columns[selector.get_support(indices=True)]]
        return data[[c for (s, c) in zip(selector.get_support(), data.columns.values) if s]]
    
    def variance_stable(self):
        from scipy.stats import skew #, kurtosis
        def perc_99(a):
            return np.percentile(a, 99)
        
        def perc_1(a):
            return np.percentile(a, 1)
        
        # a = merged.perfiles_TCGA["mirna"].profile
        percentile_99 = np.apply_along_axis(perc_99, 0, self.profile)
        percentile_1 = np.apply_along_axis(perc_1, 0, self.profile)
        for i in range(0, self.profile.shape[1]):
            if skew(self.profile[:,i]) > 2:
                self.profile[:,i][self.profile[:,i]>percentile_99[i]] = percentile_99[i]
            elif skew(self.profile[:,i]) < -2:
                self.profile[:,i][self.profile[:,i]<percentile_1[i]] = percentile_1[i]
        
    
    def __str__(self): ## imprime informacion del objeto
        return("{} {}".format("We are working with: ", self._tipoPerfil) + ". " + 
               "{} {}".format("The original shape is: ", self.profile.shape))