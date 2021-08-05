# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 16:41:20 2019

@author: da.salazarb
"""
#%% ccle
import pandas as pd
import dataTCGA

class dataCCLE(dataTCGA.dataTCGA):
    """
    
    """
    
    def __init__(self, _ruta, _tipoPerfil, _estudio, _techSeq, _preTreatedData):
        self.tipoPerfil = _tipoPerfil
        self.info = list([_ruta, _estudio, _techSeq])
        if _preTreatedData != True:
            self.profile = pd.read_csv(_ruta, index_col=0)
            if _estudio != "all":
                self.profile = self.profile.loc[self.profile["DepMap_ID"].str.contains(_estudio),:]
            self.profile.drop(columns="DepMap_ID", inplace=True)
            
            if self.tipoPerfil == "drug":
                self.profile.drop(columns="drug_CCLE_name", inplace=True)
            
            self.infoFeatures = {"rows": list(map(lambda x: x.upper(), self.profile.index)), "columns": list(map(lambda x: x.upper().replace(self.tipoPerfil.upper()+"_", self.tipoPerfil.lower()+"_"),self.profile.columns))}
            self.profile = self.profile.values
            self.profile = self.profile.astype('float32') ## agilizar computo al fijar el tipo del array
        else:
            self.profile = pd.read_csv(_ruta, index_col=0)
            self.profileSobrante = ""
            self.infoFeatures = {"columns": list(self.profile.columns), "rows": list(self.profile.index)} ## guardar barcodes y nombres moleculas
            self.profile = self.profile.values
            self.profile = self.profile.astype('float32')
