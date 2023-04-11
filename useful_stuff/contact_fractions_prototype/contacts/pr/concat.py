#!/usr/bin/env python                                                       
# coding: utf-8 

import matplotlib.pyplot as plt                                                 
import numpy as np                                                              
import pandas as pd                                                             
import scipy as sp
import seaborn as sns

# membranes = ["TailOx", "RingOx", "Mix5"]
# inhibitors = ["None", "w9K", "w9W"]
# lipids = ['AnyLipid', 'POP3', 'POP2', 'POP1', 'IPE', 'PPE', 'PNG3',
#        'PNG1', 'DBG3', 'DBG1', 'PUPS', 'PUPI', 'PUPE', 'PUPC', 'POSM', 'POPS',
#        'POPI', 'POPE', 'POPC', 'POGS', 'PNSM', 'PNGS', 'PIPI', 'PFPC', 'PBSM',
#        'PAPS', 'PAPI', 'PAPE', 'PAPC', 'PAPA', 'PAP3', 'PAP2', 'PAP1', 'PADG',
#        'OUPS', 'OUPE', 'OUPC', 'OIPE', 'OIPC', 'OAPE', 'DPSM', 'DPPS', 'DPPC',
#        'DPGS', 'DPG3', 'DPG1', 'DPCE', 'DOPC', 'DBGS', 'AnySterol', 'CLRO',
#        'OCLR', 'CHOL']

# get file list

import os
flist = os.listdir()
cflist = [f for f in flist if "contacts" in f]

# load data file by file

c_df = pd.DataFrame(columns=["Residue"])
for i in cflist:
    df = pd.read_csv(i)
    c_df = c_df.append(df)


c_df[["resname", "resid"]] = c_df["Residue"].str.split("-", expand=True)
# c_df["resid"] = pd.to_numeric(c_df["f1 resid"]).map(GLYT2_OUT_RESID) 

# pad with dummy data so there is an entry for every unique [mem, drug, rep, resid] combo in the dataset 

# reslist = list(c_df.resid.unique())

# for mem in membranes:
#     for drug in inhibitors:
#         for rep in c_df.loc[(c_df["membrane"] == mem) & (c_df["inhibitor"] == drug)]["rep"].unique():
#             for res in reslist:
#                 if len(c_df.loc[(c_df["membrane"] == mem) & (c_df["inhibitor"] == drug) & (c_df["rep"] == rep) & (c_df["resid"] == res)]) == 0 :
#                     dummy_df = pd.DataFrame({"protein" : ["GlyT2"], "conformation" : ["outward"], "membrane": [mem], "inhibitor" : [drug], "rep" : [rep], "resid" : [res]})
#                     for lipid in lipids:
#                         dummy_df[lipid] = 0
#                     dummy_df["resname"] = c_df.loc[c_df["resid"] == res]["resname"].unique()[0]
#                     dummy_df["f1 resid"] = c_df.loc[c_df["resid"] == res]["f1 resid"].unique()[0]
#                     dummy_df["Residue"] = c_df.loc[c_df["resid"] == res]["Residue"].unique()[0]
#                     c_df = pd.concat([c_df,dummy_df])


# check quality

# for res in reslist:
#     if len(c_df.loc[(c_df["resid"] == res)]) != 27: print("error:",res)

# data is now ready for processing
c_df = c_df.drop(['Unnamed: 0'], axis = 1)
c_df.to_csv("../contact_data.csv", index=False)

u_df = c_df.loc[c_df["pose"] == "upper"]
u_df.to_csv("../contact_data_upper.csv", index=False)

l_df = c_df.loc[c_df["pose"] == "lower"]
l_df.to_csv("../contact_data_lower.csv", index=False)
