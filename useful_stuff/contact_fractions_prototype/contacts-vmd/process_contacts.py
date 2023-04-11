#!/usr/bin/env python                                                       
# coding: utf-8 

import matplotlib.pyplot as plt                                                 
import numpy as np                                                              
import pandas as pd                                                             
import scipy as sp
import seaborn as sns

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("protein", type=str, help="protein name")
# parser.add_argument("conf", type=str, help="inward outward")
parser.add_argument("membrane", type=str, help="Neuronal TailOx RingOx Mix5 Mix25")
parser.add_argument("inhibitor", type=str, help="None E25R E25S E73")
parser.add_argument("site", type=str, help="None S1 or LAS")
parser.add_argument("pose", type=str, help="None S1 upper lower flex phenyl antiphenyl  ")
parser.add_argument("rep", type=int, help="replicate number as int")
args = parser.parse_args()


import os
flist = os.listdir()
cfracs = [f for f in flist if "contact-" in f]

d = {f.split("-")[1].split('.')[0]:f for f in flist if "contact-" in f}

def loadcontacts(fname, group):
    df = pd.read_csv(fname, header=4,sep=" \t\t").rename(columns={"fraction" : group})
    df["Residue"] = df["Residue"].str[1:]
    df[group] = df[group].str[:-1]
    return df

c_df = pd.DataFrame(columns=["Residue"])

for k,v in d.items():
    df = loadcontacts(v , k)
    c_df = c_df.merge(df, how="outer", on="Residue")

c_df["protein"] = args.protein
# c_df["conformation"] = args.conf
c_df["membrane"] = args.membrane
c_df["inhibitor"] = args.inhibitor
c_df["site"] = args.site
c_df["pose"] = args.pose
c_df["rep"] = args.rep

sys_str = "_".join([args.inhibitor, args.pose, args.protein, args.membrane, "r" + str(args.rep)])

c_df["system"] = sys_str

fname = "/media/uqadaqu1/GUMPTION/Atomwise2/analysis/contacts/pr/" + sys_str + "_contacts.csv"
fname

c_df.to_csv(fname)

