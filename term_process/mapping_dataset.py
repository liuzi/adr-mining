import os, sys
# import glob
from os.path import join
import numpy as np
import pandas as pd
import subprocess

sys.path.append(os.path.abspath('/home/liu/project/adr-mining'))
from term_process.icd9_conversion import decimal_to_short
# from icd9_conversion import decimal_to_short
from utils._tools import *
from utils._path import read_prefix

# import re

def load_mrconso_args():

    python_path="/home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server/S2/term_process/"
    mapping_prefix="/data/liu/mimic3/MAPPING/TERMINOLOGIES"
    mrconiso_icd9_file="MRCONSO_ICD9.RRF"
    mrconiso_icd9_title=["CUI", "LAT", "TS", "LUI", "STT" ,"SUI" ,"ISPREF", "AUI","SAUI","SCUI",\
        "SDUI","SAB","TTY","CODE","STR","SRL","SUPPRESS","CVF","NONE"]
    mrconiso_icd9_selected_title=["CUI","TTY","CODE","STR"]
    mrconiso_icd9_ab_file="MRCONSO_ICD9CM_AB"

    return python_path, mapping_prefix, mrconiso_icd9_file, \
        mrconiso_icd9_title, mrconiso_icd9_selected_title, mrconiso_icd9_ab_file


def isOneToOne(df, col1, col2):
    first = df.groupby(col1)[col2].count().max()==1
    if(not first):
        first_df=df.groupby(col1)[col2].count()\
            .reset_index(name='count').sort_values(['count'], ascending=False)
        print("First 20 %s which have 2 values in %s"%(col1, col2))
        print(first_df.head(20))
        print(len(first_df[first_df['count']==2]))
    second = df.groupby(col2)[col1].count().max()==1
    if(not second):
        second_df=df.groupby(col2)[col1].count()\
            .reset_index(name='count').sort_values(['count'], ascending=False)
        print("First 20 %s which have 2 values in %s"%(col2, col1))
        print(second_df.head(20))
        print(len(second_df[second_df['count']==2]))
    print("%s only have one value in %s : %s"%(col1, col2, "True" if first else "False"))
    print("%s only have one value in %s : %s"%(col2, col1, "True" if second else "False"))

# -------------------------------------------------------
# NOTE: ICD9->UMLS CUI
# -------------------------------------------------------
def get_mrconso_icd9_cui_df():
    python_path, mapping_prefix, mrconiso_icd9_file, \
        mrconiso_icd9_title, mrconiso_icd9_selected_title, mrconiso_icd9_ab_file = load_mrconso_args()
    if(not os.path.exists(join(mapping_prefix, mrconiso_icd9_file))):
        os.system("%s %s"%(join(python_path,"mrconso_icd9_cui.sh"), mapping_prefix))

    if(not os.path.exists(join(mapping_prefix, "%s.csv"%mrconiso_icd9_ab_file))):
        # NOTE: ignore last column of df from RRF file since this kind of file ends with "|", 
        # which will lead to additional Null column
        mrconso_icd9_df = pd.read_csv(join(mapping_prefix, mrconiso_icd9_file),\
            sep="|", header=None, names=mrconiso_icd9_title).iloc[:,:-1][mrconiso_icd9_selected_title]
        mrconso_icd9_df=mrconso_icd9_df[mrconso_icd9_df["TTY"]=="AB"] 
        # cannot directly remove decimal, 695.2 and 69.52 are different, so use decimal to short
        mrconso_icd9_df["CODE"]=mrconso_icd9_df["CODE"].apply(decimal_to_short)
        write2file(mrconso_icd9_df.drop_duplicates(),join(mapping_prefix,mrconiso_icd9_ab_file))
    else:
        mrconso_icd9_df = read_data(join(mapping_prefix, mrconiso_icd9_ab_file))
    

    print("Numner of unique CUIs: %d"%len(mrconso_icd9_df["CUI"].unique()))
    print("Numner of unique ICD9s: %d"%len(mrconso_icd9_df["CODE"].unique()))
    isOneToOne(mrconso_icd9_df,"CUI","CODE")

    ## CHECK wiht mimiciii ICD9
    mimic_icd9= read_data(join("/data/MIMIC3","D_ICD_DIAGNOSES"),dtype=str)["ICD9_CODE"]
    icd9_mimic_notin_mrconso = list(set(mimic_icd9)-set(mrconso_icd9_df["CODE"]))
    print("ICD9s from MIMIC-III which cannot be found in MRCONSO: %s"%((', ').join(icd9_mimic_notin_mrconso)))

    return mrconso_icd9_df

def load_rxnsat_args():

    # python_path="/home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server/S2/term_process/"
    mapping_prefix="/data/liu/mimic3/MAPPING/TERMINOLOGIES/RAR/rrf"
    rxnsat_ndc_file="RXNSAT.RRF"
    rxnsat_ndc_title=["RXCUI", "LUI", "SUI", "RXAUI", "STYPE" ,"CODE" ,"ATUI", "SATUI","ATN","SAB",\
        "ATV","SUPPRESS","CVF","NONE"]
    rxnsat_ndc_selected_title=["RXCUI", "RXAUI","STYPE","CODE","ATN","SAB","ATV","CVF"]


    return mapping_prefix, rxnsat_ndc_file, rxnsat_ndc_title, rxnsat_ndc_selected_title

def get_rxnsat_rxnorm_ndc_df():
    mapping_prefix, rxnsat_ndc_file, rxnsat_ndc_title, rxnsat_ndc_selected_title = load_rxnsat_args()
    rxnsat_ndc_df = pd.read_csv(join(mapping_prefix, rxnsat_ndc_file),\
        sep="|", header=None, names=rxnsat_ndc_title).iloc[:,:-1]
        # [rxnsat_ndc_selected_title]
    return rxnsat_ndc_df

# tt=get_rxnsat_rxnorm_ndc_df().query('ATN=="RXNORM" & ATV==40209')
# tt=get_rxnsat_rxnorm_ndc_df().query('RXCUI==40209')

# print(tt.head())


# mrconso_icd9_cui()

# pres_df = read_data(
#     join(read_prefix,'PRESCRIPTIONS'),
#     dtype={"NDC":str})[["NDC"]].dropna(subset=['NDC'])
# pres_df = pres_df[pres_df['NDC']!="0"].drop_duplicates()
# pres_df=pres_df.assign(NDC_LENGTH=pres_df["NDC"].apply(len))
# print(pres_df["NDC_LENGTH"].unique())