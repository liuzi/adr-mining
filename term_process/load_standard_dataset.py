import os, sys
sys.path.append(os.path.abspath('/home/liu/project/adr-mining'))
import re
from utils._path import *
from utils._tools import *
from utils._get_item_args import get_item_args
import pandas as pd
import numpy as np
from os.path import join
import glob
from term_process.mapping_dataset import isOneToOne

def reset_name_db_columns(item_name_db,item):
    item_args=get_item_args(item)
    code_field, item_field_name=item_args["CODE_FIELD"],item_args["ITEM_FIELD_NAME"]
    item_name_db.columns=[code_field,item_field_name]
    item_name_db[item_field_name]=item_name_db[item_field_name].str.title()
    return item_name_db

def __load_itemid_name_db(item):

    db_paths, key_values = get_item_args(item)["NAME_DB"].values()
    # print(db_path)
    if(isinstance(db_paths, str)):
        db_path, key_value = db_paths, key_values
        db_files=glob.glob(db_path)
        # print(db_files)
        item_name_db=pd.concat([
            read_data(db_file,dtype=str) for db_file in db_files
        ], axis=0, ignore_index=True)[key_value].sample(frac=1).drop_duplicates(key_value[0])

    else:
        item_name_db_list=[]
    # first from clamp, second from UMLS ICD9CM
        for db_path, key_value in zip(db_paths, key_values):
            db_files=glob.glob(db_path)
            # print(db_files)
            item_name_db=pd.concat([
                read_data(db_file,dtype=str) for db_file in db_files
            ], axis=0, ignore_index=True)[key_value].sample(frac=1).drop_duplicates(key_value[0])
            item_name_db.columns=key_values[0]
            item_name_db_list=item_name_db_list+[item_name_db]
        item_name_db1, item_name_db2=item_name_db_list
        item_name_db1_only=pd.merge(item_name_db1,item_name_db2, indicator=True, how='outer')\
            .query('_merge=="left_only"').drop('_merge', axis=1)
        item_name_db=pd.concat([item_name_db1,item_name_db2],axis=0, ignore_index=True)


    return reset_name_db_columns(item_name_db.dropna(),item)

def append_rxnorm_name_db(item):
    # main_db=__load_itemid_name_db(item)
    sub_name_df_keys = np.array(get_item_args(item)["NAME_DB_1"])

    sub_name_df_left,sub_name_df_right=list(map(
        lambda db_path, key_value:\
            read_data(db_path,dtype=str)[key_value].drop_duplicates(),
        sub_name_df_keys[:,0],
        sub_name_df_keys[:,1]
    ))
    intersec_field=list(set.intersection(*map(set,sub_name_df_keys[:,1])))
    sub_name_df=left_join(sub_name_df_left,sub_name_df_right,intersec_field).dropna()
    remain_column=[col for col in sub_name_df.columns if col not in intersec_field]
    
    return reset_name_db_columns(
        sub_name_df[remain_column].sample(frac=1).drop_duplicates(remain_column[0]),
        item)
    # print(sub_name_df_right.head())

def load_icd9_name_db():
    return __load_itemid_name_db("DISEASE")

def load_umlscui_name_db():
    return __load_itemid_name_db("DISEASE_CUI")

def load_rxnorm_name_db():
    item="DRUG"
    code_field=get_item_args(item)["CODE_FIELD"]
    base_name_db = __load_itemid_name_db(item)
    ## Main db is more formal, keep all codes in main, remove duplicate codes from sub db
    sub_name_db = append_rxnorm_name_db(item)
    sub_name_db = sub_name_db[~sub_name_db[code_field].isin(base_name_db[code_field])]
    return pd.concat(
        [base_name_db,
        sub_name_db],
        axis=0,ignore_index=True)

def load_sider_cid_cui_db():
    item_args=get_item_args("SIDER4")
    code_field, item_field_name, key_fields=\
        item_args["CODE_FIELD"], item_args["ITEM_FIELD_NAME"], item_args["KEY_FIELDS"]
    sider_db=read_data(join("/data/liu/mimic3/SIDER","meddra_all_se.tsv"),sep="\t",header=None).iloc[:,[0,4,5]]
    sider_db.columns=key_fields
    return sider_db




# print(load_icd9_name_db().head())
# print(load_umlscui_name_db().head())

