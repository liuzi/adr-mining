import os, sys
# sys.path.append(os.path.abspath(os.path.join(
#     '/home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server', 'S2')))

from os.path import join
import numpy as np
import pandas as pd
import re
from term_process.load_standard_dataset import load_sider_cid_cui_db
import urllib.request 
import xmltodict
# from utils._tools import *

sider_drug_file="/data/liu/mimic3/SIDER/drug_names.tsv"
umls_term_file="/data/liu/mimic3/MAPPING/TERMINOLOGIES/MRCONSO.RRF"
single_drug_term_file="/data/liu/mimic3/CLAMP_NER/single_drug_analysis/drug_rxnorm_name.csv"
discode_dict={"icd9": "MTHICD9", "sndus":"SNOMEDCT_US"}
drugcode_dict={"rxn":"RXNORM","atc":"ATC"}


def process_cui_item(grep_result):
    if(not grep_result):
        return None
    result_parts=grep_result.split("|")  
    if(len(result_parts)==19):
        return result_parts[-5]
    else:
        raise ValueError('Incomplete Data for Disease/Drug Searched by CUI/RXNORM')

# def process_rxnorm_drug(grep_result):
#     if(not grep_result):
#         return None
#     result_parts=grep_result.split("|")  
#     print(len(result_parts))
#     if(len(result_parts)==19):
#         return result_parts[-5]
#     else:
#         raise ValueError('Incomplete Data for Drug Searched by RXNORM')

def process_cui_drug(grep_result):
    return 0

def search_cui(umls_cui, itemcode):
    language_version="|ENG|"
    grep_command="grep '{}{}' {} | grep -m 1 {}".format(
        umls_cui, language_version, umls_term_file, itemcode)
    result=os.popen(grep_command)
    if itemcode in discode_dict.values():
        return process_cui_item(result.read())
    elif itemcode in drugcode_dict.values():
        return process_cui_item(result.read())
    else:
        raise ValueError('{} is not found, plz add in the code')

def search_rxnorm(rxnorm):
    grep_command="grep '{}' {}".format(
        rxnorm, single_drug_term_file)
    result=os.popen(grep_command).read()
    return result.split(",")[-1].strip()

def search_rxnorm_long(rxnorm):
    language_version="|ENG|"
    grep_command="grep '{}' {} | grep -m 1 '{}||RXNORM'".format(
        language_version, umls_term_file,rxnorm)
    result=os.popen(grep_command).read()
    return process_cui_item(result)
    # if(result):

    # return result

def get_drugname_byrxcui_api(rxcui):
    # return None
    result=search_rxnorm(rxcui)
    if(result):
        return result
    else:
        with urllib.request.urlopen("https://rxnav.nlm.nih.gov/REST/rxcui/%s"%rxcui) as url:
            data = url.read()
            data = xmltodict.parse(data)
            return data['rxnormdata']['idGroup'].get('name',rxcui)

def search_cid_for_rxnorm(rxnorm):
    drug_name=search_rxnorm(rxnorm)
    grep_command_prefix="grep -i '{}' {} "
    result=os.popen(grep_command_prefix.format(
        drug_name, sider_drug_file)).read()
    if(result):
        return result.split('\t')[0]
    else:
        result=os.popen(grep_command_prefix.format(\
            drug_name.split(' ')[0], sider_drug_file)).read()
        if(result):
            return result.split('\t')[0]
        else:      
            drug_name = (search_rxnorm_long(rxnorm)).split(' ')[0]
            result=os.popen(grep_command_prefix.format(\
                drug_name, sider_drug_file)).read()
            if(result):
                return result.split('\t')[0]
            else:
                if(drug_name.lower()=="vancomycin"):
                    return "CID100005651"

def get_side_effects_cui(rxnorm):
    sider4=load_sider_cid_cui_db()
    # sider_dict={}
    # for rxnorm_id in rxnorms:
    cid=search_cid_for_rxnorm(rxnorm)
    return sider4[sider4["STICH_ID_FLAT"]==cid].iloc[:,1:].drop_duplicates()
        # sider_dict[rxnorm_id]=sider4[sider4["STICH_ID_FLAT"]==cid].iloc[:,1:].drop_duplicates()
    # return sider_dict


# if __name__ == "__main__":
#     test_rxnomrs=["1658259","866924","966571","885257","836358","855334",\
#         "855290","1808219","1724383","1719291","1807516","213169","1659151"]
#     get_side_effects_cui(*test_rxnomrs)
    # print(search_rxnorm("212033"))
#     # search_cui("C0428977","SNOMEDCT_US")
#     print(search_cui("C001880","MTHICD9"))
#     # print("ATC" in drugcode_dict.values())

