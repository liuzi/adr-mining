import os, sys
sys.path.append(os.path.abspath(os.path.join(
    '/home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server', 'S2')))
from os.path import join
from utils._path import *


def get_item_args(item):
    args=[
        "CODE_FIELD","INDEX_FILE","INDEX_NAME_FILE","CODE_INDEX",
        "FIG_SUPTITLE","FILENAME_ALL","FILENAME_UNIQUE","FILENAME_STATS",
        "ITEM_NAME","ITEM_FIELD_NAME","NAME_DB"]
    disease_dict={
        "CODE_FIELD":"ICD9_CODE",
        "INDEX_FILE":"disease_index_icd9",
        "INDEX_NAME_FILE":"disease_index_icd9_name",
        "CODE_INDEX":"Disease Index",
        "FIG_SUPTITLE":"TOP 20 Frequent%s Diseases - Count of Episodes, Cluster: %s",
        "FILENAME_ALL":"COMBINED_%s%s_diseases_%s",
        "FILENAME_UNIQUE":"DISTINCT_%s%s_diseases_%s",
        "FILENAME_STATS":"STATS_%s%s_diseases_%s",
        "ITEM_NAME":"Disease",
        "ITEM_FIELD_NAME":"SHORT_TITLE",
        "NAME_DB":{"DB_path":join(read_prefix,"D_ICD_DIAGNOSES.csv"), "key_value":["ICD9_CODE","SHORT_TITLE",]},
        "CODE_NAME_FILE":"disease_icd9_name"

    }
    drug_dict={
        "CODE_FIELD":"RxNorm Code",
        "INDEX_FILE":"drug_index_rxnorm",
        "INDEX_NAME_FILE":"drug_index_rxnorm_name",
        "CODE_INDEX":"Drug Index",
        "FIG_SUPTITLE":"TOP 20 Frequent%s Drugs - Count of Episodes, Cluster: %s",
        "FILENAME_ALL":"COMBINED_%s%s_drugs_%s",
        "FILENAME_UNIQUE":"DISTINCT_%s%s_drugs_%s",
        "FILENAME_STATS":"STATS_%s%s_drugs_%s",
        "ITEM_NAME":"Drug",
        "ITEM_FIELD_NAME":"Drug NAME",
        "NAME_DB":{"DB_path":join(concat_clamp_prefix,"*CUI_RXNORM_NAME.csv"), "key_value":["RxNorm","NE"]},
        "NAME_DB_1":[
            [join(singledrug_featurepreprocess_prefix,"pres_rxnorm_df") ,["RxNorm","NDC"]],
            [join(read_prefix,"PRESCRIPTIONS.csv"), ["NDC","DRUG"]]],
        "CODE_NAME_FILE":"drug_rxnorm_name"



    }
    disease_cui_dict={
        "CODE_FIELD":"CUI",
        "INDEX_FILE":"disease_index_cui",
        "INDEX_NAME_FILE":"disease_index_cui_name",
        "CODE_INDEX":"Disease Index",
        "FIG_SUPTITLE":"TOP 20 Frequent%s Diseases - Count of Episodes, Cluster: %s",
        "FILENAME_ALL":"COMBINED_%s%s_diseases_%s",
        "FILENAME_UNIQUE":"DISTINCT_%s%s_diseases_%s",
        "FILENAME_STATS":"STATS_%s%s_diseases_%s",
        "ITEM_NAME":"Disease",
        "ITEM_FIELD_NAME":"SHORT_TITLE",
        "NAME_DB":{"DB_path":\
            [join(concat_clamp_prefix,"*CUI_SNOMED_NAME.csv"), join(term_path,"MRCONSO_ICD9CM_AB*")], \
                "key_value":[["CUI","NE"],["CUI","STR"]]},
        "CODE_NAME_FILE":"disease_cui_name"

    }

    sider_dict={
        "KEY_FIELDS":["STICH_ID_FLAT","CUI","SIDE_EFFECT_NAME"],
        "CODE_FIELD":"CUI",
        "ITEM_FIELD_NAME":"SIDE_EFFECT_NAME"
    }
    
    args_dict={
        "DISEASE":disease_dict,
        "DRUG":drug_dict,
        "DISEASE_CUI":disease_cui_dict,
        "SIDER4":sider_dict
    }
    return args_dict[item]