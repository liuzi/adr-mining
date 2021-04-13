import os, sys
sys.path.append(os.path.abspath('/home/liu/project/adr-mining'))
from utils._tools import read_data, write2file, append_csv_byrow, left_join
# from utils._path import read_prefix, product_ndc_file_path_prefix, processed_map_prefix
# import pubchempy as pcp
from os.path import join
import itertools
import pandas as pd
import threading 
import re
import numpy as np
import urllib.request 
import xmltodict


stitch_prefix="/data/liu/mimic3/STITCH"
kegg_prefix="/data/liu/mimic3/KEGG"

def get_most_common_value(value_list, n_values=5):
    try:
        # common_value=value_list.value_counts().argmax()
        common_value=value_list.value_counts().nlargest(n_values).index.tolist()
    except:
        return None
    return (",").join(common_value)

# NOTE: Get NDC amd drug names from MIMIC, only choose the most common used names for each NDC
def get_ndc_from_mimic():
    ndc_field, drug_field, drug_name_field="NDC", "DRUG", "DRUG_NAME_GENERIC"
    pres_df=read_data(
        join(read_prefix,"PRESCRIPTIONS"),usecols=[ndc_field,drug_field,drug_name_field],dtype=str
    ).drop_duplicates()
    pres_df=pres_df.query('%s != "0"'%ndc_field)

    pres_df_processed=pres_df.groupby([ndc_field]).agg(
        MOST_DRUG=(drug_field, get_most_common_value),
        MOST_DRUG_NAME=(drug_name_field, get_most_common_value)
    ).reset_index()
    print("# of NDCs from MIMIC: %d"%(len(pres_df_processed)))
    pres_df_processed.to_csv(join(stitch_prefix, "MIMICT_NDC_COMMON_DRUG_NAME.csv"), sep="\t",index=False)
    return pres_df_processed

def get_KEGG_code(ndc,drug_name_list,n_kegg=5):
    drug_names=[drug_name.strip() for drug_name in drug_name_list.split(",")]
    for drug_name in drug_names:
        regex_drug_name=("+").join(
            re.findall(r'(?:[^\W\d_]+\d|\d+[^\W\d_])[^\W_]*|[^\W\d_]+', drug_name.lower())
        )
        # kegg_api_url="http://rest.kegg.jp/find/drug/%s"%(regex_drug_name.replace(" ","+"))
        kegg_api_url="http://rest.kegg.jp/find/drug/%s"%(regex_drug_name)

        # with urllib.request.urlopen(kegg_api_url) as url:
            # data = url.read()
        result_df=pd.read_csv(kegg_api_url, header=None, names=["KEGG", "DETAILED_NAME"], sep=("\t"))
        if(result_df.empty):
            print("No KEGG codes found for drug: %s-%s"%(ndc,drug_name))
            continue
        else:
            # print(result_df)
            result_df["KEGG"]=result_df["KEGG"].apply(
            # NOTE: remove prefix "dr:" for KEGG codes
                lambda x: x[3:]
            )
            matched_kegg=(",").join(result_df["KEGG"][:n_kegg].tolist())
            drug_kegg_name=(",").join(result_df["DETAILED_NAME"][:n_kegg].tolist())
            print("Found KEGG codes %s for drug: %s-%s"%(matched_kegg,ndc,drug_name))
            row=[ndc, regex_drug_name, drug_kegg_name, matched_kegg]
            append_csv_byrow(row, join(stitch_prefix ,"MI_NDC_DRUG_KEGG"), sep="\t")
            break

def search_KEGG_for_smalldf(df):
    # NOTE: df columns: NDC, DRUGNAME
    return df.apply(lambda x: get_KEGG_code(*x), axis=1)

def process_name_KEGG_searching_forall():
    drug_name_field="MOST_DRUG"
    if(os.path.exists(join(stitch_prefix,"MIMICT_NDC_COMMON_DRUG_NAME.csv"))):
        mimic_ndc_name=read_data(
            join(stitch_prefix, "MIMICT_NDC_COMMON_DRUG_NAME"), sep="\t", usecols=["NDC",drug_name_field], dtype=str
        )
    else:
        mimic_ndc_name=get_ndc_from_mimic()

    if(os.path.exists(join(stitch_prefix ,"MI_NDC_DRUG_KEGG.csv"))):
        os.remove(join(stitch_prefix ,"MI_NDC_DRUG_KEGG.csv"))

    append_csv_byrow(
        ["NDC",drug_name_field,"KEGG_NAME","KEGG"], join(stitch_prefix ,"MI_NDC_DRUG_KEGG"), sep="\t"
    )

    batch_size=20
    n_groups=int(len(mimic_ndc_name)/batch_size)
    n_threads=2
    frames=[mimic_ndc_name.iloc[i*batch_size:(i+1)*batch_size].copy() for i in range(n_groups+1)]
    # check_ndcs=[list(frame["NDC"].unique()) for frame in frames]
    # merged_ndc= list(itertools.chain(*check_ndcs))
    # print(len(merged_ndc))
    for frame in frames:
        thread_frame=np.array_split(frame, n_threads)
        # [print(thread_frame_single["NDC"]) for thread_frame_single in thread_frame]
        # quit()
        thread_list=[]
        for i in range(n_threads):
            thread_list=thread_list+[
                threading.Thread(target=search_KEGG_for_smalldf, args=(thread_frame[i],))
            ]
        [thread.start() for thread in thread_list]
        [thread.join() for thread in thread_list]
        del thread_list

def fix_missed_ndc_due_tothread(mapped_file=join(stitch_prefix ,"MI_NDC_DRUG_KEGG.csv"),append_flag=True):
    mapped_ndcs=read_data(
        mapped_file, sep="\t",dtype=str,usecols=["NDC"]
    ).drop_duplicates()
    print("# of NDCs are mapped before revised: %d"%(len(mapped_ndcs)))
    mimic_ndc_name=read_data(
        join(stitch_prefix, "MIMICT_NDC_COMMON_DRUG_NAME"), sep="\t", usecols=["NDC","MOST_DRUG"], dtype=str
    ).drop_duplicates()
    unmapped_ndcs=set(mimic_ndc_name["NDC"])-set(mapped_ndcs["NDC"])
    print("# of unmmaped NDCs which are re-processed: %d"%(len(unmapped_ndcs)))
    unmapped_ndc_name=left_join(
        pd.DataFrame({"NDC":list(unmapped_ndcs)}),\
            mimic_ndc_name,"NDC")
    # print(len(unmapped_ndc_name))
    if(append_flag):
        unmapped_ndc_name.apply(
            lambda x: get_KEGG_code(*x), axis=1
        )




def get_pubchem_for_KEGGs(ndc, drug_name, kegg_list, kegg_cid_map, write_file_path):
    for kegg in kegg_list.split(","):
        cid=kegg_cid_map.get(kegg,None)
        if(cid):
            row=[ndc, kegg, cid, drug_name]
            append_csv_byrow(row, write_file_path, sep="\t")
            break
        else:
            continue

def map_KEGG_to_pubchem(
    ndc_drug_kegg_file=join(stitch_prefix,"MI_NDC_DRUG_KEGG"),
    write_file_path=join(stitch_prefix ,"MI_NDC_DRUG_KEGG_PUCHEM.csv")
):
    drug_name_field="MOST_DRUG"
    # ndc_drug_kegg_df=read_data(
    #     join(stitch_prefix,"MI_NDC_DRUG_KEGG"), dtype=str, sep="\t"
    # )
    ndc_drug_kegg_cols=["NDC", "MOST_DRUG", "KEGG"]
    ndc_drug_kegg_df=read_data(
        ndc_drug_kegg_file, dtype=str, sep="\t"
    )[ndc_drug_kegg_cols]
    print(ndc_drug_kegg_df.head())

    kegg_pubchem_df= pd.read_csv(
        join(stitch_prefix, "chemical_KEGG.tsv"), dtype=str, \
            sep="\t", header=None, names=["CID","SID","no","KEGG"]
    ).iloc[1:,[0,3]]

    kegg_pubchem_df["CID"]=kegg_pubchem_df["CID"].apply(
        lambda x: x.replace("m","1")
    )

    kegg_pubchem_map=kegg_pubchem_df.set_index("KEGG")["CID"].to_dict()

    if(os.path.exists(write_file_path)):
        os.remove(write_file_path)

    append_csv_byrow(
        ["NDC","KEGG","CID",drug_name_field], \
            write_file_path, sep="\t"
    )

    ndc_drug_kegg_df.apply(
        lambda x: get_pubchem_for_KEGGs(*x,kegg_pubchem_map,write_file_path), axis=1
    )

##--------------------------------------------------------
# HACK: stitch
# get_ndc_from_mimic()
##test
# get_KEGG_code("00002735501", "Vancomycin HCl,Vancomycin Oral Liquid,Vancomycin")
# process_name_KEGG_searching_forall()
# fix_missed_ndc_due_tothread(False)
# map_KEGG_to_pubchem()
# fix_missed_ndc_due_tothread(join(stitch_prefix,"MI_NDC_DRUG_KEGG_PUCHEM.csv"),False)
# test=["fds","fds"]
# print("%s%s%s"%(*test,"fds"))

##--------------------------------------------------------
# HACK: kegg

def ndc_format(raw_ndc):
    return str([len(part) for part in raw_ndc.split("-")])

def ndc_normalization(raw_ndc):
    ## Thus all ndc prefixes are of 9 length
    if(ndc_format(raw_ndc)=='[4, 4]'):
        return ("0"+raw_ndc).replace("-","")
    elif(ndc_format(raw_ndc)=='[5, 3]'):
        parts=raw_ndc.split("-")
        return parts[0]+"0"+("").join(parts[1:])
    elif(ndc_format(raw_ndc)=='[5, 4]'):
        parts=raw_ndc.split("-")
        return ("").join(parts[:2])

def process_ndc_kegg():
    ndc_prefix_field, kegg_field="NDC_PREFIX", "KEGG"

    ndc_kegg_df=pd.read_csv(
        join(kegg_prefix,"ndc_kegg.tsv"), dtype=str, sep="\t", names=[ndc_prefix_field, kegg_field]
    )
    # NOTE: remove prefix "ndc:"
    ndc_kegg_df[ndc_prefix_field]=ndc_kegg_df[ndc_prefix_field].apply(
        lambda x: ndc_normalization(x[4:])
    )
    # NOTE: remove prefix "dr:"
    ndc_kegg_df[kegg_field]=ndc_kegg_df[kegg_field].apply(lambda x: x[3:])
    ndc_kegg_df.to_csv(
        join(kegg_prefix, "norm_ndc_kegg.tsv"), index=False, sep="\t"
    )
    # print(ndc_kegg_df[ndc_prefix_field].apply(ndc_format).unique())

def append_kegg_to_mimic_ndc():
    ndc_field, drug_field, drug_name_field,ndc_prefix_field=\
        "NDC", "MOST_DRUG", "MOST_DRUG_NAME_GENERIC", "NDC_PREFIX", 
    mimic_ndc_name=read_data(
        join(stitch_prefix, "MIMICT_NDC_COMMON_DRUG_NAME"), \
            sep="\t", usecols=[ndc_field,drug_field], dtype=str
    )
    ndc_kegg_df=read_data(
        join(kegg_prefix, "norm_ndc_kegg.tsv"), sep="\t", dtype=str
    )
    mimic_ndc_name=mimic_ndc_name.assign(
        NDC_PREFIX=mimic_ndc_name[ndc_field].apply(
            lambda x: x[:9]
        )
    )
    mimic_ndc_name_kegg=left_join(
        mimic_ndc_name, ndc_kegg_df, ndc_prefix_field
    ).dropna(subset=["KEGG"])
    mimic_ndc_name_kegg.to_csv(
        join(kegg_prefix, "mimic_ndcprefix_name_kegg.tsv"), index=False, sep="\t"
    )
    
def combine_ndc_pubchem_maps():
    map_file_list=[
        join(stitch_prefix,"MI_NDC_DRUG_KEGG_PUCHEM.csv"),
        join(kegg_prefix, "mimic_ndc_name_kegg_pubchem.tsv")
    ]
    map_df=[
        read_data(
            map_file, dtype=str, sep="\t"
        )[["NDC","KEGG","CID"]].drop_duplicates() for map_file in map_file_list
    ]
    [print("# of Unique NDCs %d"%(len(df["NDC"].unique()))) for df in map_df]
    combined_map_df=pd.concat(
        map_df,axis=0,sort=False
    ).drop_duplicates()
    print("# of Unique NDCs if combining all mapping files: %d"%(len(combined_map_df["NDC"].unique())))
    combined_map_df.to_csv(
        join(kegg_prefix, "combined_ndc_kegg_pubchem.tsv"), index=False, sep="\t"
    )

def get_kegg_stereo_map():
    kegg_df=pd.read_csv(
        join(stitch_prefix,"chemical_KEGG.tsv"), dtype=str, sep="\t",\
            header=None, names=["CID","SID","no","KEGG"]
    ).iloc[1:,:]
    kegg_df["SID"]=kegg_df["SID"].apply(
        lambda x: x.replace("s","0")
    )

    kegg_df_map=kegg_df.set_index("KEGG")["SID"].to_dict()
    # kegg_pubchem_df["CID"]=kegg_pubchem_df["CID"].apply(
    #     lambda x: x.replace("m","1")
    # )
    # print(kegg_df.head())
    # print(kegg_df_map["C00002"])
    return kegg_df_map

def add_stereo_pubchem():
    combined_ndc_pubrxnorm_df=read_data(
        join(kegg_prefix,"combined_ndc_kegg_pubchem_rxnorm.tsv"), sep="\t", dtype=str
    )
    kegg_sid_map=get_kegg_stereo_map()
    # print(combined_ndc_pubrxnorm_df.head())
    combined_ndc_pubrxnorm_df=combined_ndc_pubrxnorm_df.assign(
        SID=combined_ndc_pubrxnorm_df["KEGG"].apply(
            lambda x: kegg_sid_map.get(x,None)
        )
    )
    combined_ndc_pubrxnorm_df.to_csv(
        join(kegg_prefix, "combined_ndc_kegg_pubchem_rxnorm_sid.tsv"), index=False, sep="\t"
    )


#-----------------------------------------------------------------

# process_ndc_kegg()
# append_kegg_to_mimic_ndc()
# map_KEGG_to_pubchem(
#     ndc_drug_kegg_file=join(kegg_prefix,"mimic_ndcprefix_name_kegg.tsv"),
#     write_file_path=join(kegg_prefix ,"mimic_ndc_name_kegg_pubchem.tsv")
# )

# combine_ndc_pubchem_maps()
# NOTE: additional process required for TWOSIDES->OFFSIDES
# get_kegg_stereo_map()
add_stereo_pubchem()

#-----------------------------------------------------------------
# NOTE: add rxnorm to combine file
def load_rxnsat_df(file_path):

    # python_path="/home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server/S2/term_process/"
    # mapping_prefix="/data/liu/mimic3/MAPPING/TERMINOLOGIES/RAR/rrf"
    # rxnsat_ndc_file="RXNSAT.RRF"
    rxnsat_ndc_title=["RXCUI", "LUI", "SUI", "RXAUI", "STYPE" ,"CODE" ,"ATUI", "SATUI","ATN","SAB",\
        "ATV","SUPPRESS","CVF","NONE"]
    # rxnsat_ndc_selected_title=["RXCUI", "RXAUI","STYPE","CODE","ATN","SAB","ATV","CVF"]
    rxnsat_ndc_df = pd.read_csv(file_path,\
        sep="|", header=None, names=rxnsat_ndc_title,dtype=str).iloc[:,:-1]
        # [rxnsat_ndc_selected_title]
    return rxnsat_ndc_df

def get_sat_ndc_rxnorm_map(file_path=join(kegg_prefix,"NDC_RXNORM","NDC_RXNORM.RRF")):
    df=load_rxnsat_df(file_path)

    if(df["ATV"].apply(len)[0]>11):
        df["ATV"]=df["ATV"].apply(lambda x: x[1:])
   
    df_map=df.set_index("ATV")["RXCUI"].to_dict()
    return df_map

def add_rxnorm_pubchem(
    precessed_file=join(kegg_prefix, "combined_ndc_kegg_pubchem.tsv"),
    ndc_rx_mapfile="NDC_RXNORM.RRF"
):
    file_prefix=join(kegg_prefix, "NDC_RXNORM")
    # processed_ori_file = join(file_prefix, precessed_file)
    map_file=join(file_prefix, ndc_rx_mapfile)

    ndc_rxnorm_map=get_sat_ndc_rxnorm_map(map_file)
    ndc_pubchem_df=read_data(
        precessed_file , dtype=str, sep="\t"
    )[["NDC","KEGG","CID"]]
    ndc_pubchem_df=ndc_pubchem_df.assign(
        RXNORM=ndc_pubchem_df["NDC"].apply(
            lambda x: ndc_rxnorm_map.get(x, None)
        )
    )
    ndc_pubchem_none_df, ndc_pubchem_rxnorm_df= [
        x for _, x in ndc_pubchem_df.groupby(ndc_pubchem_df['RXNORM'].notnull())
    ]
    ndc_pubchem_rxnorm_df.to_csv(
        join(file_prefix,"combined_ndc_kegg_pubchem_%s.tsv"%(ndc_rx_mapfile.split(".")[0])),
        sep="\t", index=False
    )
    ndc_pubchem_none_df.to_csv(
        join(file_prefix,"combined_ndc_kegg_pubchem_%s.tsv"%("NONE")),
        sep="\t", index=False
    )

    print("# of NDCs mapped to RXNORM: %d"%len(ndc_pubchem_rxnorm_df["NDC"].unique()))
    print("# of NDCs not mapped to RXNORM: %d"%len(ndc_pubchem_none_df["NDC"].unique()))

def search_name_for_ndc(ndc, ndc_name_file=join(kegg_prefix,"MIMICT_NDC_COMMON_DRUG_NAME.csv")):
    # drug_name=search_rxnorm(rxnorm)
    grep_command_prefix="grep -i '{}' {} "
    result=os.popen(grep_command_prefix.format(
        ndc, ndc_name_file)).read()
    if(result):
        return result.split('\t')[1]
    else:
        return None

def get_rxcuid_byrugname_api(ndc):
    # return None
    for name in [name.strip().replace(" ","+") for name in search_name_for_ndc(ndc).split(",")]:

        with urllib.request.urlopen("https://rxnav.nlm.nih.gov/REST/rxcui?name=%s"%name) as url:
            data = url.read()
            data = xmltodict.parse(data)
            rxcui= data['rxnormdata']['idGroup'].get('rxnormId',None)
        if(rxcui):
            return rxcui
            break
        else:
            continue
    return None

def add_rxnorm_using_restapi():
    processed_file_path=join(kegg_prefix,"NDC_RXNORM","combined_ndc_kegg_pubchem_NONE.tsv")
    none_rxnorm_df=read_data(
        processed_file_path,
        dtype=str, sep="\t")
    none_rxnorm_df["RXNORM"]=none_rxnorm_df["NDC"].apply(get_rxcuid_byrugname_api)
    none_rxnorm_df.to_csv(
        processed_file_path, sep="\t", index=None
    )




# # print(get_sat_ndc_rxnorm_map()["51927275500"])
# # NOTE: ndc_rxnorm
# add_rxnorm_pubchem()
# # NOTE: ndc_MMSL, NDC_NDDF, NDC_VANDF
# add_rxnorm_pubchem(join(kegg_prefix,"NDC_RXNORM","combined_ndc_kegg_pubchem_NONE.tsv"),"NDC_MMSL.RRF")
# add_rxnorm_pubchem(join(kegg_prefix,"NDC_RXNORM","combined_ndc_kegg_pubchem_NONE.tsv"),"NDC_NDDF.RRF")
# add_rxnorm_pubchem(join(kegg_prefix,"NDC_RXNORM","combined_ndc_kegg_pubchem_NONE.tsv"),"NDC_VANDF.RRF")
## ---------------additional process add rxnorm using name--------------
# add_rxnorm_using_restapi()


