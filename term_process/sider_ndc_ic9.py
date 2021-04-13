
# from astropy.table import Table
import os, sys
sys.path.append(os.path.abspath('/home/liu/project/adr-mining'))
from utils._tools import read_data, write2file, append_csv_byrow, left_join
from utils._path import read_prefix, product_ndc_file_path_prefix, processed_map_prefix
import pubchempy as pcp
from os.path import join
import itertools
import pandas as pd

# t=Table.read('/data/liu/mimic3/MAPPING/RAR/ndctext/product.txt')
# t.iteratorstyle = 'dict';
# resolved={}
# unresolved={}
# for row in t:
#     ndc = row['PRODUCTNDC']   
#     print(ndc,'\t',end='')
#     d={}
#     for s in [w.strip() for w in row['SUBSTANCENAME'].split(';')]:
#         print(s)
#         quit()
#         if s in resolved:
#             print(s)
#             d[s] = resolved[s]          
#         elif s in unresolved:
#             print(s)
#             print('failed to resolve substance name',s,file=sys.stderr)
#             unresolved[s] = True

def _zero_pad(x, n=9):
    if len(x) < n:
        x = "1"+(n - len(x)-1) * "0" + x
    return x

def get_cid_by_compound(compound):
    try:
        results = pcp.get_compounds(compound, 'name') 
    except:
        print("%s,"%compound)
        return None
    if results:
        cid_list= ["CID"+_zero_pad(str(compund_class.cid)) for compund_class in results]
        print("%s, %s"%(compound, cid_list))
        return cid_list
    else:
        print("%s,"%compound)
        return None

def get_cid_by_compound_list(compound_list):
    cid_results=filter(None, list(map(get_cid_by_compound, compound_list)))
    if cid_results:
        return list(set(list(itertools.chain(*cid_results))))
    else:
        return cid_results

def get_ndcprefix_pubchemid():
    # product_ndc_file_path_prefix='/data/liu/mimic3/MAPPING/RAR/ndctext'
    product_ndc_field='PRODUCTNDC'
    compound_field='NONPROPRIETARYNAME'
    # compound_field='SUBSTANCENAME'
    product_df=read_data(join(product_ndc_file_path_prefix,'product.txt'),\
        sep='\t',usecols=[product_ndc_field, compound_field],dtype=str)
    # product_df=product_df.assign(NDC_Length=product_df['PRODUCTNDC'].apply(len))
    # print(product_df['NDC_Length'].unique())
    # print(product_df.query('NDC_Length==10').head())
    # print(product_df.columns)
    # print(product_df['NONPROPRIETARYNAME'].head())
    # quit()

    # NOTE: add header to csv file
    append_csv_byrow(
        ["PRODUCTNDC","NONPROPRIETARYNAME","CID_LIST"],
        join(product_ndc_file_path_prefix ,"product_ndc_cid_rbr")
    )
    product_df.apply(
        lambda row: append_csv_byrow(
            [
                row[product_ndc_field],
                row[compound_field],
                str(get_cid_by_compound(row[compound_field]))
            ],join(product_ndc_file_path_prefix ,"product_ndc_cid_rbr")), 
        axis=1
    )  

    #--------------------------------------------------------------------------------------------------------
    # product_df=product_df.assign(CID_LIST=product_df[compound_field].apply(get_cid_by_compound))
    # write2file(product_df,join(product_ndc_file_path_prefix ,"product_ndc_cid"))
    #--------------------------------------------------------------------------------------------------------
    # for substance in product_df['NONPROPRIETARYNAME'].head():
    #     results = pcp.get_compounds(substance, 'name')
    #     print(results[0].cid)

def get_ndcprefix_pubchemid_sep_compd():
    # product_ndc_file_path_prefix='/data/liu/mimic3/MAPPING/RAR/ndctext'
    product_ndc_field='PRODUCTNDC'
    compound_field='NONPROPRIETARYNAME'

    product_df=read_data(join(product_ndc_file_path_prefix,'product.txt'),\
        sep='\t',usecols=[product_ndc_field, compound_field],dtype=str)

    product_df[compound_field]=product_df[compound_field].apply(
        lambda x: [s.strip() for s in str(x).split(",")])
    
    # NOTE: add header to csv file
    append_csv_byrow(
        ["PRODUCTNDC","NONPROPRIETARYNAME","CID_LIST"],
        join(product_ndc_file_path_prefix ,"product_ndc_cid_separate_compd_rbr")
    )

    product_df.apply(
        lambda row: append_csv_byrow(
            [
                row[product_ndc_field],
                str(row[compound_field]),
                str(get_cid_by_compound_list(row[compound_field]))
            ],join(product_ndc_file_path_prefix ,"product_ndc_cid_separate_compd_rbr")), 
        axis=1
    )

    # product_df.apply(
    #     lambda row: print(str.join(',',
    #         [
    #             row[product_ndc_field],
    #             str(row[compound_field]),
    #             str(get_cid_by_compound_list(row[compound_field]))
    #         ])),axis=1
    #         # ,join(product_ndc_file_path_prefix ,"product_ndc_cid_separate_compd_rbr"))

    # )



        # lambda row: append_csv_byrow((", ").join(
            # [row[product_ndc_field], str(row[compound_field]), str(get_cid_by_compound_list(row[compound_field]))],\
                # join(product_ndc_file_path_prefix ,"product_ndc_cid_separate_compd_rbr"))), \
            # axis=1
    # )

    # product_df=product_df.assign(CID_LIST=product_df[compound_field].apply(get_cid_by_compound_list))
         
    # write2file(product_df,join(product_ndc_file_path_prefix ,"product_ndc_cid_separate_compd"))
def ndc_format(raw_ndc):
    return str([len(part) for part in raw_ndc.split("-")])

def ndc_normalization(raw_ndc):
    if(ndc_format(raw_ndc)=='[4, 4, 2]'):
        return ("0"+raw_ndc).replace("-","")
    elif(ndc_format(raw_ndc)=='[5, 3, 2]'):
        parts=raw_ndc.split("-")
        return parts[0]+"0"+("").join(parts[1:])
    elif(ndc_format(raw_ndc)=='[5, 4, 1]'):
        parts=raw_ndc.split("-")
        return ("").join(parts[:2])+"0"+parts[-1]
        
def isOneToOne(df, col1, col2):
    first = df.groupby(col1)[col2].count().max()==1
    if(not first):
        first_df=df.groupby(col1)[col2].count()\
            .reset_index(name='count').sort_values(['count'], ascending=False)
        print("First 20 %s which have more than 2 values in %s"%(col1, col2))
        print(first_df.head(20))
        print(len(first_df[first_df['count']==2]))
    second = df.groupby(col2)[col1].count().max()==1
    if(not second):
        second_df=df.groupby(col2)[col1].count()\
            .reset_index(name='count').sort_values(['count'], ascending=False)
        print("First 20 %s which have more than 2 values in %s"%(col2, col1))
        print(second_df.head(20))
        print(len(second_df[second_df['count']==2]))
    print("%s only have one value in %s : %s"%(col1, col2, "True" if first else "False"))
    print("%s only have one value in %s : %s"%(col2, col1, "True" if second else "False"))

def save_ndc_from_mimic_2mappingdta():
    # NOTE:step1
    pres_df_ndc=read_data(join(read_prefix,"PRESCRIPTIONS"),usecols=["NDC","DRUG","DRUG_NAME_GENERIC"],dtype=str).drop_duplicates()
    pres_df_ndc=pres_df_ndc.query('NDC != "0"')
    # write2file(pres_df_ndc,join(product_ndc_file_path_prefix,"MIMIC_NDC"))

    # NOTE: get cleaned rxnorm_ndc dataframe
    # rxnorm_ndc_df=read_data(
    #     join(processed_map_prefix,"rxnorm_ndc_df"), usecols=["CUI_CODE","CODE"],dtype=str
    # ).drop_duplicates()
    # rxnorm_ndc_df.columns=["RXNORM","NDC"]
    # write2file(rxnorm_ndc_df,join(processed_map_prefix,"only_rxnorm_ndc_df"))
    # print(pres_df_ndc.head())
    #------------------------------------------------------------------------------------------

    rxnorm_ndc_df=read_data(
        join(processed_map_prefix,"only_rxnorm_ndc_df"),dtype=str
    )
    pres_df_rxnorm=left_join(pres_df_ndc[["NDC"]].dropna().drop_duplicates(), rxnorm_ndc_df, "NDC")
    isOneToOne(pres_df_rxnorm,"NDC","RXNORM")
    print(len(pres_df_rxnorm.dropna())/len(pres_df_rxnorm))
    write2file(pres_df_rxnorm, join(product_ndc_file_path_prefix, "MIMIC_NDC_RXNORM"))

def process_package_txt(filename="package.txt"):
    # NOTE:step2
    package_df = read_data(
        join(product_ndc_file_path_prefix,filename),sep="\t",dtype=str,\
            usecols=["PRODUCTNDC","NDCPACKAGECODE"]
    )
    package_df=package_df.assign(
        NDC=package_df["NDCPACKAGECODE"].apply(ndc_normalization)
    ).drop_duplicates()
    write2file(package_df,join(product_ndc_file_path_prefix,"%s_ndc_norm"%(filename.split('.')[0])))

    package_df["PACKAGE_VERSION"]=package_df["NDCPACKAGECODE"].apply(
        lambda x: str([len(part) for part in x.split("-")])
    )
    # print(package_df["PACKAGE_VERSION"].unique())

    print(package_df.query('PACKAGE_VERSION=="[2]"').head(10))

def compare_mimic_package_ndc():
    mimic_ndc=read_data(join(product_ndc_file_path_prefix,"MIMIC_NDC"),dtype=str,usecols=["NDC"]).drop_duplicates()
    package_ndc=pd.concat(
        [read_data(join(product_ndc_file_path_prefix,filename),dtype=str) for filename in
        ["package_ndc_norm", "unfinished_package_ndc_norm"]], axis=0
    ).drop_duplicates()
    # print(package_ndc.head())
    # quit()
    # TODO: only count NDC number
    print(len(mimic_ndc))
    # 4204
    mimic_package_map=left_join(mimic_ndc,package_ndc,"NDC")
    # print(mimic_package_map.query('PRODUCTNDC.isnull()').head(20))
    print(len(mimic_package_map.dropna()))
    # 1282
    # 1286(package+unfinished_package)

#  NOTE: STEP1: ndc -> ndcprefix -> pubchem
def mimic_ndc_to_pubchem():
    if(not os.path.exists(join(product_ndc_file_path_prefix,"all_package_ndc_norm"))):
        # process_package_txt()
        # process_package_txt("unfinished_package.txt")
        package_ndc=pd.concat(
            [read_data(join(product_ndc_file_path_prefix,filename),dtype=str) for filename in
            ["package_ndc_norm", "unfinished_package_ndc_norm"]], axis=0
        ).drop_duplicates()
        write2file(package_ndc, join(product_ndc_file_path_prefix, "all_package_ndc_norm"))
    else:
        package_ndc=read_data(join(product_ndc_file_path_prefix,"all_package_ndc_norm"),dtype=str)
    print("package_ndc %d"%(len(package_ndc)))
    mimic_ndc=read_data(
        join(product_ndc_file_path_prefix,"MIMIC_NDC"),dtype=str,usecols=["NDC"]
    ).dropna().drop_duplicates()
    print("mimic_ndc %d"%(len(mimic_ndc)))

    ndcpre_pubchem=read_data(
        join(product_ndc_file_path_prefix,"product_ndc_cid_separate_compd_rbr"),sep="\t",dtype=str
    )
    print("ndcpre_pubchem %d"%(len(ndcpre_pubchem)))

    mimic_package_map=left_join(mimic_ndc,package_ndc,"NDC")
    print("mimic_package_map first %d"%(len(mimic_package_map)))
    mimic_package_map=left_join(mimic_package_map,ndcpre_pubchem,"PRODUCTNDC")
    print("mimic_package_map second %d"%(len(mimic_package_map)))
    print("mimic_package_map not null %d"%len(mimic_package_map.dropna(subset=["PRODUCTNDC"])))
    print(mimic_package_map.head(20))
    write2file(mimic_package_map,join(product_ndc_file_path_prefix,"MIMIC_NDC_PREFIX_PUBCHEM"))




# process_package_txt("unfinished_package.txt")


# get_ndcprefix_pubchemid()

# get_ndcprefix_pubchemid_sep_compd()
# print((",").join([str([1,23]),"32442",str(["2433,234"])]))
# save_ndc_from_mimic_2mappingdta()
# process_package_txt()
# process_package_txt("unfinished_package.txt")
# compare_mimic_package_ndc()

# from pharmpy.atc import ATCEngine
# ae = ATCEngine(root_url="https://rxnav.nlm.nih.gov/REST")
# res = ae.get_atc("50090347201")
# print(json.dumps(res, indent=2))

# save_ndc_from_mimic_2mappingdta()
# drug_atc=pd.read_csv(join("/data/liu/mimic3/SIDER","drug_atc.tsv"),sep="\t",dtype=str,header=None,\
#     names=["CID","ATC"])
# isOneToOne(drug_atc,"CID","ATC")

# results=pcp.get_compounds("Sodium Heparin", 'name') 
# print(results)
# compare_mimic_package_ndc()

mimic_ndc_to_pubchem()