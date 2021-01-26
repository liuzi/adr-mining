
# from astropy.table import Table
import os, sys
sys.path.append(os.path.abspath('/home/liu/project/adr-mining'))
from utils._tools import read_data, write2file
import pubchempy as pcp
from os.path import join

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

def get_ndcprefix_pubchemid():
    product_ndc_file_path_prefix='/data/liu/mimic3/MAPPING/RAR/ndctext'
    product_ndc_field='PRODUCTNDC'
    compound_field='NONPROPRIETARYNAME'
    # compound_field='SUBSTANCENAME'
    product_df=read_data(join(product_ndc_file_path_prefix,'product.txt'),\
        sep='\t',usecols=[product_ndc_field, compound_field]).head(20)
    # product_df=product_df.assign(NDC_Length=product_df['PRODUCTNDC'].apply(len))
    # print(product_df['NDC_Length'].unique())
    # print(product_df.query('NDC_Length==10').head())
    # print(product_df.columns)
    # print(product_df['NONPROPRIETARYNAME'].head())
    # quit()
    product_df=product_df.assign(CID_LIST=product_df[compound_field].apply(get_cid_by_compound))
    # print(product_df)
    write2file(product_df,join(product_ndc_file_path_prefix ,"product_ndc_cid"))
    # for substance in product_df['NONPROPRIETARYNAME'].head():
    #     results = pcp.get_compounds(substance, 'name')
    #     print(results[0].cid)


get_ndcprefix_pubchemid()