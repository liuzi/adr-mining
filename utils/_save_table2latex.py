import re 
import pandas as pd
import numpy as np
from itertools import chain
from os.path import join
from utils._path import read_prefix, singledrug_prefix
from utils._tools import read_data
import itertools


def template1():
    return r'''\documentclass[preview]{{standalone}}
\usepackage{{booktabs}}
\begin{{document}}
%s
\end{{document}}'''

def template2():
    return r'''\documentclass{article}
\usepackage{amsfonts}
\usepackage{booktabs}
\usepackage{siunitx}
\usepackage[letterpaper,margin=0.5in]{geometry}
\usepackage{titlesec}
\titleformat*{\section}{\fontsize{12}{12}\bfseries}
\titleformat*{\subsection}{\fontsize{10}{10}\bfseries}
\begin{document}
%s
\end{document}'''

def get_df_stats(df,value_index=1):

    values=df.iloc[:,value_index]
    stats=values.describe().map(
        lambda x: round(x,2)
    )

    return pd.DataFrame(stats).transpose()


def each_drug_latex(dfs,drug_vars):
    section_drug_template=r'''\section{Drug rxnorm=%s, drug name=%s}'''
    ## NOTE: top num-> len(per df in dfs)
    sub_unique_section=r'''\subsection{Top %d new diseases that are unique to each cluster}
-'''%(len(dfs[0]))
    sub_section=r'''\subsection{Top %d new diseases}
-'''%(len(dfs[0]))

    section_drug=section_drug_template%(tuple(drug_vars[:2]))

    df_stats_list=list(chain.from_iterable(zip(dfs, list(map(get_df_stats,dfs))))) 

    
    dfs_latex=list(map(lambda df:df.to_latex(index=False), df_stats_list))

    dfs_latex_sum=('\n\n').join(
        [section_drug,sub_unique_section]+dfs_latex[:4]+\
            [sub_section]+dfs_latex[4:])
    # print(dfs_latex_sum)
    # quit()
    return dfs_latex_sum



def save_new_dis_df_as_latextable(rxnorm_dflists,save_path):
    # plt.subplots(nrows,2,sharey=True)
    # filename = join(save_path,'out.tex')

    template=template2()
    drug_latex_list=[]

    for (rxnorm,drug_name,dfs) in rxnorm_dflists:
        # print(dfs)
        drug_vars=[
            rxnorm,drug_name]
            # ("\"{}\"").format(rxnorm), 
            # ("\"{}\"").format(drug_name)]
        each_drug_latex(dfs,drug_vars)
        drug_latex_list.append(
            each_drug_latex(dfs,drug_vars))

    drug_latex=template%(('\n\n').join(drug_latex_list))
    with open(save_path, 'wb') as f:
        f.write(bytes(drug_latex,'UTF-8'))


def stas_orignial_mimic():
    table_files={
        "DIAGNOSES_ICD":"ICD9_CODE",
        "PRESCRIPTIONS":"NDC",
        "LABEVENTS":"ITEMID",
        "PROCEDURES_ICD":"ICD9_CODE",
        "NOTEEVENTS":None,
        "PATIENTS":None,
        "ADMISSIONS":None}
    
    check_field=["HADM_ID","SUBJECT_ID"]
    stats_field=check_field+["ITEM","ITEM_NAME"]
    stats_df_list=[]

    for key, value in table_files.items():

        iter_check_field=check_field.copy()
        if(value):
            iter_check_field.append(value)

        df=read_data(
            join(read_prefix, key),
            dtype=dict(zip(iter_check_field,
            itertools.repeat(str)))).reindex(iter_check_field,axis=1).\
                dropna(axis=1, how='all')
        stats_series=df.nunique()
        if(value):
            stats_series_df=pd.DataFrame(
                {key:stats_series.append(pd.Series({stats_field[-1]:value}))})
            stats_series_df.index=stats_field
        else:
            stats_series_df=pd.DataFrame({key:stats_series})
        stats_df_list.extend([stats_series_df])
        del df

    stats_df=pd.concat(stats_df_list,axis=1,ignore_index=False,sort=False)\
        .transpose().reset_index()
    stats_df.rename(columns={"index":"Orignial Table Name"}, inplace=True)

    with open(join(singledrug_prefix,"mimicstats"), 'wb') as f:
        f.write(bytes(stats_df.to_latex(index=False),'UTF-8'))





# if __name__ == '__main__':
# #     plot_df_as_table([],"",["12341","dfge","that are uiqune to each cluster"])
#     # print((*["12341","dfge","that are uiqune to each cluster"]))
#     stas_orignial_mimic()  


