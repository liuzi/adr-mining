import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import os
from os.path import join
from sklearn.decomposition import LatentDirichletAllocation
from scipy.spatial import distance
from scipy.stats import chi2_contingency, ks_2samp
import ast
from utils._path import singledrug_featurepreprocess_prefix, concat_clamp_prefix
from utils._tools import read_data, left_join, write2file, inner_join, append_csv_byrow, append_csv_bydf
import subprocess


import seaborn as sns
import matplotlib.pyplot as plt

# from tools import *

# NOTE: common paths for LDA analysis
# read_prefix= "/data/MIMIC3/"
# write_prefix = "/data/liu/LDA"
# res_prefix = "/data/liu/LDA/lda_result"
# res_r_prefix = "/data/liu/LDA/lda_R_result/"
joint_lda_prefix="/data/liu/mimic3/LDA_MODEL/JOINT_LDA"
charac_prefix="/data/liu/mimic3/CLAMP_NER/single_drug_analysis/FEATURE/PRE_PROCESS"
n2c2_prefix="/data/liu/mimic3/N2C2/"


### Option: Get top columns of original matrix, only run for one time
## Sparsity metric for matrixes
def get_sparsity(matrix):
    sparsity = 1.0 - np.count_nonzero(matrix) / matrix.size
    return sparsity
## Sparsity metric for columns

def col_sparse(matrix):
    width = matrix.shape[1]-1
    sparse_list = []
    for col in range(0,width):
        sparse_list = sparse_list + [get_sparsity(matrix.iloc[:,col])]
    return sparse_list

def df_top_col(matrix, sparse_list, top_drug_list,top_percentage=0.2):
    m_len = len(sparse_list)
    top_N = int(m_len*top_percentage)
    filter_list = sorted(range(m_len), key=lambda i: sparse_list[i])[:top_N]
    final_list = list(set(filter_list).union(set(remain)))
    return matrix.iloc[:, final_list]

def load_diag_and_pres_matrix():
    diag_matrix, pres_matrix = [
        read_data(join(singledrug_featurepreprocess_prefix,file_name)).head(20) for file_name\
            in ["diag_matrix","pres_rxnorm_matrix"]
    ]
    return diag_matrix, pres_matrix

def get_top_drug_list(top_num=10):
    # TODO:
        ## Only run for the first time for geting top N NDCs in SIDER4
    # top10_NDC_ade = ade_df.groupby('NDC')['ICD9_CODE'].count().reset_index(name='count').sort_values(['count'], ascending=False).head(10)
    # top10_NDCs = top10_NDC_ade['NDC']
    # write2file(top10_NDC_ade,join(write_prefix,'top10_NDC_ade'))


    ## Get top NDC, top columns should include these NDCs regardless of their frequency
    # top10_NDC_ade = read_data(join(write_prefix,'top10_NDC_ade'),dtype={'NDC':str})
    # top10_NDCs = top10_NDC_ade['NDC']
    # top10_NDC_indexes = get_top_index(pres_matrix,top10_NDCs)
    # top10_NDC_indexes

    top_drug_list=""
    return top_drug_list

def run_lda_py(pres_matrix, diag_matrix, pres_n_comp=5,diag_n_comp=5):
    lda = LatentDirichletAllocation(n_components=pres_n_comp, random_state=2019)
    pres_Z = lda.fit_transform(pres_matrix) 
    pres_Y = lda.components_


    lda_2 = LatentDirichletAllocation(n_components=diag_n_comp, random_state=2019)
#                                      ,doc_topic_prior=0.001,topic_word_prior=0.00001)
    diag_Z = lda_2.fit_transform(diag_matrix) 
    diag_Y = lda_2.components_
    
    return pres_Z, pres_Y, diag_Z, diag_Y

def get_filterd_diag_and_pres_matrix(diag_matrix, pres_matrix, top_percentage=0.2):
    ## Calculating sparsity
    print("Sparsity of diagnosis matrix: %f"%get_sparsity(diag_matrix))
    print("Sparsity of prescriptions matrix: %f"%get_sparsity(pres_matrix))
    diag_sparse = col_sparse(diag_matrix)
    pres_sparse = col_sparse(pres_matrix)

    ## Only get columns which are top N least sparse
    top_drug_list=get_top_drug_list()

    diag_matrix_filtered = df_top_col(diag_matrix,diag_sparse,[])
    pres_matrix_filtered = df_top_col(pres_matrix,pres_sparse,top_drug_list)
    return diag_matrix_filtered ,pres_matrix_filtered

def run_rlda():
    subprocess.call ("/home/liu/anaconda3/bin/Rscript --vanilla ./rlda_model.r", shell=True)


def get_ditance(pres_theta, diag_theta):   
    pres_n_comp=pres_theta.shape[1]
    diag_n_comp=diag_theta.shape[1] 

    # dis_list = []
    l = []
    for i in range(pres_n_comp):
        # for j in range(diag_n_comp):
        #     sim_i_j = distance.euclidean(pres_Z[:,i],diag_Z[:,j])
        #     l = l + [sim_i_j]
        for j in range(diag_n_comp):
            sim_i_j = distance.cosine(pres_theta.iloc[:,i],diag_theta.iloc[:,j])
            l = l + [sim_i_j]

        # for j in range(diag_n_comp):
        #     sim_i_j = distance.correlation(pres_Z[:,i],diag_Z[:,j])
        #     l = l + [sim_i_j]            
        # dis_list = dis_list + l

    dis_array = np.reshape(l, [pres_n_comp,-1])
    # dis_3array =np.array_split(dis_array,3,axis=1)
    columns = ["DIAG_CLUSTER_%d"%(n+1) for n in list(range(diag_n_comp))]
    index = ["PRES_CLUSTER_%d"%(n+1) for n in list(range(pres_n_comp))]
    
    # dis_3df = [sum_dis(pd.DataFrame(array,columns=columns,index=index)) for array in dis_3array]

    return pd.DataFrame(dis_array,columns=columns,index=index)

# TODO: Get statistics for 10 clusters:

def load_demographics():
    # --------------------------load characteristics dta-----------------------------
    patient_epis_df=read_data(
        join(singledrug_featurepreprocess_prefix,"patient_epis_map"),dtype=str
    )
    # print(patient_epis_df.head())
    demographic_df=read_data(join(charac_prefix,"demographic_matrix"),dtype={"SUBJECT_ID":str})
    epis_demo_df=left_join(demographic_df,patient_epis_df,"SUBJECT_ID")
    epis_demo_df=epis_demo_df.query("AGE >= 0")
    # epis_demo_df=epis_demo_df.drop('SUBJECT_ID', 1)
    # print(epis_demo_df.head())
    # --------------------------load characteristics dta-----------------------------
    return epis_demo_df

def plot_death_age(theta_patient_label, demographic_df, write_prefix, data_name="DIAG"):
    theta_patient_label_demographics = left_join(
        theta_patient_label, demographic_df, "HADM_ID"
    )
    theta_patient_label_demographics.rename(
        columns={'EXPIRE_FLAG': 'DEATH'}, inplace=True)

    g = sns.FacetGrid(theta_patient_label_demographics, col='DEATH')
    g.map(plt.hist, 'AGE', bins=20)
    g.fig.savefig(join(write_prefix,"%s_AGE_DEATH_ALL.png"%data_name))

    grid = sns.FacetGrid(theta_patient_label_demographics, col='DEATH', row='LABEL', size=2.2, aspect=1.6)
    grid.map(plt.hist, 'AGE', alpha=.5, bins=20)
    grid.add_legend()
    grid.fig.savefig(join(joint_lda_prefix,write_prefix,"%s_AGE_DEATH_EACH_CLUSTER.png"%data_name))


def import_theta(args="ngib700_ncomp10_gama0.01_alpha0.01"):
    diag_epis_df=read_data(
        join(singledrug_featurepreprocess_prefix,"diag_matrix"),
        usecols=["HADM_ID"],dtype=str)
    pres_epis_df=read_data(
        join(singledrug_featurepreprocess_prefix,"pres_rxnorm_matrix"),
        usecols=["HADM_ID"],dtype=str)
    # NOTE: Original diagdf and presdf (49955,) (58925,) 
    diag_theta_df=read_data(join(joint_lda_prefix,args,"diag_Full_theta.csv"))
    diag_theta_df.index=list(diag_epis_df["HADM_ID"])

    pres_theta_df=read_data(join(joint_lda_prefix,args,"pres_rxnorm_Full_theta.csv"))
    pres_theta_df.index=list(pres_epis_df["HADM_ID"])
    diag_pres_epis=list(set(diag_epis_df["HADM_ID"]).intersection(set(pres_epis_df["HADM_ID"])))

    pres_theta_df = pres_theta_df.loc[diag_pres_epis,:]
    diag_theta_df = diag_theta_df.loc[diag_pres_epis,:]
    # print(len(diag_pres_epis))
    diag_theta_patient_label= diag_theta_df\
        .idxmax(axis="columns").to_frame(name="LABEL").rename_axis('HADM_ID').reset_index()
    pres_theta_patient_label= pres_theta_df\
        .idxmax(axis="columns").to_frame(name="LABEL").rename_axis('HADM_ID').reset_index()
    
    return pres_theta_df, diag_theta_df, pres_theta_patient_label, diag_theta_patient_label

def get_stats_table(theta_demo):
    fields=["LENGTH_STAY","AGE","GENDER","EXPIRE_FLAG"]
    theta_demo_group = theta_demo.groupby("LABEL")
    # NOTE:discrete distribution
    stats_l=[]
    for field in fields[:2]:
        theta_field = theta_demo_group.apply(
            lambda subdf: (", ").join([
                '{:0.4f}'.format(metric) for metric in \
                    [subdf[field].mean(), subdf[field].std()]]))
        stats_l=stats_l+[theta_field]
    # NOTE:binary distribution
    for field in fields[2:]:
        theta_field = theta_demo_group.apply(
            lambda subdf: (", ").join([str(key_value[1])
                for key_value in sorted(subdf[field].value_counts().items())])
        )
        stats_l=stats_l+[theta_field]

    stats_df=pd.concat(stats_l, axis=1)
    stats_df.columns=["LENGTH_STAY(μ/std)","AGE(μ/std)","GENDER(F/M)","EXPIRE_FLAG(0/1)"]

    return stats_df.reset_index()


def get_stats_of_ldacluster(args="ngib700_ncomp10_gama0.01_alpha0.01"):

# ------------------------------------------------------------
    demographic_df = load_demographics()

    pres_theta_df, diag_theta_df, \
        pres_theta_patient_label, diag_theta_patient_label = import_theta()
    
    pres_theta_demo = inner_join(
        demographic_df, pres_theta_patient_label, "HADM_ID"
    )
    diag_theta_demo = inner_join(
        demographic_df, diag_theta_patient_label, "HADM_ID"
    )
    pres_theta_stats = get_stats_table(pres_theta_demo)
    write2file(pres_theta_stats, join(joint_lda_prefix,args,"pres_theta_stats.csv"))
    diag_theta_stats = get_stats_table(diag_theta_demo)
    write2file(diag_theta_stats, join(joint_lda_prefix,args,"diag_theta_stats.csv"))

#  -------------------------------plot----------------------------
    # plot_death_age(
    #     diag_theta_patient_label, demographic_df,
    #     join(joint_lda_prefix,args),data_name="DIAG")
    # plot_death_age(pres_theta_patient_label, demographic_df, join(joint_lda_prefix,args), data_name="PRES")

    # pres_diag_dis=get_ditance(pres_theta_df, diag_theta_df)
    # write2file(pres_diag_dis,join(joint_lda_prefix,args,"pres_diag_distance.csv"))
    # pres_diag_dis_link=pres_diag_dis.\
    #     idxmin(axis="columns").to_frame(name="DIAG_CLUSTER").rename_axis('PRES_CLUSTER').reset_index()
    # print(pres_diag_dis_link)
    # # TODO: save distance result get prob distribution distance of new diseases drugs
    # write2file(pres_diag_dis_link,join(joint_lda_prefix,args,"pres_diag_link.csv"))

# ------------------------------------------------------------


def get_newitems_of_ldacluster(args="ngib700_ncomp10_gama0.01_alpha0.01"):
    _, _, \
        pres_theta_patient_label, diag_theta_patient_label = import_theta()

    # new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newCUI.csv"),dtype={"HADM_ID":str})
    new_disease_matrix =pd.read_csv(join(concat_clamp_prefix,"allepis_newCUI.csv"),dtype={"HADM_ID":str})
    print(new_disease_matrix.head())
    new_disease_matrix_label = left_join(
        new_disease_matrix, pres_theta_patient_label[["HADM_ID","LABEL"]],"HADM_ID"
    ).dropna(subset=["LABEL"]).iloc[:,1:]
    new_disease_matrix_label=new_disease_matrix_label.groupby("LABEL").apply(
        lambda subdf: subdf.iloc[:,:-1].sum(axis=0)
    )
    gl= []
    pl =[]
    ksl=[]
    kspl=[]
    cluster_index=list(new_disease_matrix_label.index.values)
    for row_cluster in cluster_index:
        # for j in range(diag_n_comp):
        #     sim_i_j = distance.euclidean(pres_Z[:,i],diag_Z[:,j])
        #     l = l + [sim_i_j]
        for column_cluster in cluster_index:

            g, p, _, _ = chi2_contingency(pd.crosstab(
                new_disease_matrix_label.loc[row_cluster,:],
                new_disease_matrix_label.loc[column_cluster,:]))
            kt,kp=ks_2samp(
                new_disease_matrix_label.loc[row_cluster,:],
                new_disease_matrix_label.loc[column_cluster,:]
            )
            gl = gl + [g]
            pl= pl + [p]
            ksl=ksl + [kt]
            kspl=kspl+[kp]

    larray_list=[
            np.reshape(l, [len(cluster_index),-1]) for l in [gl,pl,ksl,kspl]]
    # pl_chi2_array=np.reshape(pl, [len(cluster_index),-1])

    gl_chi2_df, pl_chi2_df, \
        ksl_df, kspl_df=[
            pd.DataFrame(larray,columns=cluster_index,index=cluster_index).reset_index() \
                for larray in larray_list]

    # [write2file(gl_chi2_df,join(joint_lda_prefix,args]
    write2file(gl_chi2_df,join(joint_lda_prefix,args,"new_dis_chi2_g.csv"))
    write2file(pl_chi2_df,join(joint_lda_prefix,args,"new_dis_chi2_pl.csv"))
    write2file(ksl_df,join(joint_lda_prefix,args,"new_dis_ksl.csv"))
    write2file(kspl_df,join(joint_lda_prefix,args,"new_dis_kspl.csv"))

    
def get_n2c2_ades():
    try:
        os.remove(join(n2c2_prefix,"n2c2_ade_drug.csv"))
    except OSError:
        pass
    # append_csv_byrow(["HADM_ID", "ID","ADE","DRUG"],join(n2c2_prefix,"n2c2_ade_drug"))
    test_prefix=join(n2c2_prefix, "training_and_test")
    ann_files=list(filter(
        lambda file:re.match(".*\\.(ann)$",file),\
            os.listdir(test_prefix)))

    for ann_file in ann_files:
    # for ann_file in ["102027.ann"]:
        print(ann_file)
        anns = pd.read_csv(
            join(test_prefix, ann_file), sep='\t', header=None,
            names=["ID","ANNOTATION","NAME"]
        ).dropna(subset=["ANNOTATION"])

        anns_ade = anns[anns['ANNOTATION'].str.contains("ADE-Drug")]

        if(len(anns_ade)):
            anns=anns.set_index("ID")
            anns_ade = anns_ade.assign(ADE_DRUG=anns_ade["ANNOTATION"].apply(
                lambda anno: (",").join([anns["NAME"][arg.split(":")[1]] \
                    for arg in anno.split(" ")[1:]]
            )))
            anns_ade[["ADE","DRUG"]]=anns_ade["ADE_DRUG"].str.split(',', 1, expand=True)
            anns_ade.insert(0, 'HADM_ID',ann_file.split(".")[0])
            append_csv_bydf(anns_ade[["HADM_ID", "ID","ADE","DRUG"]],\
                join(n2c2_prefix,"n2c2_ade_drug"),sep="\t")

def show_n2c2_ades(file_path="/data/liu/mimic3/N2C2/n2c2_ade_drug"):
    _, _, \
        pres_theta_patient_label, diag_theta_patient_label = import_theta()

    # NOTE: N2C2 True ADE-----------------------------------------------------------------
    test_prefix=join(n2c2_prefix, "training_and_test")
    ann_files=list(filter(
        lambda file:re.match(".*\\.(ann)$",file),\
            os.listdir(test_prefix)))
    ann_files=list(map(lambda x:x.split(".")[0], ann_files))
    n2c2_ades=read_data(file_path,sep="\t",dtype=str)[["HADM_ID"]].drop_duplicates()
    n2c2_ades=n2c2_ades.assign(ADE_FLAG=1)
    n2c2_all_episodes=pd.DataFrame({"HADM_ID":ann_files})
    n2c2_all_episodes=left_join(
        n2c2_all_episodes,n2c2_ades,"HADM_ID"
    ).fillna(value={"ADE_FLAG":0})    
    # NOTE: N2C2 True ADE-----------------------------------------------------------------


    # print(n2c2_ades.T.apply(lambda x: x.nunique(), axis=1))
    # ragu_n2c2_pkl=["dict_pat_drug_ade_codes.pkl","dict_pat_drug_ade_codes_only_both_match.pkl",\
    #     "dict_pat_drug_ade_cui.pkl","dict_sider_mimic_n2c2_data_v7.pkl"]
    # for i in range(4)[:1]:
    #     object = pd.read_pickle(\
    #         '/data/ragu/from_home/cdal/dcmtf_heuristic/data/adr_mimic_n2c2/%s'%ragu_n2c2_pkl[i])
    #     # print(object.keys())
    #     print(len(object.keys()))
    # print(object["100035"])


    for theta_patient_label, dataname in zip(
        [pres_theta_patient_label, diag_theta_patient_label],["PRES","DIAG"]
    ):
    # NOTE: N2C2 True ADE-----------------------------------------------------------------
        n2c2_all_label=left_join(
            n2c2_all_episodes,theta_patient_label,"HADM_ID"
        ).dropna(subset=["LABEL"]).sort_values(['LABEL', 'ADE_FLAG'], ascending=[1, 0])
        n2c2_label_ade_count=n2c2_all_label.groupby(["LABEL","ADE_FLAG"])["HADM_ID"].count().reset_index()
        write2file(n2c2_label_ade_count,join(n2c2_prefix,"N2C2_LABEL_ADE_COUNT.csv"))
        print(n2c2_label_ade_count)
        # print(n2c2_ade_label.groupby("LABEL")["HADM_ID"].count().reset_index())
        g = sns.FacetGrid(n2c2_all_label, col='ADE_FLAG')
        g.map(plt.hist, 'LABEL')
        for ax in g.axes.flat:
            for label in ax.get_xticklabels():
                label.set_rotation(40)
        g.fig.savefig(join(n2c2_prefix,"%s_ADE_LABEL.png"%dataname))
    # NOTE: N2C2 True ADE-----------------------------------------------------------------

def drop_zeros(df):
    # NOTE: remove zero rows
    df = df.loc[~(df==0).all(axis=1)]
    # NOTE: remove zero columns
    df = df.loc[:, ~(df==0).all(axis=0)]
    return df

def show_newdis_e9(file_path="/data/liu/mimic3/N2C2/n2c2_ade_drug"):       
    _, _, \
        pres_theta_patient_label, diag_theta_patient_label = import_theta()

    # TODO: New diseases E9codes
    new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newCUI"),dtype={"HADM_ID":str})
    for theta_patient_label, dataname in zip(
        [pres_theta_patient_label, diag_theta_patient_label],["PRES","DIAG"]
    ):
        theta_patient_new_dis = left_join(pres_theta_patient_label,new_disease_matrix,"HADM_ID").dropna()
        theta_patient_new_dis = drop_zeros(theta_patient_new_dis.set_index(["HADM_ID","LABEL"]))
        print(theta_patient_new_dis.head())
    # print(new_disease_matrix.head())
    # quit()

# get_stats_of_ldacluster()
# get_newitems_of_ldacluster()
# get_n2c2_ades()
# show_n2c2_ades()
show_newdis_e9()





