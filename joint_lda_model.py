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
from term_process.umls_api import search_rxnorm, get_side_effects_cui, search_cid_for_rxnorm


import seaborn as sns
import matplotlib.pyplot as plt

# from tools import *

# NOTE: common paths for LDA analysis
# read_prefix= "/data/MIMIC3/"
# write_prefix = "/data/liu/LDA"
# res_prefix = "/data/liu/LDA/lda_result"
# res_r_prefix = "/data/liu/LDA/lda_R_result/"
def name_for_rxnorm(rxnorm):
    return search_rxnorm(rxnorm).split(" ")[0]

joint_lda_prefix="/data/liu/mimic3/LDA_MODEL/JOINT_LDA"
charac_prefix="/data/liu/mimic3/CLAMP_NER/single_drug_analysis/FEATURE/PRE_PROCESS"
n2c2_prefix="/data/liu/mimic3/N2C2/"
cluster_stats_prefix="/data/liu/mimic3/LDA_MODEL/JOINT_LDA/cluster_stats"
# top_10_drug={'854236': 'CID100000143', '1791232': 'CID100003075', '198148': 'CID100004900', '854239': 'CID100000160', '855346': 'CID105001396', '237205': 'CID100000085', '854232': 'CID100000137','854242': 'CID100000085', '855290': 'CID105001396', '854253': 'CID100000085'}
# 855290
# Coumadin
# 966571
# Hydralazine
# 1659151
# Zosyn
# 854242
# Lovenox
# 198148
# Predisone
# 1791232
# Diltiazem
# 237205
# Nitroglycerine
# top_10_drug={'1791232': 'CID100003075', '855290': 'CID105001396', '854239': 'CID100000160', '854242': 'CID100000085', '237205': 'CID100000085'}
# selected_drugs=[*top_10_drug]
# for drug in selected_drugs:
#     print(drug)
#     print(name_for_rxnorm(drug))
selected_drugs=['855290', '966571', '1659151', '854242', '198148', '1791232','237205'] 
# selected_drugs=['966571', '866924', '1807516', '1658259', '855290', \
#     '1808219', '213169', '1659151', '885257', '836358', '1719291'] 
selected_drug_names=list(map(name_for_rxnorm, selected_drugs))



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


def run_lda_py(pres_matrix, diag_matrix, pres_n_comp=5,diag_n_comp=5):
    lda = LatentDirichletAllocation(n_components=pres_n_comp, random_state=2019)
    pres_Z = lda.fit_transform(pres_matrix) 
    pres_Y = lda.components_


    lda_2 = LatentDirichletAllocation(n_components=diag_n_comp, random_state=2019)
#                                      ,doc_topic_prior=0.001,topic_word_prior=0.00001)
    diag_Z = lda_2.fit_transform(diag_matrix) 
    diag_Y = lda_2.components_
    
    return pres_Z, pres_Y, diag_Z, diag_Y


def run_rlda():
    subprocess.call ("/home/liu/anaconda3/bin/Rscript --vanilla ./rlda_model.r", shell=True)

def get_distance(pres_theta, diag_theta):   
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


def sum_presPhi(pres_Y, n_digits=4):
    
    # columns = ["drug_%d"%(r) for r in list(range(R))]
    index = ["PRES_CLUSTER_%d"%(n) for n in list(range(pres_Y.shape[0]))]
    # index_dict=dict(zip(list(range(pres_Y.shape[1])), index))
    presY_top_sum = pres_Y.loc[:,selected_drugs]
    presY_top_sum = presY_top_sum.assign(
        index=index
    ).set_index("index")
    presY_top_sum['sum'] = np.sum(presY_top_sum,axis=1).round(n_digits)
    presY_top_sum['std'] = np.std(presY_top_sum.iloc[:,:-1],axis=1).round(n_digits)

    presY_top_sum['drug_set'] = [(',').join([selected_drugs[i] for i in np.argwhere(ll>np.mean(ll)).flatten()]) 
                                 for ll in pres_Y.loc[:,selected_drugs].values]
    return presY_top_sum

def get_predict_top_diag(diag_phi, feature_index,topN):
    all_CUIs=diag_phi.columns
    diag_Y = diag_phi.values
    current_feature = diag_Y[feature_index]
    mean = np.mean(current_feature)
    top_diag_indices = current_feature.argsort()[-topN:][::-1]
    select_CUIs = [all_CUIs[i] for i in top_diag_indices]
    return select_CUIs

def import_phi(args="ngib700_ncomp10_gama0.01_alpha0.01"):
    pres_phi_df=read_data(join(joint_lda_prefix,args,"pres_rxnorm_Full_phi.csv"))
    diag_phi_df=read_data(join(joint_lda_prefix,args,"diag_Full_phi.csv"))

    pres_phi_df.columns=[col[1:] for col in pres_phi_df.columns]
    # for i,rxnorm in enumerate(selected_drugs):
    #     if(rxnorm in pres_phi_df.columns):
    #         print(i,rxnorm)
    # print(pres_phi_df.loc[:,selected_drugs])
    # quit()
    return sum_presPhi(pres_phi_df), diag_phi_df


def import_theta(args="ngib700_ncomp10_gama0.01_alpha0.01"):
    if("new" in args):
        diag_epis_df=read_data(
            join(concat_clamp_prefix,"allepis_newCUI"),
            usecols=["HADM_ID"],dtype=str)
        pres_epis_df=read_data(
            join(concat_clamp_prefix,"allepis_newRxNorm"),
            usecols=["HADM_ID"],dtype=str)
    else:
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


def get_validation(folder_path="ngib700_ncomp10_gama0.01_alpha0.01", n_top_disease=30, union=True):
    pres_theta, diag_theta, \
        _, _ = import_theta(folder_path)
    pres_diag_dis=get_distance(pres_theta, diag_theta)
    pres_diag_dis_link=pres_diag_dis.\
        idxmin(axis="columns").to_frame(name="DIAG_CLUSTER").rename_axis('PRES_CLUSTER').reset_index()
    # print(pres_diag_dis_link)
    pre_min_diag_int=[int(s.split('_')[-1]) for s in pres_diag_dis_link["DIAG_CLUSTER"]]
    # print(pre_min_diag_int)
    # quit()
    # print(pres_diag_dis_link["DIAG_CLUSTER"])
    pres_phi_sum, diag_phi = import_phi(folder_path)
    len_true_ade=[]
    l=[]
    drug_set_list=pres_phi_sum["drug_set"].to_list()

    for i in range(pres_phi_sum.shape[0]):
        drugs = drug_set_list[i].split(",")
        predict_disease = get_predict_top_diag(diag_phi,pre_min_diag_int[i],n_top_disease)
        selected_ades = [get_side_effects_cui(rxnorm)["CUI"] for rxnorm in drugs]
        # print(selected_ades)            
        if(union==False):
            #intersection of each drug in a drug set
            actual_CUIs = set.intersection(*[set(selected_ade) for selected_ade in selected_ades])
        else: 
            actual_CUIs = set.union(*[set(selected_ade) for selected_ade in selected_ades])

        len_true_ade=len_true_ade+[len(actual_CUIs)]
        matched_disease = set(predict_disease).intersection(set(actual_CUIs))
        precision = len(matched_disease)/len(predict_disease)*100
        if(len(actual_CUIs)>0):
            recall = (len(matched_disease)/len(actual_CUIs)*100)
        else: recall=0            
        l = l +[precision, recall]

    valid_df = pd.DataFrame(np.reshape(l,[pres_phi_sum.shape[0],-1]),
                            index=["PRES_CLUSTER_%d"% i for i in list(range(pres_phi_sum.shape[0]))],
                            columns=["precision%","recall%"]).round(2)

    valid_df.insert(loc=0, column='closest_diag_c', value=pres_diag_dis_link["DIAG_CLUSTER"].to_list())
    valid_df.insert(loc=1, column='drug_set', value=[\
        (", ").join([name_for_rxnorm(rxnorm) for rxnorm in dd.split(",")]) for dd in drug_set_list])
    valid_df.insert(loc=2, column='len(true_ade)', value=len_true_ade)
    write2file(valid_df.reset_index(), join(joint_lda_prefix,folder_path,"join_lda_valid.csv"))
    # print(valid_df)
    # return valid_df

# get_validation("new_ngib700_ncomp10_gama0.01_alpha0.01")

def get_diag_pres_df(new_flag=True):
    if(new_flag):
        diag_epis_df=read_data(
            join(concat_clamp_prefix,"allepis_newCUI"),
            dtype={"HADM_ID":str})
        pres_epis_df=read_data(
            join(concat_clamp_prefix,"allepis_newRxNorm"),
            dtype={"HADM_ID":str})
    else:
        diag_epis_df=read_data(
            join(singledrug_featurepreprocess_prefix,"diag_matrix"),
            dtype={"HADM_ID":str})
        pres_epis_df=read_data(
            join(singledrug_featurepreprocess_prefix,"pres_rxnorm_matrix"),
            dtype={"HADM_ID":str})
    return [df.set_index("HADM_ID") for df in [diag_epis_df, pres_epis_df]]
    # return diag_epis_df, pres_epis_df

# diag_epis_df, pres_epis_df = get_diag_pres_df(False)
# print(diag_epis_df.shape)
# print(pres_epis_df.shape)
# quit()

# NOTE: Two baseline model. Model one, frequency-based
def frequency_based_model(topn=30, union=True, n_digits=4,args="ngib700_ncomp10_gama0.01_alpha0.01"):
    diag_matrix, pres_full_df = get_diag_pres_df(False)
    pres_phi_sum,_=import_phi(args=args)
    drug_set_list = pres_phi_sum["drug_set"].to_list()

    l=[]
    for drug_string in drug_set_list:
        drugs=drug_string.split(",")
        print("drug set: {}".format(drugs))
        res=[]
        pres_df=pres_full_df.loc[:,drugs]
        pres_df['prod'] = pres_df.prod(axis=1)
        # print(pres_df[pres_df["prod"]==1])
        # print(pres_df.head())
        all_index = pres_df.index
        treat_index = pres_df[pres_df['prod']==1].index
        # print(treat_index[:30])
        # print(pres_df[pres_df['prod']==1].head())

        control_index = list(set(all_index)-set(treat_index))
        # print(control_index[:30])
        # print(treat_index)
        # print(control_index)
        treat_diag=diag_matrix.reindex(treat_index).agg('mean').fillna(0)
        control_diag=diag_matrix.reindex(control_index).agg('mean').fillna(0)
        compare_diag = treat_diag-control_diag
        # print(treat_diag)
        # print(control_diag)
        # print(compare_diag)
        predict_ade = compare_diag[compare_diag>0].sort_values(ascending=False).index
        print("length of predict ade: %d"%(len(predict_ade)))
        # selected_ades = ade_df[ade_df['NDC'].isin(ndc_list)]
        selected_ades = [get_side_effects_cui(rxnorm)["CUI"] for rxnorm in drugs]      
        if(union==False):
            #intersection of each drug in a drug set
            actual_CUIs = set.intersection(*[set(selected_ade) for selected_ade in selected_ades])
        else: 
            actual_CUIs = set.union(*[set(selected_ade) for selected_ade in selected_ades])

        true_ade_len = len(actual_CUIs)
        
        res=res+[(',').join([name_for_rxnorm(rxnorm) for rxnorm in drugs]), len(treat_index), true_ade_len]
        # for topn in topNs:
        predict_ade_top = predict_ade[:topn]
        true_positive = set(predict_ade_top).intersection(actual_CUIs)

        if(len(true_positive)==0):
            res=res+[0,0]
        else:
            precision = round(len(true_positive)/len(predict_ade_top)*100,n_digits)
            recall = round(len(true_positive)/true_ade_len*100,n_digits)
            res=res+[precision, recall]

        l = l + res
    
    l_df = pd.DataFrame(np.reshape(l,[-1,len(res)]),\
                    columns=['drug_set','num(treated_patient)','len(true_ade)']+['precision%','recall%'])
    write2file(l_df, join(joint_lda_prefix,args,"baseline1_valid.csv"))
    return l_df

# l_df = frequency_based_model(args="ngib700_ncomp10_gama0.01_alpha0.01")
# print(l_df)
# res, len_res=frequency_based_model()
# print(res)

def get_top_10_drugs(top_drug_num=20):
# 855290
# Coumadin
# 966571
# Hydralazine
# 1659151
# Zosyn
# 854242
# Lovenox
# 198148
# Predisone
    topn=800
    newpres_epis_df=read_data(
        join(concat_clamp_prefix,"allepis_newRxNorm"),
        dtype={"HADM_ID":str}).set_index("HADM_ID")
    common_drugs=set(selected_drugs).intersection(set(newpres_epis_df.columns))
    newpres_epis_df_sum = newpres_epis_df\
        .sum(axis=0).sort_values(ascending=False)[:topn]\
            .index.values.tolist()
    # print(newpres_epis_df_sum[:800][-1])
    # quit()
    pres_epis_df=read_data(
        join(singledrug_featurepreprocess_prefix,"pres_rxnorm_matrix"),
        dtype={"HADM_ID":str}).set_index("HADM_ID")
    pres_epis_df_sum = pres_epis_df.sum(axis=0)\
        .sort_values(ascending=False)[:topn]\
            .index.values.tolist()
    # print(pres_epis_df_sum[:2000][-1])
    common_drugs=set(newpres_epis_df_sum).intersection(set(pres_epis_df_sum))
    print(len(common_drugs))
    # quit()
    top_10_drug={}
    for drug in common_drugs:
        drug_cid=search_cid_for_rxnorm(drug)
        if(drug_cid):
            top_10_drug[drug]=drug_cid
        if(len(top_10_drug)>top_drug_num):
            break
    print(top_10_drug)

    # print([search_cid_for_rxnorm(drug) for drug in common_drugs])
# get_top_10_drugs()

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

def show_dis_newdis_e9():

    _, _, \
        pres_theta_patient_label, diag_theta_patient_label = import_theta()
    patient_epis_df=read_data(
        join(singledrug_featurepreprocess_prefix,"patient_epis_map"),dtype=str
    ).set_index("HADM_ID")
    # print(patient_epis_df.head())
    # quit()
    # TODO: New diseases E9codes
    # NOTE: read newdis and dis matrix
    # new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newCUI"),dtype={"HADM_ID":str})
    # disease_matrix = read_data(join(singledrug_featurepreprocess_prefix,"diag_matrix"),dtype={"HADM_ID":str})
    disease_icd9_matrix=read_data(join(
        singledrug_featurepreprocess_prefix,"original_code","diag_matrix.csv"),
        dtype={"HADM_ID":str})
    e9_cols=["HADM_ID"]+[col for col in disease_icd9_matrix.columns if "E9" in col]
    # print(disease_icd9_matrix[e9_cols])
    # quit()

    for theta_patient_label, dataname in zip(
        [pres_theta_patient_label, diag_theta_patient_label],["PRES","DIAG"]
    ):
        theta_patient_dis = left_join(
            theta_patient_label,
            disease_icd9_matrix[e9_cols],"HADM_ID").dropna()
        theta_patient_dis = drop_zeros(theta_patient_dis.set_index(["HADM_ID"]))
        theta_dis_df = theta_patient_dis.groupby("LABEL").apply(
            lambda label_df: len(drop_zeros(\
                patient_epis_df.reindex(label_df.index.values.tolist())).drop_duplicates())
        ).reset_index()
        print(theta_dis_df)
        write2file(theta_dis_df,join(cluster_stats_prefix,"%s_episodes_count"%(dataname)))


        # # NOTE: count of diseases and new diseases
        # for dis_matrix, dis_flag in zip([
        #     # new_disease_matrix, 
        #     disease_matrix],[
        #         # "NEW_DIS",
        #         "DIS"]):
        #     theta_patient_dis = left_join(theta_patient_label,dis_matrix,"HADM_ID").dropna()
        #     # NOTE: drop zero row or column
        #     # theta_patient_new_dis = drop_zeros(theta_patient_new_dis.set_index(["HADM_ID","LABEL"]))
        #     theta_patient_dis = drop_zeros(theta_patient_dis.set_index(["HADM_ID"]))
        #     theta_dis_df = theta_patient_dis.groupby("LABEL").apply(
        #         lambda label_df: len(drop_zeros(label_df).columns)
        #     ).reset_index()
        #     write2file(theta_dis_df,join(cluster_stats_prefix,"%s_%s_count"%(dataname,dis_flag)))



# run_rlda()
# get_stats_of_ldacluster()
# get_newitems_of_ldacluster()
# get_n2c2_ades()
# show_n2c2_ades()
show_dis_newdis_e9()





