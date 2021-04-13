import os, sys
import glob
from os.path import join
from sklearn.cluster import KMeans
from collections import Counter
import math
import urllib.request 
import xmltodict
import textwrap
import re
from scipy.stats import wasserstein_distance
import pandas as pd
import matplotlib.pyplot as plt
# from matplotlib import gridspec
import numpy as np
from datetime import datetime
from operator import itemgetter
import seaborn as sns


# from utils._preprocess_mimic import *
# from S2.utils._tools import *
# from S2.utils._path import*
from utils._path import*
from utils._save_table2latex import save_new_dis_df_as_latextable
from utils._tools import *
from utils._get_item_args import get_item_args
from term_process.load_standard_dataset import *

from term_process.umls_api import search_cui, search_rxnorm, get_side_effects_cui
# from pathlib import Path


class feature_creation():
    """ Generate features from different data sources using 
    
    Parameters
    ----------
    
    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    Attributs
    ---------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.
    
    See also
    --------
    
    
    Notes
    -----
    
    
    Examples
    --------

    
    """
    

    

    def __init__(self,rxnorm_id):
        """
        Parameters
        ----------
        drug_id : string, default="197380" 
            RxNorm code (e.g. Atenolol, RxNorm: 197380, NDC:51079075920) for drug. 
            Determine the study object. All analysis are based on
            the selected drug/drugs.
        """



        # self.drugid = drugid
        # self.outputpath = singledrug_prefix%drugid
        # print(self.outputpath)
        ## folders for features
        self.rxnorm_id=rxnorm_id
        # self.sample_epis_file='HADM_ID_SAMPLE_PER_PATIENT'
        self.feature_folder=join(singledrug_prefix,"FEATURE")
        self.preprocess_folder=join(self.feature_folder,"PRE_PROCESS")
        self.drug_result_path=join(singledrug_feature_prefix,("_").join(self.rxnorm_id))
        self.epis_field = "HADM_ID"
        self.item_dict={
            "_drugs_":"DRUG",
            "_new_drugs_":"DRUG",
            "_diseases_":"DISEASE",
            "_new_diseases":"DISEASE_CUI"}
        # self.hadm_sampled = None
    
    # HACK: FEATURE0
    def sampling(self, df, size, sample_group="SUBJECT_ID", sample_unit="HADM_ID"):
        """Get randomly selected samples from each subgroup across all 

        Args:
            df (DataFrame): patient-episode record
            size (int): numbers of samples wanna be selected from each subgroups
            sample_group (str, optional): column name of subgroup. Defaults to "SUBJECT_ID".
            sample_unit (str, optional): column name of sampled subject. Defaults to "HADM_ID".

        Returns:
            [type]: 
        """

        fn = lambda obj: obj.loc[np.random.choice(obj.index, size),:]
        return df.groupby(sample_group, as_index=True).apply(fn)
    
    def get_sampled_epis(self, treated_epis, size=1):
        """Sample episode for each patient and generate the file of HADM IDs. 
        If the file already exists. skip this method
        
        Notes
        -----
        A patient may enter hospital several times so they have multiple episodes (HADM ID) 
        Thus we randomly select one episode for each patient
        """
        
        if os.path.exists(join(singledrug_prefix,"%s.csv"%"pres_subject_hadmid")):
            subject_hadm = read_data(join(singledrug_prefix,"pres_subject_hadmid"),dtype=str)
            # hadm_sampled = read_data(join(singledrug_prefix,self.sample_epis_file), dtype=str)
        else:  
            subject_hadm = read_data(join(read_prefix,'PRESCRIPTIONS'),
                dtype={"SUBJECT_ID":str,"HADM_ID":str})[['SUBJECT_ID',"HADM_ID"]].drop_duplicates()
            write2file(subject_hadm,join(singledrug_prefix,"pres_subject_hadmid"))
                    ## RUN ONLY FOR THE FIRST TIME! ramdonly select one hospital stay for each patient
        # print(len(treated_epis))
        treated_patient_epis=left_join(treated_epis,subject_hadm,"HADM_ID")
        size = 1        # sample size
        ## ramdonly get a sample hadm_id from each patient's record
        hadm_sampled = self.sampling(treated_patient_epis[['SUBJECT_ID','HADM_ID']].drop_duplicates(),size)['HADM_ID']
        return hadm_sampled
#         pres_patient_hadm = pres_ade_df[['SUBJECT_ID','HADM_ID']].drop_duplicates()
#         hadm_sampled = pres_patient_hadm.groupby('SUBJECT_ID', as_index=True).apply(fn)['HADM_ID']
        # self.hadm_sampled = hadm_sampled
        

    # HACK: FEATURE 1
    def df2matrix(self,df):
        cols=df.columns
        print("Original Data:")
        print_patient_stats(df)
        ## HACK: sampling should be after creating matrix
        # df_sampled = left_join(self.hadm_sampled,df,list(self.hadm_sampled.columns))
        # print("Sampled Data:")
        # print_patient_stats(df_sampled)

        if("VALUE" not in cols):
            df["VALUE"]=1
        
        df_matrix = df.pivot_table(
            index=cols[1], columns=cols[2], values="VALUE",fill_value=0).reset_index()

        print("Number of %s: %d"%(cols[2],len(df_matrix.columns)-1))
        return df_matrix

    def mat_noZero(self,mat,axis_num=0):
        if(axis_num):
            ## 1 for drop row
            mat = mat[(mat!=0).any(axis=1)]
        else:
            ## drop column
            mat = mat.loc[:,(mat!=0).any(axis=0)]
        return mat

    ## simple process on diagnosis table

    # TODO: 1ST: PUT AFTER condition: if(PREPROCESS dir not extits):run else:read
    def create_fivedata_repre128(self):
        
        create_folder(self.preprocess_folder)
        read_dtype={"SUBJECT_ID":str,"HADM_ID":str}
        ## pre-process four input data respectively
        #NOTE:1) diagnosis: ICD9_CODE, value=1
        diag_df = read_data(join(
            read_prefix,"DIAGNOSES_ICD"),
            dtype=read_dtype).dropna(subset=["ICD9_CODE"])        
        diag_matrix = self.df2matrix(
            diag_df[[*read_dtype]+["ICD9_CODE"]].drop_duplicates())
        write2file(diag_matrix,join(self.preprocess_folder,"diag_matrix"))
    
   
        # NOTE:2) prescription: NDC, value=1
        pres_df = read_data(
            join(read_prefix,'PRESCRIPTIONS'),
            dtype={**read_dtype,**{"NDC":str}})[[*read_dtype]+["NDC"]].dropna(subset=['NDC'])
        pres_df = pres_df[pres_df['NDC']!="0"].drop_duplicates()
        write2file(self.df2matrix(pres_df) ,join(self.preprocess_folder,"pres_matrix"))

        # HACK:
        # NOTE:3) labevents: ITEMID, randomly selected VALUE
        # get_labmatrix()

        # NOTE:4) procedure: ICD9_CODE, value=1
        procedure_df=read_data(
            join(read_prefix,'PROCEDURES_ICD'),
            dtype=read_dtype).dropna(subset=["ICD9_CODE"])[[*read_dtype]+["ICD9_CODE"]].drop_duplicates() 
        write2file(self.df2matrix(procedure_df) ,join(self.preprocess_folder,"procedure_matrix"))


        # NOTE: 5) demographic: []
        # get_demographic_df()


    
    # HACK:Feature 2



    # def create_dissum_feature(self):
        # TODO: 2ND
        ## import original concat dataframe of clamp result
        # section_titles=['HOSPITAL_COURSE', 'MEDICAL_HISTORY', 'MEDICATIONS_ON_ADMISSION']

        # input_files=glob.glob(os.path.join(
        #     clamp_output_prefix,"CONCAT_Results", "%s*")%section_titles[0])
        # print(input_files)  
        # return 0

    
    # HACK: ALL FEATURES CREATION
    def create_data(self):
        """
        Notes:
        1) First Feature:

        2) Second Feature:
            Drugs: PRESCRIPTIONS.csvm only remain rows with drugs that can be found in SIDER
        3) Third Feature:
        """

        create_folder(self.feature_folder)

        # NOTE: feature0
        ## PATIENT PRESCRIPTION LOG
        pres_df=read_data(join(
            read_prefix,'PRESCRIPTIONS'),dtype={'NDC':str}).dropna(subset=['NDC'])
        ## DRUG-ADE IN SIDER4, !!SIDER HAVE DUPLICATED RECORDS
        ade_df = read_data(
            join(sideffect_prefix, 'ndc_icd9_side_effects'), 
            dtype={'NDC':str,'ICD9_CODE':str},usecols=['NDC','ICD9_CODE']).drop_duplicates()

        ## GET LIST OF DRUGS FROM SIDER4
        ade_drug=ade_df['NDC'].drop_duplicates()
        # NOTE:
        ## Remove records from Prescriptions where drugs cannot be found in SIDER
        pres_ade_df = pres_df[pres_df['NDC'].isin(ade_drug)].drop_duplicates()
        ## Sample only one episode for each patient
        self.get_sampled_epis(pres_ade_df)
        # self.create_pres_ade_feature(pres_ade_df, ade_df)
        # self.create_pres_ade_feature()

        # TODO:
        ## NOTE: feature1
        # self.create_fivedata_repre128()

        # NOTE: feature2
        # self.create_dissum_feature()

    



    
    def runKMeans(self,data,n_clusters):
        km = KMeans(n_clusters=n_clusters).fit(data)
        print(Counter(km.labels_))
        return km.labels_

    def get_ade_rxnorm_df(self):
        if(os.path.exists(join(sideffect_prefix,'rxnorm_icd9_side_effects.csv'))):
            ade_df = read_data(join(sideffect_prefix,'rxnorm_icd9_side_effects'),dtype=str)
        else:
            ade_related_prefix="/data/liu/mimic3/CLAMP_NER/ADE_Related/"
            ndc_cui_map = read_data(join(ade_related_prefix,"ndc_cui_map"),dtype={'CODE':str,'CUI_CODE':str})
            ndc_cui_map = ndc_cui_map.set_index('CODE')['CUI_CODE'].to_dict()
            ade_df = read_data(
                join(sideffect_prefix, 'ndc_icd9_side_effects'), 
                dtype={'NDC':str,'ICD9_CODE':str},usecols=['NDC','ICD9_CODE']).drop_duplicates()
            ade_df['RxNorm']=ade_df['NDC'].apply(lambda x: ndc_cui_map.get(x, None))
            write2file(
                ade_df.dropna(subset=['RxNorm'])[["RxNorm","ICD9_CODE"]],
                join(sideffect_prefix,'rxnorm_icd9_side_effects'))   
        return ade_df

    def compare_sider_newdisease(self,treated_epis,new_icd9_matrix):
        # if(os.path.exists(join(sideffect_prefix,'rxnorm_icd9_side_effects.csv'))):
        #     ade_df = read_data(join(sideffect_prefix,'rxnorm_icd9_side_effects'),dtype=str)
        # else:
        #     ade_related_prefix="/data/liu/mimic3/CLAMP_NER/ADE_Related/"
        #     ndc_cui_map = read_data(join(ade_related_prefix,"ndc_cui_map"),dtype={'CODE':str,'CUI_CODE':str})
        #     ndc_cui_map = ndc_cui_map.set_index('CODE')['CUI_CODE'].to_dict()
        #     ade_df = read_data(
        #         join(sideffect_prefix, 'ndc_icd9_side_effects'), 
        #         dtype={'NDC':str,'ICD9_CODE':str},usecols=['NDC','ICD9_CODE']).drop_duplicates()
        #     ade_df['RxNorm']=ade_df['NDC'].apply(lambda x: ndc_cui_map.get(x, None))
        #     write2file(
        #         ade_df.dropna(subset=['RxNorm'])[["RxNorm","ICD9_CODE"]],
        #         join(sideffect_prefix,'rxnorm_icd9_side_effects'))
        ade_df=get_ade_rxnorm_df()
        
        single_drug_ades=ade_df[ade_df['RxNorm']==self.rxnorm_id]['ICD9_CODE']
        
        # new_icd9_matrix =read_data(join(concat_clamp_prefix,"allepis_newICD9_CODE"),dtype={"HADM_ID":str})
        treated_new_icd9=self.mat_noZero(
            inner_join(treated_epis,new_icd9_matrix,"HADM_ID"), axis_num=0
        )
        newdurgs_sider=self.mat_noZero(
            treated_new_icd9.reindex(columns=["HADM_ID"]+list(single_drug_ades)).dropna(axis=1, how='all'),
            axis_num=1)
        print("newdurgs_sider: %s"%str(newdurgs_sider.shape))
        newdrugs_notin_sider=self.mat_noZero(
            treated_new_icd9.reindex(
                columns=list(set(treated_new_icd9.columns)-set(single_drug_ades))
            ),axis_num=1)
        print("newdurgs_notin_sider: %s"%str(newdrugs_notin_sider.shape))
        item_list_prefix=join(singledrug_prefix,self.rxnorm_id)
        write2file(newdurgs_sider,join(item_list_prefix,"newdurgs_sider"))
        write2file(newdrugs_notin_sider,join(item_list_prefix,"newdurgs_not_sider"))


        
    # def extra_dis_sider_intersection_anlz(self,new_disease_matrix):
    #   new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newCUI"),dtype={"HADM_ID":str})
  

    def run_model_plot(self):
        ##NOTE: 1)DRUGS IN PRES 
        if(not os.path.exists(join(self.preprocess_folder,"pres_rxnorm_matrix.csv"))):
            pres_rxnorm_matrix = self.df2matrix(
                read_data(join(
                    self.preprocess_folder,"pres_rxnorm_df"),dtype=str)[
                        ['SUBJECT_ID','HADM_ID','RxNorm']].drop_duplicates().dropna(subset=["RxNorm"])
            ).fillna(0)
            write2file(pres_rxnorm_matrix,join(self.preprocess_folder,"pres_rxnorm_matrix"))
        else:
            pres_rxnorm_matrix =read_data(join(self.preprocess_folder,"pres_rxnorm_matrix"),dtype={"HADM_ID":str})
        
        # TODO: ICD9->CUI
        ###NOTE: 3)Disease in diagnosie
        diag_matrix=read_data(join(self.preprocess_folder,"diag_matrix"),dtype={"HADM_ID":str})

        ##NOTE: 4)New Diseases
        # new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newICD9_CODE"),dtype={"HADM_ID":str})
        new_disease_matrix =read_data(join(concat_clamp_prefix,"allepis_newCUI"),dtype={"HADM_ID":str})
        create_folder(join(singledrug_prefix,rxnorm_id))

        ##ll randomly select one episode of treated patients
        if(not os.path.exists(join(singledrug_prefix,self.rxnorm_id,"treated_epis.csv"))):

            treated_epis = self.mat_noZero(
                pres_rxnorm_matrix[["HADM_ID",self.rxnorm_id]].set_index("HADM_ID"),1).reset_index()['HADM_ID']
            treated_epis=self.get_sampled_epis(treated_epis)
            write2file(pd.DataFrame({"HADM_ID":treated_epis}),join(singledrug_prefix,self.rxnorm_id,"treated_epis"))
        else:
            treated_epis=read_data(join(singledrug_prefix,self.rxnorm_id,"treated_epis"),dtype=str)

        # TABLE:selected_drugs
        print("\n\nRXNORM={}, DRUG_NAME={}".format(self.rxnorm_id,get_drugname_byrxcui_api(self.rxnorm_id)))
        print("#Episodes (MIMIC PRES): {:d}".format(len(treated_epis)))
        ###NOTE: 2)NEW DRUGS
        new_rxnorm_matrix =read_data(join(concat_clamp_prefix,"allepis_newRxNorm"),dtype={"HADM_ID":str})

        # self.compare_sider_newdisease(treated_epis,new_icd9_matrix)

        # #HACK: Import feature data
        # TODO: Feature1 pres diag
        # pres_diag_sider_matrix=read_data(
        #     join(self.feature_folder,"pres_diag_sider_matrix"),dtype=str).fillna(0)
        pres_diag_sider_matrix=inner_join(treated_epis,diag_matrix,"HADM_ID")
        # TODO: GET ADE CUI
        ades_df=get_side_effects_cui(self.rxnorm_id)
        ades_list=ades_df["CUI"].unique()
        # ades_df=get_ade_rxnorm_df()
        # ades_list=ades_df[ades_df['RxNorm']==self.rxnorm_id]['ICD9_CODE'].unique()

        # TABLE:selected_drugs episodes, # of Side Effects
        print("#Side Effects (SIDER): %d"%len(ades_list))

        pres_diag_sider_matrix=pres_diag_sider_matrix.reindex(columns=['HADM_ID']+list(ades_list)).dropna(axis=1, how='all')

        dissum_autoencoder=pd.concat(
            [read_data(join(self.feature_folder,"dissum_Autoencoder_EPIS"),dtype=str),
            read_data(join(self.feature_folder,"dissum_Autoencoder_128"))],axis=1,sort=False)
        five_autoencoder=pd.concat(
            [read_data(join(self.feature_folder,"five_Autoencoder_EPIS"),dtype=str),
            read_data(join(self.feature_folder,"five_Autoencoder_128"))],axis=1,sort=False)
        ##HACK: Import feature data
        feature_list=["pres_diag_sider_matrix","dissum_autoencoder","five_autoencoder"]
        feature_folder_list=[join(
            singledrug_prefix,self.rxnorm_id,folder) for folder in feature_list]
        [create_folder(folder) for folder in feature_folder_list]
        data_list=[pres_diag_sider_matrix, dissum_autoencoder,five_autoencoder]

        ## get features
        # for i in [0]:
        for i in range(len(feature_folder_list)):
            print("Feature=%s"%feature_list[i])
        # range(0,len(feature_folder_list)):
            treated_feature=inner_join(treated_epis,data_list[i],"HADM_ID")

            # NOTE: i=0, then use pres_diag_sider_matrix
            if(i==0):
                # TABLE:selected_drugs episodes, J_d
                print("# Episodes (MIMIC PRES^pres_diag_sider_matrix), J_d+1: {:d}, {:d}".format(*treated_feature.shape))
                # write2file(treated_feature,join(feature_folder_list[i],"pres_diag_sider_matrix"))
            else:
                # TABLE:selected_drugs episodes, 
                print("# Episodes (MIMIC PRES^{}): {:d}".format(feature_list[i],treated_feature.shape[0]))
                
            updated_treated_epis=treated_feature["HADM_ID"]
            # NOTE: 1. For TOP 20 FREQUENT DRUGS 2. For TOP 20 FREQUENT NEW DRUGS 3. For TOP 20 FREQUENT DISEASES  4. For TOP 20 FREQUENT NEW DISEASES
            frequent_item_treated = [self.mat_noZero(
                inner_join(updated_treated_epis,frequent_item_matrix,"HADM_ID"), axis_num=0
            ) for frequent_item_matrix in [pres_rxnorm_matrix,new_rxnorm_matrix,diag_matrix,new_disease_matrix]]
            
            # TABLE:selected_drugs episodes, # Episodes (MIMIC PRES^FEATURE^ITEMS)
            print("# Episodes (MIMIC PRES^{}^ITEM) drug/new/disease/new: {}".format(feature_list[i],[
                len(item_table) for item_table in frequent_item_treated
            ]))


            # for n_clusters in range(2,6,2):
            for n_clusters in [2]:

                if(os.path.exists(join(feature_folder_list[i],'CLUSTER_label_C%d.csv'%n_clusters))):
                    label_df=read_data(join(feature_folder_list[i],'CLUSTER_label_C%d'%n_clusters),dtype=str)
                    counter_labels=Counter(label_df['LABEL'])
                else:
                    labels=self.runKMeans(treated_feature.iloc[:,1:],n_clusters)
                    counter_labels=Counter(labels)
                    label_df = pd.DataFrame({self.epis_field:updated_treated_epis,'LABEL':labels})
                    write2file(label_df,join(feature_folder_list[i],'CLUSTER_label_C%d'%n_clusters))

                distance_list=[]
                for (item,treated_frequent_df, newitem_flag) in zip(
                    ["DRUG"]*2+["DISEASE_CUI"]*2,frequent_item_treated,[False,True,False,True]):
                    treated_frequent_label=inner_join(
                        label_df,treated_frequent_df,"HADM_ID")
                    grouped = treated_frequent_label.groupby("LABEL")
                    distance_list.append(plot_top20items(
                        item,grouped,n_clusters,counter_labels,figure_path=feature_folder_list[i],
                        feature_name=feature_list[i],new_flag=newitem_flag))
                # if()
                cluster_distance_df = pd.DataFrame(
                    distance_list,
                    columns=['N_DISTANCE','M_DISTANCE','Q_DISTANCE'])
                cluster_distance_df.insert(0, 'ITEMS', ['DRUG','NEW_DRUG','DISEASE','NEW_DISEASE'])
                write2file(cluster_distance_df,join(feature_folder_list[i],'CLUSTER_DISTANCE_C%d'%n_clusters))
        
            ## deprecated
            # range(0,len(feature_folder_list)):
                # self.model_plot(data_list[i],"kmeans",n_clusters,prescribed_patients,feature_folder_list[i])

    def debug(self):
        # NOTE: 
        pres_rxnorm_df=read_data(join(self.preprocess_folder,'pres_rxnorm_df'),dtype=str).dropna(subset=['RxNorm'])
        pres_rxnorm_df_group= pres_rxnorm_df.groupby(
            'RxNorm')['HADM_ID'].count().sort_values(ascending=False).to_frame().reset_index()
        pres_rxnorm_df_group.rename(columns={'HADM_ID':'NUM_OF_EPISODES'},inplace=True)

        ade_df=get_ade_rxnorm_df()[['RxNorm']].drop_duplicates()
        pres_rxnorm_df_group_sider=inner_join(pres_rxnorm_df_group,ade_df,"RxNorm").head(100)
        pres_rxnorm_df_group_sider.insert(
            1,"Drug Name",list(map(get_drugname_byrxcui_api,pres_rxnorm_df_group_sider['RxNorm'])))
        
        write2file(pres_rxnorm_df_group_sider,join(singledrug_prefix,"TOP100rxnorm_in_presANDsider"))

    def extra_two_drugs_anaylsis(self,rxnorm_list,not_new_dis=True):
        dis_item="DISEASE_CUI"
        result_prefix=join(
            singledrug_prefix,"CONCAT_RESULTS","DISEASE_SIDER_INTERSECTION")
            # join(result_path,"%s_%d_of_drugs"%(dis_item, len(rxnorm_id))))
        create_folder(result_prefix)
        epis_field="HADM_ID"
        dis_args=get_item_args(dis_item)
        dis_code_field=dis_args.get("CODE_FIELD")
        print("loading pres_rxnorm_matrix...")
        pres_rxnorm_matrix =read_data(join(self.preprocess_folder,"pres_rxnorm_matrix"),dtype={"HADM_ID":str})
        print("complete loading pres_rxnorm_matrix.")
        print("loading diag_matrix...")
        disease_file_path= join(self.preprocess_folder,"diag_matrix")\
            if not_new_dis else join(concat_clamp_prefix,"allepis_newCUI")
        disease_matrix=read_data(disease_file_path,dtype={epis_field:str})
        print("complete loading diag_matrix.")

        if(len(rxnorm_list[0])==1):
            result_path=join(result_prefix,"SINGLE_DRUG%s_%s"%("" if not_new_dis else "_New",dis_item))
            create_folder(result_path)
        else:
            result_path=join(result_prefix,"DOUBLE_DRUG%s_%s"%("" if not_new_dis else "_New",dis_item))
            create_folder(result_path)
            result_path=join(result_path,rxnorm_list[0][0])
            create_folder(result_path)

        for rxnorm_id in rxnorm_list:
            # print("\n")
            # print([epis_field,*rxnorm_id])

            pres_rxnorm_epis=pres_rxnorm_matrix[[epis_field,*rxnorm_id]].set_index(epis_field)
            pres_rxnorm_epis["PRESENT"]=pres_rxnorm_epis.prod(axis=1)
            treated_epis = self.mat_noZero(
                pres_rxnorm_epis[["PRESENT"]],1)\
                    .reset_index()[epis_field]

            # if(not len(treated_epis)):
                # print("%s, %d"%(("_").join(rxnorm_id), 0))
                # continue

            treated_epis=self.get_sampled_epis(treated_epis)
            # print("%s, %d"%(("_").join(rxnorm_id), len(treated_epis)))
            # continue
            
            pres_diag_sider_matrix=inner_join(treated_epis,disease_matrix,epis_field)
            if(dis_item=="DISEASE"):
                ades_df=get_ade_rxnorm_df()
                ades_list=ades_df[ades_df['RxNorm'].isin(rxnorm_id)][dis_code_field].unique()
            else:
                ades_df=pd.concat(
                    [get_side_effects_cui(rxnorm) for rxnorm in rxnorm_id],
                    axis=0,
                    ignore_index=True)
                ades_list=ades_df[dis_code_field].unique()
                
            # TABLE:selected_drugs
            print("\nRXNORM=%s, DRUG_NAME=%s"%(rxnorm_id,\
                list(map(search_rxnorm, rxnorm_id))))
            print("#Episodes (MIMIC PRES): {:d}".format(len(treated_epis)))
            # TABLE:selected_drugs episodes, # of Side Effects
            print("#Side Effects (SIDER): %d"%len(ades_list))

            pres_diag_sider_matrix=pres_diag_sider_matrix\
                .reindex(columns=[epis_field]+list(ades_list)).dropna(axis=1, how='all')
            # TABLE:selected_drugs episodes, J_d
            print("# Episodes (MIMIC PRES^pres_diag_sider_matrix), J_d+1: {:d}, {:d}".format(*pres_diag_sider_matrix.shape))
            write2file(pres_diag_sider_matrix,join(result_path,("_").join(rxnorm_id)))

    def get_dissumary_drug_feature_newdis(self,feature_flag, group_cui, unique_flag=False):
        epis_field, rank_field, dissum_field="HADM_ID", "GROUP_RANK", "TEXT"
        discharge_summary_file="/data/liu/mimic3/CLAMP_NER/input/ReRule0_Discharge summary_All.csv"
        result_prefix=join(singledrug_prefix,"NEWDIS_SIDER_DISCHARG",get_feature_name(feature_flag))
        create_folder(result_prefix)
        result_prefix=join(result_prefix,self.rxnorm_id)
        create_folder(result_prefix)
        if(unique_flag):
            result_prefix=join(result_prefix,"UNIQUE")
        else:
            result_prefix=join(result_prefix,"COMMON")
        create_folder(result_prefix)

        dissumary_df=read_data(discharge_summary_file, dtype=str,\
            usecols=[epis_field, rank_field, dissum_field])

        epis_cluster_df_dict=dict(tuple(read_data(join(singledrug_prefix,\
            self.rxnorm_id,\
                get_feature_name(feature_flag),\
                    "CLUSTER_label_C2"),dtype=str).groupby("LABEL")))

        all_sd_cuis=list(set(group_cui[0]).union(set(group_cui[1])))
        new_disease_matrix = read_data(
            join(concat_clamp_prefix,"allepis_newCUI"),dtype={"HADM_ID":str}, \
                usecols=[epis_field]+all_sd_cuis)

        for cluster_id in epis_cluster_df_dict.keys():
            result_path=join(result_prefix, "CLUSTER_%s"%cluster_id)
            create_folder(result_path)
            group_cluster=left_join(\
                epis_cluster_df_dict[cluster_id][[epis_field]],\
                    new_disease_matrix,epis_field).set_index(epis_field).loc[:,group_cui[int(cluster_id)]]
            for cui in group_cluster.columns:
                print(cui)
                result_final_path=join(result_path,cui)
                create_folder(result_final_path)
                ds=group_cluster[cui].dropna()
                # print(ds)
                treated_cui_epis=list(ds[ds>0].index)
                print(len(treated_cui_epis))
                print(len(dissumary_df.query('HADM_ID in @treated_cui_epis')))
                # print(dissumary_df.query('HADM_ID in @treated_cui_epis').head())
                treated_cui_epis_dissummary=dissumary_df.query('HADM_ID in @treated_cui_epis')
                treated_cui_epis_dissummary.apply(
                    lambda row: write2txt(row[dissum_field],\
                        join(result_final_path,"%s_%s"%(row[epis_field],row[rank_field]))), axis=1
                )

            # print(group_cluster.head())
        # quit()
        # group_0=left_join(epis_cluster_df["0"][[epis_field]],new_disease_matrix,epis_field)
        # print(group_0.head())
        # quit()
        # ds=group_0[epis_field]
        # print(ds)
        # ds=ds[ds>0]
        # print(len(ds))
        # print(ds)
        # print(len(dissumary_df.query('HADM_ID in @ds')))
        # print(dissumary_df.query('HADM_ID in @ds').head())
        # print(epis_cluster_df["1"]["C0008031"].head())


def get_feature_name(feature_flag):
    if feature_flag=="DS":
        return "pres_diag_sider_matrix"
    elif feature_flag=="EN":
        return "dissum_autoencoder"
    elif feature_flag=="ES":
        return "five_autoencoder"
    else:
        return None



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

def get_diseasename_byumlscui_api(umlscui):

    return None
    # for diseaae_code in discode_dict.values():
    #     disease_name=search_cui(umlscui,diseaae_code)
    #     if(disease_name):
    #         return disease_name
    # return None

def get_item_by_api(code_id, item):
    if(item=="DISEASE_CUI"):
        return get_diseasename_byumlscui_api(code_id)
    elif(item=="DRUG"):
        return get_drugname_byrxcui_api(code_id)


def norm_wasserstein_distance(list_vals):
    list_valsarray=list(map(np.array,list_vals))
    list_valsarray=[[i/array.sum() for i in array] for array in list_valsarray]
    w_distance=wasserstein_distance(*list_valsarray)
    return w_distance


def plot_top20items(item,grouped,n_clusters,counter_labels,figure_path,feature_name,new_flag=False):
    args=get_item_args(item)
    code_field, index_file, code_index, item_name, item_name_file, item_field_name=\
    	args["CODE_FIELD"],args["INDEX_FILE"],\
            args["CODE_INDEX"],args["ITEM_NAME"],args["CODE_NAME_FILE"],\
                args["ITEM_FIELD_NAME"]

    # DF[code_field,item_field_name]
    item_name_db_map = read_data(
        join(singledrug_prefix,item_name_file),dtype=str)\
            .set_index(code_field)[item_field_name].to_dict()
    # ##import drug list 
    # if os.path.exists(join(singledrug_prefix,"%s.csv"%index_file)):
    #     all_item_df=read_data(join(singledrug_prefix,index_file),dtype=str)
    #     all_item_dict=dict(zip(all_item_df[code_field],all_item_df[code_index]))
    # else:
    #     all_item_dict={}

    nrows = int(math.ceil(n_clusters/2.))
    plt.close('all')
    fig, axs = plt.subplots(nrows,2,sharey=True)

    fig.set_size_inches(30, 15*(nrows*0.75))
    fig.subplots_adjust(left=0.2,top=1.6, wspace=0.2,hspace=0.5)

    fig.suptitle(
        args["FIG_SUPTITLE"] % (" New"*new_flag, counter_labels),
        fontweight='bold',fontsize=16)

    ## top 20 frequent serires
    all_series_df_list=[]
    ## all serires
    combined_series_df_list=[]

    for (name, groupdf), ax in zip(grouped, axs.flatten()):      
        serires_whole=groupdf.iloc[:,2:].sum(axis=0).sort_values(ascending=False)
        serires=serires_whole.head(20)
        combined_series_df_list=combined_series_df_list+[serires_whole[serires_whole>0]]

        ## remove codes which are already recoded in "all_item_dict" from current series
        # group_item_dict=list(set(serires.index).difference(set([*all_item_dict])))
        # group_item_dict=dict(zip(
        #     group_item_dict,
        #     ["%s_%d"%(item_name,id) for id in list(range(len(all_item_dict),len(group_item_dict)+len(all_item_dict)))]))
        # all_item_dict={**all_item_dict,**group_item_dict}


        ax.bar(
            # [all_item_dict[item_code] for item_code in serires.index],
            serires.index,
            # textwrap.fill(get_drugname_byrxcui_api(rxcui)[:35],25)+"..." for rxcui in serires.index], 
            list(serires)
        )

        series_df= serires.to_frame().reset_index()
        series_df.columns=[code_field,"Count of Episodes"]
        series_df["Cluster_id"]=name
        series_df[item_field_name]=series_df[code_field].apply(
            lambda x: item_name_db_map.get(x,None) )
        # series_df[code_index]=series_df[code_field].apply(lambda x:all_item_dict[x])
        series_df["New %s"%(item_name)]=new_flag
        all_series_df_list=all_series_df_list+[series_df]

        
        ax.set_title("Cluster Label: %s"%name)
        plt.setp(
            ax.get_xticklabels(), rotation=30, 
            fontsize=10,
            horizontalalignment='right')
        
    plt.setp([a.get_yticklabels() for a in np.reshape(axs,(-1,2))[:, 1]], visible=False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(join(figure_path,'TOP20%s_%ss_C%d.png'%("_new"*new_flag,item_name,n_clusters)))
    
    # all_item_df=pd.DataFrame(all_item_dict.items(), columns=[code_field,code_index])

    # write2file(all_item_df,join(join(singledrug_prefix,index_file)))

    all_series_df=pd.concat(all_series_df_list,sort=False,ignore_index=True)
    counter_labels_string=re.sub('[^0-9a-zA-Z]+','_',str(counter_labels))
    write2file(all_series_df,join(figure_path,args["FILENAME_STATS"]%(
        feature_name,"_new"*new_flag,counter_labels_string)))

    # NOTE: find distinct items for each cluster
    common_items=list(set.intersection(
        *[set(serires_whole.index) for serires_whole in combined_series_df_list]))
    common_items_df=pd.DataFrame({code_field:common_items,"Common":True})

    combined_cluster_df_list=[]
    actual_num_clusers=len(combined_series_df_list)
    for i in range(actual_num_clusers):
        combined_cluster_df_list=combined_cluster_df_list+[
            combined_series_df_list[i].to_frame(name="Cluster_%d"%i)]
 
    combined_cluster_df=pd.concat(
        combined_cluster_df_list,sort=True,axis=1).fillna(0).reset_index().rename(
            {'index': code_field},axis=1)
    combined_cluster_df=left_join(combined_cluster_df,common_items_df, code_field).fillna({"Common":False})
    combined_cluster_df[item_field_name]=combined_cluster_df[code_field].apply(
            lambda x: item_name_db_map.get(x,None) )
    write2file(
        combined_cluster_df,
        join(figure_path,args["FILENAME_ALL"]%(
        feature_name,"_new"*new_flag,counter_labels_string)))
    distinct_cluster_df=combined_cluster_df[combined_cluster_df['Common']==False]
    write2file(
        distinct_cluster_df,
        join(figure_path,args["FILENAME_UNIQUE"]%(
        feature_name,"_new"*new_flag,counter_labels_string)))

    # NOTE: wasserstein DISTANCE: 1) ALL DISTRIBUTION n 2) TOP 20 m  3) TOP 20 unique q
    if(actual_num_clusers>1):
        distances=[]
        unique_series_df_list=[(
            distinct_cluster_df.loc[:,"Cluster_%d"%cluster_id]).sort_values(ascending=False).head(20) 
            for cluster_id in range(actual_num_clusers)]
        for distribution in [
            combined_series_df_list,
            [series.head(20) for series in combined_series_df_list],
            unique_series_df_list]:
            distances.append(norm_wasserstein_distance(set_union_index(distribution)))
        return distances    
    else:
        return [0,0,0]


def set_union_index(pd_series_list):
    union_index=pd_series_list[0].index
    for series in pd_series_list[1:]:
        union_index=union_index.union(series.index)
    return [series.reindex(union_index,fill_value=0) for series in pd_series_list]

def concat_wasserstein_distance(distance_filename_prefix="CLUSTER_DISTANCE_C2"):
    item_args=get_item_args("DRUG")
    code_field, item_field_name = itemgetter(
        "CODE_FIELD","ITEM_FIELD_NAME")(item_args)

    os.chdir(singledrug_prefix)
    dirs = [f for f in glob.glob(r"[0-9]*") if (os.path.isdir(f))]
    drug_names = [search_rxnorm(dir_num).split(" ")[0] for dir_num in  dirs]
    
    # for drug_order in range(5):
    feature_list=["pres_diag_sider_matrix","dissum_autoencoder","five_autoencoder"]
    result_path=join(singledrug_prefix,"CONCAT_RESULTS")
    # create_folder_overwrite(join(singledrug_prefix,"CONCAT_RESULTS"))
    create_folder(join(singledrug_prefix,"CONCAT_RESULTS"))


    group_item_order={"DISEASE":0, "NEW_DISEASE":1, "DRUG":2, "NEW_DRUG":3}
    plt.figure(figsize=(36, 28))
    fig, axs = plt.subplots(3,4,sharey=True,sharex=True)
    plt.suptitle("Wasserstein Distance of Distribution for 2 Clusters")
    i=0

    for feature_name in feature_list:
        all_distance_df=[]
        for drug_order in range(len(dirs)):
            drug_path=join(singledrug_prefix,dirs[drug_order])
            distance_df=read_data(join(drug_path,feature_name,distance_filename_prefix))
            distance_df.insert(0,item_field_name,drug_names[drug_order])
            distance_df.insert(0,code_field,dirs[drug_order])
            all_distance_df.append(distance_df)

        all_distance_df=pd.concat(all_distance_df, axis=0, sort=False)
        # NOTE: Uncomment
        # write2file(all_distance_df,join(
        #     result_path,"CONCAT_%s_%s"%(distance_filename_prefix,feature_name)))
        # NOTE: subplot line chart

        # for (item, groupdf), ax in zip(all_distance_df.groupby("ITEMS"),axs[i]):
        for item, groupdf in all_distance_df.groupby("ITEMS"):
            ax=axs[i,group_item_order.get(item)]
            if(i==0):
                ax.set_title(item,size='x-small')
            groupdf_resetindex=groupdf.set_index(code_field) 
            # groupdf_resetindex=groupdf.set_index(item_field_name) 

            for col in ["N_DISTANCE","M_DISTANCE","Q_DISTANCE"]:        
                ax.plot(
                    # groupdf_resetindex[item_field_name],
                    groupdf_resetindex[col],
                    label=col)
            if(i==len(feature_list)-1):
                ax.set_xticklabels(groupdf_resetindex[item_field_name], fontsize=5)
                plt.setp(ax.get_xticklabels(), horizontalalignment='center',rotation=90)

        i=i+1

    lines, labels = fig.axes[-1].get_legend_handles_labels()
    fig.legend(lines, labels, loc = 'upper right', frameon=False, fontsize='x-small')
    for ax, feature in zip(axs[:,0], feature_list):
        ax.set_ylabel(feature, rotation=90, size='xx-small')
    # fig.autofmt_xdate(rotation=30)
    # ax.xaxis.get_label().set_fontsize(20)
    plt.subplots_adjust(top=0.78, bottom=0.22, hspace=0.1, right=0.95,\
        wspace=0.15)
    # fig.tight_layout()
    plt.savefig(join(
        result_path,"%s_PLOT_ALL.png"%(
            distance_filename_prefix)),dpi=600)        

        # for item, groupdf in all_distance_df.groupby("ITEMS"):
        #     groupdf_resetindex=groupdf.set_index(code_field) 
        #     plt.figure()
        #     plt.xlabel(code_field)
        #     plt.ylabel('Wasserstein Distance')
        #     plt.title("Wasserstein Distance of %s Distribution"%item)
        #     for col in ["N_DISTANCE","M_DISTANCE","Q_DISTANCE"]:
        #         plt.plot(
        #             # groupdf_resetindex[code_field], 
        #             groupdf_resetindex[col], 
        #             label = col)
        #     plt.legend(loc='upper right', frameon=False, fontsize='x-small')
        #     plt.xticks(rotation=30, 
        #         fontsize=6,
        #         horizontalalignment='right')
        #     # plt.show()
        #     plt.savefig(join(
        #         result_path,"%s_PLOT_%s_%s.png"%(
        #             distance_filename_prefix, 
        #             item,
        #             feature_name)),
        #         bbox_inches='tight',dpi=120)
        #     plt.clf()


# def get_icd9_name_dict():
#     icd9_name_file="/data/liu/mimic3/EXTERNEL_DATASET/dexur_ICD9_DISEASENAME/CONCAT_RESULTS/ICD9_NAME_ALL_COMBINE.csv"
#     icd9_name_df=read_data(icd9_name_file,dtype=str)
#     return icd9_name_df.set_index("ICD9_CODE")["SHORT_TITLE"].to_dict()

# TODO:
def get_ade_rxnorm_df():

    if(os.path.exists(join(sideffect_prefix,'rxnorm_icd9_side_effects.csv'))):
        ade_df = read_data(join(sideffect_prefix,'rxnorm_icd9_side_effects'),dtype=str)
    else:
        ade_related_prefix="/data/liu/mimic3/CLAMP_NER/ADE_Related/"
        ndc_cui_map = read_data(join(ade_related_prefix,"ndc_cui_map"),dtype={'CODE':str,'CUI_CODE':str})
        ndc_cui_map = ndc_cui_map.set_index('CODE')['CUI_CODE'].to_dict()
        ade_df = read_data(
            join(sideffect_prefix, 'ndc_icd9_side_effects'), 
            dtype={'NDC':str,'ICD9_CODE':str},usecols=['NDC','ICD9_CODE']).drop_duplicates()
        ade_df['RxNorm']=ade_df['NDC'].apply(lambda x: ndc_cui_map.get(x, None))
        write2file(
            ade_df.dropna(subset=['RxNorm'])[["RxNorm","ICD9_CODE"]],
            join(sideffect_prefix,'rxnorm_icd9_side_effects'))   
    return ade_df


## drug_folders: None: not specifying drug combinations, list of folders: specifying drug combinations
def concat_plot_disease_sider(drug_folders=None,get_top_epis=False, dis_item="DISEASE_CUI", \
    matrix_path=dis_sider_intersection_path):

    top_num=30
    drug_code_field,drug_name_field=get_item_args("DRUG")["CODE_FIELD"],get_item_args("DRUG")["ITEM_FIELD_NAME"]
    dis_code_field=get_item_args(dis_item)["CODE_FIELD"]
    epis_field, count_epis_field="HADM_ID", "COUNT_HADM_ID"
    result_path=matrix_path
    # n_clusters=2
    pres_diag_sider_X_sum_list=[]
    pres_diag_sider_Y_sum_list=[]
    if(drug_folders is None):
        drug_folders=list(filter(
            lambda file:re.match("[0-9]{3,}",file),os.listdir(matrix_path)))
        drug_folders=list(map(lambda x: x.split('.')[0],drug_folders))

    # feature ="pres_diag_sider_matrix"
    for drug_rxnorm in drug_folders:

        # if(drug_folders):
        #     file_path=join(matrix_path,drug_rxnorm.split("_")[0],drug_rxnorm)
        # else:
        #     file_path=join(matrix_path,drug_rxnorm)
        file_path=join(matrix_path,drug_rxnorm)
        pres_diag_sider_matrix=read_data(
            file_path,dtype={epis_field:str,drug_code_field:str}).set_index(epis_field)

        # NOTE: X, row: axis=1, Y, column: axis=0
        pres_diag_sider_matrix_sum=pres_diag_sider_matrix.sum(axis=1)
        pres_diag_sider_X_sum_list.append(pd.DataFrame({
            drug_code_field: drug_rxnorm, 
            epis_field: pres_diag_sider_matrix_sum.index.values.tolist(),
            count_epis_field: pres_diag_sider_matrix_sum}))
        pres_diag_sider_Y_sum_list.append(pd.DataFrame({
            drug_code_field: drug_rxnorm, 
            dis_code_field: pres_diag_sider_matrix.sum(axis=0)}))

    pres_diag_sider_X_sum_df=pd.concat(pres_diag_sider_X_sum_list, axis=0, ignore_index=True )
    pres_diag_sider_Y_sum_df=pd.concat(pres_diag_sider_Y_sum_list, axis=0, ignore_index=True )

    
    if(get_top_epis):
        pres_diag_sider_X_topdis=pres_diag_sider_X_sum_df.copy()\
            .groupby(drug_code_field).apply(
                lambda x: x\
                    .assign(HADM_ID_COUNT=lambda x: x[[epis_field,count_epis_field]].astype(str).agg(': '.join,axis=1))
                    .sort_values([count_epis_field],ascending=False)\
                    .head(top_num)
                    .assign(ORDER=list(range(top_num)))
            ).set_index(["ORDER",drug_code_field])["HADM_ID_COUNT"].unstack("ORDER").reset_index()
        write2file(pres_diag_sider_X_topdis,join(result_path,"pres_%s_sider_top_%d_%s")%(dis_item,top_num,epis_field))


    pres_diag_sider_Y_percentage_df=pres_diag_sider_Y_sum_df.copy()
    pres_diag_sider_Y_percentage_df.iloc[:,1]=pres_diag_sider_Y_percentage_df.iloc[:,1].value_counts(normalize=True) * 100
    pres_diag_sider_Y_sqrt_df=pres_diag_sider_Y_sum_df.copy()
    pres_diag_sider_Y_sqrt_df.iloc[:,1]=pres_diag_sider_Y_sqrt_df.iloc[:,1].apply(np.sqrt)
    

    for df, (axis, field, sum_item, color) in zip(
        [pres_diag_sider_X_sum_df,\
            # pres_diag_sider_Y_k_df,
            pres_diag_sider_Y_sum_df,\
                # pres_diag_sider_Y_logsum_df],
            # pres_diag_sider_Y_percentage_df,
            pres_diag_sider_Y_sqrt_df],
        [("X",count_epis_field,"Sum of Side Effects","Blues"),
        ("Y",dis_code_field, "Sum of Patients","summer_r"),
        # ("Y_1k",dis_code_field, "Sum of Patients/1k"),
        # ("Y_Log",dis_code_field, "Log(Sum of Patients)"),
        # ("Y_Percentage",dis_code_field, "Percentage(Sum of Patients)"),
        ("Y_Sqrt",dis_code_field, "Sqrt(Sum of Patients)","summer_r")
        ]):
        # write2file(df,join(result_path,"PRES_DIAG_SIDER_SUM_%s"%axis))
        dd=pd.melt(df,id_vars=[drug_code_field],value_vars=[field],var_name=axis)
        ax=sns.boxplot(x=drug_code_field,y='value',data=dd,hue=axis,palette=color)
        # ax.legend(loc="upper left")
        ax.set(ylabel=sum_item)
        ax.set(xlabel=drug_name_field)
        if(drug_folders):
            xlabel_names = [
                (" and ").join([search_rxnorm(x).split(" ")[0] for x in t.get_text().split("_")])  \
                    for t in ax.get_xticklabels()]
        else:
            xlabel_names = [search_rxnorm(t.get_text()).split(" ")[0]  for t in ax.get_xticklabels()]
        ax.set_xticklabels(xlabel_names)
        plt.legend([],[], frameon=False)
        plt.xticks(rotation=30, 
            fontsize=8,
            horizontalalignment='right')
        plt.savefig(
            join(result_path,"%s_%s_%s.png")%(axis,field,"Box_Whisker_Plot"),
            bbox_inches='tight',dpi=160)
        plt.clf()

def concat_top_newdiseases(top_num=20):
    root_path=singledrug_prefix
    all_new_dis_items_filename="COMBINED*_new_diseases_*.csv"   
    feature_list=['pres_diag_sider_matrix','dissum_autoencoder', 'five_autoencoder']
    new_dis_item_args=get_item_args("DISEASE_CUI")
    new_dis_code_field, new_dis_item_field_name=itemgetter("CODE_FIELD","ITEM_FIELD_NAME")(new_dis_item_args)
    drug_code_field=get_item_args("DRUG")["CODE_FIELD"]
    dis_code_field=get_item_args("DISEASE")["CODE_FIELD"]
    cluster_field, epis_field="Cluster_%s", "HADM_ID"
    result_path=join(singledrug_prefix,"CONCAT_RESULTS")
    n_clusters=2
    
    # ade_rxnorm_df=get_ade_rxnorm_df()
    # icd9_name_dict=get_icd9_name_dict()

    drug_folders=list(filter(
        lambda file:re.match("[0-9]{3,}",file),os.listdir(root_path)))
    
    for feature in feature_list:
        top_df_list=[]
        for drug_rxnorm in drug_folders:
            sub_path=join(root_path,drug_rxnorm)
            sub_feature_path=join(sub_path,feature)

            
            combined_file=list(filter(
                lambda file:all([word in file for word in ["COMBINED","new_diseases"]]),
                os.listdir(sub_feature_path)
            ))[0]
            combined_df=read_data(join(sub_feature_path,combined_file))
            distinct_df=combined_df[combined_df["Common"]==False]
            unique_top_table_list=[]
            all_top_table_list=[]
            for i in range(n_clusters):
                unique_all_top_dfs=list(map(lambda df: df[[new_dis_code_field,cluster_field%i,new_dis_item_field_name]] \
                    .set_index(new_dis_code_field)\
                        .sort_values(by=[cluster_field%i],ascending=False)\
                            .head(top_num).rename(
                                columns={cluster_field%i:"NUM_EPISODES_"+cluster_field%i}
                                ).reset_index(), [distinct_df,combined_df]))
                ## match new diseases to sider
                unique_all_top_dfs=[
                    left_join(df,get_side_effects_cui(drug_rxnorm),new_dis_code_field)\
                        for df in unique_all_top_dfs]
                unique_top_table_list.append(unique_all_top_dfs[0])
                all_top_table_list.append(unique_all_top_dfs[1])
            top_df_list.append(unique_top_table_list+all_top_table_list)

        # #if(save_top20_newdis_2_latex):
        rxnormdf_lists=zip(
            drug_folders,
            list(map(get_drugname_byrxcui_api,drug_folders)),
            top_df_list
            )        
        save_new_dis_df_as_latextable(
            rxnorm_dflists,join(result_path,"{}_Top{}_NewDisease.tex".format(feature,top_num)))

  


def load_save_name_db():

    for item in ["DRUG","DISEASE_CUI"]:
        item_name_file=get_item_args(item)["CODE_NAME_FILE"]
        item_name_db=load_rxnorm_name_db() if item=="DRUG" \
            else (load_icd9_name_db() if item=="DISEASE" \
                else load_umlscui_name_db())
        
        write2file(item_name_db,join(singledrug_prefix,item_name_file))
    

import itertools
if __name__ == '__main__':
    doc_drugs=["1658259","866924","966571","885257","836358","855334",\
        "855290","1808219","1724383","1719291","1807516","213169","1659151"]
        
    # print("START: %s"%(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')))
    # ##STEP:1 load code_name files for drugs, diseases and new diseases
    # load_save_name_db()

    # ##STEP:2 cluster analysis for each drugs, save tables and plots
    # for rxnorm_id in doc_drugs:
    #     fc = feature_creation(rxnorm_id)
    #     fc.run_model_plot()

    # print("END: %s"%(datetime.today().strftime('%Y-%m-%d-%H:%M:%S')))
    
    #STEP:3 concat wasserstein distance and plot
    # concat_wasserstein_distance()
    ##STEP:4 concat new diseases 
    # concat_top_newdiseases()
    ##STEP:5 pres_diag_sider sum and plot
    # concat_plot_disease_sider()

    ##STEP:6 check disease-sider intersection on two drugs
    # -----------------------------------------------------------------
    ## one drug for disease
    # doc_drugs=[[rxnorm] for rxnorm in doc_drugs]
    # fc = feature_creation('NONE')
    # fc.extra_two_drugs_anaylsis(doc_drugs,not_new_dis=True)
    # concat_plot_disease_sider(get_top_epis=False, dis_item="DISEASE_CUI", \
    #     matrix_path=join(dis_sider_intersection_path, "SINGLE_DRUG_DISEASE_CUI"))
    # -----------------------------------------------------------------
    # top_10_pairs_list=[["1658259","866924"],["1658259","1808219"],["1658259","885257"],["866924","885257"],\
    #     ["866924","1808219"],["1658259","836358"],["866924","966571"],["1658259","966571"],["866924","836358"],["1658259","1724383"]]
    # # ## two drugs: disease
    # top_10_pairs = list(map(("_").join, top_10_pairs_list))
    # # fc = feature_creation('NONE')
    # # fc.extra_two_drugs_anaylsis(top_10_pairs_list,not_new_dis=True)
    # concat_plot_disease_sider(drug_folders=top_10_pairs, get_top_epis=True, dis_item="DISEASE_CUI", \
    #     matrix_path=join(dis_sider_intersection_path, "DOUBLE_DRUG_DISEASE_CUI/1658259/"))
    # -----------------------------------------------------------------
    # dis_item="DISEASE_CUI"
    # for i in range(0,len(doc_drugs)):
    #     doc_drug_lists=[[doc_drugs[i],rxnorm] for rxnorm in (set(doc_drugs)-set([doc_drugs[i]]))]
    # #     # print(doc_drug_lists)
    #     # fc = feature_creation('NONE')
    #     # fc.extra_two_drugs_anaylsis(doc_drug_lists,dis_item=dis_item)
    #     concat_plot_disease_sider(get_top_epis=True, dis_item=dis_item, \
    #         matrix_path=join(dis_sider_intersection_path, "DOUBLE_DRUG_%s"%dis_item,doc_drugs[i]))

    # -----------------------------------------------------------------

    # -----------------------------------------------------------------
    # # top_10_pairs_list=[["1658259","866924"],["866924","885257"],["1658259","885257"]\
    # #     ,["1658259","966571"],["1658259","213169"], ["866924","966571"], ["866924","836358"],["1658259","836358"],\
    # #         ["1658259","855334"],["866924","1808219"]]
    # dis_item="DISEASE_CUI"
    # print("%s: \n\n"%dis_item)
    # #     # print(doc_drug_lists)
    # fc = feature_creation('NONE')
    # fc.extra_two_drugs_anaylsis(top_10_pairs_list,dis_item=dis_item)
    # print("\n\n%s: \n\n"%dis_item)
    # dis_item="DISEASE"
    # fc = feature_creation('NONE')
    # fc.extra_two_drugs_anaylsis(top_10_pairs_list,dis_item=dis_item)

    # top_10_pairs = list(map(("_").join, top_10_pairs_list))

    # for dis_item in ["DISEASE_CUI","DISEASE"]:
    #     concat_plot_disease_sider(drug_folders=top_10_pairs, get_top_epis=True, dis_item=dis_item, \
    #         matrix_path=join(dis_sider_intersection_path, "CORRECT_DOUBLE_DRUG_%s/1658259"%dis_item))
    # concat_plot_disease_sider(get_top_epis=True, dis_item=dis_item, \
        # matrix_path=join(dis_sider_intersection_path, "DOUBLE_DRUG_%s"%dis_item,doc_drugs[i]))

    # top_10_pairs=["1658259_866924","866924_885257","1658259_885257"\
    #     ,"1658259_966571","1658259_213169", "866924_966571", "866924_836358","1658259_836358",\
    #         "1658259_855334","866924_1808219"]
    # for dis_item in ["DISEASE_CUI","DISEASE"]:
    #     concat_plot_disease_sider(drug_folders=top_10_pairs, get_top_epis=True, dis_item=dis_item, \
    #         matrix_path=join(dis_sider_intersection_path, "CORRECT_DOUBLE_DRUG_%s/1658259"%dis_item))
    # drug_paris=[[*pair] for pair in itertools.combinations(doc_drugs, 2)]
    # dis_item="DISEASE_CUI"
    # print(dis_item)
    # fc = feature_creation('NONE')
    # fc.extra_two_drugs_anaylsis(drug_paris,dis_item=dis_item)
    # dis_item="DISEASE"
    # print(dis_item)
    # fc = feature_creation('NONE')
    # fc.extra_two_drugs_anaylsis(drug_paris,dis_item=dis_item)
    # -----------------------------------------------------------------


    # SUP_STEP:check extracted section titles

    # ----------------------------------------------------------------- 
    # df=read_data(join("/data/liu/mimic3/CLAMP_NER/input",\
    #     "ReRule0_Discharge summary_All"))
    # print(df.TITLE.unique())

    ## test for one rxnorm
    # rxnorm_id = "197380"
    # print(get_item_args("DISEASE")['INDEX_FILE'])

    # -----------------------------------------------------------------
    # # # NOTE: Get preliminary results
    # for rxnorm_id in ["1658259","866924","966571","885257","836358","855334",
        # "855290","1808219","1724383","1719291","1807516","213169","1659151"][1:2]:

    # for rxnorm_id in [
    #     "197380","1807632","1807639","1807634","1658259",
    #     "866924","1807633","1860466","847630","866508","1361615"][:1]:
    #     fc = feature_creation(rxnorm_id)
    #     fc.run_model_plot()

        # TODO:delete function add item_name
        # fc.add_item_name_to_CATEGORY(cat_title="STATS")
        # fc.add_item_name_to_CATEGORY(cat_title="DISTINCT")
        # fc.add_item_name_to_CATEGORY(cat_title="COMBINED")
    # concat_wasserstein_distance()
    # # NOTE: Get preliminary results
    # -----------------------------------------------------------------

    # concat_top_newdiseases()
    # path=(join(singledrug_prefix,"drug_index_rxnorm_name"))

    # append_csv_byrow(["test1","test2","test3"],path)

    # # NOTE: supplement disease names
    # -----------------------------------------------------------------
    # STEP: get_dissumary_drug_feature_newdis
    # # #NOTE:metropolol
    # rxnorm_id="866924"
    # group_newdis=[["C0008031","C0039231","C0013404"],["C0008031","C0013404"]]
    # #NOTE:vancomycin
    # rxnorm_id="1807516"
    # group_newdis=[["C0015672"],[]]
    # group_newdis=[["C0002871","C0013404","C0039231"],["C0002871","C0039231"]]

    # #NOTE:furosemide
    # rxnorm_id="1719291"
    # group_newdis_unique=[[],["C0019151"]]
    # group_newdis=[["C0002871"],["C0002871"]]
    # feature_flag="DS"

    # feature_flag="ES"
    # -----------------------------------------------------------------
    # #NOTE:furosemide
    # rxnorm_id="1719291"
    # group_newdis_unique=[["C0000737","C0546884","C0009806"],[]]
    # group_newdis=[["C0002871"],[]]
    # # #NOTE:metropolol
    # rxnorm_id="866924"
    # group_newdis=[["C0008031","C0476273","C0039231"],["C0008031", "C0039231","C0013404"]]
    # #NOTE:vancomycin
    # rxnorm_id="1807516"
    # group_newdis_unique=[["C0002871","C0039231","C0015967","C0013404"],[]]
    # group_newdis=[["C0002871","C0039231"],["C0428977"]]
    # feature_flag="EN"

    # -----------------------------------------------------------------
    # #NOTE:furosemide
    # rxnorm_id="1719291"
    # group_newdis_unique=[["C0000737","C0546884","C0009806"],[]]
    # group_newdis=[["C0002871"],[]]
    # # #NOTE:metropolol
    # rxnorm_id="866924"
    # group_newdis=[["C0008031","C0039231","C0013404"],["C0008031", "C0476273","C0039231"]]
    # #NOTE:vancomycin
    rxnorm_id="1807516"
    group_newdis_unique=[["C0002871","C0039231","C0015967","C0013404"],[]]
    group_newdis=[["C0002871","C0039231"],["C0428977"]]
    feature_flag="ES"

    fc = feature_creation(rxnorm_id)
    fc.get_dissumary_drug_feature_newdis(feature_flag ,group_newdis_unique,True)
    fc.get_dissumary_drug_feature_newdis(feature_flag,group_newdis)


# def main():
#     # rxnorm_id = "197380"
#     # fc = feature_creation(rxnorm_id)
#     # fc.run_model_plot()
    
#     if len(sys.argv) != 2:
#         print("Wrong command format, please follwoing the command format below:")
#         print("python single_drug_analysis.py [RxNorm]")
#         exit(0)

#     if len(sys.argv) == 2:    
#         rxnorm_id = sys.argv[1]         
#         fc = feature_creation(rxnorm_id)
#         fc.run_model_plot()

        






        
        
        