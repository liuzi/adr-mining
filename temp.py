

def plot_top20drugs(grouped,n_clusters,counter_labels,figure_path,drug_list_path,feature_name,newdrug=False):
    drug_code_field="RxNorm Code"
    ##import drug list 
    if os.path.exists(join(singledrug_prefix,"drug_index_rxnorm.csv")):
        all_drug_df=read_data(join(singledrug_prefix,"drug_index_rxnorm"),dtype=str)
        all_drug_dict=dict(zip(all_drug_df[drug_code_field],all_drug_df["Drug Index"]))
    else:
        all_drug_dict={}

    nrows = int(math.ceil(n_clusters/2.))
    plt.close('all')
    fig, axs = plt.subplots(nrows,2,sharey=True)

    fig.set_size_inches(30, 15*(nrows*0.75))
    fig.subplots_adjust(left=0.2,top=1.6, wspace=0.2,hspace=0.5)

    fig.suptitle(
        "TOP 20 Frequent%s Drugs - Count of Episodes, Cluster: %s" % (" New"*newdrug, counter_labels),
        fontweight='bold',fontsize=16)
    ## top 20 frequent serires
    all_serires_df_list=[]
    ## all serires
    combined_serires_df_list=[]
    
    for (name, groupdf), ax in zip(grouped, axs.flatten()):      
        serires_whole=groupdf.iloc[:,2:].sum(axis=0).sort_values(ascending=False)
        serires=serires_whole.head(20)
        combined_serires_df_list=combined_serires_df_list+[serires_whole[serires_whole>0]]

        ## remove rxnorm codes which are already recoded in "all_drug_dict" from current series
        group_drug_dict=list(set(serires.index).difference(set([*all_drug_dict])))
        group_drug_dict=dict(zip(
            group_drug_dict,
            ["Drug_%d"%id for id in list(range(len(all_drug_dict),len(group_drug_dict)+len(all_drug_dict)))]))
        all_drug_dict={**all_drug_dict,**group_drug_dict}

        serires_df= serires.to_frame().reset_index()
        serires_df.columns=["RxNorm Code","Count of Episodes"]
        serires_df["Cluster_id"]=name
        serires_df["Drug Index"]=serires_df['RxNorm Code'].apply(lambda x:all_drug_dict[x])
        serires_df["New Drug"]=newdrug
        all_serires_df_list=all_serires_df_list+[serires_df]


        ax.bar(
            [all_drug_dict[rxcui] for rxcui in serires.index],

            list(serires)
        )
        
        ax.set_title("Cluster Label: %s"%name)
        plt.setp(
            ax.get_xticklabels(), rotation=30, 
            fontsize=10,
            horizontalalignment='right')
        
    plt.setp([a.get_yticklabels() for a in np.reshape(axs,(-1,2))[:, 1]], visible=False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(join(figure_path,'TOP20%s_drugs_C%d.png'%("_new"*newdrug,n_clusters)))


    all_drug_df=pd.DataFrame(all_drug_dict.items(), columns=['RxNorm Code',"Drug Index"])
    write2file(all_drug_df,join(join(singledrug_prefix,"drug_index_rxnorm")))
    
    all_serires_df=pd.concat(all_serires_df_list,sort=False,ignore_index=True)
    counter_labels_string=re.sub('[^0-9a-zA-Z]+','_',str(counter_labels))
    write2file(all_serires_df,join(figure_path,"STATS_%s%s_drugs_%s"%(
        feature_name,"_new"*newdrug,counter_labels_string)))

    # NOTE: find distinct items for each cluster
    common_items=list(set.intersection(
        *[set(serires_whole.index) for serires_whole in combined_serires_df_list]))
    common_items_df=pd.DataFrame({drug_code_field:common_items,"Common":True})

    combined_cluster_df_list=[]
    for i in range(n_clusters):
        combined_cluster_df_list=combined_cluster_df_list+[
            combined_serires_df_list[i].to_frame(name="Cluster_%d"%i)]
 
    combined_cluster_df=pd.concat(
        combined_cluster_df_list,sort=True,axis=1).fillna(0).reset_index().rename(
            {'index': drug_code_field},axis=1)
    combined_cluster_df=left_join(combined_cluster_df,common_items_df, drug_code_field).fillna({"Common":False})
    write2file(
        combined_cluster_df,
        join(figure_path,"COMBINED_%s%s_drugs_%s"%(
        feature_name,"_new"*newdrug,counter_labels_string)))

    
    distinct_cluster_df=combined_cluster_df[combined_cluster_df['Common']==False]
    write2file(
        distinct_cluster_df,
        # .reset_index().rename({'index':'RxNorm Code'},axis=1),
        join(figure_path,"DISTINCT_%s%s_drugs_%s"%(
        feature_name,"_new"*newdrug,counter_labels_string)))

    # NOTE: wasserstein DISTANCE: 1) ALL DISTRIBUTION n 2) TOP 20 m  3) TOP 20 unique q
    distances=[]
    unique_serires_df_list=[(
        distinct_cluster_df.loc[:,"Cluster_%d"%cluster_id]).sort_values(ascending=False).head(20) 
        for cluster_id in range(n_clusters)]
    for distribution in [
        combined_serires_df_list,
        [series.head(20) for series in combined_serires_df_list],
        unique_serires_df_list]:
        distances.append(norm_wasserstein_distance(set_union_index(distribution)))
    # n_distance=wasserstein_distance(*set_union_index(combined_serires_df_list))
    # m_distance=wasserstein_distance(
    #     *set_union_index([series.head(20) for series in combined_serires_df_list]))

    # q_distance=wasserstein_distance(*set_union_index(unique_serires_df_list))
    return distances
    # return [n_distance, m_distance, q_distance]


def plot_top20diseases(grouped,n_clusters,counter_labels,figure_path,feature_name,newdisease=False):
    disease_code_field="ICD9_CODE"
    ##import drug list 
    if os.path.exists(join(singledrug_prefix,"disease_index_icd9.csv")):
        all_disease_df=read_data(join(singledrug_prefix,"disease_index_icd9"),dtype=str)
        all_disease_dict=dict(zip(all_disease_df[disease_code_field],all_disease_df["Disease Index"]))
    else:
        all_disease_dict={}

    nrows = int(math.ceil(n_clusters/2.))
    plt.close('all')
    fig, axs = plt.subplots(nrows,2,sharey=True)

    fig.set_size_inches(30, 15*(nrows*0.75))
    fig.subplots_adjust(left=0.2,top=1.6, wspace=0.2,hspace=0.5)

    fig.suptitle(
        "TOP 20 Frequent%s Diseases - Count of Episodes, Cluster: %s" % (" New"*newdisease, counter_labels),
        fontweight='bold',fontsize=16)

    ## top 20 frequent serires
    all_serires_df_list=[]
    ## all serires
    combined_serires_df_list=[]

    for (name, groupdf), ax in zip(grouped, axs.flatten()):      
        serires_whole=groupdf.iloc[:,2:].sum(axis=0).sort_values(ascending=False)
        serires=serires_whole.head(20)
        combined_serires_df_list=combined_serires_df_list+[serires_whole[serires_whole>0]]

        ## remove icd9 codes which are already recoded in "all_disease_dict" from current series
        group_disease_dict=list(set(serires.index).difference(set([*all_disease_dict])))
        group_disease_dict=dict(zip(
            group_disease_dict,
            ["Disease_%d"%id for id in list(range(len(all_disease_dict),len(group_disease_dict)+len(all_disease_dict)))]))
        all_disease_dict={**all_disease_dict,**group_disease_dict}
        # print("Cluster_%s"%name)
        # print(group_drug_dict)
        ax.bar(
            [all_disease_dict[icd9] for icd9 in serires.index],
            # textwrap.fill(get_drugname_byrxcui_api(rxcui)[:35],25)+"..." for rxcui in serires.index], 
            list(serires)
        )

        serires_df= serires.to_frame().reset_index()
        serires_df.columns=["ICD9_CODE","Count of Episodes"]
        serires_df["Cluster_id"]=name
        serires_df["Disease Index"]=serires_df['ICD9_CODE'].apply(lambda x:all_disease_dict[x])
        serires_df["New Disease"]=newdisease
        all_serires_df_list=all_serires_df_list+[serires_df]
        
        ax.set_title("Cluster Label: %s"%name)
        plt.setp(
            ax.get_xticklabels(), rotation=30, 
            fontsize=10,
            horizontalalignment='right')
        
    plt.setp([a.get_yticklabels() for a in np.reshape(axs,(-1,2))[:, 1]], visible=False)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(join(figure_path,'TOP20%s_diseases_C%d.png'%("_new"*newdisease,n_clusters)))
    
    all_disease_df=pd.DataFrame(all_disease_dict.items(), columns=['ICD9_CODE',"Disease Index"])

    write2file(all_disease_df,join(join(singledrug_prefix,"disease_index_icd9")))

    all_serires_df=pd.concat(all_serires_df_list,sort=False,ignore_index=True)
    counter_labels_string=re.sub('[^0-9a-zA-Z]+','_',str(counter_labels))
    write2file(all_serires_df,join(figure_path,"STATS_%s%s_diseases_%s"%(
        feature_name,"_new"*newdisease,counter_labels_string)))

    # NOTE: find distinct items for each cluster
    common_items=list(set.intersection(
        *[set(serires_whole.index) for serires_whole in combined_serires_df_list]))
    common_items_df=pd.DataFrame({disease_code_field:common_items,"Common":True})
    combined_cluster_df_list=[]
    for i in range(n_clusters):
        combined_cluster_df_list=combined_cluster_df_list+[
            combined_serires_df_list[i].to_frame(name="Cluster_%d"%i)]
 
    combined_cluster_df=pd.concat(
        combined_cluster_df_list,sort=True,axis=1).fillna(0).reset_index().rename(
            {'index': disease_code_field},axis=1)
    combined_cluster_df=left_join(combined_cluster_df,common_items_df, disease_code_field).fillna({"Common":False})
    write2file(
        combined_cluster_df,
        join(figure_path,"COMBINED_%s%s_diseases_%s"%(
        feature_name,"_new"*newdisease,counter_labels_string)))
    distinct_cluster_df=combined_cluster_df[combined_cluster_df['Common']==False]
    write2file(
        distinct_cluster_df,
        join(figure_path,"DISTINCT_%s%s_diseases_%s"%(
        feature_name,"_new"*newdisease,counter_labels_string)))

    # NOTE: wasserstein DISTANCE: 1) ALL DISTRIBUTION n 2) TOP 20 m  3) TOP 20 unique q
    distances=[]
    unique_serires_df_list=[(
        distinct_cluster_df.loc[:,"Cluster_%d"%cluster_id]).sort_values(ascending=False).head(20) 
        for cluster_id in range(n_clusters)]
    for distribution in [
        combined_serires_df_list,
        [series.head(20) for series in combined_serires_df_list],
        unique_serires_df_list]:
        distances.append(norm_wasserstein_distance(set_union_index(distribution)))
    return distances    


    # n_distance=wasserstein_distance(*set_union_index(combined_serires_df_list))
    # m_distance=wasserstein_distance(
    #     *set_union_index([series.head(20) for series in combined_serires_df_list]))
    # unique_serires_df_list=[(
    #     distinct_cluster_df.loc[:,"Cluster_%d"%cluster_id]).sort_values(ascending=False).head(20) 
    #     for cluster_id in range(n_clusters)]
    # q_distance=wasserstein_distance(*set_union_index(unique_serires_df_list))
    # # print(combined_serires_df_list)
    # return [n_distance, m_distance, q_distance]


# for (treated_frequent_df, newdrug_flag) in zip(
#     frequent_item_treated[:2],[False, True]):
#     treated_frequent_label=inner_join(
#         label_df,treated_frequent_df,"HADM_ID")
#     grouped = treated_frequent_label.groupby("LABEL")
#     distance_list.append(plot_top20drugs(
#         grouped,n_clusters,counter_labels,figure_path=folder_list[i],
#         drug_list_path=join(singledrug_prefix,self.rxnorm_id),
#         newdrug=newdrug_flag,feature_name=feature_list[i]))

# for (treated_frequent_df, newdisease_flag) in zip(
#     frequent_item_treated[2:],[False,True]):
#     treated_frequent_label=inner_join(
#         label_df,treated_frequent_df,"HADM_ID")
#     grouped = treated_frequent_label.groupby("LABEL")
#     distance_list.append(plot_top20diseases(
#         grouped,n_clusters,counter_labels,figure_path=folder_list[i],
#         drug_list_path=join(singledrug_prefix,self.rxnorm_id),
#         newdisease=newdisease_flag,feature_name=feature_list[i]))

