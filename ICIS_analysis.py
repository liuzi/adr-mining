import pandas as pd
import matplotlib.pyplot as plt
from utils._path import singledrug_prefix
from os.path import join
import re
from adjustText import adjust_text
# from cluster_analysis import get_all_distance_ade
from term_process.umls_api import search_rxnorm
import seaborn as sns

pd.options.display.max_rows = 4000
result_path = join(singledrug_prefix, "ICIS_PLOT")

def corrFilter(x: pd.DataFrame, bound: float):
    xCorr = x.corr(method='spearman')
    xFiltered = xCorr[((xCorr >= bound) | (xCorr <= -bound)) & (xCorr !=1.000)]
    xFlattened = xFiltered.unstack().sort_values().drop_duplicates().dropna()
    return xFlattened

def diff_corr_analysis(df):
    two_sum_distance="2sum_distance"
    two_sum_diff_distance="2sum_diff_in_dis"
    # diff of inter-cluster Q dist between new drugs and all drugs distributions when clustering is done using ES feature
    # NOTE: want signals of ades, patients have similar diseases but different new diseases
    df["Diff_ND-D_ES_Q"] = (df["ND_ES_Q"] - df["D_ES_Q"])
    # diff of inter-cluster Q dist between new dis and all dis distributions when clustering is done using ES feature
    df["Diff_NS-S_ES_Q"] = (df["NS_ES_Q"] - df["S_ES_Q"])
    #sum of Q dist new dis distributions and Q dist new drug distributions when clustering is done using ES feature
    # NOTE: both are high
    df[two_sum_distance] = (df["NS_ES_Q"] + df["ND_ES_Q"])

    #sum of the two differences, measure singals of ade
    df[two_sum_diff_distance] = df["Diff_NS-S_ES_Q"] + df["Diff_ND-D_ES_Q"]

    x_labels_dict={two_sum_distance: "SUM-Q",\
         two_sum_diff_distance: "SUM-DIFF-Q"}
    for x_field in [two_sum_distance ,two_sum_diff_distance]:
        z=df[x_field].to_list()
        y=df["norm_ES_CM_CBOTH"].to_list()
        plt.scatter(z,y)
            # df[x_field], df["norm_ES_CM_CBOTH"])
        texts = []    
        for xx, yy, tt in zip(z, y, df["Name"].to_list()):
            texts.append(plt.text(xx, yy, tt))
        # for i, name in enumerate(df["Name"].to_list()):
        #     plt.annotate(name, (z[i], y[i]),fontsize=6)
        # plt.title("Correlation Plots")
        plt.xlabel(x_labels_dict[x_field])
        plt.ylabel("Normalized Number of New Diseases Listed in SIDER")
        plt.tight_layout()
        adjust_text(texts)
        # plt.show()
        plt.savefig(
            join(result_path,"corr_%s.png"%x_field),
            bbox_inches='tight')
        plt.clf()


def get_box_plot(df,cols=["ES_S_CM_CBOTH", "ES_NS_CM_CBOTH"]):
    ##please add box-and-whiskers plot comparing ES_S_CM_CBOTH and ES_NS_CM_CBOTH
    plot_df = df[cols]
    plot_df.columns=["Disease","New Disease"]
    
    plot_df = pd.melt(plot_df)
    plot_columns=["Entity", "Number of Diseases Listed in SIDER"]
    plot_df.columns = plot_columns
    # print(plot_df)
    # quit()
    sns.boxplot(x=plot_columns[0], y=plot_columns[1], data=plot_df)
    plt.savefig(
        join(result_path,"ES_S_NS_boxplot.png"),
        bbox_inches='tight',dpi=160)
    plt.clf()


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

            for col in ["N_DISTANCE","Q_DISTANCE"]:        
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

def plot_distance():
    # entity_dict={"D":"Drug", "ND":"New Drug", "S":"Disease", "NS":"New Disease"}
    distance_ade=pd.read_csv(
        join(singledrug_prefix,"distance_num_sde.csv"),dtype={"RXNORM":str}
    ).drop([9,12])
    distance_ade=distance_ade.assign(
        DRUG_NAME=distance_ade["RXNORM"].apply(
            lambda rxnorm: search_rxnorm(rxnorm).split(" ")[0]
        )
    )
    # print(distance_ade.head())
    # quit()
    # drug_names = [search_rxnorm(dir_num).split(" ")[0] for dir_num in  dirs]
    ES_N_col=[col for col in distance_ade.columns if 'ES_N' in col][:4]
    ES_Q_col=[col for col in distance_ade.columns if 'ES_Q' in col]
    # rename_col_dict=dict(zip(
        # ES_N_col+ES_Q_col, [col.replace("ES_","") for col in (ES_N_col+ES_Q_col)]
    # ))
    # ES_N_col=list(map(lambda col: col.replace("ES_",""),ES_N_col))
    # ES_Q_col=list(map(lambda col: col.replace("ES_",""),ES_Q_col))
    # distance_ade=distance_ade.rename(rename_col_dict,inplac)
    # print(distance_ade.columns)
    # quit()
    
    # num_legends=len(ES_N_col)+len(ES_Q_col)
    # distance_ade=distance_ade.loc[:,["RXNORM"]+ES_N_col+ES_Q_col]
    for col, symbol, color, size in zip(
        ES_N_col, ["X","p","P","D"], \
            ["lightcoral","limegreen","yellow","teal"],
            # ["green","limegreen","turquoise","teal"],
            [18,14,10,6]):
        plt.plot( 'DRUG_NAME', col, label=col.replace("ES_",""),data=distance_ade, \
            marker=symbol, color=color,linestyle=':',markersize=size)
    
    for col, symbol, color in zip(
        ES_Q_col,["o","^","*","d"],["orange","olive","red","blue"]
    ):
        plt.scatter(\
            'DRUG_NAME', col, data=distance_ade,label=col.replace("ES_",""),
            marker=symbol, color=color, s=50)
    
    plt.xticks(rotation=30, 
        fontsize=8,
        horizontalalignment='right')
    plt.legend()
    plt.tight_layout()

    # print(distance_ade)
    # print(ES_N_col)
    # print(ES_Q_col)
    plt.savefig(join(
        result_path,"wasserstein_distance.png"),dpi=600) 

    # distance_ade
    # distance_ade=get_all_distance_ade()
    # print(distance_ade)


def norm_distance_numsde():
    fname=join(singledrug_prefix ,"distance_num_sde.csv")
    df = pd.read_csv(fname,dtype={"RXNORM":str})
    ## TODO: modiefy patientcounts in get_all_distance_ade() 
    df["patientcounts"] = [7774,15459,4509,15869,3980,4865,5392,3719,11387,4480,7161,2154,4865]
    # df["Name"] = ["Hydralazine","Metoprolol","Vancomycin","Magnesium","Coumadin","Propofol","Plavix","Zosyn","Pantoprazole","Hydromorphone","Ipratropium","Furosemide","Coumadin"]
    df["Name"]=df["RXNORM"].apply(
        lambda rxnorm: search_rxnorm(rxnorm).split(" ")[0] 
    )

    # df["norm_SIDER"] = df["SIDER"]/df["patientcounts"]
    df['norm_ES_CM_CBOTH'] = df['ES_NS_CM_CBOTH']/df["patientcounts"]
    df['norm_ES_UN_CBOTH'] = df['ES_NS_UN_CBOTH']/df["patientcounts"]
    df['norm_ES_S_CM_CBOTH'] = df['ES_S_CM_CBOTH']/df["patientcounts"]
    df['norm_ES_S_UN_CBOTH'] = df['ES_S_UN_CBOTH']/df["patientcounts"]
    df = df.drop([9,12]) #some problem with Hydromorphone and didn't want Coumadin twice
    # print(df["RXNORM"].tolist())
    # quit()
    # NOTE: plot1
    # get_box_plot(df)
    # NOTE: plot2
    diff_corr_analysis(df)
    quit()
    print("Processed columns of distance_num_sde: {}".format(df.columns))

    for col in ['ES_S_CM_CBOTH', 'ES_NS_CM_CBOTH']:
        print("Description of Dataframe with column %s"%col)
        print(df[col].describe())



norm_distance_numsde()
# plot_distance()