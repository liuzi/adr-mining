from bs4 import BeautifulSoup, SoupStrainer
from utils._tools import *
import os, sys
import glob
from os.path import join
import requests
from urllib.parse import urlparse, urljoin
from itertools import chain
import numpy as np
import pandas as pd
import re


def filter_icd9link(link_list,filter_pattern):
    return list(filter(lambda x: re.match(filter_pattern,x),link_list))

def crawl_dexur(
    url="https://dexur.com/icd9/",
    slash_num=2,
    filter_pattern="/icd9/[A-Z0-9]*-[A-Z0-9]*"):

    url_prefix=url.rsplit('/',slash_num)[0]
    ## Make a GET request to fetch the raw HTML content
    html_content = requests.get(url).text
    ## Parse the html content
    soup = BeautifulSoup(html_content, "lxml")
    # print(soup.prettify()) # print the parsed data of html
    links_text= [
        "{}".format(link.get("href")) 
        for link in soup.find_all('a')]

    icd9_cat_links= [
        "{}{}".format(url_prefix,icd9_cat) 
        for icd9_cat in filter_icd9link(links_text,filter_pattern)]
    
    return icd9_cat_links

## TODO:digui
# def search_end_link(sup_link,icd9_len):

# TODO: better to use DI GUI to iterrally search all end links of tables firstly
def save_dexur_df(
    url="https://dexur.com/icd9/",
    result_path="/data/liu/mimic3/EXTERNEL_DATASET/dexur_ICD9_DISEASENAME"):
    create_folder_overwrite(result_path)

    for sup_link in crawl_dexur(url):
        code_name_all_df=[]
        sup_icd9_range=sup_link.rsplit("/",2)[-2].replace("-","_")
        # e.g, 001-139/
        print("Processing for SUP link: {}...".format(sup_link))
        sub1_link_list=crawl_dexur(sup_link,3)

        # NOTE: [1:] must be include for sub1 link list
        for sub1_link in sub1_link_list[1:]:
            ## e.g., icd9/001-009/
            # print("SUB1 link: {}".format(sub1_link))

            sub2_link_list=crawl_dexur(sub1_link,3,"/icd9/[A-Z0-9]{3,}/")
            # print(sub2_link_list[:1])
            for sub2_link in sub2_link_list:
                # e.g., icd9/001/
                # print("SUB2 link: {}".format(sub2_link))
                html_content = requests.get(sub2_link).text
                ## Parse the html content
                soup = BeautifulSoup(html_content, "lxml")

                if(soup.find_all('table')):
                    table=soup.find_all('table')[0]
                    table_title=table.find("tr").text.strip()
                    if("ICD Code" in table_title):
                        sub3_link_list=crawl_dexur(sub2_link,3,"/icd9/[A-Z0-9]{4,}/")
                        # sub3_link_list=["https://dexur.com/icd9/4282/"]
                        # e.g., icd9/0010/
                        # print("SUB3 link: {}".format(sub3_link_list[0]))
                        code_name_ndarray_list=[]
                        for sub3_link in list(set(sub3_link_list)):
                            html_content = requests.get(sub3_link).text
                            soup = BeautifulSoup(html_content, "lxml")
                            if(soup.find_all('table')):
                                sub_table=soup.find_all('table')[0]
                                sub_table_title=sub_table.find("tr").text.strip()
                                sub_rows=[i.text.strip() for i in sub_table.find_all('td')]
                                if(len(crawl_dexur(sub3_link,3,"/icd9/[A-Z0-9]{5,}/"))):
                                     ## if sub3 icd9 has subsequent icd9 codes
                                    code_name_ndarray_list.append(np.reshape(sub_rows,(-1,2)))
                                else:
                                    # NOTE:
                                    # code_name_ndarray = np.reshape(rows,(-1,2))             
                                    code_name_ndarray_list.append(np.array([list(
                                        map(lambda x: x.strip(),sub_table_title.split('-')[:2]))]))
                            else:
                                title=soup.find_all('b')[0].text
                                if("ICD 9 Diagnosis Code:" in title):
                                    titles=title.split(":",1)[-1].split("-")[:2]
                                    code_name_ndarray_list.append(np.array([list(
                                            map(lambda x: x.strip(),titles))]))
                        if(len(code_name_ndarray_list)):
                            code_name_ndarray=np.vstack(code_name_ndarray_list)
                        else:
                            break
                    else:
                        code_name_ndarray=[list(
                            map(lambda x: x.strip(),table_title.split('-')[:2]))]
                    code_name_df=pd.DataFrame(
                        code_name_ndarray,columns=["ICD9_CODE","SHORT_TITLE"])
                    code_name_all_df.append(code_name_df)
                else:
                    title=soup.find_all('b')[0].text
                    if("ICD 9 Diagnosis Code:" in title):
                        titles=title.split(":",1)[-1].split("-")[:2]
                        code_name_ndarray_list.append(np.array([list(
                                map(lambda x: x.strip(),titles))]))

        code_name_all_df=pd.concat(code_name_all_df,axis=0)
        write2file(code_name_all_df,join(result_path,"ICD9_NAME%s"%sup_icd9_range))
        # print(code_name_all_df.head(20))

def concat_dexur_df(result_path="/data/liu/mimic3/EXTERNEL_DATASET/dexur_ICD9_DISEASENAME"):

    files = [f for f in glob.glob(join(result_path,"*[0-9]*.csv"))]
    sub_icd9_name_dfs=list(map(lambda file: read_data(file,dtype=str), files))
    all_icd9_name_df=pd.concat(sub_icd9_name_dfs,axis=0)
    create_folder_overwrite(join(result_path,"CONCAT_RESULTS"))
    write2file(all_icd9_name_df,join(result_path,"CONCAT_RESULTS","ICD9_NAME_ALL"))

def combine_idc9_name_all(result_path="/data/liu/mimic3/EXTERNEL_DATASET/dexur_ICD9_DISEASENAME"):
    v0_resultpath="%s_v0"%result_path
    icd9_all_dflist=list(map(
        lambda file_path:read_data(join(file_path,"CONCAT_RESULTS","ICD9_NAME_ALL"),dtype=str),
        [v0_resultpath,result_path]))
    # icd9_all_v0=read_data(join(v0_resultpath,"CONCAT_RESULTS","ICD9_NAME_ALL"))
    icd9_alls=pd.concat(icd9_all_dflist,axis=0).drop_duplicates()
    write2file(icd9_alls,join(result_path,"CONCAT_RESULTS","ICD9_NAME_ALL_COMBINE"))


if __name__ == '__main__':
    result_path="/data/liu/mimic3/EXTERNEL_DATASET/dexur_ICD9_DISEASENAME"
    # save_dexur_df()
    # concat_dexur_df()
    combine_idc9_name_all()


    # # # NOTE:test for one icd9 cat or icd9 code
    # # # test_url="https://dexur.com/icd9/024/"
    # # test_url="https://dexur.com/icd9/4282/"
    # test_url="https://dexur.com/icd9/0011/"

    # html_content = requests.get(test_url).text
    # soup = BeautifulSoup(html_content, "lxml")
    # print(soup.find_all('b')[0].text.split(":",1)[-1].split("-")[:2][1].strip())
    # print(soup.prettify())
    # if(len(crawl_dexur(test_url,3,"/icd9/[A-Z0-9]{5,}/"))):
    #     print("yes")
    # table=soup.find_all('table')[0]
    # table_title=table.find("tr").text.strip()
    # if("ICD Code" in table_title):
    #     print(table_title)
    # else:
    #     code_name_ndarray=[list(
    #         map(lambda x: x.strip(),table_title.split('-')[:2]))]
    #     code_name_df=pd.DataFrame(code_name_ndarray,columns=["ICD9_CODE","SHORT_TITLE"])
    #     print(code_name_df)
    # NOTE:test for one icd9 cat or icd9 code








# print(link_text)
    # print("Inner Text: {}".format(link.text))
    # print("Title: {}".format(link.get("title")))
    # print("href: {}".format(link.get("href")))
    # links_text=list(link.get("href"))
    # print(links_text[0)
    # icd9_links= filter(lambda x: "")
    # [link for link in links_text if ("/icd9/" in link)]
    # print(chain(*links_text))
    # for i in tt:
    #     print(i)
    # print(len(tt))
    # links=[link for link in "{}".format(link.get("href")) if ("/icd9/" in link)]
    # print(link)


