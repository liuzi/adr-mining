import numpy as np
import pandas as pd
import re
import io

from os.path import *
# from S2.utils.tools import *
from utils._tools import *

read_prefix = "/data/MIMIC3/"
dis_sum_prefix = "/data/liu/mimic3/NOTESEVENT_EXTRACTION/ALL"
write_prefix="/data/liu/mimic3/CLAMP_NER/input"
rule_dissum_file = "ReRule%d_Discharge summary"

def clean_string(string):
    if(string):
        return string.strip('\n').lower()
    else: 
        return None

## Revised Section Extraction Script
def re_split(text, regex = r"\n\n[A-Z][A-Za-z/ ]+:"):
    """Split Discharge Summary to mulitple sections according to respective titles using regex rules

    Return:
    ------
    {titles:sections}: dict
        titiles : list. 
            1) Remove enter character
            2) Transfer all titles to lower case, making it easier to check and match
        line_number : int.
            line number of each title
        sections: list. 
            len(sections)=len(titles)+1. The first section doesn't have title as it is not beginning 
            with "\n"       
    """
    titles = re.findall(regex, text)
    sections = re.split(regex, text)[1:]
    sections_lines = [len(section.splitlines()) for section in sections]
    # split_lines=text.splitlines()
    # TODO:
    ## BUG: some line number mismatch to title
    title_search_index=0
    title_length=len(titles)
    title_linenums=[None]*title_length

    if(titles):
        with io.StringIO(text) as f:
            for num, line in enumerate(f, start=1):
                # print('{} : {}'.format(i, line.strip()))
                # print("line_num: %d"%num)
                if(titles[title_search_index].strip('\n') in line):
                    # print("find line number %d for title %s"%(num,titles[title_search_index].strip('\n')))
                    title_linenums[title_search_index]=num
                    title_search_index = title_search_index+1   
                    if(title_search_index==title_length): break;
                else: continue;        

    # # Prvious version
    # title_linenums = [num+1 for num, line in zip(range(len(split_lines)),split_lines) 
    #                     if any(title.strip('\n') in line for title in titles)]

    next_titleline = title_linenums[1:]+[None]
    next_title = titles[1:]+[None]
    section_endline = [x+y for x, y in zip (title_linenums, sections_lines)]

    ## 1) transform extracted titles to lower case  2)zip line numbers and sections
    # if len(titles):
    #     return dict(zip(map(clean_string,titles), zip(title_linenums, sections)))
    # else:
    #     return None
    

    return dict(zip(map(clean_string,titles), zip(
        title_linenums, 
        sections, 
        section_endline,
        next_titleline,
        map(clean_string,next_title))))

    
    
def extract_sections(text,rule_index = 0, rerules = [
    r"\n\n[A-Z][A-Z/ ]+:",
    r"\n[A-Z][A-Z/ ]+:",
    r"\n\n[A-Z][A-Za-z/ ]+:",
    r"\n[A-Z][A-Za-z/ ]+:"
], selected_titles = [["hospital course"],["medical history"],["medication","admission"]]):
    """

    Arguments:
    ---------
    selected_titles: list
        1)hospital course: brief hc
        2)medical history:
        3)medications on admission: admission medications
    
    Step:
    -----
    1) split original text to multiple sections with titles: dict{title -> section}
    2) select key sections by matching their titles to selected titles
    """
    section_dict = re_split(text, rerules[rule_index])  
    # return section_dict

    # if(section_dict):
    #     ## selected title -> matched title
    #     titles_dict =dict([(selected_title, title) for title in section_dict.keys() for selected_title in selected_titles if selected_title in title])
    #     ## get actually extracted titles that are matched to selected titles
    #     section_titles = [titles_dict.get(title, None) for title in selected_titles]
    #     ## unzip line numbers and sections       
    #     line_numbers, sections = list(zip(*[section_dict.get(title, (None,None)) for title in section_titles]))
    #     # return section_titles
    #     return tuple([(',').join(filter(None, section_titles)),
    #         (',').join(filter(None, map(str,line_numbers)))]+list(sections))

    # else:
    #     return tuple([None]*5)

    ## selected title -> matched title
    # BUG: medication -> [medication, admission]
    ## get actually extracted titles that are matched to selected titles
    len_selected_titles = len(selected_titles)
    section_titles=[None]*len_selected_titles
    # dict(zip(map((";").join,test),zip(range(3),test)))
    search_selected_titles = dict(
        zip(map((" ").join, selected_titles), zip(
            range(len_selected_titles), selected_titles)))

    for extract_title in section_dict.keys():
        for key in list(search_selected_titles):
            value = search_selected_titles[key]
            if all(word in extract_title for word in value[1]):
                # print("find extracted title %s for %d selected title"%(
                #     extract_title,value[0]))
                section_titles[value[0]] = extract_title
                search_selected_titles.pop(key)
    # section_titles=dict(zip(map((";").join, selected_titles),extract_selected_titles))
    # # HACK: previous version
    # titles_dict =dict([
    #     (selected_title, title) 
    #     for title in section_dict.keys() 
    #     for selected_title in selected_titles 
    #     if selected_title in title])

    # # HACK: previous version
    # section_titles = [titles_dict.get(title, None) for title in selected_titles]
    ## unzip line numbers and sections       
    title_line, sections, section_endline, next_titleline, next_title= list(
        zip(*[section_dict.get(title, (None,None,None,None,None)) for title in section_titles]))
    # return section_titles
    return tuple([(',').join(filter(None, section_titles)),
        (',').join(filter(None, map(str,title_line))),
        (',').join(filter(None, map(str,section_endline))),
        (',').join(filter(None, map(str,next_titleline))),
        (',').join(filter(None, next_title))]+list(sections))


def save_keysessions_byfile(
    notes_discharge_df,
    rule_index,
    file="%s_%s_%d",
    folder_titles = ['HOSPITAL COURSE', 'MEDICAL HISTORY', 'MEDICATIONS ON ADMISSION']):
    
    write_rule_prefix = join(write_prefix, "Rule%d"%rule_index)
    create_folder_overwrite(write_rule_prefix)
    tilfolder_path = [join(write_rule_prefix, title.replace(" ","_")) for title in folder_titles]
    [create_folder_overwrite(path) for path in tilfolder_path]
    title_path_dict = dict(zip(folder_titles, tilfolder_path))
    
    for subject_id, group in notes_discharge_df.groupby(['SUBJECT_ID']):
        hadm_ids = group['HADM_ID'].unique()
        ## create folder ~/subject_id/hadm_id
        cols = ['GROUP_RANK','HOSPITAL COURSE', 'MEDICAL HISTORY', 'MEDICATIONS ON ADMISSION']
        
        for hadm_id, subgroup in group.groupby(['HADM_ID']):    
            
            for _, row in subgroup[cols].iterrows():
                file_name = file%(str(subject_id), hadm_id, row['GROUP_RANK'])
                [write2txt(row[title], join(
                    title_path_dict[title],file_name.replace(" ","_"))) for title in folder_titles if row[title]]
                # write2txt(row["HOSPITAL COURSE"],join(course_pre,file_name))
                # write2txt(row["MEDICAL HISTORY"],join(pastdrug_pre,file_name))
                # write2txt(row["MEDICATIONS ON ADMISSION"],join(admi_medi_pre,file_name))        

def run_rule(notes_discharge_df, rule_index):
    # notes_discharge_df=read_data(join(dis_sum_prefix,"Discharge summary_All"),dtype={"HADM_ID":str})
    add_cols = ["TITLE", "TTITLE_LINE", "SECTION_END_LINE","NEXT_TITLE_LINE","NEXT_TITLE","HOSPITAL COURSE", "MEDICAL HISTORY", "MEDICATIONS ON ADMISSION"]
    notes_discharge_df[add_cols] = pd.DataFrame(notes_discharge_df['TEXT'].apply(lambda x: 
                                        extract_sections(x,rule_index)).tolist())
    # print(notes_discharge_df)
    if(notes_discharge_df['TITLE'].isnull().values.all()):
        print("Rule %d is useless"%rule_index)

    #TODO: uncomment   
    write2file(notes_discharge_df,join(write_prefix,("%s_All"%rule_dissum_file)%rule_index))
    save_keysessions_byfile(notes_discharge_df,rule_index)
    # return notes_discharge_df

def get_rule_result(rule_index):
    ## check whetehr the input file (output by last rule) exists
    # if isfile(join(write_prefix, ("%s_All.csv"%rule_dissum_file)%(rule_index))):
    #     notes_discharge_df=read_data(join(write_prefix,("%s_All"%rule_dissum_file)%(rule_index)),dtype={"HADM_ID":str})

    if(rule_index==0):
        run_rule(read_data(join(dis_sum_prefix,"Discharge summary_All"),dtype={"HADM_ID":str}),rule_index)
        
    else:
        run_rule(read_data(join(write_prefix,("%s_Remain"%rule_dissum_file)%(rule_index-1)),dtype={"HADM_ID":str}),rule_index)
    
    notes_discharge_df = read_data(join(write_prefix,("%s_All"%rule_dissum_file)%rule_index),dtype={"HADM_ID":str})
    rule_only_df, remain_df= [x for _, x in notes_discharge_df.groupby(notes_discharge_df['TITLE'].isna())]
    print("Rule %d :\nrule only df row num %d \n episode num %d \n\nremain df row num %d \n episode num %d \n\n"%(
        rule_index,len(rule_only_df),len(rule_only_df['HADM_ID'].unique()),len(remain_df),len(remain_df['HADM_ID'].unique())))
        
    write2file(rule_only_df,join(write_prefix,("%s_Only"%rule_dissum_file)%(rule_index)))
    write2file(remain_df,join(write_prefix,("%s_Remain"%rule_dissum_file)%(rule_index)))


"""
    test code: remain, 4 rules
"""
def main():
    for rule_index in range(0,4):
        # notes_discharge_df=read_data(join(write_prefix,("%s_Remain"%rule_dissum_file)%(rule_index-1)),dtype={"HADM_ID":str})
        # notes_discharge_df_rule1 = run_rule(notes_discharge_df, rule_index=1)
        get_rule_result(rule_index)

# main()


"""
    test added columns
"""
# notes_discharge_df=read_data(join(dis_sum_prefix,"Discharge summary_All"),dtype={"HADM_ID":str})
# # print(len(notes_discharge_df))
# notes_discharge_df.head()
# notes_df=notes_discharge_df
# # test_df = notes_df[notes_df['HADM_ID']=="196600"]
# test_df = notes_df.head()

# # test_df_text = test_df.loc[36180,"TEXT"]
# # print(len(test_df))
# # print(extract_sections(test_df.loc[36180,"TEXT"],rule_index=2))
# run_rule(test_df,rule_index=2)



"""
tt
"""
# notes_discharge_df.head()
# notes_df=notes_discharge_df
# test_df = notes_df[notes_df['HADM_ID']=="137810"]
# test_df_text = test_df.loc[36180,"TEXT"]
# titles = extract_sections(test_df_text,3)
# print(len(titles))


"""
    test code: one text, rule 0
"""
# notes_discharge_df=read_data(join(dis_sum_prefix,"Discharge summary_All"),dtype={"HADM_ID":str})
# print(len(notes_discharge_df))
# notes_discharge_df.head()
# notes_df=notes_discharge_df
# test_df = notes_df[notes_df['HADM_ID']=="196600"]
# test_df_text = test_df.loc[36180,"TEXT"]
# titles = extract_sections(test_df_text,3)
# print(len(titles))
