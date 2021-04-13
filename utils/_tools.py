import pandas as pd
import os
# import errno
import json
import shutil
from pathlib import Path
from csv import writer

def csv_suffix(file_path,suffix = '.csv'):
    if('.' in file_path):
        final_file_path=file_path
    else:
        final_file_path = file_path+suffix
    return final_file_path

def read_data(file_path, dtype=None, usecols=None, sep=',', header = 'infer', suffix = '.csv', pre=''):
    return pd.read_csv(csv_suffix(file_path), dtype=dtype, usecols=usecols,sep=sep, header = header, encoding='latin1')

def write2txt(string, file_path):
    textfile = open("%s.txt"%file_path, 'w')
    textfile.write(string)
    textfile.close()
    
def write2json(jdata,file_path):   
    with open("%s.json"%file_path, 'w') as fp:
        json.dump(jdata, fp)
    
# def create_folder(dir):

#     if not os.path.exists(dir):
#         try:
#             os.makedirs(dir)
#         except OSError as exc: # Guard against race condition
#             if exc.errno != errno.EEXIST:
#                 raise

def create_folder(path):
    try:
        os.mkdir(path)
    except OSError:
        print ("Directory %s already exists" % path)
    else:
        print ("Successfully create the directory %s" % path)

def create_folder_overwrite(file_path):
    if os.path.exists(file_path):
        print("%s already exists. Clear and re-create this folder."%file_path)
        shutil.rmtree(file_path)
    os.makedirs(file_path)

def write2file(df, file_path):
    df.to_csv(csv_suffix(file_path), index=False)

def write2file_nooverwrite(df, file_path):
    if os.path.exists((csv_suffix(file_path))):
        print("%s already exists."%csv_suffix(file_path))
    else:
        df.to_csv(csv_suffix(file_path), index=False)
        print("%s is successfully saved."%csv_suffix(file_path))


def left_join(left_df, right_df, joined_field):
    return pd.merge(left_df, right_df, how='left', on=joined_field)

def inner_join(left_df, right_df, joined_field):
    return pd.merge(left_df, right_df, how='inner', on=joined_field)

def print_patient_stats(df):
    
    print("# of rows: %d"%len(df))
    
    for col in df.columns:
        print("# of %s: %d"%(col,len(df[col].unique())))
        # print("# of %s: %d"%len(df['HADM_ID'].unique()))

def append_csv_byrow(row, file_name,sep="\t"):
    # NOTE: row:[]
    # Open file in append mode
    with open(csv_suffix(file_name), 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = writer(write_obj,delimiter=sep )
        # Add contents of list as last row in the csv file
        csv_writer.writerow(row)
        # print("Insert row:\"%s\" into file %s"%(
        #     row,df.shape[0],file_name))


def append_csv_bydf(df, file_name, sep="\t"):
    
    # Open file in append mode
    with open(csv_suffix(file_name), 'a') as write_obj:
        df.to_csv(write_obj, mode='a', header=write_obj.tell()==0,sep=sep,index=False)
        print("Insert %d rows into file %s"%(
            df.shape[0],file_name))
