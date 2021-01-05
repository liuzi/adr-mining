import os 
from utils._tools import *
from os.path import *
from utils._path import *

def run_autoencoder():



    nohup_cmd = "cd /home/liu/project/Clinic-Analysis/Scripts/Data_Process_Server/S2/utils && nohup python _five_dimension_reduction.py >> /home/liu/nohup/fiveauto_nohup.out 2>&1 &"
    os.system(nohup_cmd)
        # print([sub_dirs])
        
# HACK: uncomment code in next line
# run_autoencoder()






