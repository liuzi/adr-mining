import numpy as np
import pandas as pd
from os.path import join
from utils._path import singledrug_featurepreprocess_prefix, concat_clamp_prefix
from utils._tools import read_data, left_join, write2file, inner_join, append_csv_byrow, append_csv_bydf
from term_process.umls_api import search_rxnorm, get_side_effects_cui, search_cid_for_rxnorm
import seaborn as sns
import matplotlib.pyplot as plt