import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from scipy.stats import pearsonr
from pathlib import Path
from livivtra_functions import *
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Run different ML models')
parser.add_argument('--X_invitro', action='store',help='file to load pre-processed in vitro molecular features')
parser.add_argument('--X_invivo', action='store', help='file to load pre-processed in vivo molecular features')
parser.add_argument('--Y_invivo', action='store', help='file to load in vivo phenotypes')
parser.add_argument('--res_dir', action='store', help='folder to save results of the analysis')
parser.add_argument('--num_LVs', action='store', type=int,help='number of latent variables in the PLSR model',default=8)
parser.add_argument('--num_folds', action='store', type=int,help='number of of splits for CV',default=10)
args = parser.parse_args()
num_LVs = int(args.num_LVs)
X_invitro = args.X_invitro
X_invivo = args.X_invivo
Y_invivo = args.Y_invivo
res_dir = args.res_dir
num_folds = int(args.num_folds)

### Load the data

### Fit PLSR model
### SPLIT WIT SKLEARN RANDOMLY
model = PLSRegression(n_components=num_LVs,scale=False)


