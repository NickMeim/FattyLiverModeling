import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from pathlib import Path
from livivtra_functions import *
import pickle
import argparse

import decoupler as dc


### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Run different ML models')
parser.add_argument('--X_invitro_path', action='store',help='file to load pre-processed (replicate merged) in vitro molecular features')
parser.add_argument('--X_invivo_path', action='store', help='file to load pre-processed in vivo molecular features')
parser.add_argument('--Y_invivo_path', action='store', help='file to load in vivo phenotypes')
parser.add_argument('--num_LVs', action='store', type=int,help='number of latent variables in the PLSR model')
parser.add_argument('--res_dir', action='store', help='folder to save results of the analysis')
args = parser.parse_args()
num_LVs = int(args.num_LVs)
X_invitro_path = args.X_invitro_path
X_invivo_path = args.X_invivo_path
Y_invivo_path = args.Y_invivo_path
res_dir = args.res_dir

Path(res_dir).mkdir(parents=True, exist_ok=True)

### Load the data
Y_invivo = pd.read_csv(Y_invivo_path,index_col=0).to_numpy()
X_invivo = pd.read_csv(X_invivo_path,index_col=0).to_numpy()
X_invitro = pd.read_csv(X_invitro_path,index_col=0)
genes = X_invitro.columns
X_invitro = X_invitro.to_numpy()

### Get PCA matrix of the grouped in vitro data
pca_invitro = PCA()
pca_invitro.fit(X_invitro)
Wm = pca_invitro.components_.T
Wm = Wm[:,:-1]

# Manually center Xh and Yh
Xinvivo_centered = X_invivo - np.mean(X_invivo, axis=0)
Yinvivo_centered = Y_invivo - np.mean(Y_invivo, axis=0)

### Fit PLSR model
model = PLSRegression(n_components=num_LVs,scale=False)
model.fit(Xinvivo_centered,Yinvivo_centered)

### Save the model
# Save the model to a file
with open(res_dir+'model.pkl', 'wb') as f:
    pickle.dump(model, f)

### Get extra basis
phi = model.coef_.T # in PLSR in python this directly saved in the model
Wm_opt = extra_basis_analytical_solution(y=Y_invivo,
                                         W_invitro = Wm,
                                         phi = phi)
Wm_opt = pd.DataFrame(Wm_opt)
Wm_opt.index = genes
Wm_opt.columns=['V1','V2']
Wm = pd.DataFrame(Wm)
Wm.index = genes
Wm.columns= ['PC'+str(i) for i in range(1,Wm.shape[1]+1)]

W_total = pd.concat([Wm,Wm_opt],axis=1)
### Save the extra basis and the PC loadings
W_total.to_csv(res_dir+'W_total.csv')

### Interpret extra basis pathways
progeny = dc.get_progeny(organism='human', top=500)
progeny = progeny.loc[:,['source','target','weight']]
# pathway_scores= dc.decouple(
#     methods= ['mlm', 
#               'ulm', 
#               'viper', 
#               'wsum'],
#     mat=Wm_opt.T,
#     net=progeny,
#     source='source',
#     target='target',
#     weight='weight',
#     min_n = 5,
#     verbose=True
# )
# pathway_scores = pathway_scores['consensus_estimate']
pathway_scores,_ = dc.run_viper(
    mat=Wm_opt.T,
    net=progeny,
    source='source',
    target='target',
    weight=None,
    min_n = 5,
    verbose=True
)

pathway_scores = pathway_scores.reset_index().rename(columns={'index':'condition'}).melt(id_vars='condition',var_name='Pathway', value_name='activity')

## Plot and save barplots
fig_extra_basis_1 = visualize_paths_activity(pathway_scores,extra_basis = 1) 
fig_extra_basis_1.savefig(res_dir+'pathways_extra_basis_1.pdf')

fig_extra_basis_2 = visualize_paths_activity(pathway_scores,extra_basis = 2) 
fig_extra_basis_2.savefig(res_dir+'pathways_extra_basis_2.pdf')