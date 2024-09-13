import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold
from scipy.stats import pearsonr
from pathlib import Path
from liv2trans_functions import *
from matplotlib import pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import argparse

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Run different ML models')
parser.add_argument('--X_invitro_path', action='store',help='file to load pre-processed (replicate merged) in vitro molecular features')
parser.add_argument('--X_invivo_path', action='store', help='file to load pre-processed in vivo molecular features')
parser.add_argument('--Y_invivo_path', action='store', help='file to load in vivo phenotypes')
parser.add_argument('--res_dir', action='store', help='folder to save results of the analysis')
parser.add_argument('--num_LVs', action='store', type=int,help='number of latent variables in the PLSR model')
parser.add_argument('--num_folds', action='store', type=int,help='number of of splits for CV',default=10)
args = parser.parse_args()
num_LVs = int(args.num_LVs)
X_invitro_path = args.X_invitro_path
X_invivo_path = args.X_invivo_path
Y_invivo_path = args.Y_invivo_path
res_dir = args.res_dir
num_folds = int(args.num_folds)

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

# Create folds
kf = KFold(n_splits=num_folds)
train_r = []
val_r = []
train_r_backproj = []
val_r_backproj = []
train_r_extra = []
val_r_extra = []
print('Biging cross validation')
for i, (train_index, test_index) in enumerate(kf.split(Xinvivo_centered)):
    # select train and test data
    x_train = Xinvivo_centered[train_index,:]
    y_train = Yinvivo_centered[train_index,:]
    x_val = Xinvivo_centered[test_index,:]
    y_val = Yinvivo_centered[test_index,:]

    # fit model and evaluate in validation set
    model = PLSRegression(n_components=num_LVs,scale=False)
    model.fit(x_train,y_train)
    yhat_train = model.predict(x_train)
    yhat_val = model.predict(x_val)
    train_r.append(pair_pearsonr(y_train, yhat_train, axis=0).mean())
    val_r.append(pair_pearsonr(y_val, yhat_val, axis=0).mean())

    ### Get extra basis
    phi = model.coef_.T # in PLSR in python this directly saved in the model
    Wm_opt = extra_basis_analytical_solution(y=y_train,
                                             W_invitro = Wm,
                                             phi = phi)
    Wm_total = np.concatenate([Wm,Wm_opt],axis=1)

    # prediction after backprojection
    #first for training data
    Xback_train = x_train @  Wm @ Wm.T
    yhat_train_backproj = model.predict(Xback_train)
    train_r_backproj.append(pair_pearsonr(y_train, yhat_train_backproj, axis=0).mean())
    # repeat for valdiation
    Xback_val = x_val @  Wm @ Wm.T
    yhat_val_backproj = model.predict(Xback_val)
    val_r_backproj.append(pair_pearsonr(y_val, yhat_val_backproj, axis=0).mean())

    # prediction after backprojection with Wm_total
    #first for training data
    Xextra_train = x_train @  Wm_total @ Wm_total.T
    yhat_train_extra = model.predict(Xextra_train)
    train_r_extra.append(pair_pearsonr(y_train, yhat_train_extra, axis=0).mean())
    # repeat for valdiation
    Xextra_val = x_val @  Wm_total @ Wm_total.T
    yhat_val_extra = model.predict(Xextra_val)
    val_r_extra.append(pair_pearsonr(y_val, yhat_val_extra, axis=0).mean())
    print('Finished fold '+str(i))

res_val = pd.DataFrame({'human genes':val_r,'back-projected':val_r_backproj,'optimized MPS':val_r_extra})
res_val['fold'] = [xx for xx in range(num_folds)]
res_val['set'] = 'test'
res_train = pd.DataFrame({'human genes':train_r,'back-projected':train_r_backproj,'optimized MPS':train_r_extra})
res_train['fold'] = [xx for xx in range(num_folds)]
res_train['set'] = 'train'
df_res_cv = pd.concat([res_val,res_train])
df_res_cv.to_csv(res_dir+'df_res_cv.csv')

## Visualize
df_res_cv = df_res_cv.melt(id_vars=['fold','set'],var_name='input', value_name='r')

# Create the plot
sns.set_context("talk", font_scale=1.2)
# Create a boxplot split into panels by the 'set' column
g = sns.catplot(
    data=df_res_cv,
    x='input',
    y='r',
    hue='input',
    col='set', 
    kind='box',
    height=6,
    aspect=1.2
)
g.set_titles("Set: {col_name}")
g.set(ylim=(np.min([0,df_res_cv['r'].min()]), 1))  # Example: setting y-axis limits from 0 to 1
g.set(yticks=np.arange(0, 1.1, 0.1))  # Example: custom ticks
# Add horizontal grid lines to each subplot
for ax in g.axes.flatten():
    ax.grid(True, which='major', axis='y', linestyle='--')

g.set_axis_labels("Input features", "Pearson`s r")

# Save the figure
g.savefig(res_dir + 'CV_performance.pdf')

print('Saved results and done!')

