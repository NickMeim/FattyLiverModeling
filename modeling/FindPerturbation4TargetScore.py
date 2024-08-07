import torch
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import joblib
from scipy.stats import pearsonr
import pyreadr
from pathlib import Path
import argparse

def pearson_r(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = torch.mean(x, dim=0)
    my = torch.mean(y, dim=0)
    xm, ym = x - mx, y - my
    r_num = torch.sum(xm * ym,dim=0)
    x_square_sum = torch.sum(xm * xm,dim=0)
    y_square_sum = torch.sum(ym * ym,dim=0)
    r_den = torch.sqrt(x_square_sum * y_square_sum)
    r = r_num / r_den
    return r #torch.mean(r)

def pair_pearsonr(x, y, axis=0):
    mx = np.mean(x, axis=axis, keepdims=True)
    my = np.mean(y, axis=axis, keepdims=True)
    xm, ym = x-mx, y-my
    r_num = np.add.reduce(xm * ym, axis=axis)
    r_den = np.sqrt((xm*xm).sum(axis=axis) * (ym*ym).sum(axis=axis))
    r = r_num / r_den
    return r

def getSamples(N, batchSize):
    order = np.random.permutation(N)
    outList = []
    while len(order)>0:
        outList.append(order[0:batchSize])
        order = order[batchSize:]
    return outList

def L2Regularization(deepLearningModel, L2):
    weightLoss = 0.
    biasLoss = 0.
    for layer in deepLearningModel:
        if isinstance(layer, torch.nn.Linear):
            weightLoss = weightLoss + L2 * torch.sum((layer.weight)**2)
            biasLoss = biasLoss + L2 * torch.sum((layer.bias)**2)
    L2Loss = biasLoss + weightLoss
    return(L2Loss)

### Initialize the parsed arguments
parser = argparse.ArgumentParser(description='Move in vitro samples along the extra basis')
parser.add_argument('--initial_condition', metavar='N', type=str, nargs='*', help='stimuli name in the MPS to use as intialization',default=['Lean_TGF-B'])
parser.add_argument('--meta_data_location', action='store',help='location of files used in CV training of originalPLSR model',default='../preprocessing/')
parser.add_argument('--centered_data_location', action='store',help='location of external clinical files',default='../preprocessing/')
parser.add_argument('--in_vitro_matrices_loc', action='store',help='location of gene weights/loadings/extra basis matrices of the in-vitro data',default='../results/')
parser.add_argument('--in_vitro_name', action='store',help='name of the in-vitro data',default='Kostrzewski')
parser.add_argument('--in_vivo_name', action='store',help='name of the human data',default='govaere')
parser.add_argument('--phenotype', action='store',help='phenotype of interest',default='NAS')
parser.add_argument('--iters', action='store', type=int,help='number of iterations for optimization',default=100)
parser.add_argument('--learning_rate', action='store',help='learning rate for the optimization process',default=0.1)
parser.add_argument('--weight_decay', action='store',help='weight decay for the optimization process',default=0)
parser.add_argument('--l1_reg', action='store',help='L1 regularization coefficient',default=1e-6)
args = parser.parse_args()
initial_condition = args.initial_condition
meta_data_location= args.meta_data_location
centered_data_location = args.centered_data_location
in_vitro_matrices_loc = args.in_vitro_matrices_loc
in_vitro_name = args.in_vitro_name
in_vivo_name = args.in_vivo_name
phenotype = args.phenotype
iters = args.iters
learning_rate = args.learning_rate
wd = args.weight_decay
l1_reg = args.l1_reg

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print('Using device:', device)

### Load the data
Wm_total = pyreadr.read_r(in_vitro_matrices_loc+'Wm_'+in_vitro_name+'_total.rds')
Wm_total = Wm_total[None]
Wm = Wm_total.iloc[:,:-2]
Wm_opt = Wm_total.iloc[:,-2:]
## conver to numpy
Wm = Wm.to_numpy()
Wm_total = Wm_total.to_numpy()
Wm_opt = Wm_opt.to_numpy()

Xm = pyreadr.read_r(centered_data_location+in_vitro_name+'_centered_data.rds')
Xm = Xm[None]
meta_data = pyreadr.read_r(meta_data_location+in_vitro_name+'_meta_data.rds')
meta_data = meta_data[None]
treatments = np.unique(meta_data.treatment)

### Load human data and calculate the centroids in the extra basis
Xh = pyreadr.read_r(centered_data_location+in_vivo_name+'_centered_data.rds')
Xh = Xh[None]
Zh = Xh.to_numpy() @ Wm_opt
Zh = pd.DataFrame(Zh,index=Xh.index,columns=['Z'+str(i) for i in range(Zh.shape[1])])
Yh = pyreadr.read_r(centered_data_location+in_vivo_name+'_phenotype_data.rds')
Yh = Yh[None]
human_projected_data = pd.concat([Zh.reset_index(drop=True),Yh.loc[:,phenotype].reset_index(drop=True)],axis=1)
human_projected_data.index = Zh.index
centroids = human_projected_data.groupby(phenotype).mean()
if (phenotype == 'NAS'):
    centroids.iloc[:,1] = 0.
else:
    centroids.iloc[:,0] = 0.

Wm = torch.tensor(Wm).to(device)
Wm_total = torch.tensor(Wm_total).to(device)
Wm_opt = torch.tensor(Wm_opt).to(device)

res_all = pd.DataFrame()
Z_pert = pd.DataFrame()
for treatment in treatments:
    # plt.figure()
    ind = meta_data[meta_data.treatment == treatment].sampleName
    X_init = Xm.loc[ind,:]
    X_init = torch.tensor(X_init.to_numpy())
    if (ind.shape[0]==1):
        X_init = X_init.unsqueeze(0)
    X_init = X_init.to(device)
    for j in range(centroids.shape[0]):
        centroid = centroids.iloc[j,:]
        centroid = torch.tensor(centroid.to_numpy())
        centroid = centroid.to(device)
        dX = torch.nn.Parameter(torch.tensor(torch.randn(1,X_init.shape[1])).float().to(device), requires_grad = True)
        dX.retain_grad()
        optimizer = torch.optim.Adam([dX],lr= learning_rate,weight_decay=wd)
        all_loss = []  
        for i in range(iters):
            optimizer.zero_grad()
            X = X_init + dX
            Z = X @ Wm_opt
            loss = torch.sum(torch.square(Z - centroid),0).mean() + l1_reg * torch.sum(torch.abs(dX))
            loss.backward()
            optimizer.step()
            all_loss.append(loss.item())
            # if (i%50 == 0):
            #     print('Treatment: {}, {} = {} ,Iteration {}/{} : Loss = {}'.format(treatment,phenotype,centroids.index[j],i+1,iters,loss.item()))
        all_dx = dX.detach().cpu().numpy()
        X = X_init + dX
        Z = X @ Wm_opt
        dist = torch.sum(torch.square(Z - centroid),0).mean()
        dist = dist.item()
        res = pd.DataFrame(all_dx)
        res.columns = Xm.columns
        res['treatment'] = treatment
        res[phenotype] = centroids.index[j]
        res['mean_distance'] = dist
        res_all = res_all.append(res)
        tmp = pd.DataFrame(Z.detach().cpu().numpy(),
                                           index=Xm.loc[ind,:].index,
                                           columns=['Z'+str(i) for i in range(Z.shape[1])])
        tmp['treatment'] = treatment
        tmp[phenotype]=centroids.index[j]
        Z_pert= Z_pert.append(tmp)
        # plt.plot(range(iters),np.log10(all_loss))
        print('Finished: Treatment = {}, {} = {}'.format(treatment,phenotype,centroids.index[j]))
    # plt.savefig('FindPerturbation4Target'+phenotype+'_'+in_vitro_name+'_for_'+in_vivo_name+'_'+treatment+'.png')
res_all = res_all.reset_index(drop=True)
res_all.to_csv('../results/FindPerturbation4Target'+phenotype+'_'+in_vitro_name+'_for_'+in_vivo_name+'.csv')
Z_pert.to_csv('../results/Z_pert'+phenotype+'_'+in_vitro_name+'_for_'+in_vivo_name+'.csv')