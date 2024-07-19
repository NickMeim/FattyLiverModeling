import torch
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.svm import SVR, LinearSVR
from sklearn.linear_model import Lasso,Ridge,ElasticNet
# from sklearn.gaussian_process import GaussianProcessRegressor
# from sklearn.gaussian_process.kernels import DotProduct, RBF, RationalQuadratic
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor as KNN
from sklearn.multioutput import MultiOutputRegressor
from sklearn.model_selection import GridSearchCV
import joblib
import xgboost as XGB
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
parser = argparse.ArgumentParser(description='Run different ML models')
parser.add_argument('model_types', metavar='N', type=str, nargs='*', help='models to train',default=['knn' ,'svmLinear','svmRBF','svmPoly','lasso','ridge','elasticNet','neuralNet','xgboost','rf'])
parser.add_argument('--cv_files_location', action='store',help='location of files used in CV training of originalPLSR model',default='../preprocessing/TrainingValidationData/WholePipeline/crossfoldPLSR/')
parser.add_argument('--clinical_files_location', action='store',help='location of external clinical files',default='../preprocessing/TrainingValidationData/external_clinical_data/')
parser.add_argument('--clinical_datasets', metavar='N', type=str, nargs='*', help='names of the external clinical datasets',default=['Hoang','Pantano'])
parser.add_argument('--in_vitro_matrices_loc', action='store',help='location of gene weights/loadings/extra basis matrices of the in-vitro data',default='../results/')
parser.add_argument('--in_vitro_name', action='store',help='name of the in-vitro data',default='Kostrzewski')
parser.add_argument('--num_folds', action='store', type=int,help='number of folds',default=10)
args = parser.parse_args()
model_types = args.model_types
cv_files_location= args.cv_files_location
clinical_files_location = args.clinical_files_location
clinical_datasets = args.clinical_datasets
in_vitro_matrices_loc = args.in_vitro_matrices_loc
in_vitro_name = args.in_vitro_name
num_folds = args.num_folds

### Load the data
Wm_total = pyreadr.read_r(in_vitro_matrices_loc+'Wm_'+in_vitro_name+'_total.rds')
Wm_total = Wm_total[None]
Wm = Wm_total.iloc[:,:-2]
## conver to numpy
Wm = Wm.to_numpy()
Wm_total = Wm_total.to_numpy()

### Initialize all the models I am going to use
models = []
for mdl in model_types:
    if mdl == 'knn':
        params = {
            'n_neighbors': [2,3,5,10,15,20],
            }
        model = GridSearchCV(estimator=KNN(),param_grid = params, cv=5, n_jobs=-1)
        # model = KNN(n_neighbors=5)
    elif mdl == 'rf':
        # params = {
        #     'n_estimators': [10,25,50,75,100],
        #     'max_depth':[10, 25, 50, 75,100, None],
        #     'min_samples_leaf': [1, 2, 4],
        #     'min_samples_split': [2, 5, 10]
        #     }
        # model = GridSearchCV(estimator=RandomForestRegressor(criterion='absolute_error'),param_grid = params, cv=5, n_jobs=-1)
        model = RandomForestRegressor(n_estimators=100,n_jobs=-1)
    elif mdl == 'xgboost':
        model = XGB.XGBRegressor(n_estimators=100,n_jobs = -1)
    elif mdl == 'svmLinear':
        params = {
            'epsilon': [0,1e-4,1e-3,1e-2,1e-1,1],
            'C': [1e-1,1, 10, 100, 1000]
            }
        model = GridSearchCV(estimator=LinearSVR(),param_grid = params, cv=5, n_jobs=-1)
        model = MultiOutputRegressor(model)
    elif mdl == 'svmRBF':
        params = {
            'gamma': [1e-5,1e-4,1e-3,1e-2,1e-1,1,1.5,2],
            'C': [1e-1,1, 10, 100, 1000]
            }
        model = GridSearchCV(estimator=SVR(kernel='rbf'),param_grid = params, cv=5, n_jobs=-1)
        model = MultiOutputRegressor(model)
    elif mdl == 'svmPoly':
        params = {
            'gamma': [1e-5,1e-4,1e-3,1e-2,1e-1,1,1.5,2],
            'C': [1e-1,1, 10, 100, 1000],
            'degree':[2,3,4,5],
            'coef0':[0,0.1,0.5,1.,1.2,2.]
            }
        model = GridSearchCV(estimator=SVR(kernel='poly'),param_grid = params, cv=5, n_jobs=-1)
        model = MultiOutputRegressor(model)
    elif mdl == 'lasso':
        params = {
            'alpha': [0.001,0.01,0.1,1,10,100,1000]
            }
        model = GridSearchCV(estimator=Lasso(),param_grid = params, cv=5, n_jobs=-1)
    elif mdl == 'ridge':
        params = {
            'alpha': [0.001,0.01,0.1,1,10,100,1000],
            }
        model = GridSearchCV(estimator=Ridge(),param_grid = params, cv=5, n_jobs=-1)
    elif mdl == 'elasticNet':
        params = {
            'alpha': [0.001,0.01,0.1,1,10,100,1000],
            'l1_ratio': [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
            }
        model = GridSearchCV(estimator=ElasticNet(),param_grid = params, cv=5, n_jobs=-1)
    elif mdl == 'neuralNet':
        device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        model = 'define the ANN just before it is trained'
        epochs = 100
        l2_reg  = 0.01
        bs = 15
        criterion = torch.nn.MSELoss(reduction='mean')
    models.append(model)


### Initialize lists and matrices to save the results
k=0
for model in models:
    print('Begun fitting and evaluation for model: %s'%model_types[k])
    folder = '../results/MLresults/'+model_types[k]+'/'
    Path(folder).mkdir(parents=True, exist_ok=True)
    val_r_extrl = np.zeros((num_folds,len(clinical_datasets)))
    val_r_extrl_backproj = np.zeros((num_folds,len(clinical_datasets)))
    val_r_extrl_extra = np.zeros((num_folds,len(clinical_datasets)))
    val_r = []
    train_r = []
    val_r_backproj = []
    train_r_backproj = []
    val_r_extra = []
    train_r_extra = [] 
    for i in range(num_folds):
        x_train = pyreadr.read_r(cv_files_location+'Xh_train'+str(i+1)+'.rds')
        x_train = x_train[None]
        y_train = pyreadr.read_r(cv_files_location+'Yh_train'+str(i+1)+'.rds')
        y_train = y_train[None]
        x_val= pyreadr.read_r(cv_files_location+'Xh_val'+str(i+1)+'.rds')
        x_val = x_val[None]
        y_val = pyreadr.read_r(cv_files_location+'Yh_val'+str(i+1)+'.rds')
        y_val = y_val[None]
        if (model_types[k] != 'neuralNet'):
            # fit model and evaluate in validation set
            model.fit(x_train,y_train)
            yhat_train = model.predict(x_train)
            yhat_val = model.predict(x_val)
            joblib.dump(model, folder + 'model' + str(i) + '.pkl')
            train_r.append(pair_pearsonr(y_train.to_numpy(), yhat_train, axis=0).mean())
            val_r.append(pair_pearsonr(y_val.to_numpy(), yhat_val, axis=0).mean())

            # prediction after backprojection
            #first for training data
            Xback_train = x_train.to_numpy() @  Wm @ Wm.T
            Xback_train = pd.DataFrame(Xback_train, index=x_train.index, columns=x_train.columns)
            yhat_train_backproj = model.predict(Xback_train)
            train_r_backproj.append(pair_pearsonr(y_train.to_numpy(), yhat_train_backproj, axis=0).mean())
            # repeat for valdiation
            Xback_val = x_val.to_numpy() @  Wm @ Wm.T
            Xback_val = pd.DataFrame(Xback_val, index=x_val.index, columns=x_val.columns)
            yhat_val_backproj = model.predict(Xback_val)
            val_r_backproj.append(pair_pearsonr(y_val.to_numpy(), yhat_val_backproj, axis=0).mean())

            # prediction after backprojection with Wm_total
            #first for training data
            Xextra_train = x_train.to_numpy() @  Wm_total @ Wm_total.T
            Xextra_train = pd.DataFrame(Xextra_train, index=x_train.index, columns=x_train.columns)
            yhat_train_extra = model.predict(Xextra_train)
            train_r_extra.append(pair_pearsonr(y_train.to_numpy(), yhat_train_extra, axis=0).mean())
            # repeat for valdiation
            Xextra_val = x_val.to_numpy() @  Wm_total @ Wm_total.T
            Xextra_val = pd.DataFrame(Xextra_val, index=x_val.index, columns=x_val.columns)
            yhat_val_extra = model.predict(Xextra_val)
            val_r_extra.append(pair_pearsonr(y_val.to_numpy(), yhat_val_extra, axis=0).mean())

            # evaluate in external clinical datasets
            for j in range(len(clinical_datasets)):
                X = pd.read_csv(clinical_files_location+'/'+clinical_datasets[j]+'/'+'X.csv',index_col=0)
                Y = pd.read_csv(clinical_files_location+'/'+clinical_datasets[j]+'/'+'Y.csv',index_col=0)
                yhat_val_ext = model.predict(X)
                val_r_extrl[i][j] = pair_pearsonr(Y.to_numpy(), yhat_val_ext, axis=0).mean()

                # prediction after backprojection
                X2 = X.to_numpy() @  Wm @ Wm.T
                X2 = pd.DataFrame(X2, index=X.index, columns=X.columns)
                yhat_val_backproj = model.predict(X2)
                val_r_extrl_backproj[i][j] = pair_pearsonr(Y.to_numpy(), yhat_val_backproj, axis=0).mean()

                # prediction after backprojection with Wm_total
                X2 = X.to_numpy() @  Wm_total @ Wm_total.T
                X2 = pd.DataFrame(X2, index=X.index, columns=X.columns)
                yhat_val_extra = model.predict(X2)
                val_r_extrl_extra[i][j] = pair_pearsonr(Y.to_numpy(), yhat_val_extra, axis=0).mean()
        else:
            Wm_gpu = torch.tensor(Wm,dtype=torch.double).to(device)
            Wm_total_gpu = torch.tensor(Wm_total_gpu,dtype=torch.double).to(device)
            model = torch.nn.Sequential(torch.nn.Dropout(0.5),
                                        torch.nn.Linear(Wm.shape[0],1024,bias=True,dtype=torch.double),
                                        torch.nn.BatchNorm1d(num_features=1024, momentum=0.2,dtype = torch.double),
                                        torch.nn.ELU(),
                                        torch.nn.Dropout(0.2),
                                        torch.nn.Linear(1024, 256,bias=True,dtype=torch.double),
                                        torch.nn.BatchNorm1d(num_features=256, momentum=0.2,dtype = torch.double),
                                        torch.nn.ELU(),
                                        torch.nn.Dropout(0.2),
                                        torch.nn.Linear(256, 32,bias=True,dtype=torch.double),
                                        torch.nn.BatchNorm1d(num_features=32, momentum=0.2,dtype = torch.double),
                                        torch.nn.ELU(),
                                        torch.nn.Dropout(0.2),
                                        torch.nn.Linear(32, 2,bias=True,dtype=torch.double)).to(device)
            optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
            # loss_mu = []
            # loss_std = []
            # r_mu = []
            # r_std = []
            for e in range(epochs):
                trainloader = getSamples(x_train.shape[0], bs)
                all_losses=[]
                all_r = []
                for dataIndex in trainloader:
                    model.train()
                    optimizer.zero_grad()
                    dataIn = torch.tensor(x_train.iloc[dataIndex, :].values).to(device)
                    dataOut = torch.tensor(y_train.iloc[dataIndex, :].values,dtype=torch.double).to(device)
                    Yhat = model(dataIn)
                    fitLoss = criterion(dataOut, Yhat)
                    RegLoss = L2Regularization(model,L2=l2_reg)
                    loss = fitLoss + RegLoss
                    loss.backward()
                    optimizer.step()
                    all_losses.append(loss.item())
                    r = torch.mean(pearson_r(dataOut, Yhat))
                    all_r.append(r.item())
                if(e%10==0 or e==0 or e==epochs-1):
                    print('Fold {}, Epoch {}/{} : Loss = {}, r = {}'.format(i,e+1,epochs,np.mean(all_losses),np.mean(all_r)))
                # loss_mu.append(np.mean(all_losses))
                # loss_std.append(np.std(all_losses))
                # r_mu.append(np.mean(all_r))
                # r_std.append(np.std(all_r))

            model.eval()
            yhat_train = model(torch.tensor(x_train.values).to(device))
            yhat_val = model(torch.tensor(x_val.values).to(device))
            torch.save(model, folder+'model'+str(i)+'.pt')
            train_r.append(pearson_r(torch.tensor(y_train.values,dtype=torch.double).to(device), yhat_train).mean().item())
            val_r.append(pearson_r(torch.tensor(y_val.values,dtype=torch.double).to(device), yhat_val).mean().item())

            # prediction after backprojection
            #first for training data
            Xback_train = x_train.to_numpy() @  Wm @ Wm.T
            yhat_train_backproj = model(torch.tensor(Xback_train).to(device))
            train_r_backproj.append(pearson_r(torch.tensor(y_train.to_numpy(),dtype=torch.double).to(device), yhat_train_backproj).mean().item())
            # repeat for valdiation
            Xback_val = x_val.to_numpy() @  Wm @ Wm.T
            yhat_val_backproj = model(torch.tensor(Xback_val).to(device))
            val_r_backproj.append(pearson_r(torch.tensor(y_val.to_numpy(),dtype=torch.double).to(device), yhat_val_backproj).mean().item())

            # prediction after backprojection with Wm_total
            #first for training data
            Xextra_train = x_train.to_numpy() @  Wm_total @ Wm_total.T
            yhat_train_extra = model(torch.tensor(Xextra_train).to(device))
            train_r_extra.append(pearson_r(torch.tensor(y_train.to_numpy(),dtype=torch.double).to(device), yhat_train_extra).mean().item())
            # repeat for valdiation
            Xextra_val = x_val.to_numpy() @  Wm_total @ Wm_total.T
            yhat_val_extra = model(torch.tensor(Xextra_val).to(device))
            val_r_extra.append(pearson_r(torch.tensor(y_val.to_numpy(),dtype=torch.double).to(device), yhat_val_extra).mean().item())

            # evaluate in external clinical datasets
            for j in range(len(clinical_datasets)):
                X = torch.tensor(pd.read_csv(clinical_files_location+'/'+clinical_datasets[j]+'/'+'X.csv',index_col=0).values,dtype=torch.double).to(device)
                Y = torch.tensor(pd.read_csv(clinical_files_location+'/'+clinical_datasets[j]+'/'+'Y.csv',index_col=0).values,dtype=torch.double).to(device)
                yhat_val_ext = model(X)
                val_r_extrl[i][j] = pearson_r(Y, yhat_val_ext).mean().item()

                # prediction after backprojection
                X2 = X @  Wm_gpu @ Wm_gpu.T
                yhat_val_backproj = model.predict(X2)
                val_r_extrl_backproj[i][j] = pearson_r(Y, yhat_val_backproj).mean().item()
                
                # prediction after backprojection with Wm_total
                X2 = X @  Wm_total_gpu @ Wm_total_gpu.T
                yhat_val_extra = model(X2)
                val_r_extrl_extra[i][j] = pearson_r(Y, yhat_val_extra).mean().item()
        print('Finished fold %s'%i)
    res_val = pd.DataFrame({'human genes':val_r,'back-projected':val_r_backproj,'optimized MPS':val_r_extra})
    res_val['fold'] = [xx for xx in range(num_folds)]
    res_val['model'] = model_types[k]
    res_val['set'] = 'test'
    res_train = pd.DataFrame({'human genes':train_r,'back-projected':train_r_backproj,'optimized MPS':train_r_extra})
    res_train['fold'] = [xx for xx in range(num_folds)]
    res_train['model'] = model_types[k]
    res_train['set'] = 'train'
    df_res_cv = pd.concat([res_val,res_train])
    df_res_cv.to_csv(folder+'df_res_cv.csv')

    val_r_extrl = pd.DataFrame(val_r_extrl)
    val_r_extrl = val_r_extrl.reset_index()
    val_r_extrl.columns =  ['fold']+clinical_datasets
    val_r_extrl['model'] = model_types[k]
    val_r_extrl = pd.melt(val_r_extrl, id_vars=['model','fold'], var_name='dataset', value_name='r')
    val_r_extrl['set'] = 'external dataset'
    val_r_extrl['input'] = 'human genes'

    val_r_extrl_backproj = pd.DataFrame(val_r_extrl_backproj)
    val_r_extrl_backproj = val_r_extrl_backproj.reset_index()
    val_r_extrl_backproj.columns = ['fold']+clinical_datasets
    val_r_extrl_backproj['model'] = model_types[k]
    val_r_extrl_backproj = pd.melt(val_r_extrl_backproj, id_vars=['model','fold'], var_name='dataset', value_name='r')
    val_r_extrl_backproj['set'] = 'external dataset'
    val_r_extrl_backproj['input'] = 'back-projected'

    val_r_extrl_extra = pd.DataFrame(val_r_extrl_extra)
    val_r_extrl_extra = val_r_extrl_extra.reset_index()
    val_r_extrl_extra.columns =  ['fold']+clinical_datasets
    val_r_extrl_extra['model'] = model_types[k]
    val_r_extrl_extra = pd.melt(val_r_extrl_extra, id_vars=['model','fold'], var_name='dataset', value_name='r')
    val_r_extrl_extra['set'] = 'external dataset'
    val_r_extrl_extra['input'] = 'optimized MPS'
    
    df_res_external  = pd.concat([val_r_extrl,val_r_extrl_backproj,val_r_extrl_extra])
    df_res_external.to_csv(folder+'df_res_external.csv')

    k+=1
