import numpy as np
import pandas as pd

from scipy.stats import percentileofscore

def calculate_p_values(df1, df2):
    p_values = pd.DataFrame(index=df1.index, columns=df1.columns)
    
    for feature in df1.columns:
        for sample in df1.index:
            # Absolute value of the feature in the current sample
            value = abs(df1.loc[sample, feature])
            # Distribution of the feature in the second matrix
            distribution = df2[feature].abs()
            # Calculate the percentile of the current value within the distribution
            percentile = percentileofscore(distribution, value, kind='rank')
            # p-value is the probability of finding a higher or equal absolute value
            p_value = 1 - percentile / 100
            p_values.loc[sample, feature] = p_value
    
    return p_values

def p_value_label(p_value):
    if p_value <= 0.0001:
        return "****"
    elif p_value <= 0.001:
        return "***"
    elif p_value <= 0.01:
        return "**"
    elif p_value <= 0.05:
        return "*"
    elif p_value <= 0.1:
        return "\u2219"
    else:
        return "ns"

def extra_basis_analytical_solution(y, W_invitro, phi):
    Wopt = np.zeros((W_invitro.shape[0], y.shape[1]))  # Initialize Wopt with zeros
    for i in range(y.shape[1]):
        if i == 0:
            alpha = W_invitro.T @ phi[:, i]
            Wopt[:, i] = (phi[:, i] - W_invitro @ alpha) / np.sqrt(np.sum(phi[:, i]**2) - np.sum(alpha**2))
        else:
            Wnew = np.hstack((W_invitro, Wopt[:, :i]))  # Combine W_invitro and Wopt[:, :i]
            alpha = Wnew.T @ phi[:, i]
            Wopt[:, i] = (phi[:, i] - Wnew @ alpha) / np.sqrt(np.sum(phi[:, i]**2) - np.sum(alpha**2))
    
    return Wopt

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

