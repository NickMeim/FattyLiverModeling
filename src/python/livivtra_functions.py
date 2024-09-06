import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

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

def pair_pearsonr(x, y, axis=0):
    mx = np.mean(x, axis=axis, keepdims=True)
    my = np.mean(y, axis=axis, keepdims=True)
    xm, ym = x-mx, y-my
    r_num = np.add.reduce(xm * ym, axis=axis)
    r_den = np.sqrt((xm*xm).sum(axis=axis) * (ym*ym).sum(axis=axis))
    r = r_num / r_den
    return r

def visualize_paths_activity(scores,extra_basis,lim=8.0, cmap='vlag', figsize=(12, 8)):
    i = extra_basis
    # Filter data based on condition
    filtered_data = scores[scores['condition'] == f'V{i}']
    # Sort the data based on 'activity' from highest to lowest
    filtered_data = filtered_data.sort_values(by='activity', ascending=False)

    # Create the plot
    plt.figure(figsize=figsize)

    # Barplot
    barplot = sns.barplot(
        data=filtered_data,
        x='activity',
        y='Pathway',
        palette=cmap
    )
    # Titles and labels
    plt.title(f'LV extra {i}', fontsize=24, family='Arial')
    plt.xlabel('Activity', fontsize=24, family='Arial')
    plt.ylabel('Pathway', fontsize=24, family='Arial')
    # Customize the font size and legend position
    plt.xticks(fontsize=18, family='Arial')
    plt.yticks(fontsize=18, family='Arial')
    # Get the current figure
    fig = plt.gcf()
    return fig