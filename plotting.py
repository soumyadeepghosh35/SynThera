# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 13:23:53 2021

Contains plotting functions for MACAW project.

@author: Vincent
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error, root_mean_squared_error
from sklearn.metrics import auc

from matplotlib import rc
matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.sans-serif'] = ['DejaVu Sans']
matplotlib.rcParams['font.size']        = 10
matplotlib.rcParams['figure.dpi'] = 96
matplotlib.rcParams['savefig.dpi'] = 300
matplotlib.rc('figure', figsize=(5.0, 4.0))
#matplotlib.rcParams['font.family']      = 'serif'
#matplotlib.rcParams['font.serif']       = ['Times New Roman']
#matplotlib.rcParams['mathtext.fontset'] = 'stix'




# ----- Plotting functions -----


def parity_plot(
    x,
    y,
    x_test=None,
    y_test=None,
    y_train_std=None,
    y_test_std=None,
    xlabel='True value',
    ylabel='Predicted',
    title=None,
    savetitle=None,
    save_formats=['svg', 'png'],  # NEW: Accept multiple formats
):
    x = np.array(x)
    y = np.array(y)

    # Plot Figures
    plt.figure(figsize=(4.4, 4.0))
    
    # Find the boundaries of X and Y values
    if x_test is not None:
        x_test = np.array(x_test)
        y_test = np.array(y_test)
        min1 = min(x.min(), y.min(), x_test.min(), y_test.min())
        max1 = max(x.max(), y.max(), x_test.max(), y_test.max())
    else:
        min1 = min(x.min(), y.min())
        max1 = max(x.max(), y.max())
    rng1 = max1 - min1
    bounds = (min1 - 0.05 * rng1, max1 + 0.05 * rng1)

    # Reset the limits
    ax = plt.gca()
    ax.set_xlim(bounds)
    ax.set_ylim(bounds)

    # plot the diagonal
    ax.plot([0, 1], [0, 1], 'k-', lw=0.5, transform=ax.transAxes)

    # Plot the data
    plt.errorbar(x=x, y=y, yerr=y_train_std, color='blue', fmt='.', elinewidth=.5, alpha=0.5)

    if x_test is not None:
        plt.errorbar(x=x_test, y=y_test, yerr=y_test_std, color='red', fmt='.', elinewidth=.5, alpha=0.5)
        
        # Calculate train metrics
        train_r2 = r2_score(x, y)
        train_mae = mean_absolute_error(x, y)
        train_rmse = np.sqrt(mean_squared_error(x, y))
        
        # Calculate test metrics
        test_r2 = r2_score(x_test, y_test)
        test_mae = mean_absolute_error(x_test, y_test)
        test_rmse = np.sqrt(mean_squared_error(x_test, y_test))
        
        # Train metrics (top-left, blue)
        ax.text(0.02, 0.98, 
                f'$R^2_{{\\mathrm{{Train}}}} = {train_r2:.2f}$\n'
                f'$\\mathrm{{MAE}}_{{\\mathrm{{Train}}}} = {train_mae:.2f}$\n'
                f'$\\mathrm{{RMSE}}_{{\\mathrm{{Train}}}} = {train_rmse:.2f}$',
                transform=ax.transAxes, fontsize=10, va='top', ha='left',
                color='blue')
        
        # Test metrics (bottom-right, red)
        ax.text(0.98, 0.02, 
                f'$R^2_{{\\mathrm{{Test}}}} = {test_r2:.2f}$\n'
                f'$\\mathrm{{MAE}}_{{\\mathrm{{Test}}}} = {test_mae:.2f}$\n'
                f'$\\mathrm{{RMSE}}_{{\\mathrm{{Test}}}} = {test_rmse:.2f}$',
                transform=ax.transAxes, fontsize=10, va='bottom', ha='right',
                color='red')
    else:
        r2_text = f"$R^2 = {r2_score(x, y):0.2f}$"
        rmse_text = f"RMSE = {np.sqrt(mean_squared_error(x, y)):0.2f}"
        mae_text = f"MAE = {mean_absolute_error(x, y):0.2f}"
        plt.gca().text(0.05, 0.93, r2_text, transform=plt.gca().transAxes, fontsize=10.)
        plt.gca().text(0.05, 0.87, mae_text, transform=plt.gca().transAxes, fontsize=10.)
        plt.gca().text(0.05, 0.81, rmse_text, transform=plt.gca().transAxes, fontsize=10.)

    # Title and labels
    if title:
        plt.title(title)
    plt.xlabel(xlabel, fontsize=10.)
    plt.ylabel(ylabel, fontsize=10.)

    # Save the figure in multiple formats
    if savetitle:
        # NEW: Save in multiple formats
        import os
        base_path = os.path.splitext(savetitle)[0]  # Remove extension
        
        for fmt in save_formats:
            output_path = f"{base_path}.{fmt}"
            plt.savefig(output_path, format=fmt, dpi=300, bbox_inches='tight', transparent=False)
            print(f"Saved: {output_path}")
    
    plt.show()
    
    plt.close()  # Close figure to free memory


def plot_precision_vs_recall(precisions, recalls, precisions_test=None, recalls_test=None, title=None, savetitle=None):
    
    plt.figure(figsize=(4.4, 4.0))
    plt.plot(recalls, precisions, 'b-', linewidth=1.)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.axis([0, 1, 0, 1])

    if precisions_test is None:
        text1 = f"$AUPRC = {auc(recalls, precisions):0.3f}$"
    else:
        text1 = f"$AUPRC_{{train}} = {auc(recalls, precisions):0.3f}$"
        
        plt.plot(recalls_test, precisions_test, 'r-', linewidth=2)
        text2 = f"$AUPRC_{{test}} = {auc(recalls_test, precisions_test):0.3f}$"
        plt.gca().text(0.05, 0.13, text2, transform=plt.gca().transAxes, fontsize=10., c='red')
    
    plt.gca().text(0.05, 0.05, text1, transform=plt.gca().transAxes, fontsize=10., c='blue')
    
    if title:
        plt.title(title)
    
    if savetitle:
        plt.savefig(savetitle, format='svg', dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()


def plot_histogram(
    Y, xlabel="Property value", ylabel="No. of compounds", title='', savetitle=None
):
    plt.figure(figsize=(5.0, 3.5))
    n, bins, patches = plt.hist(x=Y, bins=20, alpha=1., rwidth=1., edgecolor='black', linewidth=.5)
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 5)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if title is not None:
        if len(title) == 0:
            title = f"{len(Y)} compounds"
        plt.title(title)
    if savetitle:
        plt.savefig(savetitle, format='svg', dpi=300, bbox_inches='tight', transparent=False)
    else:
        plt.show()
