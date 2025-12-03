import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from matplotlib.figure import Figure
from pathlib import Path

def save_plot(figure,
              filepath:Path):
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    figure.savefig(filepath,
                   bbox_inches = 'tight',
                   )

    return


def plot_recist_accuracy(recist_accuracy,
                         save_path: Path | None = None
                         ) -> Figure:
    fig = plt.figure(figsize=(8, 6))
    plt.axvline(5, color='red', linestyle='--', label='RECIST v1.1')
    plt.scatter(range(1, 11), 1 - np.array(recist_accuracy) / 100, marker='o', s=75, label = 'Observation')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=2)
    sns.despine()
    plt.xlabel('Target Lesions')
    plt.ylabel('Misclassified Patients')

    if save_path:
        save_file = save_path / "recist_accuracy.png"
        save_plot(fig, save_file)

    return fig


def plot_pd_sensitivity(pd_sensitivity,
                        save_path: Path | None = None
                        ) -> Figure:
    fig = plt.figure(figsize=(8, 6))
    plt.axvline(5, color='red', linestyle='--', label='RECIST v1.1')
    plt.scatter(range(1, 11), np.array(pd_sensitivity) / 100, marker='o', s=75, label = 'Observation')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=2)
    sns.despine()
    plt.xlabel('Target Lesions')
    plt.ylabel('Sensitivity')

    if save_path:
        save_file = save_path / "PD_sensitivity.png"
        save_plot(fig, save_file)

    return fig


def plot_acc_and_sens(accuracy,
                      sensitivity,
                      save_path: Path | None = None
                      ) -> Figure:
    # Plot sensitivity and misclassification rate on the same plot
    fig, ax1 = plt.subplots(figsize=(8, 6))

    color = 'tab:blue'
    ax1.set_xlabel('Target Lesions')
    ax1.set_ylabel('Misclassified Patients', color=color)
    ax1.scatter(range(1, 11), 1 - np.array(accuracy) / 100, marker='o', s=75, label='Misclassified Patients', color=color)
    ax1.tick_params(axis='y', labelcolor=color)


    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:brown'
    ax2.set_ylabel('Sensitivity', color=color)  # we already handled the x-label with ax1
    ax2.scatter(range(1, 11), np.array(sensitivity) / 100, marker='o', s=75, label='Sensitivity', color=color)
    ax2.tick_params(axis='y', labelcolor=color, length=0)
    ax1.axvline(5, color='red', linestyle='--', label='RECIST v1.1')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    # Combine handles and labels from both axes
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles1 + handles2
    labels = labels1 + labels2

    # Place legend outside, 1 row x 3 columns
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=3)
    sns.despine()

    if save_path:
        save_file = save_path / "recist_accuracy_and_sensitivity.png"
        save_plot(fig, save_file)

    return fig


def plot_vol_vs_diameter(lesion_data: pd.DataFrame,
                        save_path: Path | None = None):

    diameters = lesion_data['diameter_pre'].values / 10
    volumes = lesion_data['volume_cc_contoured'].values

    expected_volume = 4/3 * np.pi * (np.linspace(0,25,100)/2)**3

    # remove any points where the volume is >100 if the diameter is < 2
    idx_to_keep = ~np.logical_and(diameters <= 3, volumes >= 100)
    diameters = diameters[idx_to_keep]
    volumes = volumes[idx_to_keep]

    volume_variation = volumes - 4/3 * np.pi * (diameters/2)**3
    
    # Scatter plot of volume versus diameter
    vol_v_diam_fig = plt.figure(figsize=(8, 7))
    plt.scatter(diameters, volumes, alpha=0.5, label='Observed volume')
    plt.plot(np.linspace(0, 25, 100), expected_volume, color='red', label='Expected volume')
    plt.ylim(0, 2000)
    plt.xlim(0, 20)
    plt.xlabel('Diameter (cm)')
    plt.ylabel(r'Volume ($cm^3$)')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=2)
    sns.despine(trim=True, offset=5)

    # Residuals plot of volume variation
    vol_var_v_diam_fig = plt.figure(figsize=(8, 7))
    plt.scatter(diameters, volume_variation, alpha=0.5, label='Volume Variation')
    plt.axhline(0, color='red', linestyle='--', label='Expected Volume')
    plt.ylim(-1000, 1000)
    plt.xlim(0, 20)
    plt.xlabel('Diameter (cm)')
    plt.ylabel('Volume variation ($cm^3$)')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.09), ncol=2)
    sns.despine()

    if save_path:
        vol_v_diam_file = save_path / "expected_vs_observed_volume.png"
        vol_var_v_diam_file = save_path / "expected_vs_variation_volume.png"
        
        save_plot(vol_v_diam_fig, vol_v_diam_file)
        save_plot(vol_var_v_diam_fig, vol_var_v_diam_file)

    return vol_v_diam_fig, vol_var_v_diam_fig
    