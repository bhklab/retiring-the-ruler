from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from damply import dirs
from pathlib import Path

def save_plot(figure,
              filepath:Path):
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)

    figure.savefig(filepath)

    return

def plot_recist_accuracy(recist_accuracy,
                         save_path: Path | None = None
                         ) -> Figure:
    fig = plt.figure(figsize=(8, 6))
    plt.axvline(5, color='red', linestyle='--', label='RECIST v1.1')
    plt.scatter(range(1, 11), 1 - np.array(recist_accuracy) / 100, marker='o', s=75, label = 'Observation')
    plt.legend(bbox_to_anchor=(0.90, 1), loc='upper left')
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
    plt.legend(bbox_to_anchor=(0.9, 1), loc='upper left')
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
    fig, ax1 = plt.subplots(figsize=(8, 4))

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
    fig.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left', ncol=3)
    sns.despine()

    if save_path:
        save_file = save_path / "recist_accuracy_and_sensitivity.png"
        save_plot(fig, save_file)

    return fig