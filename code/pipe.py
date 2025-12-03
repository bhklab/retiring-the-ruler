
import pandas as pd
import numpy as np
from damply import dirs
from pathlib import Path
from synthetic_gen import generate_synthetic_patients
from recist import recist_assess, select_target_lesions, recist_metrics_by_target_count
from plot import (
    plot_recist_accuracy, 
    plot_pd_sensitivity,
    plot_acc_and_sens,
    plot_vol_vs_diameter
)


def pipe(radiomic_features_filepath: str,
         num_sim_patients: int = 10000,
         expected_num_lesions: int = 10,
         location_label: str = "LABEL",
         save_out: bool = False,
         random_seed: int | None = None
         ):
    dataset_name = Path(radiomic_features_filepath).parent.stem
    # Load radiomics data
    rad_data = pd.read_csv(radiomic_features_filepath)

    # Extract the slice thickness from the radiomics data
    rad_data['slice_thickness'] = [float(i.split(', ')[-1].strip('()')) for i in rad_data['diagnostics_Image-interpolated_Spacing']]

    # Calculate the true volume of the segmentation
    rad_data['volume_cc_contoured'] = rad_data['original_shape_VoxelVolume'] * rad_data['slice_thickness'] / 1000

    lesion_selection_rng = np.random.default_rng(random_seed)

    # Generate lesion measurements for N synthetic patients with 1 to M synthetic lesions based on real lesion data
    synth_lesions = generate_synthetic_patients(num_sim_patients=num_sim_patients,
                                                base_radiomic_data=rad_data,
                                                expected_num_lesions=expected_num_lesions,
                                                location_label=location_label,
                                                lesion_selection_rng=lesion_selection_rng)

    # Assess the RECIST response category of each synthetic patient using all lesions
    synth_response = recist_assess(synth_lesions)

    # Reassess RECIST iteratively using 1-10 target lesions per patient (max 2 per location)
    for num_targets in range(1, 11):

        target_lesions = select_target_lesions(num_lesions=num_targets,
                                               lesion_data=synth_lesions,
                                               lesion_selection_rng=lesion_selection_rng
                                               )
        # Reassess RECIST categorization with subset of lesions
        select_target_response = recist_assess(target_lesions)
        # Add RECIST category for this number of target lesions
        synth_response[f"RECIST ({num_targets} targets)"] = select_target_response['RECIST (all)']

    if save_out:
        out_path = dirs.PROCDATA / dataset_name / f"sim_{num_sim_patients}_pats"
        out_path.mkdir(parents=True, exist_ok=True)

        synth_lesions.to_csv(out_path / f"{dataset_name}_synthetic_lesions.csv", index_label="index")
        synth_response.to_csv(out_path / f"{dataset_name}_synthetic_patient_response.csv", index_label="index")

        plot_path = dirs.RESULTS / dataset_name / f"sim_{num_sim_patients}_pats"
    else:
        plot_path = None

    # Calculate the classification accuracy for RECIST as a function of the number of target lesions
    recist_accuracy, pd_sensitivity = recist_metrics_by_target_count(patient_response=synth_response,
                                                                     max_targets=11)
    acc_plot = plot_recist_accuracy(recist_accuracy, plot_path)
    pd_sense_plot = plot_pd_sensitivity(pd_sensitivity, plot_path)
    acc_sense_plot = plot_acc_and_sens(recist_accuracy, pd_sensitivity, plot_path)
    vol_v_diam_plot, vol_var_v_diam_plot = plot_vol_vs_diameter(synth_lesions, plot_path)
    
    return


if __name__ == '__main__':
    pipe("data/rawdata/SARC021/SARC021_radiomics.csv",
         num_sim_patients = 100,
         expected_num_lesions = 10,
         save_out=True,
         random_seed=165)