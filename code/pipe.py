
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from damply import dirs
from plot import (
    plot_acc_and_sens,
    plot_pd_sensitivity,
    plot_recist_accuracy,
    plot_vol_vs_diameter,
)
from recist import recist_assess, recist_metrics_by_target_count, select_target_lesions
from synthetic_gen import generate_synthetic_patients

logger = logging.getLogger(__name__)

def pipe(radiomic_features_filepath: str,
         num_sim_patients: int = 10000,
         expected_num_lesions: int = 10,
         location_label: str = "LABEL",
         save_out: bool = False,
         random_seed: int | None = None
         ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retiring the Ruler Simulation pipeline
    
       Generate synthetic lesion data with pre- and post-treatment diameters and volumes to perform RECIST assessment.

    Parameters
    ----------
    radiomic_features_filepath: str
        Real radiomic feature data with diameter and volume data to generate simulation from.
    num_sim_patients: int = 10000
        Number of patients to generate simulations for
    expected_num_lesions: int = 10
        Expected number of lesions per patient. Used to set up a Poisson distribution to select from.
    location_label: str = "LABEL"
        Label in the radiomic feature data identifying where the ground truth tumor is located.
    save_out: bool = False
        Whether to save out the simulated lesion data and plots for analysis.
    random_seed: int | None = None
        Random seed to use for reproducible results. Used to initialize the random number generators.
    
    Returns
    -------
    synth_lesions: pd.DataFrame
        Measurment data for each synthetic lesion, including diameters, volumes, and location.
        If save_out is set to True, will be saved out to the data/procdata directory. 
    synth_response: pd.DataFrame
        Response data for each synthetic patient, including sum of longest diameters (SLD) and RECIST response classification with different numbers of target lesions.
        If save_out is set to True, will be saved out to the data/procdata directory.
    """
    dataset_name = Path(radiomic_features_filepath).parent.stem

    logging.basicConfig(filename=f'logs/pipe_{dataset_name}.log', level=logging.INFO, format='%(asctime)s %(message)s')
    logger.info('\n')
    logger.info(f'Retiring the Ruler pipeline started for {dataset_name}')

    # Load radiomics data
    rad_data = pd.read_csv(radiomic_features_filepath)

    # Extract the slice thickness from the radiomics data
    rad_data['slice_thickness'] = [float(i.split(', ')[-1].strip('()')) for i in rad_data['diagnostics_Image-interpolated_Spacing']]

    # Calculate the true volume of the segmentation
    rad_data['volume_cc_contoured'] = rad_data['original_shape_VoxelVolume'] * rad_data['slice_thickness'] / 1000

    lesion_selection_rng = np.random.default_rng(random_seed)

    logger.info('Starting synthetic lesion generation.')
    # Generate lesion measurements for N synthetic patients with 1 to M synthetic lesions based on real lesion data
    synth_lesions = generate_synthetic_patients(num_sim_patients=num_sim_patients,
                                                base_radiomic_data=rad_data,
                                                expected_num_lesions=expected_num_lesions,
                                                location_label=location_label,
                                                lesion_selection_rng=lesion_selection_rng)
    logger.info('Synthetic lesion generation finished.')
    logger.info('Performing RECIST assessment of synthetic lesions.')
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

    logger.info('Finished RECIST assessments.')
    logger.info('Calculating accuracy and PD sensitivity.')
    # Calculate the classification accuracy for RECIST as a function of the number of target lesions
    recist_accuracy, pd_sensitivity = recist_metrics_by_target_count(patient_response=synth_response,
                                                                     max_targets=11)
    

    if save_out:
        logger.info('Saving out the synthetic lesion and response data.')
        out_path = dirs.PROCDATA / dataset_name / f"sim_{num_sim_patients}_pats"
        out_path.mkdir(parents=True, exist_ok=True)

        synth_lesions.to_csv(out_path / f"{dataset_name}_synthetic_lesions.csv", index_label="index")
        synth_response.to_csv(out_path / f"{dataset_name}_synthetic_patient_response.csv", index_label="index")

        logger.info('Plotting analysis results.')
        plot_path = dirs.RESULTS / dataset_name / f"sim_{num_sim_patients}_pats"
        plot_recist_accuracy(recist_accuracy, plot_path)
        plot_pd_sensitivity(pd_sensitivity, plot_path)
        plot_acc_and_sens(recist_accuracy, pd_sensitivity, plot_path)
        plot_vol_vs_diameter(synth_lesions, plot_path)
    
    return synth_lesions, synth_response


if __name__ == '__main__':
    pipe("data/rawdata/SARC021/SARC021_radiomics.csv",
         num_sim_patients = 100,
         expected_num_lesions = 10,
         save_out=True,
         random_seed=165)