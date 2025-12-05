import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tqdm import tqdm

# RECIST response category thresholds
PD_THRESHOLD = 20
SD_THRESHOLD = -30
PR_THRESHOLD = -100

def recist_thresholding(sld_chg):
    if sld_chg > PD_THRESHOLD:
        return 'PD'
    elif SD_THRESHOLD < sld_chg <= PD_THRESHOLD:
        return 'SD'
    elif PR_THRESHOLD < sld_chg <= SD_THRESHOLD:
        return 'PR'
    elif sld_chg == PR_THRESHOLD:
        return 'CR'


def recist_assess(lesion_data:pd.DataFrame,
                  parallel:bool = False
                  ) -> pd.DataFrame:
    """Calculate RECIST response category from lesion diameter data."""
    patients, lesion_counts = np.unique(lesion_data['patient_id'], return_counts=True)

    # Calculate the sum of longest diameters (sld) for pre, post
    sld_pre = [np.sum(lesion_data['diameter_pre'][lesion_data['patient_id'] == p]) for p in patients]
    sld_post = [np.sum(lesion_data['diameter_post'][lesion_data['patient_id'] == p]) for p in patients]

    # Calculte the percentage different in the sum of longest diameters (sld) pre and post
    sld_chg = [(sld_post[i] - sld_pre[i])/sld_pre[i] * 100 for i in range(len(sld_pre))]

    # Initialize RECIST response vector
    recist_response = []

    # Classify each patient's response based on the percent change in sld
    if parallel:
        recist_response = Parallel()(
            delayed(recist_thresholding)(
                sld_chg=sld_chg[pat_idx]
            )
            for pat_idx in tqdm(range(len(patients)), 
                                desc="Performing RECIST assessment...", 
                                total=len(patients))
        )
    else:
        recist_response = [
            recist_thresholding(
                sld_chg = sld_chg[pat_idx]
            ) 
            for pat_idx in tqdm(range(len(patients)), 
                                desc="Performing RECIST assessment...", 
                                total=len(patients))
        ]

    patient_response = pd.DataFrame({'patient_id': patients, 
                                'num_lesions': lesion_counts, 
                                'sld_pre': sld_pre, 
                                'sld_post': sld_post, 
                                'sld_chg': sld_chg, 
                                'RECIST (all)': recist_response})
    patient_response = patient_response.sort_values(by='patient_id', ignore_index=True)
    return patient_response



def select_target_lesions(num_lesions:int,
                          lesion_data:pd.DataFrame,
                          lesion_selection_rng: np.random.Generator,
                          parallel: bool = False) -> pd.DataFrame:
    """Randomly select specified number of target lesions from lesion data, ensuring no more than 2 are selected for each location following RECIST specifications.
       
       Returns the selected target lesion rows from lesion data.
       """
    patients = np.unique(lesion_data['patient_id'])

    def selection_process(patient):
        # ensure that there are not more than 2 target lesions per location
        # get the locations for the patient
        locations = np.unique(lesion_data['location'][lesion_data['patient_id'] == patient])
        selected_lesion_idx = []

        # for each location, select up to 2 lesions
        for location in locations:
            # get the indices of the lesions at the location
            targets_at_loc_idx = lesion_data[(lesion_data['patient_id'] == patient) & (lesion_data['location'] == location)].index
            
            # select up to 2 lesions
            if len(targets_at_loc_idx) > 2:
                # select 2 random indices
                targets_at_loc_idx = lesion_selection_rng.choice(list(targets_at_loc_idx), 2, replace=False)
                selected_lesion_idx.extend(targets_at_loc_idx)
            else:
                selected_lesion_idx.extend(targets_at_loc_idx)
        
        # If there are more selected lesions from this location than the required number of lesions, randomly select from these
        if len(selected_lesion_idx) > num_lesions:
            return lesion_selection_rng.choice(selected_lesion_idx, num_lesions)
        else:
            return selected_lesion_idx
    
    if parallel:
        target_lesions = Parallel()(
            delayed(selection_process)(patient=patient)
            for patient in patients
        )
    else:
        target_lesions = [selection_process(patient=patient) for patient in patients]

    lesion_idxs = [idx for patient in target_lesions for idx in patient]


    # Select out the rows of the lesion data for the selected target lesions
    selected_target_lesions = lesion_data.copy().iloc[lesion_idxs]
    # Sort the selected lesion data by index and return
    return selected_target_lesions.sort_index()



def recist_metrics_by_target_count(patient_response: pd.DataFrame,
                                    max_targets: int = 11
                                    ) -> tuple[list, list]:
    """Calculate RECIST response categorization accuracy and sensitivity to progressive disease (PD) when using different numbers of target lesions in the range 1 to max_targets. If there are no patients with a target value, 0 will be saved for that target number.
    
    Parameters
    ----------
    patient_response: pd.DataFrame,
        Response data for a set of patients, like that output by recist_assess.
    max_targets: int = 11
        Maximum number of target lesions to test. 

    Returns
    -------
    accuracy:list
        Accuracy values for each number of targets selected.
    pd_sensitivity: list
        Sensitivity values for each number of targets selected.
    """
    accuracy = []
    pd_sensitivity = []
    for num_targets in range(1, max_targets):
        # Get patients with more than num_targets lesions
        acc_idxs = patient_response['num_lesions'] > num_targets
        # Check that there are patients with enough lesions
        if acc_idxs.any():
            # Get percentage of patients that have the same response category for RECIST (all) as RECIST (num_targets targets) 
            accuracy.append(np.sum(patient_response['RECIST (all)'][acc_idxs] == patient_response[f"RECIST ({num_targets} targets)"][acc_idxs]) / np.sum(acc_idxs) * 100)
        else:
            accuracy.append(0)

        # Get number of patients with more than num_targets and RECIST category is PD
        sens_idxs = np.logical_and(patient_response['num_lesions'] > num_targets, patient_response['RECIST (all)'] == 'PD')
        # Check that there are patients with enough targets with category PD
        if sens_idxs.any():
            # Calculate sensitivity for detecting progressive disease (PD)
            pd_sensitivity.append(np.sum(patient_response[f"RECIST ({num_targets} targets)"][sens_idxs] == 'PD') / np.sum(sens_idxs) * 100)
        else:
            pd_sensitivity.append(0)

    return accuracy, pd_sensitivity
