import numpy as np
import pandas as pd


def recist_assess(lesion_data:pd.DataFrame) -> pd.DataFrame:
    """Calculate RECIST response category from lesion diameter data."""
    patients, lesion_counts = np.unique(lesion_data['patient_id'], return_counts=True)

    # Calculate the sum of longest diameters (SLD) for pre, post
    SLD_pre = [np.sum(lesion_data['diameter_pre'][lesion_data['patient_id'] == p]) for p in patients]
    SLD_post = [np.sum(lesion_data['diameter_post'][lesion_data['patient_id'] == p]) for p in patients]

    # Calculte the percentage different in the sum of longest diameters (SLD) pre and post
    SLD_chg = [(SLD_post[i] - SLD_pre[i])/SLD_pre[i] * 100 for i in range(len(SLD_pre))]

    # Initialize RECIST response vector
    RECIST_response = []

    # Classify each patient's response based on the percent change in SLD
    for idx, _patient in enumerate(patients):
        if SLD_chg[idx] > 20:
            RECIST_response.append('PD')
        elif -30 < SLD_chg[idx] <= 20:
            RECIST_response.append('SD')
        elif -100 < SLD_chg[idx] <= -30:
            RECIST_response.append('PR')
        elif SLD_chg[idx] == -100:
            RECIST_response.append('CR')

    patient_response = pd.DataFrame({'patient_id': patients, 
                                'num_lesions': lesion_counts, 
                                'SLD_pre': SLD_pre, 
                                'SLD_post': SLD_post, 
                                'SLD_chg': SLD_chg, 
                                'RECIST (all)': RECIST_response})
    
    return patient_response



def select_target_lesions(num_lesions:int,
                          lesion_data:pd.DataFrame,
                          lesion_selection_rng: np.random.Generator) -> pd.DataFrame:
    """Randomly select specified number of target lesions from lesion data, ensuring no more than 2 are selected for each location
       following RECIST specifications.
       
       Returns the selected target lesion rows from lesion data.
       """
    patients = np.unique(lesion_data['patient_id'])

    # ensure that there is not more than 2 target lesions per location
    target_lesions = []
    for patient in patients:
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
            target_lesions.extend(lesion_selection_rng.choice(selected_lesion_idx, num_lesions))
        else:
            target_lesions.extend(selected_lesion_idx)
    
    # Select out the rows of the lesion data for the selected target lesions
    selected_target_lesions = lesion_data.copy().iloc[target_lesions]
    # Sort the selected lesion data by index and return
    return selected_target_lesions.sort_index()



def recist_metrics_by_target_count(patient_response: pd.DataFrame,
                                    max_targets: int = 11
                                    ) -> tuple[list, list]:
    accuracy = []
    pd_sensitivity = []
    for num_targets in range(1, max_targets):
        # Get patients with more than num_targets lesions
        acc_idxs = patient_response['num_lesions'] > num_targets
        # Get percentage of patients that have the same response category for RECIST (all) as RECIST (num_targets targets) 
        accuracy.append(np.sum(patient_response['RECIST (all)'][acc_idxs] == patient_response[f"RECIST ({num_targets} targets)"][acc_idxs]) / np.sum(acc_idxs) * 100)

        # Get number of patients with more than num_targets and RECIST category is PD
        sens_idxs = np.logical_and(patient_response['num_lesions'] > num_targets, patient_response['RECIST (all)'] == 'PD')
        # Calculate sensitivity for detecting progressive disease (PD)
        pd_sensitivity.append(np.sum(patient_response[f"RECIST ({num_targets} targets)"][sens_idxs] == 'PD') / np.sum(sens_idxs) * 100)

    return accuracy, pd_sensitivity
