import numpy as np
import pandas as pd


def recist_assess(lesion_data:pd.DataFrame) -> list:

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

    synth_response = pd.DataFrame({'patient_id': patients, 
                                'num_lesions': lesion_counts, 
                                'SLD_pre': SLD_pre, 
                                'SLD_post': SLD_post, 
                                'SLD_chg': SLD_chg, 
                                'RECIST (all)': RECIST_response})
    
    return synth_response
