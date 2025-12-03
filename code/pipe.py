import numpy as np
import pandas as pd

from utils import (
    generate_synthetic_patients,
)


def pipe(radiomic_features_filepath: str,
         num_sim_patients: int = 10000,
         expected_num_lesions: int = 10,
         location_label: str = "LABEL",
         random_seed: int | None = None
         ):
    # Load radiomics data
    rad_data = pd.read_csv(radiomic_features_filepath)

    # Extract the slice thickness from the radiomics data
    rad_data['slice_thickness'] = [float(i.split(', ')[-1].strip('()')) for i in rad_data['diagnostics_Image-interpolated_Spacing']]

    # Calculate the true volume of the segmentation
    rad_data['volume_cc_contoured'] = rad_data['original_shape_VoxelVolume'] * rad_data['slice_thickness'] / 1000

    synth_lesions = generate_synthetic_patients(num_sim_patients=num_sim_patients,
                                                base_radiomic_data=rad_data,
                                                expected_num_lesions=expected_num_lesions,
                                                location_label=location_label,
                                                random_seed=random_seed)
    
    patient_ids, lesion_counts = np.unique(synth_lesions['patient_id'], return_counts=True)

    print(patient_ids, lesion_counts)
    
    return None


if __name__ == '__main__':
    pipe("data/rawdata/SARC021/SARC021_radiomics.csv",
         num_sim_patients = 2,
         expected_num_lesions = 3,
         random_seed=10)