import pandas as pd
from joblib import Parallel, delayed
from utils import (
    generate_synthetic_lesions,
    truncate_normal_distribution,
    
)


def pipe(radiomic_features_filepath: str,
         num_sim_patients: int = 10000,
         expected_num_lesions: int = 10,
         location_label: str = "LABEL",
         parallel: bool = False,
         n_jobs: int = -1,
         random_seed: int | None = None
         ):
    # Load radiomics data
    rad_data = pd.read_csv(radiomic_features_filepath)

    # Extract the slice thickness from the radiomics data
    rad_data['slice_thickness'] = [float(i.split(', ')[-1].strip('()')) for i in rad_data['diagnostics_Image-interpolated_Spacing']]

    # Calculate the true volume of the segmentation
    rad_data['volume_cc_contoured'] = rad_data['original_shape_VoxelVolume'] * rad_data['slice_thickness'] / 1000

    # Truncated normal distribution setup - used for diameter change
    trunc_normal_dist = truncate_normal_distribution(min=-1,
                                                     max=3,
                                                     mean=0,
                                                     std_dev=0.3)
    
    if parallel:
        synthetic_lesions = Parallel(n_jobs=n_jobs)(
            delayed(generate_synthetic_lesions)(
                diameter_change_dist=trunc_normal_dist,
                base_radiomic_data=rad_data,
                expected_num_lesions=expected_num_lesions,
                patient_id=pat_idx,
                location_label=location_label,
                random_seed=random_seed
            )
            for pat_idx in range(num_sim_patients)
        )
    else:
        synthetic_lesions = [
            generate_synthetic_lesions(
                diameter_change_dist=trunc_normal_dist,
                base_radiomic_data=rad_data,
                expected_num_lesions=expected_num_lesions,
                patient_id=pat_idx,
                location_label=location_label,
                random_seed=random_seed
            )
            for pat_idx in range(num_sim_patients)
        ]
    
    print(synthetic_lesions)

    return None


if __name__ == '__main__':
    pipe("data/rawdata/SARC021/SARC021_radiomics.csv",
         num_sim_patients = 1,
         expected_num_lesions = 3,
         random_seed=10)