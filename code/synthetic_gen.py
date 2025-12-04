import numpy as np
import pandas as pd
from scipy.stats import rv_continuous, truncnorm


def truncate_normal_distribution(min: int = -1,
                                 max: int = 3,
                                 mean: float = 0,
                                 std_dev: float = 0.3
                                 ) -> rv_continuous:
    """Generate a truncated normal distribution.

    Parameters
    ----------
    min: int
    max:
    mean:
    std_dev:

    Returns
    -------
    rv_continuous instance
        A truncated normal continuous random variable.
    """
    return truncnorm((min - mean) / std_dev, (max - mean) / std_dev, loc=mean, scale=std_dev)


def volume_calc(diameter):
    """Calculate volume of lesion based on the diameter and assumption it is a sphere"""
    return (4/3 * np.pi * (diameter/2)**3) / 1000


def generate_synthetic_lesions(base_radiomic_data: pd.DataFrame,
                               expected_num_lesions: int = 10,
                               patient_id: int | str = 0,
                               location_label: str = "LABEL",
                               lesion_selection_rng: int | None = None,
                               ) -> pd.DataFrame:
    """
    Generate synthetic lesion measurement data. 
    
    Parameters
    ----------
    base_radiomic_data: pd.DataFrame
        Ground truth radiomic data to base simulations on. Must contain the following columns:
            * original_shape_Maximum2DDiameterSlice
            * volume_cc_contoured (original_shape_VoxelVolume * slice_thickness / 1000)
            * location_label matching input argument
            * original_shape_Maximum3DDiameter
            * original_shape_MajorAxisLength
            * original_shape_MinorAxisLength
    expected_num_lesions: int
        Expected number of lesions to use for poisson distribution
    patient_id: int | str
        Patient ID to label all these lesions with
    location_label: str
        Name of tumour location label in the 
    """
    if location_label not in base_radiomic_data.columns:
        raise ValueError(f"{location_label} is not a column name in the provided radiomic data.")

    n_lesions = 0
    # Select a number of lesions until the value is between 1 and 30
    while n_lesions < 1 or n_lesions > 30:
        n_lesions = lesion_selection_rng.poisson(expected_num_lesions)

    # Generated truncated normal distribution setup - used for diameter change
    diameter_change_dist = truncate_normal_distribution(min=-1,
                                                        max=3,
                                                        mean=0,
                                                        std_dev=0.3)

    # Initialize list to hold lesion data
    lesion_data = {}
    
    for lesion_idx in range(n_lesions):
        # Select a random index
        ind = lesion_selection_rng.choice(base_radiomic_data.index)

        # Generate diameter (pre-treatment)
        diameter_pre = base_radiomic_data['original_shape_Maximum2DDiameterSlice'][ind]

        # Generate diameter change
        diameter_change = diameter_change_dist.rvs(size=1)[0]

        # Generate diameter (post-treatment)
        diameter_post = diameter_pre + diameter_pre * diameter_change

        # Generate location tag
        location = base_radiomic_data[location_label][ind]

        # Volume (contoured)
        volume = base_radiomic_data['volume_cc_contoured'][ind]
        # Other data
        diameter_3D_max = base_radiomic_data['original_shape_Maximum3DDiameter'][ind]
        diameter_major_ax = base_radiomic_data['original_shape_MajorAxisLength'][ind]
        diameter_minor_ax = base_radiomic_data['original_shape_MinorAxisLength'][ind]

        # Append metadata for this lesion to main list
        lesion_data[lesion_idx] = {'patient_id': str(patient_id),
                                   'lesion_idx': lesion_idx,
                                   'diameter_pre': diameter_pre,
                                   'diameter_change': diameter_change,
                                   'diameter_post': diameter_post,
                                   'location': location,
                                   'volume_cc_contoured': volume,
                                   'diameter_3D_max': diameter_3D_max,
                                   'diameter_major_ax': diameter_major_ax,
                                   'diameter_minor_ax': diameter_minor_ax}
        
    
    # Convert to DataFrame
    synthetic_lesions = pd.DataFrame.from_dict(lesion_data, orient='index')
    
    return synthetic_lesions



def generate_synthetic_patients(num_sim_patients: int, 
                                base_radiomic_data: pd.DataFrame,
                                lesion_selection_rng: np.random.Generator,
                                expected_num_lesions: int = 10,
                                location_label: str = "LABEL",
                                ) -> pd.DataFrame:
    """Generate synthetic patient lesion data from an existing dataset"""
    synth_lesion_list = [
        generate_synthetic_lesions(
            base_radiomic_data=base_radiomic_data,
            expected_num_lesions=expected_num_lesions,
            patient_id=pat_idx,
            location_label=location_label,
            lesion_selection_rng=lesion_selection_rng
        )
        for pat_idx in range(num_sim_patients)
    ]

    synth_lesion_pd = pd.concat(synth_lesion_list, ignore_index=True)

    # Add volumes
    synth_lesion_pd['volume_cc_pre'] = volume_calc(synth_lesion_pd['diameter_pre'])
    synth_lesion_pd['volume_cc_post'] = volume_calc(synth_lesion_pd['diameter_post'])
    synth_lesion_pd['volume_cc_3Dmax'] = volume_calc(synth_lesion_pd['diameter_3D_max'])
    synth_lesion_pd['volume_cc_majorAx'] = volume_calc(synth_lesion_pd['diameter_major_ax'])
    synth_lesion_pd['volume_cc_minorAx'] = volume_calc(synth_lesion_pd['diameter_minor_ax'])

    return synth_lesion_pd