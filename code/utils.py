import numpy as np
import pandas as pd
import random
from scipy.stats import truncnorm, rv_continuous


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


def generate_synthetic_lesions(diameter_change_dist: rv_continuous,
                               base_radiomic_data: pd.DataFrame,
                               expected_num_lesions: int = 10,
                               patient_id: int | str = 0,
                               location_label: str = "LABEL",
                               random_seed: int | None = None,
                               ) -> pd.DataFrame:
    """
    Generate synthetic lesion measurement data. 
    
    Parameters
    ----------
    diameter_change_dist: rv_continuous
        A truncated normal distribution to generate the diameter change for the simulated lesion from.
        Can be made by running `truncate_normal_distribution` function.
    base_radiomic_data: pd.DataFrame
        Ground truth radiomic data to base simulations on. Must contain the following columns:
            * original_shape_Maximum2DDiameterSlice
            * volume_cc_contoured (original_shape_VoxelVolume * slice_thickness / 1000)
            * location_label matching input argument
            * original_shape_Maximum3DDiameter
            * original_shape_MajorAxisLength
            * original_shape_MinorAxisLength
    mu: int
        Expected number of events
    """
    # Generate random number of lesions with poisson distribution
    rng = np.random.default_rng(random_seed)
    n_lesions = 0
    # Select a number of lesions until the value is between 1 and 30
    while n_lesions < 1 or n_lesions > 30:
        n_lesions = rng.poisson(expected_num_lesions)

    # Initialize list to hold lesion data
    lesion_data = {}

    for lesion_idx in range(n_lesions):
        # Select a random index
        ind = rng.choice(base_radiomic_data.index)

        # Generate diameter (pre-treatment)
        diameter_pre = base_radiomic_data['original_shape_Maximum2DDiameterSlice'][ind]

        # Generate diameter change
        diameter_change = diameter_change_dist.rvs(size=1, random_state=random_seed)[0]

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
        lesion_data[lesion_idx] = {'patient_id': patient_id,
                                   'lesion_idx': lesion_idx,
                                   'diameter_pre': diameter_pre,
                                   'diameter_change': diameter_change,
                                   'diameter_post': diameter_post,
                                   'location': location,
                                   'volume_cc_contoured': volume,
                                   'diameter_3D_max': diameter_3D_max,
                                   'diameter_3D_max': diameter_major_ax,
                                   'diameter_minor_ax': diameter_minor_ax}
        
    
    # Convert to DataFrame
    synthetic_lesions = pd.DataFrame.from_dict(lesion_data, orient='index')
    
    
    return synthetic_lesions