import pandas as pd
from sksurv.nonparametric import kaplan_meier_estimator
from damply import dirs
from plot import plot_survival_curve
from pathlib import Path
import click

def calculate_survival(response_df: pd.DataFrame):
    """Calculate survival probabilities and confidence intervals using the Kaplan-Meier estimator for a given response DataFrame."""
    time, survival_prob, conf_int = kaplan_meier_estimator(
        response_df['E_OS'].astype(bool),
        response_df['T_OS'],
        conf_type="log-log"
    )

    survival_df = pd.DataFrame({
        'time': time,
        'survival_prob': survival_prob,
        'conf_int_lower': conf_int[0],
        'conf_int_upper': conf_int[1]
    })

    return survival_df

@click.command()
@click.option('--clinical_path', type=click.Path(exists=True), help='Path to the clinical CSV file.')
@click.option('--recist_path', type=click.Path(exists=True), help='Path to the RECIST response CSV file.')
@click.option('--response_categories', type=click.STRING, default='PD,SD', help='Comma-separated list of response categories to include in the analysis (e.g., PD,SD)')
def main(clinical_path: Path | None = None,
         recist_path: Path | None = None,
         response_categories: list[str] | str = ['PD', 'SD'],):
    """Main function to calculate survival probabilities and plot Kaplan-Meier survival curves for different response categories."""
    if clinical_path is None:
        raise ValueError("Clinical data could not be loaded. Please check the clinical_path.")
    if recist_path is None:
        raise ValueError("RECIST data could not be loaded. Please check the recist_path.")
    
    if isinstance(response_categories, str):
        response_categories = [cat.strip() for cat in response_categories.split(',')]

    # Load clinical and RECIST response data
    clinical_df = pd.read_csv(clinical_path)
    recist_df = pd.read_csv(recist_path)

    # For each patient_id, get the row with the lowest value in the "Study day of response assessment" column
    earliest_recist_df = recist_df.loc[recist_df.groupby("USUBJID")["Study day of response assessment"].idxmin()]

    # Select out only the response categories we want to include in the analysis
    earliest_recist_df = earliest_recist_df[(earliest_recist_df['RECIST Overall Response Assessment'].isin(response_categories))]

    # Merge recist data with clinical data on patient ID
    merged_df = pd.merge(earliest_recist_df, clinical_df, left_on="USUBJID", right_on="USUBJID")

    # Calculate the KM survival curve for each response category and store in a dictionary
    recist_survival = {}
    for response_category, group in merged_df.groupby("RECIST Overall Response Assessment"):
        survival_df = calculate_survival(group)
        survival_df['response_category'] = response_category
        recist_survival[response_category] = survival_df

    # Combine all the survival data into a single DataFrame for plotting
    recist_survival_df = pd.concat(recist_survival.values(), ignore_index=True)

    # Plot the survival curve for each response category
    fig = plot_survival_curve(recist_survival_df, category_col='response_category', save_path=dirs.PROCDATA / 'SARC021')

    return fig


if __name__ == "__main__":
    main()

