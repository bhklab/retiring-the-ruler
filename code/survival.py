import pandas as pd
from sksurv.nonparametric import kaplan_meier_estimator
from damply import dirs


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



clinical_path = dirs.RAWDATA / "SARC021" / "SARC021_clinical.csv"
clinical_df = pd.read_csv(clinical_path)

recist_path = dirs.RAWDATA / "SARC021" / "SARC021_RECIST.xlsx"
recist_df = pd.read_excel(recist_path, sheet_name="OVERALL")


# For each patient_id, get the row with the lowest value in the "Study day of response assessment" column
earliest_recist_df = recist_df.loc[recist_df.groupby("USUBJID")["Study day of response assessment"].idxmin()]

# Drop rows with 'NE' values in the "RECIST Overall Response Assessment" column
earliest_recist_df = earliest_recist_df[earliest_recist_df['RECIST Overall Response Assessment'] != 'NE']

# Merge recist data with clinical data on patient ID
merged_df = pd.merge(earliest_recist_df, clinical_df, left_on="USUBJID", right_on="USUBJID")

# don't think this is right
recist_survival = merged_df.groupby("RECIST Overall Response Assessment").apply(calculate_survival, include_groups=False).reset_index()

print(recist_survival)