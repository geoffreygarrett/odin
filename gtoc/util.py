import pandas as pd


def load_and_preprocess_data(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True)
    return df
