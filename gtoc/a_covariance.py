import seaborn as sns
import matplotlib.pyplot as plt
from util import load_and_preprocess_data
import numpy as np

# Define column names for the dataset
col_names = ['ID', 'epoch(MJD)', 'a(AU)', 'e', 'i(deg)', 'LAN(deg)', 'argperi(deg)', 'M(deg)']

# Import the data
df = load_and_preprocess_data('GTOC12_Asteroids_Data.txt')

# rename columns and convert to SI (we leave a as km)
df['id'] = df['ID']
df['epoch_mjd'] = df['epoch(MJD)']
df['a_km'] = df['a(AU)'] * 149597870.7
df['e'] = df['e']
df['i'] = df['i(deg)'] * (np.pi / 180)
df['RAAN'] = df['LAN(deg)'] * (np.pi / 180)
df['argp'] = df['argperi(deg)'] * (np.pi / 180)
df['M'] = df['M(deg)'] * (np.pi / 180)
df = df.drop(columns=['ID', 'epoch(MJD)', 'a(AU)', 'i(deg)', 'LAN(deg)', 'argperi(deg)', 'M(deg)'])

# Plot distribution for each selected column
for col in df.columns.drop(['id', 'epoch_mjd']):
    sns.displot(df[col])
    plt.title(f'Distribution for {col}')
    plt.show()
