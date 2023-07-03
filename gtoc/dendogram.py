import pandas as pd
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram
from matplotlib import pyplot as plt

# Import the data
df = pd.read_csv('GTOC12_Asteroids_Data.txt', delim_whitespace=True)

# Take a sample of 1000 points
# df = df.sample(n=500)

# Drop non-orbital parameters
df_orbital = df.drop(['ID', 'epoch(MJD)'], axis=1)

# Normalize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df_orbital)

# Perform hierarchical clustering
Z = linkage(df_scaled, 'ward')

# Convert the linkage matrix to a DataFrame for easier manipulation
df_Z = pd.DataFrame(Z, columns=['cluster_1', 'cluster_2', 'distance', 'new_cluster_size'])

# Save to CSV
df_Z.to_csv('linkage_matrix.csv', index=False)

# # Or save to JSON
# df_Z.to_json('linkage_matrix.json')

# # Plot the dendrogram
# fig = plt.figure(figsize=(10, 10))
# dn = dendrogram(Z)
#
# plt.show()
