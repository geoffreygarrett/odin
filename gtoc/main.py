import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Import the data
df = pd.read_csv('GTOC12_Asteroids_Data.txt', delim_whitespace=True)

# Drop non-orbital parameters
df_orbital = df.drop(['ID', 'epoch(MJD)', 'M(deg)'], axis=1)

# Normalize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df_orbital)

# Apply KMeans clustering
n_clusters = 10  # Specify the number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(df_scaled)
labels = kmeans.labels_

# Now we read the file with positions
df_positions = pd.read_csv("output.csv")
df_positions['cluster'] = labels  # Add labels to the data frame

fig = plt.figure(figsize=(15, 15))

# Create a grid of subplots
n_rows = int(np.ceil(np.sqrt(n_clusters)))
n_cols = int(np.ceil(n_clusters / n_rows))

# Determine the global min and max for all x, y, and z
global_min = df_positions[['r_x', 'r_y', 'r_z']].min().min()
global_max = df_positions[['r_x', 'r_y', 'r_z']].max().max()

for i in range(n_clusters):
    ax = fig.add_subplot(n_rows, n_cols, i+1, projection='3d')
    cluster_data = df_positions[df_positions['cluster'] == i]
    ax.scatter(cluster_data['r_x'], cluster_data['r_y'], cluster_data['r_z'], s=0.1)
    ax.set_title('Cluster ' + str(i+1))

    # Set same limits for all subplots
    ax.set_xlim([global_min, global_max])
    ax.set_ylim([global_min, global_max])
    ax.set_zlim([global_min, global_max])

    # Set the aspect ratio to be equal
    ax.set_box_aspect([1,1,1])

plt.tight_layout()
plt.show()
