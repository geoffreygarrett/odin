import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Import the data
df = pd.read_csv('GTOC12_Asteroids_Data.txt', delim_whitespace=True)

# Drop non-orbital parameters
df_orbital = df.drop(['ID', 'epoch(MJD)'], axis=1)

# Normalize the data
scaler = StandardScaler()
df_scaled = scaler.fit_transform(df_orbital)

# Apply KMeans clustering
n_clusters = 10  # Specify the number of clusters
kmeans = KMeans(n_clusters=n_clusters)
kmeans.fit(df_scaled)
labels = kmeans.labels_
centroids = scaler.inverse_transform(kmeans.cluster_centers_)  # Calculate centroids

# Now we read the file with positions
df_positions = pd.read_csv("output.csv")
df_positions['cluster'] = labels  # Add labels to the data frame

# Create a 2D plot
fig, ax = plt.subplots(figsize=(10, 10))

# Determine the global min and max for all x, y
global_min = df_positions[['r_x', 'r_y']].min().min()
global_max = df_positions[['r_x', 'r_y']].max().max()

# Plot each cluster with a different color and marker
markers = ['o', 'v', 's', 'p', '*', '+', 'x', 'd', '1', '2', '3', '4']  # Define a list of markers

for i in range(n_clusters):
    cluster_data = df_positions[df_positions['cluster'] == i]
    ax.scatter(cluster_data['r_x'], cluster_data['r_y'], s=0.3, label='Cluster ' + str(i+1) + ' (' + str(len(cluster_data)) + ' asteroids)', marker=markers[i % len(markers)])

# Plot centroids
ax.scatter(centroids[:, 0], centroids[:, 1], s=50, color='black', label='Centroids', marker='x')

# Set same limits for all subplots
ax.set_xlim([global_min, global_max])
ax.set_ylim([global_min, global_max])
ax.set_aspect('equal', 'box')  # Set the aspect ratio to be equal

# Show the legend
ax.legend(loc='upper right')

plt.show()
