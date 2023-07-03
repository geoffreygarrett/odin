import pandas as pd
from scipy.cluster.hierarchy import fcluster
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the linkage matrix
df_Z = pd.read_csv('linkage_matrix.csv')

# Convert back to a numpy array
Z = df_Z.to_numpy()

# Decide the number of clusters you want
num_clusters = 10

# Perform the cut
cluster_assignments = fcluster(Z, t=num_clusters, criterion='maxclust')

# Now let's plot the results
df_positions = pd.read_csv("output.csv")  # Load positions

# Assign the cluster labels to the data points
df_positions['cluster'] = cluster_assignments

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set global min and max for all axes
global_min = df_positions[['r_x', 'r_y', 'r_z']].min().min()
global_max = df_positions[['r_x', 'r_y', 'r_z']].max().max()

# Plot each cluster with a different color
for i in range(1, num_clusters+1):  # Clusters are numbered from 1
    cluster_data = df_positions[df_positions['cluster'] == i]
    ax.scatter(cluster_data['r_x'], cluster_data['r_y'], cluster_data['r_z'], s=0.1, label='Cluster ' + str(i))

# Set same limits for all axes
ax.set_xlim([global_min, global_max])
ax.set_ylim([global_min, global_max])
ax.set_zlim([global_min, global_max])

plt.legend()
plt.show()
