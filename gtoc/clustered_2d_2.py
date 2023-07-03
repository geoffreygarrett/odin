import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def load_and_preprocess_data(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True)
    return df


def scale_data(df, columns_to_cluster):
    scaler = StandardScaler()
    df_scaled = df.copy()
    df_scaled[columns_to_cluster] = scaler.fit_transform(df[columns_to_cluster])
    return df_scaled, scaler


def perform_kmeans_clustering(df_scaled, columns_to_cluster, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, n_init=10)
    kmeans.fit(df_scaled[columns_to_cluster])
    labels = kmeans.labels_
    return labels, kmeans.cluster_centers_


def plot_clusters(df_scaled, labels, centroids, scaler, columns_to_cluster):
    df_scaled['labels'] = labels





# ID    epoch(MJD)       a(AU)            e             i(deg)         LAN(deg)      argperi(deg)
# M(deg)
if __name__ == "__main__":
    file_path = 'GTOC12_Asteroids_Data.txt'
    n_clusters = 30
    # labels, centroids = perform_kmeans_clustering(df_scaled, columns_to_cluster, n_clusters)
    # plot_clusters(df_scaled, labels, centroids, scaler, columns_to_cluster)

    col_id = "ID"
    col_epoch = "epoch(MJD)"
    col_a = "a(AU)"
    col_e = "e"
    col_i = "i(deg)"
    col_lan = "LAN(deg)"
    col_argp = "argperi(deg)"
    col_m_anomaly = "M(deg)"

    df = load_and_preprocess_data(file_path)
    columns_to_cluster = ['a(AU)', 'e']
    df_scaled, scaler = load_and_preprocess_data(file_path, columns_to_cluster)



    main(columns_to_cluster)
