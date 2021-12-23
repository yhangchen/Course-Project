import numpy as np
import matplotlib.pyplot as plt


# Lloyds algorithm for k-means
def kmeans(X, k, max_iter=100, tolerance=10**(-3)):
    n_samples = X.shape[0]
    n_features = X.shape[1]
    classifications = np.zeros(n_samples, dtype=np.int64)

    # Choose initial cluster centroids randomly
    I = np.random.choice(n_samples, k)
    centroids = X[I, :]

    loss = 0
    for m in range(0, max_iter):
        # Compute the classifications
        for i in range(0, n_samples):
            distances = np.zeros(k)
            for j in range(0, k):
                distances[j] = np.sqrt(
                    np.sum(np.power(X[i, :] - centroids[j], 2)))
            classifications[i] = np.argmin(distances)

        # Compute the new centroids and new loss
        new_centroids = np.zeros((k, n_features))
        new_loss = 0
        for j in range(0, k):
            # compute centroids
            J = np.where(classifications == j)
            X_C = X[J]
            new_centroids[j] = X_C.mean(axis=0)

            # Compute loss
            for i in range(0, X_C.shape[0]):
                new_loss += np.sum(np.power(X_C[i, :] - centroids[j], 2))

        # Stopping criterion
        if np.abs(loss - new_loss) < tolerance:
            return new_centroids, classifications, new_loss

        centroids = new_centroids
        loss = new_loss

    print("Failed to converge!")
    return centroids, classifications, loss
