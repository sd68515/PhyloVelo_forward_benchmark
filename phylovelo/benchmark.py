"""Embedding-only velocity graph for comparison with low-dimensional methods.

This is a benchmark-derived cosine graph, not a native PhyloVelo transition
graph or a gene-expression-space scVelo graph.
"""

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def add_embedding_velocity_graph(adata, *, vkey, velocity_key, n_neighbors=30):
    positions = np.asarray(adata.obsm["X_umap"], dtype=np.float64)
    velocity = np.asarray(adata.obsm[velocity_key], dtype=np.float64)
    if positions.shape != velocity.shape or positions.shape != (adata.n_obs, 2):
        raise ValueError("X_umap and embedded velocity must be cells x 2")
    if not np.isfinite(positions).all() or not np.isfinite(velocity).all():
        raise ValueError("UMAP positions and embedded velocity must be finite")
    if adata.n_obs < 2 or n_neighbors < 1:
        raise ValueError("Need at least two cells and one neighbour")

    k = min(n_neighbors + 1, adata.n_obs)
    indices = np.asarray(cKDTree(positions).query(positions, k=k)[1]).reshape(adata.n_obs, k)
    rows = np.repeat(np.arange(adata.n_obs), k)
    cols = indices.ravel()
    keep = rows != cols
    rows, cols = rows[keep], cols[keep]
    displacement = positions[cols] - positions[rows]
    numerator = np.einsum("ij,ij->i", velocity[rows], displacement)
    denominator = np.linalg.norm(velocity[rows], axis=1) * np.linalg.norm(displacement, axis=1)
    cosine = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    shape = (adata.n_obs, adata.n_obs)
    positive = csr_matrix((np.clip(cosine, 0, 1), (rows, cols)), shape=shape)
    negative = csr_matrix((np.clip(cosine, -1, 0), (rows, cols)), shape=shape)
    positive.eliminate_zeros()
    negative.eliminate_zeros()
    adata.uns[f"{vkey}_graph"] = positive
    adata.uns[f"{vkey}_graph_neg"] = negative
    adata.uns[f"{vkey}_embedding_graph_params"] = {
        "basis": "X_umap",
        "velocity_key": velocity_key,
        "n_neighbors": int(min(n_neighbors, adata.n_obs - 1)),
        "edge_score": "cosine(2D velocity_i, X_umap_j - X_umap_i)",
        "graph_kind": "benchmark-derived 2D graph; not a native model graph",
    }
