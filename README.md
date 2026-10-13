# PhyloVelo_forward_benchmark

This repository is a research disclosure fork of [PhyloVelo](https://github.com/kunwang34/PhyloVelo) for cross-method RNA velocity benchmarking. It starts from upstream commit `30b93689a4a63f9363d796343cfb70a2bb238530` and retains the original BSD 3-Clause license. The [PhyloVelo paper](https://doi.org/10.1038/s41587-023-01887-5) defines projected trajectories in the backward, ancestral direction. Many downstream velocity comparisons instead evaluate transitions in the forward direction, so this fork makes that projection direction explicit.

## Changes in this fork

- `phylovelo.velocity_embedding(..., direction="backward")` keeps the upstream behavior and remains the default.
- `direction="forward"` projects the inferred gene-wise slope along increasing phylogenetic time. It recalculates neighborhood transition weights and the embedding vector with the forward slope; it does not negate the completed backward embedding.
- `phylovelo.benchmark.add_embedding_velocity_graph(...)` provides an optional cosine graph from two-dimensional UMAP velocities for embedding-based benchmark metrics. This graph is a derived comparison artifact, not a native gene-expression-space velocity graph.
- Packaging dependencies were declared in `setup.py`. The repository tree was reduced to the method source, core documentation, license, and dependency files. Benchmark execution and plotting scripts are maintained separately.

The MEG selection, gene-wise slope inference, and original pseudotime implementation are otherwise unchanged. PhyloVelo estimates one slope per gene; this fork does not claim a cell-by-gene velocity layer. Forward embedding results should be described as an adaptation of the published backward projection when reporting benchmark comparisons.

## Installation and use

```bash
python -m pip install -e .
```

After constructing a `phylovelo.scData` object with expression, embedding coordinates, and a lineage-time proxy as described upstream:

```python
import phylovelo as pv

pv.velocity_inference(sd, time=sd.cell_generation, target="count")
pv.velocity_embedding(sd, target="count", n_neigh=15, direction="backward")
backward_embedding = sd.velocity_embeded.copy()
pv.calc_phylo_pseudotime(sd, n_neighbors=15)
pv.velocity_embedding(sd, target="count", n_neigh=15, direction="forward")
forward_embedding = sd.velocity_embeded.copy()
```

The original graph pseudotime calculation uses the backward embedding and is therefore run before switching to the forward projection in this example.

## Reproducibility and scope

[`requirements.lock.txt`](requirements.lock.txt) records the historical local PhyloVelo inference environment and can help reconstruct it. That environment lacked scVelo and the complete plotting stack required by our plotting workflow, so benchmark figures were rendered in a separate environment. The Python version was not captured in the lock file, and its packages have not been reinstalled and validated together for this fork.

The forward projection and graph were checked on four local lineage benchmark datasets (92–95); the local scripts, input datasets, generated h5ad files, figures, and metric scores are outside this source repository. The graph uses UMAP geometry and does not supply the cell-by-gene velocities needed by gene-space metrics. For datasets whose sampling day or generation was used as the inference clock, evaluating against the same label would not provide an independent direction test.

## Attribution

Original method: Kun Wang *et al.*, “[PhyloVelo enhances transcriptomic velocity field mapping using monotonically expressed genes](https://doi.org/10.1038/s41587-023-01887-5),” *Nature Biotechnology* (2023). Original software: [kunwang34/PhyloVelo](https://github.com/kunwang34/PhyloVelo). See [LICENSE](LICENSE) for the retained BSD 3-Clause terms.
