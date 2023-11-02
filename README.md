# **Pilos: Scalable Subgraph Matching by Online Spectral Filtering.**

## **Introduction**
Subgraph matching seeks all the occurrences of a query graph inside
another graph. As it reduces to subgraph isomorphism, it is NP-
complete. Current methods reduce the computation by filtering
the set of candidates on which they run subgraph isomorphism.
Nevertheless, when the query is large, the number of candidates
grows rapidly, rendering current methods largely ineffective in
pruning and fail to return any answer even within one hour. A
primary reason for this ineffectiveness is their inability to effectively
consider the query graph structure in the computation.
In this paper, we propose Pilos, a novel matching algorithm
that combines a state-of-the-art backtracking process with (i) an
offline index-based phase, which leverages the top graph Laplacian
eigenvalues of query and data node neighborhoods to reduce can-
didates via neighborhood filtering and (ii) an online phase, which
further prunes candidates stored in an auxiliary data structure by
applying the interlacing theorem on their graph Laplacian spectra.
Our thorough experimental study shows that, on average, Pilos
resolves queries in 20% less time and leaves 12.5% fewer unresolved
queries after a lapse of 10 minutes than the state of the art.