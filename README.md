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

## **Code**
We build our algorithm on top of SIGMOD'2020 paper "In-Memory Subgraph Matching: an In-depth Study" by Dr. Shixuan Sun and Prof. Qiong Luo.
We kept all the functionalities of the framework and for mospecific details we refer to [Code](https://github.com/RapidsAtHKUST/SubgraphMatching).

## Additional Techniques Supported
|Pilos| the filtering method of Pilos |
|PilosMT| the filtering method of Pilos multi-thread |

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Execute

-dataset the dataset, -qsize query size, -qnumber, -qprop query name , -filter filter algorithm

```zsh
timeout 3600s ./SubgraphMatching.out -dataset $t -qsize $j -qnumber $i -qprop G -filter KFR2
```
### Query Generation
To create queries use -qnumber for the number of queries, -qprop the name of the queries you want -dataset the name of the dataset and -filter GQ.
For example, to generate 400 queries for dblp
```zsh
./SubgraphMatching.out -dataset dblp -qnumber 400 -qprop G -filter QG
```
### Index Generation
To create Eigenvalue index for a dataset use -filter GQ.
For example, to generate 400 queries for dblp
```zsh
./SubgraphMatching.out -dataset dblp -filter EC
```
## Reference

Please cite our work in your publications if it helps your research:

```
Paper under submission
```
## 
## Note
The code is under refactor. Expected to be available 15-11-23.

