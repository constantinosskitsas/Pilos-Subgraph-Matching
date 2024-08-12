# **Pilos: Scalable Subgraph Matching by Online Spectral Filtering.**

## **Introduction**
Subgraph matching seeks all the occurrences of a
query graph inside another graph. As it reduces to subgraph
isomorphism, it is NP-hard. Current methods reduce the computation by filtering the set of candidates on which they run
subgraph isomorphism. Nevertheless, when the query is large, the
number of candidates grows rapidly, rendering current methods
largely ineffective in pruning and fail to answer even within one
hour. A primary reason for this ineffectiveness is their inability to
effectively consider the query graph structure in the computation.
In this paper, we propose PILOS, a novel matching algorithm
that substantially improves the filtering phase of a typical
matching algorithm and computes up to 60% fewer candidates
for verification. PILOS uses (i) an offline light-weight index-based
phase, which leverages the top graph Laplacian eigenvalues of
query and data node neighborhoods to reduce candidates via
neighborhood filtering and (ii) an online phase, which further
prunes candidates stored in an auxiliary data structure by applying the interlacing theorem on their graph Laplacian spectra.
Our thorough experimental study shows that, on average, PILOS
resolves queries in 16% less time and leaves 21% fewer unresolved
queries after a lapse of 10 minutes than the state of the art.

## **Code**
We build our algorithm using SIGMOD'2020 paper "In-Memory Subgraph Matching: an In-depth Study" by Dr. Shixuan Sun and Prof. Qiong Luo.
We kept all the functionalities of the framework and for mospecific details we refer to [Code](https://github.com/RapidsAtHKUST/SubgraphMatching).

## Additional Techniques Supported
|Algorithm|Description|Execution code
|:--------:|:------------:|:------------:
|Pilos | the filtering method of Pilos | PL

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```

## Execute

```zsh
timeout 600s ./SubgraphMatching.out -dataset $t -qsize $j -qnumber $i -qprop G -filter PL -alpha $alpha -beta 0 -n 5 -num 100000
```
Example
```zsh
timeout 600s ./SubgraphMatching.out -dataset dblp -qsize 32 -qnumber 1 -qprop G -filter PL -alpha 125 -SF results
```
## General Parameters
|Execution code|Description|
|:--------:|:------------:|
|-dataset | dataset name|
|-qsize | query size|
|-qnumber | query number|
|-filter | filter algorithm ,GQL,CFL,DPiso,PL|
|-num | number of matchings|
|-SF | Save file name|

## PL Parameters
|Execution code|Description|
|:--------:|:------------:|
|-alpha | alpha parameter |
|-beta | beta parameter|

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
## Note

## Reproducability
To repeat the experiments on the paper we propose to use the scrips we provide. build/matching/*.sh

