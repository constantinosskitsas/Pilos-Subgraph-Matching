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
For more details of how to use the framework we refer to [Code](https://github.com/RapidsAtHKUST/SubgraphMatching).

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
After compiling the source code, you can find the binary file 'SubgraphMatching.out' under the 'build/matching' directory. 
Execute the binary with the following command ./SubgraphMatching.out -d data_graphs -q query_graphs
-filter method_of_filtering_candidate_vertices -order method_of_ordering_query_vertices -engine method_of_enumerating_partial_results -num number_of_embeddings,
in which -d specifies the input of the data graphs and -q specifies the input of the query graphs.
The -filter parameter gives the filtering method, the -order specifies the ordering method, and the -engine
sets the enumeration method. The -num parameter sets the maximum number of embeddings that you would like to find.
If the number of embeddings enumerated reaches the limit or all the results have been found, then the program will terminate.
Set -num as 'MAX' to find all results.

```zsh
./SubgraphMatching.out -d ../../test/sample_dataset/test_case_1.graph -q ../../test/sample_dataset/query1_positive.graph -filter SF -order GQL -engine LFTJ -num MAX
```
## Reference

Please cite our work in your publications if it helps your research:

```
Paper under submission
```


