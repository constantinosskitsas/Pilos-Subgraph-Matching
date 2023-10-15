

#ifndef SUBGRAPHMATCHING_FILTERVERTICES_H
#define SUBGRAPHMATCHING_FILTERVERTICES_H

#include "graph/graph.h"
#include <vector>
#include <Eigen/Core>

using namespace Eigen;

class FilterVertices {
public:
    static bool EFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,int top_s);
    static bool LDFFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,bool isEigenCheck,int top_s);
    static bool NLFFilter(Graph* data_graph, Graph* query_graph, ui** &candidates, ui* &candidates_count,bool isEigenCheck,int top_s);
    static bool GQLFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,bool isEigenCheck,int top_s);
    static bool TSOFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,
                          ui *&order,TreeNode *&tree,bool isEigenCheck,int top_s);
    static bool CFLFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,
                              ui *&order, TreeNode *&tree, bool isEigenCheck,int top_s);
    static bool DPisoFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,
                            ui *&order, TreeNode *&tree, bool isEigenCheck,int top_s);

    static bool CECIFilter(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count,
                           ui *&order, TreeNode *&tree,   std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                           std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,bool isEigenCheck,int top_s);

    // static bool VCFilter(const Graph* data_graph, const Graph* query_graph, ui **&candidates, ui *&candidates_count);

    static void computeCandidateWithNLF(Graph *data_graph, Graph *query_graph, VertexID query_vertex, ui &count, ui *buffer,
                                 MatrixXd &datagraph_eigen, MatrixXd &querygraph_eigen, bool isEigenCheck, int top_s);


    static void computeCandidateWithLDF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                        ui &count, ui *buffer = NULL);

    static void generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                      VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                      ui *candidates_count, ui *flag, ui *updated_flag);

    static void pruneCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                   VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                   ui *candidates_count, ui *flag, ui *updated_flag);


    static void
    printCandidatesInfo(const Graph *query_graph, ui *candidates_count, std::vector<ui> &optimal_candidates_count);

    static void sortCandidates(ui** candidates, ui* candidates_count, ui num);

    static double computeCandidatesFalsePositiveRatio(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                                          ui *candidates_count, std::vector<ui> &optimal_candidates_count);

    static bool
    LDFWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
               TreeNode *&tree,
               std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
               bool isEigenCheck, int top_s);
    static bool
    NLFWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
               TreeNode *&tree,
               std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
               bool isEigenCheck, int top_s);

    static bool
    GQLWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
               TreeNode *&tree,
               std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
               bool isEigenCheck, int top_s);

    static bool
    TSOWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
               TreeNode *&tree,
               std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
               bool isEigenCheck, int top_s);

    static bool
    CFLWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
               TreeNode *&tree,
               std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
               bool isEigenCheck, int top_s);

    static bool
    DPisoWrapper(Graph *data_graph, Graph *query_graph, ui **&candidates, ui *&candidates_count, ui *&order,
                      TreeNode *&tree, std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                      bool isEigenCheck, int top_s);
private:
    static void allocateBuffer(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);
    static bool verifyExactTwigIso(const Graph *data_graph, const Graph *query_graph, ui data_vertex, ui query_vertex,
                                   bool **valid_candidates, int *left_to_right_offset, int *left_to_right_edges,
                                   int *left_to_right_match, int *right_to_left_match, int* match_visited,
                                   int* match_queue, int* match_previous);
    static void compactCandidates(ui** &candidates, ui* &candidates_count, ui query_vertex_num);
    static bool isCandidateSetValid(ui** &candidates, ui* &candidates_count, ui query_vertex_num);



};


#endif //SUBGRAPHMATCHING_FILTERVERTICES_H
