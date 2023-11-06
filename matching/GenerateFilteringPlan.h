
#ifndef SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H
#define SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H

#include "graph/graph.h"
#include "configuration/types.h"

class GenerateFilteringPlan
{
public:
    static void generateTSOFilterPlan(Graph *data_graph, Graph *query_graph, TreeNode *&tree,
                                      VertexID *&order, int top_s);
    static void generateCFLFilterPlan(Graph *data_graph, Graph *query_graph, TreeNode *&tree,
                                      VertexID *&order, int &level_count, ui *&level_offset, bool isEigenCheck, int top_s);
    static void generateDPisoFilterPlan(Graph *data_graph, Graph *query_graph, TreeNode *&tree,
                                        VertexID *&order);
    static void generateCECIFilterPlan(Graph *data_graph, Graph *query_graph, TreeNode *&tree,
                                       VertexID *&order);

private:
    static VertexID selectTSOFilterStartVertex(Graph *data_graph, Graph *query_graph, int top_s);
    static VertexID selectCFLFilterStartVertex(Graph *data_graph, Graph *query_graph, bool isEigenCheck, int top_s);
    static VertexID selectDPisoStartVertex(Graph *data_graph, Graph *query_graph);
    static VertexID selectCECIStartVertex(Graph *data_graph, Graph *query_graph);
};

#endif // SUBGRAPHMATCHING_GENERATEFILTERINGPLAN_H
