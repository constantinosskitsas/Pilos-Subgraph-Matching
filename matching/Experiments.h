#include <string>
#include <iostream>
#include "graph/graph.h"
#include <vector>
#include <set>
#include "StudyPerformance.h"
// #include "EvaluateQuery.h"

using namespace std;
#ifndef SUBGRAPHMATCHING_EXPERIMENTS_H
#define SUBGRAPHMATCHING_EXPERIMENTS_H

struct queryMeta
{
    std::string dataset;
    std::string query_property;
    std::string query_size;
    string query_number;
    string query_path;
    string data_graph_path;
};

class Experiments
{
public:
    static string datagraphEigenMatrix;
    static vector<std::set<ui>> ground_truth_interpreter(string path);
    static bool candidate_set_correctness_check(vector<set<ui>> candidate, vector<set<ui>> candidate_true, ui query_size);
    static string experiment1(Graph *data_graph, Graph *query_graph);
    static void experiment2(string data_graph, string query_graph, string eigen);
    static matching_algo_outputs experiment3(const string data_graph_path, const std::string query_graph_path, const std::string filter, const std::string eigen, ui *order_pointer, string alpha, string beta, string thnum, string embeddingcount);
    static void experiment4(const string eigen, queryMeta meta);
};

#endif // SUBGRAPHMATCHING_EXPERIMENTS_H
