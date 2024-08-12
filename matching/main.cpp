#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "FilterVertices.h"
#include "IO.h"
#include "eigenHelper.h"
#include "Experiments.h"
#include "StudyPerformance.h"
#include "GenerateQuery.h"

pair<matching_algo_outputs, matching_algo_outputs> fakeMatchingWrapper(queryMeta meta, string filter)
{
    matching_algo_outputs original = Experiments::experiment3(meta.data_graph_path, meta.query_path, filter, "0", NULL, "alpha", "beta", "thnum", "10000");
    ui *fake_pointer = new ui[stoi(meta.query_size)];
    for (int i = 0; i < stoi(meta.query_size); i++)
    {
        ui order = original.matching_order[i];
        *&fake_pointer[i] = order;
    }
    matching_algo_outputs enhanced = Experiments::experiment3(meta.data_graph_path, meta.query_path, filter, "1", fake_pointer, "alpha", "beta", "thnum", "10000");
    delete[] fake_pointer;
    return pair(original, enhanced);
}

pair<matching_algo_outputs, matching_algo_outputs> MatchingWrapper(string datagraph, string querygraph, string filter)
{
    matching_algo_outputs original = Experiments::experiment3(datagraph, querygraph, filter, "0", NULL, "alpha", "beta", "thnum", "10000");
    matching_algo_outputs enhanced = Experiments::experiment3(datagraph, querygraph, filter, "1", NULL, "alpha", "beta", "thnum", "10000");
    return pair(original, enhanced);
}
pair<matching_algo_outputs, matching_algo_outputs> MatchingWrapperN(string datagraph, string querygraph, string filter)
{
    matching_algo_outputs original = Experiments::experiment3(datagraph, querygraph, filter, "0", NULL, "alpha", "beta", "thnum", "10000");
    matching_algo_outputs enhanced = Experiments::experiment3(datagraph, querygraph, filter, "1", NULL, "alpha", "beta", "thnum", "10000");
    return pair(original, enhanced);
}

void exact_eval(string dataset, string querysize, string querynumber, string property)
{
    Experiments::datagraphEigenMatrix = dataset + ".csv";
    queryMeta meta;
    meta.dataset = dataset;
    meta.query_property = property;
    meta.query_size = querysize;
    meta.query_number = querynumber;

    std::ostringstream data;
    data << "../../test/reallife_dataset/" << meta.dataset << "/data_graph/" << meta.dataset << ".graph";
    meta.data_graph_path = data.str();

    std::ostringstream query;
    //    query << "../../test/reallife_dataset/" << meta.dataset << "/query_graph/query_" << meta.query_property << "_"
    //          << meta.query_size << "_" << meta.query_number << ".graph";
    query << "../../test/large_query/test.graph";
    meta.query_path = query.str();

    pair<matching_algo_outputs, matching_algo_outputs> LDF = fakeMatchingWrapper(meta, "LDF");
    pair<matching_algo_outputs, matching_algo_outputs> NLF = fakeMatchingWrapper(meta, "NLF");
    pair<matching_algo_outputs, matching_algo_outputs> GQL = fakeMatchingWrapper(meta, "GQL");
    pair<matching_algo_outputs, matching_algo_outputs> TSOF = fakeMatchingWrapper(meta, "TSO");
    pair<matching_algo_outputs, matching_algo_outputs> CFL = fakeMatchingWrapper(meta, "CFL");
    pair<matching_algo_outputs, matching_algo_outputs> DPiso = fakeMatchingWrapper(meta, "DPiso");
    matching_algo_outputs KF = Experiments::experiment3(meta.data_graph_path, meta.query_path, "KF", "0", NULL, "25", "500", "5", "10000");

    vector<pair<matching_algo_outputs, matching_algo_outputs>> evaluations;
    evaluations.push_back(LDF);
    evaluations.push_back(NLF);
    evaluations.push_back(GQL);
    evaluations.push_back(TSOF);
    evaluations.push_back(CFL);
    evaluations.push_back(DPiso);

    std::ostringstream oss;
    oss << meta.query_property << "_" << meta.query_size << "_" << meta.query_number;

    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.call_count << "," << eval.second.call_count;
    }
    oss << "," << KF.call_count
        << "," << LDF.first.enumOutput.embedding_cnt;

    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.total_time << "," << eval.second.total_time;
    }
    oss << "," << KF.total_time << "," << LDF.first.enumOutput.embedding_cnt;

    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.candidate_count_sum << "," << eval.second.candidate_count_sum;
    }
    oss << "," << KF.candidate_count_sum;
    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.matching_order_string << "," << eval.second.matching_order_string;
    }
    oss << "," << KF.matching_order_string;
    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.preprocessing_time << "," << eval.second.preprocessing_time;
    }
    oss << "," << KF.preprocessing_time;
    for (auto &eval : evaluations)
    {
        oss << "," << eval.first.enumeration_time << "," << eval.second.enumeration_time;
    }
    oss << "," << KF.enumeration_time;

    std::string var = oss.str();

    cout << var << endl;

    string file_path = "";
    if (property == "sparse")
    {
        file_path = "performance_experiment/" + dataset + "_" + "s" + querysize + ".csv";
    }
    if (property == "dense")
    {
        file_path = "performance_experiment/" + dataset + "_" + "d" + querysize + ".csv";
    }

    std::ofstream myfile;
    myfile.open(file_path, std::ios_base::app);
    myfile << var << "\n";
    myfile.close();
}

void generate_datagraph_eigenvector(string data_graph_path, string csvfilename, int size)
{
    Graph *data_graph = new Graph(true);
    data_graph->loadGraphFromFile(data_graph_path);

    MatrixXd datagraph_eigenvalue(data_graph->getVerticesCount(), size);
    cout << "Start compute eigen value" << endl;

    MTcalc12A(data_graph, data_graph->getGraphMaxDegree(), datagraph_eigenvalue, true, size, 250);
    saveData(csvfilename, datagraph_eigenvalue);
    // saveData("test1", datagraph_eigenvalue);
}

void fixed_order_experiment(int argc, char **argv)
{

    //    generate_datagraph_eigenvector("../../test/reallife_dataset/eu2005/data_graph/eu2005.graph","eu2005.csv");

    // yeast,hprd,uspatents,youtube,dblp,eu2005,
    vector<string> query_size_dense_1 = {"4", "8", "16", "24", "32"};
    vector<string> query_size_sparse_1 = {"8", "16", "24", "32"};

    // human,wordnet
    vector<string> query_size_dense_2 = {"4", "8", "12", "16", "20"};
    vector<string> query_size_sparse_2 = {"8", "12", "16", "20"};

    MatchingCommand command(argc, argv);
    string dataset_name = command.getDatasetName();
    string query_size = command.getQuerySize();
    string query_number = command.getQueryNumber();
    string query_property = command.getQueryProperty();

    cout << dataset_name << endl;
    cout << query_size << endl;
    cout << query_number << endl;
    cout << query_property << endl;

    exact_eval(dataset_name, query_size, query_number, query_property);
}

int main(int argc, char **argv)
{
    MatchingCommand command(argc, argv);
    string dataset_name = command.getDatasetName();
    string query_size = command.getQuerySize();
    string query_number = command.getQueryNumber();
    string query_property = command.getQueryProperty();
    string query_filter = command.getFilterType();
    string alpha = command.getalpha();
    string beta = command.getbeta();
    string thnum = command.getThreadCount();
    string embeddingcount = command.getMaximumEmbeddingNum1();
    string datagraph = "../../dataset/" + dataset_name + "/data_graph/" + dataset_name + ".graph";
    auto start = std::chrono::high_resolution_clock::now();
    Graph *data_graph = new Graph(true);
    auto end = std::chrono::high_resolution_clock::now();
    double PT = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    string StoreFile = command.getStoreFile();
    ui aa[5] = {32, 64, 96, 128, 256};
    srand(time(NULL));
    if (query_filter == "GQ")
    {
        bool ck = false;
        int qn = stoi(query_number);
        for (int da = 0; da <= 4; da++)
            for (int di = 1; di <= qn; di++)
            {
                ui kk = aa[da];
                while (!GenerateQueryD(kk, data_graph, 1, di, dataset_name))
                    ;
                cout << "da " << kk << " di " << di << endl;
                ck = false;
                ui cos = 0;
            }
        return 0;
    }
    if (query_filter == "EC")
    {
        int numeig = stoi(query_number);
        generate_datagraph_eigenvector("../../dataset/" + dataset_name + "/data_graph/" + dataset_name + ".graph", dataset_name + ".csv", numeig);
        return 0;
    }
    Experiments::datagraphEigenMatrix = "../../" + dataset_name + ".csv";
    // string datagraph = "../../dataset/wordnet/data_graph/wordnet.graph";
    datagraph = "../../dataset/" + dataset_name + "/data_graph/" + dataset_name + ".graph";
    string querygraph = "../../dataset/" + dataset_name + "/query_graph/query_" + query_property + "_" + query_size + "_" + query_number + ".graph";
    // pair <matching_algo_outputs,matching_algo_outputs> KFE = MatchingWrapper(datagraph,querygraph,"KFE");
    matching_algo_outputs KF = Experiments::experiment3(datagraph, querygraph, query_filter, "0", NULL, alpha, beta, thnum, embeddingcount);
    vector<pair<matching_algo_outputs, matching_algo_outputs>> evaluations;
    std::ostringstream oss;
    int row = 0;
    oss << query_number << " " << query_size << " " << KF.call_count << " " << KF.enumOutput.embedding_cnt;
    oss << " " << KF.total_time << " " << KF.candidate_count_sum; //<<LDF.first.enumOutput.embedding_cnt;
    oss << " " << KF.preprocessing_time;
    oss << " " << KF.enumeration_time;
    std::string var = oss.str();
    cout << var << endl;
    string file_path = "";
    file_path = "performance_experiment/" + StoreFile + "_" + query_filter + "_" + embeddingcount + "_" + dataset_name + "_" + query_property + query_size + ".csv";
    std::ofstream myfile;
    myfile.open(file_path, std::ios_base::app);
    myfile << var << "\n";
    myfile.close();
    return 0;
}
