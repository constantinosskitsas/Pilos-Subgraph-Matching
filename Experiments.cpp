
#include "Experiments.h"
#include "IO.h"
#include "StudyPerformance.h"
#include "FilterVertices.h"

#include <iostream>
#include <numeric>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>
#include <sstream>

string Experiments::datagraphEigenMatrix;


// Interpret the ground truth store in the file;
vector<set<ui>> Experiments::ground_truth_interpreter(string path) {
    vector<set<ui>> ground_truth;
    ifstream file_in(path);
    if(!file_in){
        throw invalid_argument("ground truth file doesn't exist");
    }

    string line;
    while (getline(file_in,line)){
        istringstream ss(line);
        ground_truth.push_back(set<ui>());

        int x;
        while (ss>>x) ground_truth.back().insert(x);
    }

#ifdef PRINT_GROUND_TRUTH
    for (auto elem : ground_truth) {
        for(auto candidate : elem){
            cout << candidate << " ";
        }
        cout<<" "<<endl;
    }
#endif

    return ground_truth;
}


//Return true if candidate_true is a subset of candidate.
bool Experiments::candidate_set_correctness_check(vector<set<ui>> candidate, vector<set<ui>> candidate_true, ui query_size) {
    if(candidate.size()!=candidate_true.size()){
        return false;
    }
    for(int i =0; i< query_size;i++){
        if(!std::includes(candidate[i].begin(),candidate[i].end(),
                          candidate_true[i].begin(),candidate_true[i].end())){
            cout<<"Pruning process of vertex "<<i<<" contains false negative"<<endl;
            cout<<"candidate set"<<endl;
            for (ui const&cand: candidate[i]){
                std::cout << cand << ' ';
            }
            cout<<" "<<endl;
            cout<<"candidate set true"<<endl;
            for (ui const&cand: candidate_true[i]){
                std::cout << cand << ' ';
            }
            cout<<" "<<endl;
            return false;
        }
    }
    cout<<"The candidate set passed the correctness check, the ground truth is a subset of candidate set."<<endl;
    return true;
}

// This experiment computes the candidate number for each filters with and without EF.
string Experiments::experiment1(Graph *data_graph, Graph *query_graph){
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    ui* tso_order = NULL;
    TreeNode* tso_tree = NULL;
    ui* cfl_order = NULL;
    TreeNode* cfl_tree = NULL;
    ui* dpiso_order = NULL;
    TreeNode* dpiso_tree = NULL;
    TreeNode* ceci_tree = NULL;
    ui* ceci_order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;

    bool (*functptr[7])(Graph*, Graph* , ui**& , ui*&, ui*&,TreeNode *&,std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &,std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,bool,int);
    functptr[0] = FilterVertices::LDFWrapper;
    functptr[1] = FilterVertices::NLFWrapper;
    functptr[2] = FilterVertices::GQLWrapper;
    functptr[3] = FilterVertices::TSOWrapper;
    functptr[4] = FilterVertices::CFLWrapper;
    functptr[5] = FilterVertices::DPisoWrapper;
    functptr[6] = FilterVertices::CECIFilter;

    int sum=0;
    int top_s = 10;
    if(query_graph->getVerticesCount()==4) top_s = 4;
    if(query_graph->getVerticesCount()==8) top_s = 8;

    cout<<"top_s: "<< top_s<<endl;


    static int results[14];
    int counter = 0;

    for(int i =0; i<7;i++){
        functptr[i](data_graph, query_graph, candidates, candidates_count,ceci_order,ceci_tree,TE_Candidates,NTE_Candidates,true,top_s);
        results[counter] = accumulate(candidates_count, candidates_count+query_graph->getVerticesCount(), sum);
        counter++;
        functptr[i](data_graph, query_graph, candidates, candidates_count,ceci_order,ceci_tree,TE_Candidates,NTE_Candidates,false,top_s);
        results[counter] = accumulate(candidates_count, candidates_count+query_graph->getVerticesCount(), sum);
        counter++;
    }

    string results_string = "";
    for (int i=0; i<14; i++){
        results_string.append(std::to_string(results[i])+",");
    }
    results_string.pop_back();
    results_string.append("\n");

    return results_string;
}


//This experiment computes candidates for each filter and check if it will generate correct ground truth.
void Experiments::experiment2(string data_graph_path,string query_graph_path,string eigen){

    // -d ../../test/reallife_dataset/wordnet/data_graph/wordnet.graph -q ../../test/reallife_dataset/wordnet/query_graph/query_dense_12_2.graph -filter GQL -order GQL -engine LFTJ -num MAX -eigen 1 -tops 10
    //TODO Missing CECI
    string filters[6] = {"LDF","NLF","GQL","TSO","CFL","DPiso"};

    int candidate_counts[6];
    int results[6];
    int embeddings[6];


    for(int i = 0; i<6;i++){
        matching_algo_inputs inputs;
        inputs.dgraph_path = data_graph_path;
        inputs.qgraph_path = query_graph_path;
        inputs.filter = filters[2];
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
        inputs.eigen = eigen;
        matching_algo_outputs outputs = StudyPerformance::solveGraphQuery(inputs);
        candidate_counts[i] = outputs.candidate_count_sum;
        results[i]= outputs.enumOutput.candidate_true_count_sum;
        embeddings[i] = outputs.enumOutput.embedding_cnt;

    }

    cout<<"--------------------------"<<endl;
    cout<<"Candidates"<<endl;
    for(int i=0; i<6;i++){
        cout<<candidate_counts[i]<<",";
    }
    cout<<" "<<endl;
    cout<<"Candidates_true"<<endl;
    for(int i=0; i<6;i++){
        cout<<results[i]<<",";
    }
    cout<<" "<<endl;
    cout<<"Embedding count"<<endl;
    for(int i=0; i<6;i++){
        cout<<embeddings[i]<<",";
    }
    cout<<" "<<endl;
    cout<<"--------------------------"<<endl;
}

//This experiment compare the overall performance of each algorithm with and without eigen value enhanced filter.
matching_algo_outputs Experiments::experiment3(const string data_graph_path,const string query_graph_path,const string filter,const string eigen,ui* order_pointer) {

//    string filters[6] = {"LDF","NLF","GQL","TSO","CFL","DPiso"};

    matching_algo_inputs inputs;
    inputs.dgraph_path = data_graph_path;
    inputs.qgraph_path = query_graph_path;
    inputs.filter = filter;
    inputs.order = "GQL";
    if(inputs.filter!="CECI")
    inputs.engine = "LFTJ";
    else{
        inputs.engine = "CECI";
    }
    if(inputs.filter=="KFA"||inputs.filter=="KFB"){
        //inputs.order = inputs.filter;
        inputs.order = "GQL";
        //inputs.filter="KFE";
        if (stoi(eigen))
            inputs.engine = "LFTJA";
    }else if (inputs.filter=="KFC"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
    }
    else if (inputs.filter=="KFR"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
    }
    else if (inputs.filter=="KFE"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
    }
    else if (inputs.filter=="KFD"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
    }
    else if (inputs.filter=="KFL"){
        inputs.filter="KFE";
        inputs.order = "KF";
        inputs.engine = "LFTJ";    
    }    else if (inputs.filter=="KFLL"){
        inputs.filter="KFE";
        inputs.order = "KF";
        inputs.engine = "LFTJA";    
    }
    else if (inputs.filter=="KFEL"){
        inputs.filter="KFE";
        inputs.order = "GQL";
        inputs.engine = "LFTJA";    
    }
    else if (inputs.filter=="KF"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";
    }else if (inputs.filter=="EOL"){
        //inputs.filter = "KF";
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJA";

    }else if (inputs.filter=="EOH"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJB";
        
    }else if (inputs.filter=="EOCR"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJR";        
    }
    else if (inputs.filter=="EOCR1"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJR";        
    }
    else if (inputs.filter=="EOCR2"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJR";        
    }
    else if (inputs.filter=="EOCR3"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJR";        
    }
    else if (inputs.filter=="EOCR4"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJR";        
    }
    else if (inputs.filter=="EO"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }
    else if (inputs.filter=="EOF"){
        inputs.filter = "KFE";
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }else if (inputs.filter=="EOF1"){
        inputs.filter = "KFR";
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }
    else if (inputs.filter=="KFR2"){
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }
    else if (inputs.filter=="EOF3"){
        inputs.filter = "KF";
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }
    else if (inputs.filter=="EOFF"){
        inputs.filter = "KFD";
        inputs.order = "GQL";
        inputs.engine = "LFTJ";     
    }

    
    
    inputs.eigen = eigen;
    inputs.order_pointer = order_pointer;


    matching_algo_outputs outputs = StudyPerformance::solveGraphQuery(inputs);

    if(stoi(eigen)){
        cout<<"Matches: "<<outputs.enumOutput.embedding_cnt<<filter<<"call_count "<<outputs.call_count <<" With eigen filter total time "<<outputs.total_time<<" enumeration time: "
        <<outputs.enumeration_time<<" preprocessing time: "<<outputs.preprocessing_time<<" candidate_sum: "<<outputs.candidate_count_sum<<" matching order: ";
        for (int i=0; i<outputs.query_size;i++){
            cout<<outputs.matching_order[i]<<" ";
        }
        cout<<" "<<endl;
    } else{
        cout<<"Matches: "<<outputs.enumOutput.embedding_cnt<<filter<<"call_count "<<outputs.call_count<<" No eigen filter total time "<<outputs.total_time<<" enumeration time: "
        <<outputs.enumeration_time<<" preprocessing time: "<<outputs.preprocessing_time<<" candidate_sum: "<<outputs.candidate_count_sum<<" matching order: ";
        for (int i=0; i<outputs.query_size;i++){
            cout<<outputs.matching_order[i]<<" ";
        }
        cout<<" "<<endl;
    }

    return outputs;

}

//This experiment will generate the ground truth of each query and store them in a file.
void Experiments::experiment4(const std::string eigen, queryMeta meta) {
    matching_algo_inputs inputs;
    inputs.dgraph_path = meta.data_graph_path;
    inputs.qgraph_path = meta.query_path;
    inputs.filter = "GQL";
    inputs.order = "GQL";
    inputs.engine = "LFTJ";
    inputs.eigen = eigen;

    matching_algo_outputs outputs = StudyPerformance::solveGraphQuery(inputs);

    fstream file;

    std::ostringstream file_name;
    file_name << "ground_truth/"<<meta.dataset<<"/"<<meta.dataset<<"_"<<meta.query_property<<"_"<<meta.query_size<<"_"<<meta.query_number<<".txt";
    std::string file_path = file_name.str();

    file.open(file_path,ios_base::out);
    for (int i =0; i< outputs.query_size; i++){
        for (set<ui>::iterator it=outputs.enumOutput.candidate_true[i].begin(); it!=outputs.enumOutput.candidate_true[i].end(); ++it){
            file<<*it<<" ";
        }
        file<<"\n";
    }
    file.close();

}








// Compute the ground truth for each query.
//void Experiments::experiment4(const string data_graph_path, const string query_graph_path, const string eigen) {
//}


//int main(int argc, char** argv){
//    //    MatrixXd datagraph_eigenvalue(data_graph->getVerticesCount(), 35);
//    //    MTcalc12(data_graph,data_graph->getGraphMaxDegree(),datagraph_eigenvalue,true,35);
//    //    saveData("wordnet.csv", datagraph_eigenvalue);
//
//    //  -d ../../test/reallife_dataset/wordnet/data_graph/wordnet.graph -q ../../test/reallife_dataset/wordnet/query_graph/query_dense_8_1.graph -filter GQL -order GQL -engine LFTJ -num MAX -eigen 1 -tops 4
//    int query_number[200];
//    for (int i=0;i<200;i++){
//        query_number[i] = i+1;
//    }
//    //Yeast
//    int query_vertices_dense[5]={4,8,16,24,32};
//    int query_vertices_sparse[4]={8,16,24,32};
//
//    //Wordnet
//    int query_vertices_dense_v2[5]={4,8,12,16,20};
//    int query_vertices_sparse_v2[4]={8,12,16,20};
//
//
//    string querypath[1800];
//    int counter = 0;
//
//    for(int i=0; i<5;i++){
//        for (int j=0;j<200;j++){
//            std::ostringstream oss;
//            oss << "../../test/reallife_dataset/wordnet/query_graph/query_dense_" << query_vertices_dense_v2[i] <<"_" << query_number[j]<< ".graph";
//            std::string var = oss.str();
//            querypath[counter] = var;
//            counter++;
//        }
//    }
//
//    for(int i=0; i<4;i++){
//        for (int j=0;j<200;j++){
//            std::ostringstream oss;
//            oss << "../../test/reallife_dataset/wordnet/query_graph/query_dense_" << query_vertices_sparse_v2[i] <<"_" << query_number[j]<< ".graph";
//            std::string var = oss.str();
//            querypath[counter] = var;
//            counter++;
//        }
//    }
//
//    std::string input_data_graph_file = "../../test/reallife_dataset/wordnet/data_graph/wordnet.graph";
//    Graph* data_graph = new Graph(true);
//    data_graph->loadGraphFromFile(input_data_graph_file);
//
//    std::ofstream myfile;
//    myfile.open ("wordnetexperiment.csv",std::ios_base::app);
//    myfile << "Query_NUmber,LDF+EF,LDF,NLF+EF,NLF,GQL+EF,GQL,TSOF+EF,TSOF,CFL+EF,CFL,DPiso+EF,Dpiso,CECIF+EF,CECIF\n";
//    myfile.close();
//
//    for (int i=771; i<1800;i++){
//        Graph* query_graph = new Graph(true);
//        query_graph->loadGraphFromFile(querypath[i]);
//        query_graph->buildCoreTable();
//        string candidate_true = to_string(i)+",";
//        string results2 = candidate_true.append(experiment(data_graph,query_graph));
//        myfile.open ("wordnetexperiment.csv",std::ios_base::app);
//        myfile<<results2;
//        myfile.close();
//        cout << i << endl;
//    }
//
//
//
//
//
//}
