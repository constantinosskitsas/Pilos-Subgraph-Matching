

#include <functional>
#include <map>
#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <numeric>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <time.h>

#include "matchingcommand.h"
#include "graph/graph.h"
#include "GenerateFilteringPlan.h"
#include "FilterVertices.h"
#include "BuildTable.h"
#include "GenerateQueryPlan.h"
#include "EvaluateQuery.h"
#include "IO.h"
#include "eigenHelper.h"
#include "Experiments.h"
#include "StudyPerformance.h"
#include "KF/spectra.h"



#define NANOSECTOSEC(elapsed_time) ((elapsed_time)/(double)1000000000)
#define BYTESTOMB(memory_cost) ((memory_cost)/(double)(1024 * 1024))
//#define PRINT;
//#define ONLYCOUNTS;



//#define PRINT1 1
size_t StudyPerformance::enumerate(Graph* data_graph, Graph* query_graph, Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                ui* matching_order, size_t output_limit) {
    static ui order_id = 0;

    order_id += 1;

    auto start = std::chrono::high_resolution_clock::now();



    size_t call_count = 0;
    size_t embedding_count = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                               matching_order, output_limit, call_count).embedding_cnt;


    auto end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
#ifdef SPECTRUM
    if (EvaluateQuery::exit_) {
        printf("Spectrum Order %u status: Timeout\n", order_id);
    }
    else {
        printf("Spectrum Order %u status: Complete\n", order_id);
    }
#endif
    printf("Spectrum Order %u Enumerate time (seconds): %.4lf\n", order_id, NANOSECTOSEC(enumeration_time_in_ns));
    printf("Spectrum Order %u #Embeddings: %zu\n", order_id, embedding_count);
    printf("Spectrum Order %u Call Count: %zu\n", order_id, call_count);
    printf("Spectrum Order %u Per Call Count Time (nanoseconds): %.4lf\n", order_id, enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));


    return embedding_count;
}

void StudyPerformance::spectrum_analysis(Graph* data_graph, Graph* query_graph, Edges*** edge_matrix, ui** candidates, ui* candidates_count,
                       size_t output_limit, std::vector<std::vector<ui>>& spectrum, size_t time_limit_in_sec) {

    for (auto& order : spectrum) {
        std::cout << "----------------------------" << std::endl;
        ui* matching_order = order.data();
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);

        std::future<size_t> future = std::async(std::launch::async, [data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                                                     matching_order, output_limit](){
            return enumerate(data_graph, query_graph, edge_matrix, candidates, candidates_count, matching_order, output_limit);
        });

        std::cout << "execute...\n";
        std::future_status status;
        do {
            status = future.wait_for(std::chrono::seconds(time_limit_in_sec));
            if (status == std::future_status::deferred) {
                std::cout << "Deferred\n";
                exit(-1);
            } else if (status == std::future_status::timeout) {
#ifdef SPECTRUM
                EvaluateQuery::exit_ = true;
#endif
            }
        } while (status != std::future_status::ready);
    }
}


matching_algo_outputs StudyPerformance::solveGraphQuery(matching_algo_inputs inputs){
    int argc = 0;
    char** argv;
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = inputs.qgraph_path;
    std::string input_data_graph_file = inputs.dgraph_path;
    std::string input_filter_type = inputs.filter;
    std::string input_order_type = inputs.order;
    std::string input_engine_type = inputs.engine;
    std::string input_max_embedding_num = "1000";
    std::string input_time_limit = "300";
    std::string input_order_num = command.getOrderNum();
    std::string input_distribution_file_path = command.getDistributionFilePath();
    std::string input_csr_file_path = command.getCSRFilePath();
    std::string input_iseigen = inputs.eigen;
    std::string input_tops;
    string alpha=inputs.alpha;
    string beta=inputs.beta;
    string thnum=inputs.thnum;
    string embdcount=inputs.embcount;
    
    matching_algo_outputs outputs;
    /**
     * Output the command line information.
     */
#ifdef PRINT1
    std::cout << "Command Line:" << std::endl;
    std::cout << "\tData Graph CSR: " << input_csr_file_path << std::endl;
    std::cout << "\tData Graph: " << input_data_graph_file << std::endl;
    std::cout << "\tQuery Graph: " << input_query_graph_file << std::endl;
    std::cout << "\tFilter Type: " << input_filter_type << std::endl;
    std::cout << "\tOrder Type: " << input_order_type << std::endl;
    std::cout << "\tEngine Type: " << input_engine_type << std::endl;
    std::cout << "\tOutput Limit: " << input_max_embedding_num << std::endl;
    std::cout << "\tTime Limit (seconds): " << input_time_limit << std::endl;
    std::cout << "\tOrder Num: " << input_order_num << std::endl;
    std::cout << "\tDistribution File Path: " << input_distribution_file_path << std::endl;
    std::cout << "\tWith eigen filter?: " << input_iseigen << std::endl;
    std::cout << "\tTop-s eigen value?: " << input_tops << std::endl;

    std::cout << "--------------------------------------------------------------------" << std::endl;
#endif
    /**
     * Load input graphs.
     */
#ifdef PRINT1
    std::cout << "Load graphs..." << std::endl;
#endif


    auto start = std::chrono::high_resolution_clock::now();

    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);
    query_graph->buildCoreTable();
    outputs.query_size = query_graph->getVerticesCount();
    
    Graph* data_graph = new Graph(true);
    float **eigenVD1 = NULL;
    ui dsiz=0;
    if (input_csr_file_path.empty()) {
        if(inputs.filter=="KF"||inputs.filter=="KFE"||inputs.filter=="KFD"||inputs.filter=="PLMT"||inputs.filter=="PL"){
            data_graph->loadGraphFromFile(input_data_graph_file);
            data_graph->BuildLabelOffset();
            dsiz= data_graph->getVerticesCount();
            
            eigenVD1 = new float *[dsiz];

            for (ui i = 0; i < dsiz; ++i)
                {
                eigenVD1[i] = new float[35];
            }
    
            openData1(Experiments::datagraphEigenMatrix, eigenVD1);
        }
        else{

            data_graph->loadGraphFromFile(input_data_graph_file);
        }

    }
    else {
        std::string degree_file_path = input_csr_file_path + "_deg.bin";
        std::string edge_file_path = input_csr_file_path + "_adj.bin";
        std::string label_file_path = input_csr_file_path + "_label.bin";
        data_graph->loadGraphFromFileCompressed(degree_file_path, edge_file_path, label_file_path);
    }

    auto end = std::chrono::high_resolution_clock::now();

    double load_graphs_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    
    
    input_tops = "10";
    if(query_graph->getVerticesCount()==4) input_tops = "4";
    if(query_graph->getVerticesCount()==8) input_tops = "8";


#ifdef PRINT1
    std::cout << "-----" << std::endl;
    std::cout << "Query Graph Meta Information" << std::endl;
    query_graph->printGraphMetaData();
    std::cout << "-----" << std::endl;
    data_graph->printGraphMetaData();

    std::cout << "--------------------------------------------------------------------" << std::endl;
#endif
    /**
     * Start queries.
     */


#ifdef PRINT1
    std::cout << "Start queries..." << std::endl;
    std::cout << "-----" << std::endl;
    std::cout << "Filter candidates..." << std::endl;
#endif

    start = std::chrono::high_resolution_clock::now();

    bool isEigenCheck;
    int top_s = std::stoi( input_tops);
    istringstream(input_iseigen) >> isEigenCheck;
    //top_s=10;


    ui** candidates = NULL;
    float **EWeight=NULL;
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
    if (input_filter_type == "LDF") {
        FilterVertices::LDFFilter(data_graph, query_graph, candidates, candidates_count,isEigenCheck,top_s);
    } else if(input_filter_type=="PL"){
        int alpha1=stoi(alpha);
        int beta1=stoi(beta);
        SpectralMatching(query_graph->getVerticesCount(), data_graph, input_query_graph_file, 0,candidates,candidates_count,EWeight,eigenVD1,alpha1,beta1);
    }
    else if(input_filter_type=="PLMT"){
        int alpha1=stoi(alpha);
        int thnum1=stoi(thnum);
        int beta1=stoi(beta);
        SpectralMatchingMT(query_graph->getVerticesCount(), data_graph, input_query_graph_file, 0,candidates,candidates_count,EWeight,eigenVD1,alpha1,thnum1,beta1);
    }
    else if (input_filter_type == "NLF") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count,false,top_s);
    }else if (input_filter_type == "NLFE") {
        FilterVertices::NLFFilter(data_graph, query_graph, candidates, candidates_count,true,top_s);
    } else if (input_filter_type == "GQL") {
        FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count,isEigenCheck,top_s);
    } else if (input_filter_type == "TSO") {
        FilterVertices::TSOFilter(data_graph, query_graph, candidates, candidates_count, tso_order, tso_tree,isEigenCheck,top_s);
    } else if (input_filter_type == "CFL") {
        
        FilterVertices::CFLFilter(data_graph, query_graph, candidates, candidates_count, cfl_order, cfl_tree,isEigenCheck,top_s);
        
     }else if (input_filter_type == "DPiso") {

        FilterVertices::DPisoFilter(data_graph, query_graph, candidates, candidates_count, dpiso_order, dpiso_tree,false,top_s);

    } else if (input_filter_type == "CECI") {

        FilterVertices::CECIFilter(data_graph, query_graph, candidates, candidates_count, ceci_order, ceci_tree, TE_Candidates, NTE_Candidates,isEigenCheck,top_s);

    }  else {
        std::cout << "The specified filter type '" << input_filter_type << "' is not supported." << std::endl;
        exit(-1);
    }

    // Sort the candidates to support the set intersections
    // TODO figure out why CECI doesn't work, read the paper.
    if (input_filter_type != "CECI")
        FilterVertices::sortCandidates(candidates, candidates_count, query_graph->getVerticesCount());

    end = std::chrono::high_resolution_clock::now();
    double filter_vertices_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    int sum = 0;
    outputs.candidate_count_sum = accumulate(candidates_count, candidates_count + query_graph->getVerticesCount(), sum);

#ifdef ONLYCOUNTS
    return outputs;
#endif

    for (int i=0; i<query_graph->getVerticesCount();i++){
        outputs.candidate.push_back(set<ui>());
    }

    for (int i=0; i<query_graph->getVerticesCount();i++){
        for(int j=0;j<candidates_count[i];j++){
            outputs.candidate[i].insert(candidates[i][j]);
        }
    }


    // Compute the candidates false positive ratio.
#ifdef OPTIMAL_CANDIDATES
    std::vector<ui> optimal_candidates_count;
    double avg_false_positive_ratio = FilterVertices::computeCandidatesFalsePositiveRatio(data_graph, query_graph, candidates,
                                                                                          candidates_count, optimal_candidates_count);
    FilterVertices::printCandidatesInfo(query_graph, candidates_count, optimal_candidates_count);
#endif

#ifdef PRINT1
    std::cout << "-----" << std::endl;
    std::cout << "Build indices..." << std::endl;
#endif

    start = std::chrono::high_resolution_clock::now();

    Edges ***edge_matrix = NULL;
    if (input_filter_type != "CECI") {
        edge_matrix = new Edges **[query_graph->getVerticesCount()];
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            edge_matrix[i] = new Edges *[query_graph->getVerticesCount()];
        }
        BuildTable::buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    }

    end = std::chrono::high_resolution_clock::now();
    double build_table_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    size_t memory_cost_in_bytes = 0;
    if (input_filter_type != "CECI") {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, edge_matrix);
#ifdef PRINT1
        BuildTable::printTableCardinality(query_graph, edge_matrix);
#endif
    }
    else {
        memory_cost_in_bytes = BuildTable::computeMemoryCostInBytes(query_graph, candidates_count, ceci_order, ceci_tree,
                                                                    TE_Candidates, NTE_Candidates);
#ifdef PRINT1
        BuildTable::printTableCardinality(query_graph, ceci_tree, ceci_order, TE_Candidates, NTE_Candidates);
#endif
    }

#ifdef PRINT1
    std::cout << "-----" << std::endl;
    std::cout << "Generate a matching order..." << std::endl;
#endif

    start = std::chrono::high_resolution_clock::now();

    ui* matching_order = NULL;
    ui* pivots = NULL;
    ui** weight_array = NULL;

    size_t order_num = 0;
    sscanf(input_order_num.c_str(), "%zu", &order_num);

    std::vector<std::vector<ui>> spectrum;
    if (input_order_type == "QSI") {
        GenerateQueryPlan::generateQSIQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots);
    }
    
    else if (input_order_type=="KF"){
        ui qsiz=query_graph->getVerticesCount();
        float ED[qsiz];
        ui TOO=0;
        //ui count1=0;
        ui count2=0;
        ui count1[qsiz];
        float maxED=0;
        float maxEDW=0;
        ui maxCan=0;
        for (int d=0;d<qsiz;d++){
            ED[d]=0;
            count1[d]=0;
            for( int k=0;k<candidates_count[d];k++)
                if (EWeight[d][k]==100000)
                    count1[d]++;
                else{
                    ED[d]=ED[d]+EWeight[d][k];
                    if(EWeight[d][k]>maxED)
                        maxED=EWeight[d][k];
                }             
        }
        int aa=0.5;
        if(maxED==0)
        maxED=1;
        if(input_filter_type=="KF")
            aa=1;
        for (int d=0;d<qsiz;d++){
            ED[d]=ED[d]/maxED+aa*count1[d]+candidates_count[d];
            //cout<<ED[d]<<endl;
        }
        GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots);
    }
    else if (input_order_type=="KFX"){
       // else if (input_filter_type=="KF"){
        ui qsiz=query_graph->getVerticesCount();
        float ED[qsiz];
        ui TOO=0;
        ui count1=0;
        ui count2=0;
        float maxED=0;
        float maxEDW=0;
        ui maxCan=0;
        for (int d=0;d<qsiz;d++){
            for( int k=0;k<candidates_count[d];k++)
                if (EWeight[d][k]==10000)
                    count1++;
                else{
                    ED[d]=ED[d]+EWeight[d][k];
                    if(EWeight[d][k]>maxED)
                        maxED=EWeight[d][k];
                }
            ED[d]=ED[d]+count1*maxED;
            if(ED[d]>maxEDW)
                maxEDW=ED[d];
            if(maxCan<candidates_count[d])
                maxCan=candidates_count[d];
            maxED=0;
            count1=0;                
        }
        //TOO =0 -> 0.9 *Gives high value to size-> Number of Candidates/Max Candidates
                    //0.1 *Gives low value to EW-> EW/MEW
                    //1-p2 means low weight gets high value(looking for min)
                    //Prioritizes min candidates with max weight
        if (input_order_type=="KFA"){
        for (int i = 0; i < qsiz; i++){
            float p1=candidates_count[i]/maxCan;
            float p2=ED[i]/maxEDW;
            p1=p1*0.9;
            p2= (1-p2)*(0.1);
            ED[i]=p1+(p2);
    }
    GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots);
    }//TOO =1 -> 0.9 *Gives high value to size-> Number of Candidates/Max Candidates
                    //0.1 *Gives low value to EW-> EW/MEW
                    //p2 means low weight gets low value(looking for min)
                    //Prioritizes min candidates with min weight
     else if (input_order_type=="KFB"){
        for (int i = 0; i < qsiz; i++){
            float p1=candidates_count[i]/maxCan;
            float p2=ED[i]/maxEDW;
            p1=p1*0.9;
            p2= (p2)*(0.1);
            ED[i]=p1+(p2);
    }
    GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots);
    }
    //TOO=1 NEW-> EigenWeight/Num of Candidates
    else if (input_order_type=="KFC"){
        cout<<"malista"<<endl;
       for (int i = 0; i < qsiz; i++){
            ED[i]=ED[i]/candidates_count[i] ; 
            cout<<ED[i]<<endl;       
    }
    GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots); 
    }
        
    else if(input_order_type=="KFD")
GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots);
    else
        GenerateQueryPlan::generateGQLQueryPlanN(data_graph, query_graph,ED, matching_order, pivots);
    //}*/
    
    }
   
     else if (input_order_type == "GQL") {

        if (inputs.order_pointer == NULL){
            
            GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);

        }
        else {
            GenerateQueryPlan::GQLorderfake(data_graph, query_graph, candidates_count, inputs.order_pointer, pivots);
        }



    } else if (input_order_type == "TSO") {
        if (tso_tree == NULL) {
            GenerateFilteringPlan::generateTSOFilterPlan(data_graph, query_graph, tso_tree, tso_order,top_s);
        }
        GenerateQueryPlan::generateTSOQueryPlan(query_graph, edge_matrix, matching_order, pivots, tso_tree, tso_order);
    } else if (input_order_type == "CFL") {
        if (cfl_tree == NULL) {
            int level_count;
            ui* level_offset;
            GenerateFilteringPlan::generateCFLFilterPlan(data_graph, query_graph, cfl_tree, cfl_order, level_count, level_offset,isEigenCheck,top_s);
            delete[] level_offset;
        }
        GenerateQueryPlan::generateCFLQueryPlan(data_graph, query_graph, edge_matrix, matching_order, pivots, cfl_tree, cfl_order, candidates_count);
    } else if (input_order_type == "DPiso") {
        if (dpiso_tree == NULL) {
            GenerateFilteringPlan::generateDPisoFilterPlan(data_graph, query_graph, dpiso_tree, dpiso_order);
        }

        GenerateQueryPlan::generateDSPisoQueryPlan(query_graph, edge_matrix, matching_order, pivots, dpiso_tree, dpiso_order,
                                                   candidates_count, weight_array);
    }
    else if (input_order_type == "CECI") {
        GenerateQueryPlan::generateCECIQueryPlan(query_graph, ceci_tree, ceci_order, matching_order, pivots);
    }
    else if (input_order_type == "RI") {
        GenerateQueryPlan::generateRIQueryPlan(data_graph, query_graph, matching_order, pivots);
    }
    else if (input_order_type == "VF2PP") {
        GenerateQueryPlan::generateVF2PPQueryPlan(data_graph, query_graph, matching_order, pivots);
    }
    else if (input_order_type == "Spectrum") {
        GenerateQueryPlan::generateOrderSpectrum(query_graph, spectrum, order_num);
    }
    else {
        std::cout << "The specified order type '" << input_order_type << "' is not supported." << std::endl;
    }

    end = std::chrono::high_resolution_clock::now();
    double generate_query_plan_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();

    if (input_order_type != "Spectrum") {

        if (inputs.order_pointer == NULL){
            GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
            GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);

        }
        else {
            GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, inputs.order_pointer, pivots);
        }
    
#ifdef PRINT1
        GenerateQueryPlan::printSimplifiedQueryPlan(query_graph, matching_order);
#endif
    }
    else {
#ifdef PRINT1
        std::cout << "Generate " << spectrum.size() << " matching orders." << std::endl;
#endif

    }

#ifdef PRINT1
    std::cout << "-----" << std::endl;
    std::cout << "Enumerate..." << std::endl;
#endif

    //Add the matching order into return struct
    if(inputs.order_pointer == NULL){
        for(int i =0; i<outputs.query_size;i++){
            outputs.matching_order.push_back(matching_order[i]);
            outputs.matching_order_string.append(to_string(matching_order[i])+"-");
        }
        outputs.matching_order_string.pop_back();
    }
    else {
        for(int i =0; i<outputs.query_size;i++){
            outputs.matching_order.push_back(inputs.order_pointer[i]);
            outputs.matching_order_string.append(to_string(inputs.order_pointer[i])+"-");
        }
        outputs.matching_order_string.pop_back();
    }
    size_t output_limit = 0;
    size_t embedding_count = 0;
    if (input_max_embedding_num == "MAX") {
        output_limit = std::numeric_limits<size_t>::max();
    }
    else {
        sscanf(input_max_embedding_num.c_str(), "%zu", &output_limit);
    }
    //int embdcountaa=stoi(inputs.embcount);
    output_limit=stoi(inputs.embcount);;
#if ENABLE_QFLITER == 1
    EvaluateQuery::qfliter_bsr_graph_ = BuildTable::qfliter_bsr_graph_;
#endif

    size_t call_count = 0;
    size_t time_limit = 0;
     if(input_filter_type == "NLF"){
        output_limit=0;
    }
    sscanf(input_time_limit.c_str(), "%zu", &time_limit);

    start = std::chrono::high_resolution_clock::now();
    enumResult s;

    if (input_engine_type == "EXPLORE") {
        embedding_count = EvaluateQuery::exploreGraph(data_graph, query_graph, edge_matrix, candidates,
                                                      candidates_count, matching_order, pivots, output_limit, call_count);
    }
         else if (input_engine_type == "LFTJ"&& output_limit==0) {
        embedding_count = 0;
        outputs.call_count = 0;
        s.embedding_cnt=0;

    }
     else if (input_engine_type == "LFTJ") {
        
        
        if(inputs.order_pointer == NULL){
            s = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                    matching_order, output_limit, call_count);
           //embedding_count=0;
        } 
            
        else{
            cout<<"hi3"<<endl;
            s = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                    inputs.order_pointer, output_limit, call_count);
        }       
        embedding_count = s.embedding_cnt;
        outputs.call_count = call_count;
    }
    else if(input_engine_type=="LFTJA"){
        s = EvaluateQuery::LFTJA(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                    matching_order, output_limit, call_count,EWeight);
        embedding_count = s.embedding_cnt;
        
        outputs.call_count = call_count;
        }
        else if(input_engine_type=="LFTJB"){
        s = EvaluateQuery::LFTJB(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                    matching_order, output_limit, call_count,EWeight);
        embedding_count = s.embedding_cnt;
        outputs.call_count = call_count;
        }
        else if(input_engine_type=="LFTJR"){
        s = EvaluateQuery::LFTJR(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                    matching_order, output_limit, call_count);
        embedding_count = s.embedding_cnt;
        outputs.call_count = call_count;
        }
        else if(input_engine_type=="LFTCR"){
                    ui CC=0;
                for (int i=0; i<query_graph->getVerticesCount();i++){
                    CC+=candidates_count[i];
        }CC=CC-candidates_count[matching_order[query_graph->getVerticesCount()-1]];
           if (CC<2000)
           s = EvaluateQuery::LFTJB(data_graph, query_graph, edge_matrix, candidates, candidates_count,matching_order, output_limit, call_count,EWeight);
        else if(CC>20000)
                    s = EvaluateQuery::LFTJA(data_graph, query_graph, edge_matrix, candidates, candidates_count,matching_order, output_limit, call_count,EWeight);
        else
        s = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,matching_order, output_limit, call_count);
        embedding_count = s.embedding_cnt;
        outputs.call_count = call_count;
        }
     else if (input_engine_type == "GQL") {
        s = EvaluateQuery::exploreGraphQLStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, output_limit, call_count);

        embedding_count = s.embedding_cnt;

    } else if (input_engine_type == "QSI") {
        embedding_count = EvaluateQuery::exploreQuickSIStyle(data_graph, query_graph, candidates, candidates_count,
                                                             matching_order, pivots, output_limit, call_count);
    }
    else if (input_engine_type == "DPiso") {
        embedding_count = EvaluateQuery::exploreDPisoStyle(data_graph, query_graph, dpiso_tree,
                                                           edge_matrix, candidates, candidates_count,
                                                           weight_array, dpiso_order, output_limit,
                                                           call_count);
//        embedding_count = EvaluateQuery::exploreDPisoRecursiveStyle(data_graph, query_graph, dpiso_tree,
//                                                           edge_matrix, candidates, candidates_count,
//                                                           weight_array, dpiso_order, output_limit,
//                                                           call_count);
    }
    else if (input_engine_type == "Spectrum") {
        spectrum_analysis(data_graph, query_graph, edge_matrix, candidates, candidates_count, output_limit, spectrum, time_limit);
    }
    else if (input_engine_type == "CECI") {
        embedding_count = EvaluateQuery::exploreCECIStyle(data_graph, query_graph, ceci_tree, candidates, candidates_count, TE_Candidates,
                                                          NTE_Candidates, ceci_order, output_limit, call_count);
    }
    else {
        std::cout << "The specified engine type '" << input_engine_type << "' is not supported." << std::endl;
        exit(-1);
    }

    end = std::chrono::high_resolution_clock::now();
    double enumeration_time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();


    outputs.enumOutput = s;
//    bool isvalid_candidates = Experiments::candidate_set_correctness_check(outputs.candidate,s.candidate_true,query_graph->getVerticesCount());
//
//    if(!isvalid_candidates) throw invalid_argument("Invalid candidate set, false negative occurs.");

#ifdef DISTRIBUTION
    std::ofstream outfile (input_distribution_file_path , std::ofstream::binary);
    outfile.write((char*)EvaluateQuery::distribution_count_, sizeof(size_t) * data_graph->getVerticesCount());
    delete[] EvaluateQuery::distribution_count_;
#endif

#ifdef PRINT1
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "Release memories..." << std::endl;
#endif
    /**
     * Release the allocated memories.
     */
    delete[] candidates_count;
    delete[] tso_order;
    delete[] tso_tree;
    delete[] cfl_order;
    delete[] cfl_tree;
    delete[] dpiso_order;
    delete[] dpiso_tree;
    delete[] ceci_order;
    delete[] ceci_tree;
    delete[] matching_order;
    delete[] pivots;
    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;

    if (edge_matrix != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    if (weight_array != NULL) {
        for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
            delete[] weight_array[i];
        }
        delete[] weight_array;
    }

    delete query_graph;
    delete data_graph;

    /**
     * End.
     */

    double preprocessing_time_in_ns = filter_vertices_time_in_ns + build_table_time_in_ns + generate_query_plan_time_in_ns;
    double total_time_in_ns = preprocessing_time_in_ns + enumeration_time_in_ns;
    outputs.total_time = NANOSECTOSEC(total_time_in_ns);
    outputs.preprocessing_time = NANOSECTOSEC(preprocessing_time_in_ns);
    outputs.enumeration_time = NANOSECTOSEC(enumeration_time_in_ns);
#ifdef PRINT
    std::cout << "--------------------------------------------------------------------" << std::endl;
    printf("Load graphs time (seconds): %.4lf\n", NANOSECTOSEC(load_graphs_time_in_ns));
    printf("Filter vertices time (seconds): %.4lf\n", NANOSECTOSEC(filter_vertices_time_in_ns));
    printf("Build table time (seconds): %.4lf\n", NANOSECTOSEC(build_table_time_in_ns));
    printf("Generate query plan time (seconds): %.4lf\n", NANOSECTOSEC(generate_query_plan_time_in_ns));
    printf("Enumerate time (seconds): %.4lf\n", NANOSECTOSEC(enumeration_time_in_ns));
    printf("Preprocessing time (seconds): %.4lf\n", NANOSECTOSEC(preprocessing_time_in_ns));
    printf("Total time (seconds): %.4lf\n", NANOSECTOSEC(total_time_in_ns));
    printf("Memory cost (MB): %.4lf\n", BYTESTOMB(memory_cost_in_bytes));
    printf("#Embeddings: %zu\n", embedding_count);
    printf("Call Count: %zu\n", call_count);
    printf("Candidate count sum %zu\n",outputs.candidate_count_sum);
    printf("Candidate true sum so far %zu\n",outputs.enumOutput.candidate_true_count_sum);
    printf("Per Call Count Time (nanoseconds): %.4lf\n", enumeration_time_in_ns / (call_count == 0 ? 1 : call_count));
    std::cout << "End." << std::endl;
#endif

    return outputs;
}



