
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
static bool GenerateQueryD(ui QS,const Graph *data_graph,ui MinDen,int dd,string Dataset);
    static void randomRemoval(int **a,int *countS,int NE,ui QS);
    static bool GenerateQueryS(ui QS,const Graph *data_graph,ui MinDen,int dd,bool stop,string Dataset);

void randomRemoval(int **a,int *countS,int NE,ui QS){
    ui Edges=0;
    ui pos1;
    ui pos;
    int count=-1;
    while(NE>0){
        pos=rand() % QS;
        if(countS[pos]>1){
            pos1=rand()%countS[pos];
        for (int i=0;i<QS;i++)
            if(a[pos][i]==1){
                count++;
                if (count==pos1){
                    if(countS[i]>1){
                       countS[pos]--;
                       countS[i]--;
                       a[pos][i]=0;
                       a[i][pos]=0;
                       NE--; 
                    }
                else{ 
                    i=QS;
                }
                }

            }

        }
    }
}

bool GenerateQueryD(ui QS,const Graph *data_graph,ui MinDen,int dd,string Dataset){
    unordered_map<ui, ui> SID;
    vector<ui> ID;
    ui ED=0;
    ui dsize=data_graph->getVerticesCount();
    ui vertex=rand() % dsize;
    //ui vertex = getRandomNumber(0, dsize);
    ui nv;
    ui pos;
    ui vl;
    int countS[QS];
    ui count=0;
    int a[QS][QS];
    SID.insert({vertex,count});
    ID.emplace_back(vertex);
    count++;
    ui u_nbrs_count;
    for (int i=0;i<QS;i++)
    countS[i]=0;
    ui pp=0;
    ui aa=0;
    while(SID.size()<=QS){
        //cout<<SID.size()<<","<<aa<<" aa "<<vertex<<endl;

    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    if(u_nbrs_count==0)
    return false;
    pos=rand() % u_nbrs_count;
    //pos=pickRandomNumber(u_nbrs_count);
    //pos=getRandomNumber(0, u_nbrs_count);

    vertex=u_nbrs[pos];
    auto [it1, success] = SID.insert({vertex, count});
    
    if (success){
        
        count++;
        ID.emplace_back(vertex);
        aa=0;
        pp=0;
    }
    else{
        
    aa=aa+1;
    pp=pp+1;
    
    if(pp==10){
    pp=0;
    pos=rand() % ID.size();
    //pos=pickRandomNumber(ID.size());
    vertex=ID[pos];}
    if (aa==30){
        aa=0;
        
    bool sp=false;
    ui asap=0;
    while(!sp){
        
        vertex=ID[asap];
        const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
        for (int k=0;k<u_nbrs_count;k++){
            vertex=u_nbrs[k];
            auto [it, success] = SID.insert({vertex, count});
             if (success){
                    count++;
            ID.emplace_back(vertex);
            sp=true;
            break;
            }
            if (sp==true)
            break;
        }
        if(asap>=ID.size()-1)
        return false;
        if (!sp){
            asap=asap+1;
           //pos=rand() % ID.size();
            vertex=ID[asap]; 
        }
        
        
    }    
    }
    }
    }
    
    for (int i=0;i<QS;i++){
        vertex=ID[i];
        for (int j=i+1;j<QS;j++){
            nv=ID[j];
            if(data_graph->checkEdgeExistence(vertex,nv)){
                //if(true){
              a[SID[vertex]][SID[nv]]=1;
              a[SID[nv]][SID[vertex]]=1;
              countS[SID[vertex]]++;
              countS[SID[nv]]++;  
              ED++;
               
            }
            else{
              a[SID[vertex]][SID[nv]]=0;
              a[SID[nv]][SID[vertex]]=0;
            }
            
        }
    }
    
    if(((ED*2)/QS)>0){
        //std::ofstream outputFile("dataset\\wordnet\\query_graph\\query_dense_"+to_string(QS)+"_"+to_string(dd)+".graph");
        std::ofstream outputFile("../../dataset/"+Dataset+"/query_graph/query_G1_"+to_string(QS)+"_"+to_string(dd)+".graph");
        outputFile <<"t "<<QS<<" "<<ED<< std::endl;
        for (int i=0;i<QS;i++){
        vertex=SID[ID[i]];
        vl=data_graph->getVertexLabel(ID[i]);
        outputFile <<"v "<<vertex<<" "<<vl<<" "<<countS[vertex]<< std::endl;
        }
    for (int i=0;i<QS;i++){
        for (int j=i+1;j<QS;j++){
            if(a[i][j]==1)
            outputFile <<"e "<<i<<" "<<j<<" 0"<< std::endl;
        }
    }
    // Close the file
    ID.clear();
    SID.clear();
    outputFile.close();
        return true;
    }
    return false;
    
}

bool GenerateQueryS(ui QS,const Graph *data_graph,ui MinDen,int dd,bool stop,string Dataset){
    unordered_map<ui, ui> SID;
    vector<ui> ID;
    ui ED=0;
    ui dsize=data_graph->getVerticesCount();
    ui vertex=rand() % dsize;
    //ui vertex = getRandomNumber(0, dsize);
    ui nv;
    ui pos;
    ui vl;
    int countS[QS];
    ui count=0;
    int **a=NULL;
    
    a = new int *[QS];
    for (ui i = 0; i < QS; ++i)
    {
        a[i] = new int[QS];
    }
    SID.insert({vertex,count});
    ID.emplace_back(vertex);
    count++;
    ui u_nbrs_count;
    for (int i=0;i<QS;i++)
    countS[i]=0;
    ui pp=0;
    ui aa=0;
    while(SID.size()<=QS){
        //cout<<SID.size()<<","<<aa<<" aa "<<vertex<<endl;

    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    if(u_nbrs_count==0)
    return false;
    pos=rand() % u_nbrs_count;
    //pos=pickRandomNumber(u_nbrs_count);
    //pos=getRandomNumber(0, u_nbrs_count);

    vertex=u_nbrs[pos];
    auto [it1, success] = SID.insert({vertex, count});
    
    if (success){
        
        count++;
        ID.emplace_back(vertex);
        aa=0;
        pp=0;
    }
    else{
        
    aa=aa+1;
    pp=pp+1;
    
    if(pp==10){
    pp=0;
    pos=rand() % ID.size();
    //pos=pickRandomNumber(ID.size());
    vertex=ID[pos];}
    if (aa==30){
        aa=0;
        
    bool sp=false;
    ui asap=0;
    while(!sp){
        
        vertex=ID[asap];
        const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
        for (int k=0;k<u_nbrs_count;k++){
            vertex=u_nbrs[k];
            auto [it, success] = SID.insert({vertex, count});
             if (success){
                    count++;
            ID.emplace_back(vertex);
            sp=true;
            break;
            }
            if (sp==true)
            break;
        }
        if(asap>=ID.size()-1)
        return false;
        if (!sp){
            asap=asap+1;
           //pos=rand() % ID.size();
            vertex=ID[asap]; 
        }
        
        
    }    
    }
    }
    }
    
    for (int i=0;i<QS;i++){
        vertex=ID[i];
        for (int j=i+1;j<QS;j++){
            nv=ID[j];
            if(data_graph->checkEdgeExistence(vertex,nv)){
                //if(true){
              a[SID[vertex]][SID[nv]]=1;
              a[SID[nv]][SID[vertex]]=1;
              countS[SID[vertex]]++;
              countS[SID[nv]]++;  
              ED++;
               
            }
            else{
              a[SID[vertex]][SID[nv]]=0;
              a[SID[nv]][SID[vertex]]=0;
            }
            
        }
    }
    if (((ED*2)/QS)>MinDen&&stop){
        cout<<"hi"<<endl;
        randomRemoval(a,countS,(2*ED-(MinDen*QS))/2,QS);
        ED=MinDen*QS;
    }
    if(((ED*2)/QS)<=MinDen){
        cout<<"../../dataset/"+Dataset+"/query_graph/query_sparse_"+to_string(QS)+"_"+to_string(dd)+".graph"<<endl;
        std::ofstream outputFile("../../dataset/"+Dataset+"/query_graph/query_sparse_"+to_string(QS)+"_"+to_string(dd)+".graph");
        //std::ofstream outputFile("dataset\\wordnet\\query_graph\\query_sparse_"+to_string(QS)+"_"+to_string(dd)+".graph");
        outputFile <<"t "<<QS<<" "<<ED<< std::endl;
        for (int i=0;i<QS;i++){
        vertex=SID[ID[i]];
        vl=data_graph->getVertexLabel(ID[i]);
        outputFile <<"v "<<vertex<<" "<<vl<<" "<<countS[vertex]<< std::endl;
        }
    for (int i=0;i<QS;i++){
        for (int j=i+1;j<QS;j++){
            if(a[i][j]==1)
            outputFile <<"e "<<i<<" "<<j<<" 0"<< std::endl;
        }
    }
    // Close the file
    ID.clear();
    SID.clear();
    outputFile.close();
        return true;
    }
    return false;
}

double round_to(double value, double precision = 1.0)
{
    return std::round(value / precision) * precision;
}



pair<matching_algo_outputs,matching_algo_outputs> fakeMatchingWrapper(queryMeta meta,string filter){
    matching_algo_outputs original = Experiments::experiment3(meta.data_graph_path,meta.query_path,filter,"0",NULL);
    ui* fake_pointer = new ui[stoi(meta.query_size)];
    for (int i =0; i<stoi(meta.query_size);i++){
        ui order = original.matching_order[i];
        *&fake_pointer[i] = order;
    }
    matching_algo_outputs enhanced = Experiments::experiment3(meta.data_graph_path,meta.query_path,filter,"1",fake_pointer);
    delete[] fake_pointer;
    return pair(original,enhanced);
}

pair<matching_algo_outputs,matching_algo_outputs> MatchingWrapper(string datagraph,string querygraph,string filter){
    matching_algo_outputs original = Experiments::experiment3(datagraph,querygraph,filter,"0",NULL);
    matching_algo_outputs enhanced = Experiments::experiment3(datagraph,querygraph,filter,"1",NULL);
    return pair(original,enhanced);
}
pair<matching_algo_outputs,matching_algo_outputs> MatchingWrapperN(string datagraph,string querygraph,string filter){
    matching_algo_outputs original = Experiments::experiment3(datagraph,querygraph,filter,"0",NULL);
    matching_algo_outputs enhanced = Experiments::experiment3(datagraph,querygraph,filter,"1",NULL);
    return pair(original,enhanced);
}




void exact_eval(string dataset,string querysize,string querynumber,string property){
    Experiments::datagraphEigenMatrix = dataset+".csv";
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
    query<<"../../test/large_query/test.graph";
    meta.query_path = query.str();

    pair <matching_algo_outputs,matching_algo_outputs> LDF = fakeMatchingWrapper(meta,"LDF");
    pair <matching_algo_outputs,matching_algo_outputs> NLF = fakeMatchingWrapper(meta,"NLF");
    pair <matching_algo_outputs,matching_algo_outputs> GQL = fakeMatchingWrapper(meta,"GQL");
    pair <matching_algo_outputs,matching_algo_outputs> TSOF = fakeMatchingWrapper(meta,"TSO");
    pair <matching_algo_outputs,matching_algo_outputs> CFL = fakeMatchingWrapper(meta,"CFL");
    pair <matching_algo_outputs,matching_algo_outputs> DPiso = fakeMatchingWrapper(meta,"DPiso");
    matching_algo_outputs KF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"KF","0",NULL);

    vector<pair<matching_algo_outputs,matching_algo_outputs>> evaluations;
    evaluations.push_back(LDF);
    evaluations.push_back(NLF);
    evaluations.push_back(GQL);
    evaluations.push_back(TSOF);
    evaluations.push_back(CFL);
    evaluations.push_back(DPiso);


    std::ostringstream oss;
    oss <<meta.query_property<<"_"<<meta.query_size<<"_"<<meta.query_number;

    for(auto &eval : evaluations){
        oss<<","<<eval.first.call_count<<","<<eval.second.call_count;
    }
    oss<<","<<KF.call_count
    <<","<<LDF.first.enumOutput.embedding_cnt;

    for(auto &eval : evaluations){
        oss<<","<<eval.first.total_time<<","<<eval.second.total_time;
    }
    oss<<","<<KF.total_time<<","<<LDF.first.enumOutput.embedding_cnt;

    for(auto &eval : evaluations){
        oss<<","<<eval.first.candidate_count_sum<<","<<eval.second.candidate_count_sum;
    }
    oss<<","<<KF.candidate_count_sum;
    for(auto &eval : evaluations){
        oss<<","<<eval.first.matching_order_string<<","<<eval.second.matching_order_string;
    }
    oss<<","<<KF.matching_order_string;
    for(auto &eval : evaluations){
        oss<<","<<eval.first.preprocessing_time<<","<<eval.second.preprocessing_time;
    }
    oss<<","<<KF.preprocessing_time;
    for(auto &eval : evaluations){
        oss<<","<<eval.first.enumeration_time<<","<<eval.second.enumeration_time;
    }
    oss<<","<<KF.enumeration_time;


    std::string var = oss.str();

    cout<<var<<endl;

    string file_path = "";
    if(property=="sparse"){
        file_path = "performance_experiment/"+dataset+"_"+"s"+querysize+".csv";
    }
    if(property=="dense"){
        file_path = "performance_experiment/"+dataset+"_"+"d"+querysize+".csv";
    }

    std::ofstream myfile;
    myfile.open (file_path,std::ios_base::app);
    myfile<<var<<"\n";
    myfile.close();

}

void generate_datagraph_eigenvector(string data_graph_path,string csvfilename,int size){
    Graph* data_graph = new Graph(true);
    data_graph->loadGraphFromFile(data_graph_path);

    MatrixXd datagraph_eigenvalue(data_graph->getVerticesCount(), size);
    cout<<"Start compute eigen value"<<endl;
    MTcalc12(data_graph,data_graph->getGraphMaxDegree(),datagraph_eigenvalue,true,size,1100);
    saveData(csvfilename, datagraph_eigenvalue);

}

void fixed_order_experiment(int argc, char** argv){

    //    generate_datagraph_eigenvector("../../test/reallife_dataset/eu2005/data_graph/eu2005.graph","eu2005.csv");

    //yeast,hprd,uspatents,youtube,dblp,eu2005,
    vector<string> query_size_dense_1 = {"4","8","16","24","32"};
    vector<string> query_size_sparse_1 = {"8","16","24","32"};

    //human,wordnet
    vector<string> query_size_dense_2 = {"4","8","12","16","20"};
    vector<string> query_size_sparse_2 = {"8","12","16","20"};

    MatchingCommand command(argc,argv);
    string dataset_name = command.getDatasetName();
    string query_size = command.getQuerySize();
    string query_number = command.getQueryNumber();
    string query_property = command.getQueryProperty();

    cout<<dataset_name<<endl;
    cout<<query_size<<endl;
    cout<<query_number<<endl;
    cout<<query_property<<endl;

    exact_eval(dataset_name,query_size,query_number,query_property);
}

int main(int argc, char** argv) {
    //generate_datagraph_eigenvector("../../dataset/wordnet/data_graph/wordnet.graph","wordnet.csv",35);
    //generate_datagraph_eigenvector("../../dataset/dblp/data_graph/dblp.graph","dblp.csv",35);
    //generate_datagraph_eigenvector("../../dataset/eu2005/data_graph/eu2005.graph","eu2005.csv",35);
    //generate_datagraph_eigenvector("../../dataset//youtube/data_graph/youtube.graph","youtube.csv",35);
   // generate_datagraph_eigenvector("../../dataset/hprd/data_graph/hprd.graph","hprd.csv",35);
   //generate_datagraph_eigenvector("../../dataset/human/data_graph/human.graph","human.csv",35);
   //generate_datagraph_eigenvector("../../dataset/patents/data_graph/patents.graph","patents.csv",35);
   //generate_datagraph_eigenvector("../../dataset/yeast/data_graph/yeast.graph","yeast.csv",35);
   //return 0;
   // ExcelFormat::BasicExcel workbook;

    // Create sheets
   // workbook.New(2);
   // workbook.RenameWorksheet("Sheet1", "IT");
   // workbook.RenameWorksheet("Sheet2", "Time");

    // Get the data sheet
    
    //Query_Name,LDF,LDF+EF,NLF,NLF+EF,GQL,GQL+EF,TSOF,TSOF+EF,CFL,CFL+EF,DPiso,DPiso+EF,KF
//    string querygraph = "../../test/reallife_dataset/wordnet/query_graph/query_dense_16_1.graph";

    MatchingCommand command(argc,argv);
    string dataset_name = command.getDatasetName();
    string query_size = command.getQuerySize();
    string query_number = command.getQueryNumber();
    string query_property = command.getQueryProperty();
    string query_filter=command.getFilterType();
     string datagraph ="../../dataset/"+dataset_name+"/data_graph/"+dataset_name+".graph";
     Graph* data_graph = new Graph(true);
     data_graph->loadGraphFromFile(datagraph);
     int aa[5]={32,64,96,128,256};
    if(query_filter=="ALE"){
        cout<<"testAA"<<endl;
        bool ck=false;
    for (int da=0;da<=4;da++)
    for (int di=1;di<=400;di++){
        ui kk=aa[da];
        cout<<"kk, "<<kk<<" di "<<di<<endl;
    while(!GenerateQueryD(kk,data_graph,1,di,dataset_name));
    ck=false;
    ui cos=0;/*
    while(!ck){
        cout<<"Sparse kk, "<<kk<<" di "<<di<<endl;
        if(cos==30)
        ck=GenerateQueryS(kk,data_graph,3,di,true,dataset_name);
        else ck=GenerateQueryS(kk,data_graph,3,di,false,dataset_name);
        cos++;    
    }**/
     }
     return 0;
    }
    //query_number="1";
//    cout<<dataset_name<<endl;
//    cout<<query_size<<endl;
//    cout<<query_number<<endl;
//    cout<<query_property<<endl;

    //string datagraph = "../../test/reallife_dataset/wordnet/data_graph/wordnet.graph";
    //string querygraph = "../../test/randomwalk_queries/wordnet/50/randomwalk_50_"+query_number+".graph";
    Experiments::datagraphEigenMatrix = "../../"+dataset_name+".csv";
    //string datagraph = "../../dataset/wordnet/data_graph/wordnet.graph";
     datagraph ="../../dataset/"+dataset_name+"/data_graph/"+dataset_name+".graph";
    string querygraph = "../../dataset/"+dataset_name+"/query_graph/query_"+query_property+"_"+query_size+"_"+query_number+".graph";
     
    //string querygraph = "../../dataset/"+dataset_name+"/query_graph/query_"+query_property+"_"+query_size+"_"+query_number+".graph";
    //string querygraph = "../../dataset/wordnet/query_graph/query_dense_16_"+query_number+".graph";
//    string datagraph = "../../test/reallife_dataset/wordnet/data_graph/wordnet.graph";
//    string querygraph = "../../test/large_query/"+dataset_name+"/"+query_property+"/"+query_property+"_"+query_number+".graph";

//    string datagraph = "../../test/mydataset/youtube/data_graph/25-0/youtube.graph";
//    string querygraph = "../../test/mydataset/youtube/query_graph/25-0/query_"+query_property+"_"+query_size+"_"+query_number+".graph";


    //pair <matching_algo_outputs,matching_algo_outputs> LDF = MatchingWrapper(datagraph,querygraph,"LDF");
    //return 0;

    //pair <matching_algo_outputs,matching_algo_outputs> NLF = MatchingWrapper(datagraph,querygraph,"NLF");
    //pair <matching_algo_outputs,matching_algo_outputs> GQL = MatchingWrapper(datagraph,querygraph,"GQL");
    //pair <matching_algo_outputs,matching_algo_outputs> TSOF = MatchingWrapper(datagraph,querygraph,"TSO");
    //pair <matching_algo_outputs,matching_algo_outputs> CFL = MatchingWrapper(datagraph,querygraph,"CFL");
    //pair <matching_algo_outputs,matching_algo_outputs> DPiso = MatchingWrapper(datagraph,querygraph,"DPiso");
    //pair <matching_algo_outputs,matching_algo_outputs> KFA = MatchingWrapper(datagraph,querygraph,"KFA");
    //pair <matching_algo_outputs,matching_algo_outputs> KFB = MatchingWrapper(datagraph,querygraph,"KFB");
    //pair <matching_algo_outputs,matching_algo_outputs> KFC = MatchingWrapper(datagraph,querygraph,"KFC");
    //pair <matching_algo_outputs,matching_algo_outputs> KFD = MatchingWrapper(datagraph,querygraph,"KFD");
    //pair <matching_algo_outputs,matching_algo_outputs> KFE = MatchingWrapper(datagraph,querygraph,"KFE");
    matching_algo_outputs KF = Experiments::experiment3(datagraph,querygraph,query_filter,"0",NULL);

    vector<pair<matching_algo_outputs,matching_algo_outputs>> evaluations;
    /*
    evaluations.push_back(LDF);
    evaluations.push_back(NLF);
    evaluations.push_back(GQL);
    evaluations.push_back(TSOF);
    evaluations.push_back(CFL);
    evaluations.push_back(DPiso);
    */
    //evaluations.push_back(KFA);
    //evaluations.push_back(KFB);
    //evaluations.push_back(KFC);
    //evaluations.push_back(KFD);
    //evaluations.push_back(KFE);
    //evaluations.push_back(DPiso);

    std::ostringstream oss;
    oss <<query_property<<"_"<<query_size<<"_"<<query_number<<endl;
    //ExcelFormat::BasicExcelWorksheet* dataSheet = workbook.GetWorksheet("IT");
    int row=0;
    //for(auto &eval : evaluations){
    //    oss<<","<<eval.first.call_count<<","<<eval.second.call_count;
        //dataSheet->Cell(row, 0)->SetDouble(eval.first.call_count);
        //dataSheet->Cell(row, 1)->SetDouble(eval.second.call_count);
    //    row++;
        

    //}
    oss<<","<<KF.call_count<<","<<KF.enumOutput.embedding_cnt;
       //<<","<<LDF.first.enumOutput.embedding_cnt;
       //dataSheet->Cell(row, 0)->SetDouble(KF.call_count);
    //ExcelFormat::BasicExcelWorksheet* summarySheet = workbook.GetWorksheet("Time");
    
    //row=0;
    //for(auto &eval : evaluations){
    //    oss<<","<<eval.first.total_time<<","<<eval.second.total_time;
        //dataSheet->Cell(row, 0)->SetDouble(eval.first.total_time);
        //dataSheet->Cell(row, 1)->SetDouble(eval.second.total_time);
    //    row++;
   // }
   // dataSheet->Cell(row, 0)->SetDouble(KF.total_time);
    oss<<","<<KF.total_time<<","<<KF.candidate_count_sum;//<<LDF.first.enumOutput.embedding_cnt;

   // workbook.SaveAs("magkas.xlsx");
    /*
    for(auto &eval : evaluations){
        oss<<","<<eval.first.candidate_count_sum<<","<<eval.second.candidate_count_sum;
    }
    oss<<","<<KF.candidate_count_sum;
    for(auto &eval : evaluations){
        oss<<","<<eval.first.matching_order_string<<","<<eval.second.matching_order_string<<endl;
    }
    oss<<","<<KF.matching_order_string<<endl;*/
    //for(auto &eval : evaluations){
    //    oss<<","<<eval.first.preprocessing_time<<","<<eval.second.preprocessing_time;
   // }
    oss<<","<<KF.preprocessing_time;
    //for(auto &eval : evaluations){
    //    oss<<","<<eval.first.enumeration_time<<","<<eval.second.enumeration_time;
   // }
    oss<<","<<KF.enumeration_time;
    

    std::string var = oss.str();

    cout<<var<<endl;

    string file_path = "";

   // file_path = "performance_experiment/LArgeQ_"+query_filter+"_emb100000_"+dataset_name+"_"+query_property+query_size+".csv";
   
//file_path = "performance_experiment/LastTotalTime"+query_filter+"_emb100000_"+dataset_name+"_"+query_property+query_size+".csv";
    file_path = "performance_experiment/MultiThreading"+query_filter+"_emb100000_"+dataset_name+"_"+query_property+query_size+".csv";
    std::ofstream myfile;
    myfile.open (file_path,std::ios_base::app);
    myfile<<var<<"\n";
    myfile.close();

    return 0;

}


























//for(int i=2; i<5;i++){
//    for(int j=10;j<11;j++){
//        //yeast,hprd,uspatents,youtube,dblp,eu2005,
//        vector<string> query_size_dense_1 = {"4","8","16","24","32"};
//        vector<string> query_size_sparse_1 = {"8","16","24","32"};
//
//        //human,wordnet
//        vector<string> query_size_dense_2 = {"4","8","12","16","20"};
//        vector<string> query_size_sparse_2 = {"8","12","16","20"};
//
//
//        Experiments::datagraphEigenMatrix = "youtube.csv";
//        queryMeta meta;
//        meta.dataset = "youtube";
//        meta.query_property = "dense";
//        meta.query_size = query_size_dense_1[i];
//        meta.query_number = to_string(j);
//
//        std::ostringstream data;
//        data << "../../test/reallife_dataset/" << meta.dataset << "/data_graph/" << meta.dataset << ".graph";
//        meta.data_graph_path = data.str();
//
//        std::ostringstream query;
//        query << "../../test/reallife_dataset/" << meta.dataset << "/query_graph/query_" << meta.query_property << "_"
//        << meta.query_size << "_" << meta.query_number << ".graph";
//        meta.query_path = query.str();
//
//
//        //Query_Name,LDF,LDF+EF,NLF,NLF+EF,GQL,GQL+EF,TSOF,TSOF+EF,CFL,CFL+EF,DPiso,DPiso+EF,KF
//
//        matching_algo_outputs LDF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"LDF","0");
//        matching_algo_outputs LDF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"LDF","1");
//        matching_algo_outputs NLF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"NLF","0");
//        matching_algo_outputs NLF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"NLF","1");
//        matching_algo_outputs GQL = Experiments::experiment3(meta.data_graph_path,meta.query_path,"GQL","0");
//        matching_algo_outputs GQL_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"GQL","1");
//        matching_algo_outputs TSOF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"TSO","0");
//        matching_algo_outputs TSOF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"TSO","1");
//        matching_algo_outputs CFL = Experiments::experiment3(meta.data_graph_path,meta.query_path,"CFL","0");
//        matching_algo_outputs CFL_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"CFL","1");
//        matching_algo_outputs DPiso = Experiments::experiment3(meta.data_graph_path,meta.query_path,"DPiso","0");
//        matching_algo_outputs DPiso_EF =Experiments::experiment3(meta.data_graph_path,meta.query_path,"DPiso","1");
//        matching_algo_outputs KF =Experiments::experiment3(meta.data_graph_path,meta.query_path,"KF","0");
//
//        std::ostringstream oss;
//        oss <<meta.query_property<<"_"<<meta.query_size<<"_"<<meta.query_number<<"," <<LDF.total_time<<"," <<LDF_EF.total_time
//        <<","<<NLF.total_time<<"," <<NLF_EF.total_time<<","<<GQL.total_time<<"," <<GQL_EF.total_time<<"," <<TSOF.total_time<<"," <<TSOF_EF.total_time<<"," <<CFL.total_time
//        <<"," <<CFL_EF.total_time<<","<<DPiso.total_time<<"," <<DPiso_EF.total_time<<"," <<KF.total_time<<","<<LDF.enumOutput.embedding_cnt<<","<<LDF_EF.enumOutput.embedding_cnt
//        <<","<<NLF.enumOutput.embedding_cnt<<","<<NLF_EF.enumOutput.embedding_cnt<<","<<GQL.enumOutput.embedding_cnt
//        <<","<<GQL_EF.enumOutput.embedding_cnt<<","<<TSOF.enumOutput.embedding_cnt<<","<<TSOF_EF.enumOutput.embedding_cnt<<","<<CFL.enumOutput.embedding_cnt
//        <<","<<CFL_EF.enumOutput.embedding_cnt<<","<<DPiso.enumOutput.embedding_cnt<<","<<DPiso_EF.enumOutput.embedding_cnt<<","<<KF.enumOutput.embedding_cnt
//        <<","<<LDF.candidate_count_sum<<","<<LDF.enumOutput.candidate_true_count_sum<<","<<LDF_EF.candidate_count_sum<<","<<LDF_EF.enumOutput.candidate_true_count_sum
//        <<","<<NLF.candidate_count_sum<<","<<NLF.enumOutput.candidate_true_count_sum<<","<<NLF_EF.candidate_count_sum<<","<<NLF_EF.enumOutput.candidate_true_count_sum
//        <<","<<GQL.candidate_count_sum<<","<<GQL.enumOutput.candidate_true_count_sum<<","<<GQL_EF.candidate_count_sum<<","<<GQL_EF.enumOutput.candidate_true_count_sum
//        <<","<<TSOF.candidate_count_sum<<","<<TSOF.enumOutput.candidate_true_count_sum<<","<<TSOF_EF.candidate_count_sum<<","<<TSOF_EF.enumOutput.candidate_true_count_sum
//        <<","<<CFL.candidate_count_sum<<","<<CFL.enumOutput.candidate_true_count_sum<<","<<CFL_EF.candidate_count_sum<<","<<CFL_EF.enumOutput.candidate_true_count_sum
//        <<","<<DPiso.candidate_count_sum<<","<<DPiso.enumOutput.candidate_true_count_sum<<","<<DPiso_EF.candidate_count_sum<<","<<DPiso_EF.enumOutput.candidate_true_count_sum
//        <<","<<KF.candidate_count_sum<<","<<KF.enumOutput.candidate_true_count_sum;
//
//        std::string var = oss.str();
//
//        cout<<var<<endl;
//
//        std::ofstream myfile;
//        myfile.open ("performance_experiment/youtube.csv",std::ios_base::app);
//        myfile<<var<<"\n";
//        myfile.close();
//    }
//}

//
//for(int i=4; i<5;i++){
//for(int j=1;j<36;j++){
////yeast,hprd,uspatents,youtube,dblp,eu2005,
//vector<string> query_size_dense_1 = {"4","8","16","24","32"};
//vector<string> query_size_sparse_1 = {"8","16","24","32"};
//
////human,wordnet
//vector<string> query_size_dense_2 = {"4","8","12","16","20"};
//vector<string> query_size_sparse_2 = {"8","12","16","20"};
//
//
//Experiments::datagraphEigenMatrix = "wordnet.csv";
//queryMeta meta;
//meta.dataset = "wordnet";
//meta.query_property = "dense";
//meta.query_size = query_size_dense_2[i];
//meta.query_number = to_string(j);
//
//std::ostringstream data;
//data << "../../test/reallife_dataset/" << meta.dataset << "/data_graph/" << meta.dataset << ".graph";
//meta.data_graph_path = data.str();
//
//std::ostringstream query;
//query << "../../test/reallife_dataset/" << meta.dataset << "/query_graph/query_" << meta.query_property << "_"
//<< meta.query_size << "_" << meta.query_number << ".graph";
//meta.query_path = query.str();
//
//
////Query_Name,LDF,LDF+EF,NLF,NLF+EF,GQL,GQL+EF,TSOF,TSOF+EF,CFL,CFL+EF,DPiso,DPiso+EF,KF
//
//matching_algo_outputs LDF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"LDF","0");
//matching_algo_outputs LDF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"LDF","1");
//matching_algo_outputs NLF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"NLF","0");
//matching_algo_outputs NLF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"NLF","1");
//matching_algo_outputs GQL = Experiments::experiment3(meta.data_graph_path,meta.query_path,"GQL","0");
//matching_algo_outputs GQL_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"GQL","1");
//matching_algo_outputs TSOF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"TSO","0");
//matching_algo_outputs TSOF_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"TSO","1");
//matching_algo_outputs CFL = Experiments::experiment3(meta.data_graph_path,meta.query_path,"CFL","0");
//matching_algo_outputs CFL_EF = Experiments::experiment3(meta.data_graph_path,meta.query_path,"CFL","1");
//matching_algo_outputs DPiso = Experiments::experiment3(meta.data_graph_path,meta.query_path,"DPiso","0");
//matching_algo_outputs DPiso_EF =Experiments::experiment3(meta.data_graph_path,meta.query_path,"DPiso","1");
//matching_algo_outputs KF =Experiments::experiment3(meta.data_graph_path,meta.query_path,"KF","0");
//
//std::ostringstream oss;
//oss <<meta.query_property<<"_"<<meta.query_size<<"_"<<meta.query_number
//<<","<<LDF.candidate_count_sum<<","<<LDF_EF.candidate_count_sum
//<<","<<NLF.candidate_count_sum<<","<<NLF_EF.candidate_count_sum
//<<","<<GQL.candidate_count_sum<<","<<GQL_EF.candidate_count_sum
//<<","<<TSOF.candidate_count_sum<<","<<TSOF_EF.candidate_count_sum
//<<","<<CFL.candidate_count_sum<<","<<CFL_EF.candidate_count_sum
//<<","<<DPiso.candidate_count_sum<<","<<DPiso_EF.candidate_count_sum;
//
//std::string var = oss.str();
//
//cout<<var<<endl;
//
//std::ofstream myfile;
//myfile.open ("pruning_power_experiment/wordnet.csv",std::ios_base::app);
//myfile<<var<<"\n";
//myfile.close();
//}
//}