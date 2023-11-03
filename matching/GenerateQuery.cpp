
#include "GenerateQuery.h"

double round_to(double value, double precision = 1.0)
{
    return std::round(value / precision) * precision;
}

// RI -> flaged to improve the code quality
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