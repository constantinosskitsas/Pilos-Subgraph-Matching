#include <iostream>
#include <vector>
#include <Eigen/Core>
#include "../GrM.h"
using namespace std;
using  namespace Eigen;



 class CSV{
    public:  

VertexID ID;
vector<pair<VertexID,VertexID>> edges;
int* Nedge;
bool NedgeC=false;
bool change=true;
bool Ichange=true;
bool IPchange=true;
bool deleted=false;
float ED=100000;
public: 
CSV(int eigens,VertexID IDV,ui maxDeg){

ID=IDV;
deleted=false;
change=true;
Ichange=true;
}

CSV(ui IDV){
ID=IDV;
deleted=false;
change=true;
Ichange=true;
}

CSV(){
    ID=-1;
}

};

