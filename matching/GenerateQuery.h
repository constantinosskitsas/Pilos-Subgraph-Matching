

#include <future>
#include <thread>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

#include "graph/graph.h"
#include "IO.h"
 bool GenerateQueryD(ui QS,const Graph *data_graph,ui MinDen,int dd,string Dataset);
 void randomRemoval(int **a,int *countS,int NE,ui QS);
 bool GenerateQueryS(ui QS,const Graph *data_graph,ui MinDen,int dd,bool stop,string Dataset);