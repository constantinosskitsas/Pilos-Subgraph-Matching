
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