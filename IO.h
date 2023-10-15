#include <iostream>
#include<vector>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <utility> 
using namespace std;
using  namespace Eigen;
void openData1(string fileToOpen,float **&eigsData);
MatrixXd openData(string fileToOpen);
void saveData(string fileName, MatrixXd  matrix);
void addEigenstoCVS();
