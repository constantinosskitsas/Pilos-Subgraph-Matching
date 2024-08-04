#include "cs.h"
#include <chrono>
#include <thread>
#include <Eigen/StdVector>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <utility>
#include <cmath>
#include <Eigen/SparseCore>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <list>
void ExtractUI2h(vector<ui> &Deg, vector<map<ui, int>> &QueryNlabel2, Graph *query_graph, int qsiz, int **&VS);
typedef Eigen::Triplet<double> T;
inline bool OneHopEigenVGBM(CSV &cvertex, int qid, map<ui, int> EvalNeigb, Graph *query_graph, ui *&flag, ui *&updated_flag, int *&flag1, int *offsetQ, int *left_to_right_offset, int *left_to_right_edges,
                            int *left_to_right_match, int *right_to_left_match, int *match_visited,
                            int *match_queue, int *match_previous);
void EdgesCSBasicRL(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, ui *&flag, ui *&updated_flag);
inline bool OneHopEigenVG(CSV &cvertex, map<ui, int> EvalNeigb, Graph *query_graph, ui *&flag, ui *&updated_flag);

int PILOS(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix);
inline bool SHEigenb(vector<pair<VertexID, VertexID>> &q_curr, ui *&flag, ui *&flagq, ui *&updated_flag, ui *&updated_flagq, int *&Qindex2,
                     int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex, int beta);
bool RFNV(vector<map<ui, int>> NLabel, vector<map<ui, int>> NLabel2, vector<vector<CSV>> &FCS,
          int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> DM, int twohop, int alpha, int beta, ui *&flag, ui *&updated_flag, int **&VS);
inline bool OHEPM(CSV &cvertex, ui *&flag, ui *&flagq, ui *updated_flag, ui *updated_flagq, int *IDDLC,
                  float **&LM, map<ui, int> &EvalNeigb2, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
bool RefinementNV(vector<map<ui, int>> NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, Graph *data_graph, ui GDegree, ui *&flag, ui *&updated_flag);
void VerticesNF(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVq1, float **&eigenVD1);
void ExtractSecHopIndex(vector<ui> &Deg, int **&Qindex2, Graph *query_graph, int qsiz);
inline bool OneHopEigenVV(CSV &cvertex, map<ui, int> EvalNeigb, Graph *query_graph, ui *&flag, ui *&updated_flag);
void EdgesCSBasicSetVV(vector<vector<ui>> &verticesID, vector<ui> &NumOfverticesID, vector<vector<vector<ui>>> &EQID,
                       vector<vector<vector<ui>>> &EVID, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph);
void VerticesVe(vector<vector<ui>> &verticesID, vector<ui> &NumOfverticesID, vector<vector<bool>> &ValidverticesID, vector<vector<vector<ui>>> &EQID, vector<vector<vector<ui>>> &EVID, int qsiz,
                int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVq1, vector<map<ui, int>> &QueryNlabel, float **&eigenVD1);
int CSInitV(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix, size_t **&candidatesHC, unordered_map<size_t, vector<ui>> *&idToValues2);
int SpectralMatchingV(int sizd, Graph *data_graph, Graph *query_graph, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix, size_t **&candidatesHC, unordered_map<size_t, vector<ui>> *&idToValues2);
bool k_con(vector<vector<CSV>> &FCS, int qid, int vid, int qsiz, Graph *query_graph, int ksiz_);
void naiveCalcP2(int qsiz, ui **&candidates, ui *&candidates_count, int **&candidatesP, unordered_map<int, vector<int>> &idToValues, size_t **&candidatesHC, vector<vector<CSV>> &FCS);
bool check_edgesP(vector<pair<VertexID, VertexID>> &v1, vector<pair<VertexID, VertexID>> &v2);
void naiveCalcP(vector<vector<CSV>> &FCS, ui **&candidatesP, unordered_map<int, vector<int>> &idToValues, size_t **&candidatesHC);
inline bool SecHopEigenLMbeta(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, ui *&SID, map<ui, int> &EvalNeigb2,
                              int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex, int beta);
int SpectralMatchingMT(int sizd, Graph *data_graph, string input_query_graph_file, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int thnum, int beta);
int SpectralMatching(int sizd, Graph *data_graph, Graph *query_graph, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix, float *&eigenQS);
int CSInit(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix);
int CSInitMT(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int thnum, int beta);
inline void removeVertexAndEgjesFKNPMT(vector<vector<CSV>> &FCS, int i, int deli);
inline bool OneHopEigenMapMT(CSV &cvertex, map<ui, int> EvalNeigb, Graph *query_graph);
bool ReverseRefinementNOTESNMT(vector<map<ui, int>> &NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, ui GDegree, int thnum);
void EdgesCSBasicSetMT(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, int thnum);
bool RefinementEigenMT2(vector<map<ui, int>> &NLabel, vector<map<ui, int>> &NLabel2, vector<vector<CSV>> &FCS,
                        int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> &DM, int twohop, int alpha, int thnum, int beta);
void VerticesMT2(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVq1, vector<map<ui, int>> &QueryNlabel, float **&eigenVD1, int thnum);
bool ReverseRefinementNOTESN1(vector<map<ui, int>> NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, ui GDegree);
inline bool OneHopEigenPM(CSV &cvertex, unordered_map<ui, ui> &ID, ui *&SID, int *IDDLC,
                          float **&LM, map<ui, int> &EvalNeigb2, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
inline bool SecHopEigenLM(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, ui *&SID, map<ui, int> &EvalNeigb2,
                          int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex);
inline bool OneHopEigenPC(CSV &cvertex, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, int *IDDLC,
                          vector<T> &tripletList, map<ui, int> &EvalNeigb2, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax, ui Omax2);
inline bool SecHopEigenCP(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, map<ui, int> &EvalNeigb2,
                          int *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, Graph *query_graph, ui Omax, ui Omax2, int qID, ui DM);
void allocateBufferFCS1(vector<vector<CSV>> &FCS, const Graph *query_graph, ui **&candidates,
                        ui *&candidates_count, float **&EWeight);
inline bool OneHopEigen(CSV &cvertex, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, int *IDDLC,
                        vector<T> &tripletList, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
inline bool SecHopEigenPNSSJE(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, unordered_set<ui> &SID,
                              ui *IDDLC, vector<vector<CSV>> &FCS, vector<T> &tripletList, Graph *query_graph, ui Omax, int qID);
inline bool OneHopEigen(CSV &cvertex, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, ui *IDDLC,
                        vector<T> &tripletList, map<ui, int> EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
bool RefinementEigen(vector<map<ui, int>> NLabel, vector<map<ui, int>> NLabel2, vector<vector<CSV>> &FCS,
                     int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> DM, int twohop, int alpha, int beta);

bool RefinementEigen(vector<map<ui, int>> NLabel, vector<map<ui, int>> NLabel2, vector<vector<CSV>> &FCS,
                     int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> DM, int twohop);
bool RefinementEigenMT(vector<map<ui, int>> NLabel, vector<map<ui, int>> NLabel2, vector<vector<CSV>> &FCS,
                       int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> DM, int twohop, int alpha, int thnum);
inline bool SHEigen(vector<pair<VertexID, VertexID>> &q_curr, ui *&flag, ui *&flagq, ui *&updated_flag, ui *&updated_flagq, map<ui, int> &EvalNeigb2,
                    int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex);
int CSSizeReal(vector<vector<CSV>> &FCS, int qsiz);
inline void removeVertexAndEgjesFKtest(vector<vector<CSV>> &FCS, int i, int deli);
bool ReverseRefinementNOTES(vector<map<ui, int>> NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, ui GDegree);
bool ReverseRefinementNOTESN(vector<map<ui, int>> NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, ui GDegree);
int CSInit(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1);
inline VertexID findIndBS(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ);
inline void removeVertexAndEgjesFK(vector<vector<CSV>> &FCS, int i, int deli);
void clearWrong(vector<vector<CSV>> &FCS);
void allocateBufferFCS(vector<vector<CSV>> &FCS, const Graph *query_graph, ui **&candidates,
                       ui *&candidates_count);
void ExtractUI(vector<ui> &Deg, Graph *query_graph, int qsiz);
void ExtractNImap(vector<map<ui, int>> &QueryNlabel, Graph *query_graph, int qsiz);
bool InitPrunTCSR(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph);
void EdgesCSBasicSet(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph);
void Vertices(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVD1, float **&eigenVq1, vector<map<ui, int>> &QueryNlabel);
inline void removeVertexAndEgjesFKNP(vector<vector<CSV>> &FCS, int i, int deli);
int SpectralMatching(int sizd, Graph *data_graph, string input_query_graph_file, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, float *&eigenQS);
void fillEN(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph);
inline bool OneHopEigen(CSV &cvertex, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, int *IDDLC,
                        vector<T> &tripletList, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
static inline bool OneHopEigenP(CSV &cvertex, unordered_map<ui, ui> &ID, unordered_set<ui> &SID, int *IDDLC,
                                vector<T> &tripletList, map<ui, int> &EvalNeigb2, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax);
