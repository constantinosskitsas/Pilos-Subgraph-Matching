
#include "spectra.h"
#include "../../graph/graph.h"
#include <unordered_map>
#include "../eigenHelper.h"
#include "../GrM.h"
#include "../IO.h"
#include "../Experiments.h"
#include <utility/graphoperations.h>
#include <numeric>
#include <mutex>
#include <functional>
#include <algorithm>
#include <cstring>
#include "../EvaluateQuery.h"
using namespace std;
using namespace Eigen;
using namespace std::chrono;
std::mutex fcsMutex;
std::mutex Mymutex;
std::mutex Mymutex1;
std::mutex mtx;

inline bool OHEPM(CSV &cvertex, ui *&flag, ui *&flagq, ui *updated_flag, ui *updated_flagq, int *IDDLC,
                  float **&LM, map<ui, int> &EvalNeigb2, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax)
{

    // store the flag as the IDDLC+1
    // so when we look the position in array will be
    // x=flag[x]-1
    flag[cvertex.ID] = 1;         // value in LM matrix -1 same as count
    updated_flag[0] = cvertex.ID; //
    IDDLC[0]++;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    ui count2 = 0;
    ui labela = 0;
    // nedd IDDLC[1] reverse flag counter
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;
        flagq[cvertex.edges[k].first] = 1;
        if (flag[cvertex.edges[k].second] == 0)
        {
            flag[cvertex.edges[k].second] = IDDLC[0] + 1;
            updated_flag[IDDLC[0]] = cvertex.edges[k].second;
            IDDLC[0]++;
            if (IDDLC[1] > 0)
            {
                labela = query_graph->getVertexLabel(cvertex.edges[k].first);
                EvalNeigb2[labela]--;
                EvalNeigb[labela]--;
                if (EvalNeigb2[labela] == 0)
                    IDDLC[1]--;
                if (EvalNeigb[labela] == 0)
                    IDDLC[2]--;
            }
            if ((IDDLC[0]) > Omax)
            {
                return true;
            }
            else
            {
                LM[0][IDDLC[0] - 1] = -1;
                LM[IDDLC[0] - 1][0] = -1;
            }
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    return true;
}

inline bool OneHopEigenPM(CSV &cvertex, unordered_map<ui, ui> &ID, ui *&SID, int *IDDLC,
                          float **&LM, map<ui, int> &EvalNeigb2, map<ui, int> &EvalNeigb, Graph *query_graph, vector<pair<VertexID, VertexID>> &q_curr, ui Omax)
{

    ID.insert({cvertex.ID, 0});
    IDDLC[0] = 1;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    ui count2 = 0;
    ui stop = EvalNeigb2.size();
    ui labela = 0;
    for (int k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;

        SID[cvertex.edges[k].first] = 1;
        auto result = ID.insert({cvertex.edges[k].second, IDDLC[0]});

        if (result.second)
        {
            if (IDDLC[1] > 0)
            {
                labela = query_graph->getVertexLabel(cvertex.edges[k].first);
                EvalNeigb2[labela]--;
                EvalNeigb[labela]--;
                if (EvalNeigb2[labela] == 0)
                    IDDLC[1]--;
                if (EvalNeigb[labela] == 0)
                    IDDLC[2]--;
            }
            if (IDDLC[0] > Omax)
            {
                return true;
            }
            else
            {
                LM[0][IDDLC[0]] = -1;
                LM[IDDLC[0]][0] = -1;
                IDDLC[0]++;
            }
        }
        q_curr.emplace_back(cvertex.edges[k]);
    }
    return true;
}
/*Second Hop Evaluation
 *Omax size limit to calculate eigenvalues and Omax2 limit of 2-hop
 * Extension of degree and Neighborhood safety pruning
 *Add edges in triplet format.
 */

inline bool SecHopEigenLM(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, ui *&SID, map<ui, int> &EvalNeigb2,
                          int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex)
{
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID DN = 0;
    unordered_set<ui> EdgeF;
    int kk = 0;
    ui labela = 0;
    VertexID vtemp = 0;
    vector<ui> tempx;
    ui counter = 0;
    kk = 0;
    while (kk < q_curr.size())
    {
        temp1 = q_curr[kk];
        kk++;
        if (temp1.first == 1000000)
        {
            continue;
        }
        tempxx = findIndBS(FCS, temp1.second, temp1.first);
        vx1 = ID[FCS[temp1.first][tempxx].ID];
        for (int i = 0; i < FCS[temp1.first][tempxx].edges.size(); i++)
        {

            if (FCS[temp1.first][tempxx].edges[i].first == 1000000 || FCS[temp1.first][tempxx].edges[i].first == qID || FCS[temp1.first][tempxx].edges[i].second == cvertex.ID)
                continue;
            SID[FCS[temp1.first][tempxx].edges[i].first] = 1;
            vtemp = FCS[temp1.first][tempxx].edges[i].second;
            auto [it1, success] = ID.try_emplace(vtemp, IDDLC[0]);

            if (success)
            {
                vx2 = IDDLC[0];
                IDDLC[0]++;
                if (IDDLC[1] > 0)
                {
                    labela = query_graph->getVertexLabel(FCS[temp1.first][tempxx].edges[i].first);
                    EvalNeigb2[labela]--;
                    if (EvalNeigb2[labela] == 0)
                        IDDLC[1]--;
                }
                if (IDDLC[0] > Omax)

                {
                    return true;
                }
            }
            else
            {
                vx2 = ID[vtemp];
            }

            LM[vx1][vx2] = -1;
            LM[vx2][vx1] = -1;
        }
    }
    q_curr.clear();
    return (true);
}

inline bool SHEigenb(vector<pair<VertexID, VertexID>> &q_curr, ui *&flag, ui *&flagq, ui *&updated_flag, ui *&updated_flagq, map<ui, int> &EvalNeigb2,
                     int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex, int beta)
{
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID DN = 0;
    unordered_set<ui> EdgeF;
    int kk = 0;
    ui labela = 0;
    VertexID vtemp = 0;
    vector<ui> tempx;
    ui counter = 0;
    kk = 0;
    ui omaxUP = beta;
    while (kk < q_curr.size())
    {
        temp1 = q_curr[kk];
        kk++;

        if (temp1.first == 1000000)
        {

            continue;
        }

        tempxx = findIndBS(FCS, temp1.second, temp1.first);
        vx1 = flag[FCS[temp1.first][tempxx].ID] - 1;
        for (int i = 0; i < FCS[temp1.first][tempxx].edges.size(); i++)
        {
            counter++;
            if (FCS[temp1.first][tempxx].edges[i].first == 1000000 || FCS[temp1.first][tempxx].edges[i].first == qID || FCS[temp1.first][tempxx].edges[i].second == cvertex.ID)
                continue;
            flagq[FCS[temp1.first][tempxx].edges[i].first] = 1;
            // with possible reverse flag but really  needed?
            vtemp = FCS[temp1.first][tempxx].edges[i].second;

            if (flag[vtemp] == 0)
            {
                if ((IDDLC[0]) > Omax)
                {
                    return true;
                }
                flag[vtemp] = IDDLC[0] + 1;
                updated_flag[IDDLC[0]] = vtemp;
                vx2 = IDDLC[0];
                IDDLC[0]++;
                if (IDDLC[1] > 0)
                {
                    labela = query_graph->getVertexLabel(FCS[temp1.first][tempxx].edges[i].first);
                    EvalNeigb2[labela]--;
                    if (EvalNeigb2[labela] == 0)
                        IDDLC[1]--;
                }
            }
            else
            {
                vx2 = flag[vtemp] - 1;
            }
            LM[vx1][vx2] = -1;

            LM[vx2][vx1] = -1;

            if (counter > omaxUP)
            {
                IDDLC[0] = Omax + 1;
                return true;
            }
        }
    }
    q_curr.clear();
    return (true);
}

inline bool SHEigen(vector<pair<VertexID, VertexID>> &q_curr, ui *&flag, ui *&flagq, ui *&updated_flag, ui *&updated_flagq, map<ui, int> &EvalNeigb2,
                    int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex)
{
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID DN = 0;
    unordered_set<ui> EdgeF;
    int kk = 0;
    ui labela = 0;
    VertexID vtemp = 0;
    vector<ui> tempx;
    ui counter = 0;
    kk = 0;
    ui tlqb = 0;
    while (kk < q_curr.size())
    {
        temp1 = q_curr[kk];
        kk++;

        if (temp1.first == 1000000)
        {
            continue;
        }

        tempxx = findIndBS(FCS, temp1.second, temp1.first);
        vx1 = flag[FCS[temp1.first][tempxx].ID] - 1;
        for (int i = 0; i < FCS[temp1.first][tempxx].edges.size(); i++)
        {
            tlqb = FCS[temp1.first][tempxx].edges[i].first;
            vtemp = FCS[temp1.first][tempxx].edges[i].second;
            if (tlqb == 1000000 || tlqb == qID || vtemp == cvertex.ID)
                continue;
            flagq[tlqb] = 1;
            // with possible reverse flag but really  needed?

            if (flag[vtemp] == 0)
            {

                flag[vtemp] = (IDDLC[0] + 1);
                updated_flag[IDDLC[0]] = vtemp;
                vx2 = IDDLC[0];
                IDDLC[0]++;
                if (IDDLC[1] > 0)
                {
                    labela = query_graph->getVertexLabel(tlqb);
                    EvalNeigb2[labela]--;
                    if (EvalNeigb2[labela] == 0)
                        IDDLC[1]--;
                }
                if (IDDLC[0] > Omax)
                {
                    return true;
                }
            }
            else
            {
                vx2 = flag[vtemp] - 1;
            }
            LM[vx1][vx2] = -1;
            LM[vx2][vx1] = -1;
        }
    }
    q_curr.clear();
    return (true);
}

inline bool SecHopEigenLMbeta(vector<pair<VertexID, VertexID>> &q_curr, unordered_map<ui, ui> &ID, ui *&SID, map<ui, int> &EvalNeigb2,
                              int *IDDLC, vector<vector<CSV>> &FCS, float **&LM, Graph *query_graph, ui Omax, int qID, CSV &cvertex, int beta)
{
    pair<VertexID, VertexID> temp1;
    VertexID tempxx = 0;
    vector<VertexID> temp2;
    VertexID vx1 = 0;
    VertexID vx2 = 0;
    VertexID DN = 0;
    unordered_set<ui> EdgeF;
    int kk = 0;
    ui labela = 0;
    VertexID vtemp = 0;
    vector<ui> tempx;
    ui counter = 0;
    kk = 0;
    ui omaxUP = beta;
    while (kk < q_curr.size())
    {
        temp1 = q_curr[kk];
        kk++;
        if (temp1.first == 1000000)
        {

            continue;
        }

        tempxx = findIndBS(FCS, temp1.second, temp1.first);
        vx1 = ID[FCS[temp1.first][tempxx].ID];
        for (int i = 0; i < FCS[temp1.first][tempxx].edges.size(); i++)
        {
            counter++;
            if (FCS[temp1.first][tempxx].edges[i].first == 1000000 || FCS[temp1.first][tempxx].edges[i].first == qID || FCS[temp1.first][tempxx].edges[i].first == qID || FCS[temp1.first][tempxx].edges[i].second == cvertex.ID)
                continue;
            SID[FCS[temp1.first][tempxx].edges[i].first] = 1;
            vtemp = FCS[temp1.first][tempxx].edges[i].second;
            auto [it1, success] = ID.try_emplace(vtemp, IDDLC[0]);

            if (success)
            {
                vx2 = IDDLC[0];
                IDDLC[0]++;
                if (IDDLC[1] > 0)
                {
                    labela = query_graph->getVertexLabel(FCS[temp1.first][tempxx].edges[i].first);
                    EvalNeigb2[labela]--;
                    if (EvalNeigb2[labela] == 0)
                        IDDLC[1]--;
                }
                if (IDDLC[0] > Omax)
                {
                    return true;
                }
            }
            else
            {
                vx2 = ID[vtemp];
            }
            LM[vx1][vx2] = -1;
            LM[vx2][vx1] = -1;
            if (counter > omaxUP)
            {
                IDDLC[0] = Omax + 1;
                return true;
            }
        }
    }
    q_curr.clear();
    return (true);
}

void ExtractSecHopIndex(vector<ui> &Deg, int **&Qindex2, Graph *query_graph, int qsiz)
{
    const VertexID *u_nbrs;
    ui u_nbrs_count;
    const VertexID *u_nbrsN;
    ui u_nbrs_countN;
    set<ui> QueryVec;
    ui labela = 0;
    for (int i = 0; i < qsiz; i++)
    {
        QueryVec.insert(i);
        for (int j = 0; j < u_nbrs_count; j++)
        { // First hop
            Qindex2[i][query_graph->getVertexLabel(u_nbrs[j])]++;
            QueryVec.insert(u_nbrs[j]).second;
            // second Hop
            u_nbrsN = query_graph->getVertexNeighbors(u_nbrs[j], u_nbrs_countN);
            for (int k = 0; k < u_nbrs_countN; k++)
            {

                if (QueryVec.insert(u_nbrsN[k]).second)
                {
                    Qindex2[i][query_graph->getVertexLabel(u_nbrsN[k])]++;
                }
            }
        }
        Deg.emplace_back(QueryVec.size());
        QueryVec.clear();
    }
}

/*Extract Label Information for the query for 2hops
 *Createas a map with Label ID and count that we can easily compare
 *with a candidate node. Self Included!
 */
void ExtractUI2h(vector<ui> &Deg, vector<map<ui, int>> &QueryNlabel2, Graph *query_graph, int qsiz, int **&VS)
{
    const VertexID *u_nbrs;
    ui u_nbrs_count;
    const VertexID *u_nbrsN;
    ui u_nbrs_countN;
    set<ui> QueryVec;
    map<ui, int> QueryVec1;
    ui labela = 0;
    int matrix[qsiz][qsiz] = {-1};
    for (int i = 0; i < qsiz; i++)
        for (int j = 0; j < qsiz; j++)
            matrix[i][j] = -1;
    for (int i = 0; i < qsiz; i++)
    {
        QueryVec.insert(i);
        u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);

        for (int j = 0; j < u_nbrs_count; j++)
        { // First hop
            if (QueryVec.insert(u_nbrs[j]).second)
            {
                matrix[i][u_nbrs[j]] = i;
                matrix[u_nbrs[j]][i] = i;
                labela = query_graph->getVertexLabel(u_nbrs[j]);
                if (QueryVec1.count(labela) == 0)
                {
                    // Key does not exist, add it with a value of 1
                    QueryVec1[labela] = 1;
                }
                else
                {
                    // Key exists, increment its value
                    QueryVec1[labela]++;
                }
            }
            // second Hop
            u_nbrsN = query_graph->getVertexNeighbors(u_nbrs[j], u_nbrs_countN);
            for (int k = 0; k < u_nbrs_countN; k++)
            {
                matrix[u_nbrsN[k]][u_nbrs[j]] = i;
                matrix[u_nbrs[j]][u_nbrsN[k]] = i;
                if (u_nbrsN[k] == i)
                    continue;
                labela = query_graph->getVertexLabel(u_nbrsN[k]);
                if (QueryVec.insert(u_nbrsN[k]).second)
                {
                    if (QueryVec1.count(labela) == 0)
                    {
                        // Key does not exist, add it with a value of 1
                        QueryVec1[labela] = 1;
                    }
                    else
                    {
                        // Key exists, increment its value
                        QueryVec1[labela]++;
                    }
                }
            }
        }
        int sumr = 0;
        for (int dd = 0; dd < qsiz; dd++)
        {
            sumr = 0;
            for (int kk = 0; kk < qsiz; kk++)
            {
                if (matrix[dd][kk] == i)
                    sumr++;
            }
            VS[i][dd] = sumr;
        }
        Deg.emplace_back(QueryVec.size());
        QueryVec.clear();
        QueryNlabel2.emplace_back(QueryVec1);
        QueryVec1.clear();

        std::sort(VS[i], VS[i] + qsiz, std::greater<int>());

    } 
}

inline bool OneHopEigenVG(CSV &cvertex, map<ui, int> EvalNeigb, Graph *query_graph, ui *&flag, ui *&updated_flag)
{
    ui count2 = 0;
    ui labela = 0;
    ui k;
    ui updated_flag_count = 0;
    for (k = 0; k < cvertex.edges.size(); k++)
    {
        if (cvertex.edges[k].first == 1000000)
            continue;
        labela = query_graph->getVertexLabel(cvertex.edges[k].first);
        // check label only if it didnt pass the check yet
        if (EvalNeigb[labela] <= 0)
            continue;
        // Count only unique ID
        if (flag[cvertex.edges[k].second] != 1)
        {
            flag[cvertex.edges[k].second] = 1;
            updated_flag[updated_flag_count++] = cvertex.edges[k].second;
            EvalNeigb[labela]--;
            if (EvalNeigb[labela] == 0)
            {
                count2++;
                if (count2 == EvalNeigb.size())
                {
                    for (ui aa = 0; aa <= (updated_flag_count + 1); aa++)
                    {
                        flag[updated_flag[aa]] = 0;
                    }
                    return true;
                }
            }
        }
    }
    for (ui aa = 0; aa <= (updated_flag_count + 1); aa++)
    {
        flag[updated_flag[aa]] = 0;
    }
    return false;
}

void allocateBufferFCS1(vector<vector<CSV>> &FCS, const Graph *query_graph, ui **&candidates,
                        ui *&candidates_count, float **&EWeight)
{
    ui query_vertex_num = query_graph->getVerticesCount();
    candidates_count = new ui[query_vertex_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui *[query_vertex_num];
    EWeight = new float *[query_vertex_num];
    for (ui i = 0; i < query_vertex_num; ++i)
    {
        candidates[i] = new ui[FCS[i].size()];
        EWeight[i] = new float[FCS[i].size()];
    }
}
/*Extract Label Information for the query
 *Createas a map with Label ID and count that we can easily compare
 *with a candidate node
 */
void ExtractNImap(vector<map<ui, int>> &QueryNlabel, Graph *query_graph, int qsiz)
{
    const VertexID *u_nbrs;

    ui u_nbrs_count;
    for (int i = 0; i < qsiz; i++)
    {
        u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
        map<ui, int> QueryVec;
        for (int j = 0; j < u_nbrs_count; j++)
        {
            ui labela = query_graph->getVertexLabel(u_nbrs[j]);
            if (QueryVec.count(labela) == 0)
            {
                // Key does not exist, add it with a value of 1
                QueryVec[labela] = 1;
            }
            else
            {
                // Key exists, increment its value
                QueryVec[labela]++;
            }
        }
        QueryNlabel.emplace_back(QueryVec);
        QueryVec.clear();
    }
}


/*Initial Pruning. After the CS creation we remove nodes that have
 *less edges that their CS.
 */
bool InitPrunTCSR(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph)
{
    int jj = 0;
    ui VDP;
    VertexID i = 0;
    VertexID rev;
    bool ret = false;
    for (VertexID kk = 0; kk < qsiz; kk++)
    {
        jj = FCS[kk].size();
        VDP = query_graph->getVertexDegree(kk);
        while (jj > 0)
        {
            jj--;
            if (FCS[kk][jj].Ichange == true)
            {
                if (FCS[kk][jj].edges.size() == 0)
                {
                    ret = true;
                    FCS[kk][jj].deleted = true;
                }
                // pruning rule
                else if (FCS[kk][jj].edges.size() < VDP)
                {
                    i = 0;

                    while (i < FCS[kk][jj].edges.size())
                    {
                        if (FCS[kk][jj].edges[i].first == 1000000)
                        {
                            i++;
                            continue;
                        }
                        rev = findIndBS(FCS, FCS[kk][jj].edges[i].second, FCS[kk][jj].edges[i].first); // vertex to remove ID?
                        FCS[FCS[kk][jj].edges[i].first][rev].Ichange = true;
                        for (int dd = 0; dd < FCS[FCS[kk][jj].edges[i].first][rev].edges.size(); dd++)
                            if (FCS[FCS[kk][jj].edges[i].first][rev].edges[dd].first == kk && FCS[FCS[kk][jj].edges[i].first][rev].edges[dd].second == FCS[kk][jj].ID)
                            {
                                FCS[FCS[kk][jj].edges[i].first][rev].edges[dd].first = 1000000;
                                i++;
                                break;
                            }
                    }
                    FCS[kk][jj].deleted = true;
                    ret = true;
                }

                FCS[kk][jj].Ichange = false;
            }
        }
    }

    return ret;
}


/*Add edges to the candidate space based on the paper rules.
 */

void EdgesCSBasicRL(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, ui *&flag, ui *&updated_flag)
{

    ui u_nbrs_count = 0;
    ui u_nbrs_countD = 0;
    int sizA = 0;
    int sizC;
    VertexID VID = 0;
    VertexID de = 0;
    VertexID cne = 0;
    VertexID labela = 0;
    // for CS. for every node of the CS(i)
    for (VertexID a = 0; a < qsiz; a++)
    { // take the neiboors of the query node FCS[a]
        const VertexID *u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        // VertexID* u_nbrs = query_graph->getVertexNeighbors(a, u_nbrs_count);
        sizA = FCS[a].size();
        // Now for every node of the neigboors of FCS[a]---
        for (VertexID c = 0; c < u_nbrs_count; c++)
        { // we start checking query nodes with smaller id to higher
            // thus is the query neigboor has smaller ID we already
            // added the edge to the CS.
            ui updated_flag_count = 0;
            if (u_nbrs[c] < a)
                continue;
            for (int nn = 0; nn < FCS[u_nbrs[c]].size(); nn++)
            {
                flag[FCS[u_nbrs[c]][nn].ID] = nn + 1;
                updated_flag[updated_flag_count++] = FCS[u_nbrs[c]][nn].ID;
            }

            cne = u_nbrs[c];
            // For every node of the CS(i)-> the candidates query node we evaluate
            // candidate vertex for a query a is FCS[a][b].ID
            for (VertexID b = 0; b < sizA; b++)
            {
                VID = FCS[a][b].ID;
                // const VertexID* u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); //real neigboors of the candidate vertex
                // get all the neigbors of the FCS[a][b].ID in the data graph
                const VertexID *u_nbrsD = data_graph->getVertexNeighbors(FCS[a][b].ID, u_nbrs_countD); // real neigboors of the candidate vertex
                                                                                                       // for every neigboor of the candidate vertex of the real graph
                                                                                                       // check if the node exists in the set of neigboors
                for (VertexID e = 0; e < u_nbrs_countD; e++)
                {
                    ui NID = u_nbrsD[e];
                    if (flag[NID] != 0)
                    {
                        FCS[a][b].edges.emplace_back(make_pair(cne, NID));
                        FCS[cne][flag[NID] - 1].edges.emplace_back(make_pair(a, FCS[a][b].ID));
                    }
                }
            }
            for (ui aa = 0; aa < updated_flag_count; ++aa)
            {
                flag[updated_flag[aa]] = 0;
            }
        }
    }
}

/* Add vertices that Pass NLF and Eigen rule to candidate Space.
 * OpenData1 need to be removed from here and pass the eigenVD1 as parameter(index)
 **310 is hardcodes max number of label ID-> can be just extracted from data graph->number of Labels.
 **I preload for every label all the possible candidate nodes, so I will not have to reload them if
 **a query node has the same label.
 **Important to note by construction the nodes are ordered by ID so FCS[0][0].ID<FCS[0][1].ID
 **1)Pruning 1 l(q)=l(v) So labelsNum[reverseLab[label]] already has nodes that have the same label
 **2)Pruning 2 d(q)<=d(v) Degree
 **3)Pruning 3 LE(q)<=LE(v) top S laplacian
 **4)Pruning 4 N(q)<=N(v) Neigboorhood labels
 */

void Vertices(vector<vector<CSV>> &FCS, int qsiz, int dsiz, Graph *data_graph, Graph *query_graph, float **&eigenVq1, vector<map<ui, int>> &QueryNlabel, float **&eigenVD1)
{
    VectorXd devalues;
    VectorXd qevalues;
    bool con = true;
    ui com = data_graph->getGraphMaxLabelFrequency();

    ui copies = query_graph->getLabelsCount() + 1;
    ui labelsNum[copies];
    ui kk;
    ui i;
    ui j;
    ui reverseLab[310];
    ui u_nbrs_countD = 0;
    const ui *labelData[copies];
    LabelID label = 0;
    for (i = 0; i < 310; i++)
        reverseLab[i] = 310;
    int pos = 0;
    for (i = 0; i < qsiz; i++)
    {
        label = query_graph->getVertexLabel(i);
        ui vdata_vertex_num = 0;
        if (reverseLab[label] == 310)
        {
            reverseLab[label] = pos;
            labelData[pos] = data_graph->getVerticesByLabel(label, vdata_vertex_num);
            labelsNum[pos] = vdata_vertex_num;
            pos++;
        }
    }
    ui reserveS;
    vector<CSV> CS;
    CS.reserve(com);
    FCS.reserve(qsiz);
    ui k = 0;
    ui degree = 0;
    ui data_vertex_num;

    int prunES = qsiz - 3;
    prunES = 30;
    // for every C(q)
    ui vdata_vertex_num = 0;
    for (i = 0; i < qsiz; i++)
    { 
        
        label = query_graph->getVertexLabel(i);
        degree = query_graph->getVertexDegree(i);
        data_vertex_num = 0;
        for (j = 0; j < labelsNum[reverseLab[label]]; ++j)
        {
            // get Vertex ID
            VertexID data_vertex = labelData[reverseLab[label]][j];
            if (data_graph->getVertexDegree(data_vertex) >= degree)
            {
                con = true;

                // Eigen Value Pruning up to pruneEs value
                for (kk = 0; kk < prunES; kk++)
                {
                    if (eigenVq1[i][kk] <= -1)
                        break;
                    if (eigenVD1[data_vertex][kk] < eigenVq1[i][kk])
                        // Rounding errors for eigenvalue
                        if ((eigenVq1[i][kk] - eigenVD1[data_vertex][kk]) > 0.0001)
                        {
                            con = false;
                            break;
                        }
                }

                if (con)
                {
                    // Neigborhood Check
                    for (auto it = QueryNlabel[i].begin(); it != QueryNlabel[i].end(); ++it)
                    {
                        data_graph->getNeighborsByLabelCount(data_vertex, it->first, u_nbrs_countD);
                        if (u_nbrs_countD < it->second)
                        {
                            con = false;
                            break;
                        }
                    }
                    // If all rules true -> add to candidate space
                    if (con)
                    {
                        CSV cat(data_vertex);
                        CS.emplace_back(cat);
                    }
                }
            }
        } 
        FCS.emplace_back(CS);
        CS.clear();
    }
}


/*Removed all nodes from the Candidate space that are set to be pruned.
**While iterating we set edges to 1000000 as a max value
**and notes to csv deleted. Then we gather all the nodes and edges to
** remove them all together to avoid reallocations and reorderings
*/
inline void clearWrong(vector<vector<CSV>> &FCS)
{
    for (auto &row : FCS)
    {
        row.erase(remove_if(row.begin(), row.end(), [&](CSV &csv)
                            {
        if (csv.deleted) {
            return true;
        }
        auto newEnd = remove_if(csv.edges.begin(), csv.edges.end(), [](const pair<VertexID, VertexID> &edge) {
            return edge.first == 1000000;
        });
        csv.edges.erase(newEnd, csv.edges.end());
        return false; }),
                  row.end());
    }
}

/*Prune Node and edges, set node to deleted and edge to 1000000.
**Also find all neigboors of node and set the edges to the pruned node to 1000000
**Checks also if the remove edge makes the node to be pruned by the Nedge rule.
*/

inline void removeVertexAndEgjesFK(vector<vector<CSV>> &FCS, int i, int deli)
{
    VertexID vx1;
    for (int j = 0; j < FCS[i][deli].edges.size(); j++)
    {
        if (FCS[i][deli].edges[j].first == 1000000)
            continue;
        vx1 = findIndBS(FCS, FCS[i][deli].edges[j].second, FCS[i][deli].edges[j].first);
        FCS[FCS[i][deli].edges[j].first][vx1].change = true;
        FCS[FCS[i][deli].edges[j].first][vx1].IPchange = true;
        FCS[FCS[i][deli].edges[j].first][vx1].Nedge[i]--;

        if (FCS[FCS[i][deli].edges[j].first][vx1].Nedge[i] == 0)
            FCS[FCS[i][deli].edges[j].first][vx1].NedgeC = true;
        for (int k = 0; k < FCS[FCS[i][deli].edges[j].first][vx1].edges.size(); k++)
        {
            if (FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first == 1000000)
                continue;
            if (FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first == i && FCS[FCS[i][deli].edges[j].first][vx1].edges[k].second == FCS[i][deli].ID)
            {
                FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first = 1000000;
                break;
            }
        }
    }

    FCS[i][deli].edges.clear();
    FCS[i][deli].deleted = true;
}
/*Prune Node and edges, set node to deleted and edge to 1000000.
**Also find all neigboors of node and set the edges to the pruned node to 1000000
**Checks also if the remove edge makes the node to be pruned by the Nedge rule.
*/
inline void removeVertexAndEgjesFKNP(vector<vector<CSV>> &FCS, int i, int deli)
{
    VertexID vx1;
    ui j;
    ui k;
    for (j = 0; j < FCS[i][deli].edges.size(); j++)
    {
        // BSCHange
        if (FCS[i][deli].edges[j].first == 1000000)
            continue;
        vx1 = findIndBS(FCS, FCS[i][deli].edges[j].second, FCS[i][deli].edges[j].first);
        FCS[FCS[i][deli].edges[j].first][vx1].IPchange = true;
        FCS[FCS[i][deli].edges[j].first][vx1].change = true;
        Mymutex.lock();
        FCS[FCS[i][deli].edges[j].first][vx1].Nedge[i]--;

        if (FCS[FCS[i][deli].edges[j].first][vx1].Nedge[i] == 0)
            FCS[FCS[i][deli].edges[j].first][vx1].NedgeC = true;
        Mymutex.unlock();
        for (k = 0; k < FCS[FCS[i][deli].edges[j].first][vx1].edges.size(); k++)
        {
            if (FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first == 1000000)
                continue;
            if (FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first == i && FCS[FCS[i][deli].edges[j].first][vx1].edges[k].second == FCS[i][deli].ID)
            {
                FCS[FCS[i][deli].edges[j].first][vx1].edges[k].first = 1000000;
                break;
            }
        }
    }

    FCS[i][deli].edges.clear();
    FCS[i][deli].deleted = true;
}

int SpectralMatching(int sizd, Graph *data_graph, Graph *query_graph, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix, float *&eigenQS)
{
    auto start = std::chrono::high_resolution_clock::now();
    ui **candidates1 = NULL;
    ui *candidates_count1 = NULL;
    int sizq = query_graph->getVerticesCount();
    ui Eprun = sizq - 3;
    Eprun = 30;
    MatrixXd eigenVq1(sizq, Eprun);
    int oMax = sizq * 3;
    oMax = 300;
    MTcalc12(query_graph, query_graph->getGraphMaxDegree(), eigenVq1, true, Eprun, oMax);
    float **eigenQ = NULL;
    eigenQ = new float *[sizq];
    for (ui i = 0; i < sizq; ++i)
    {
        eigenQ[i] = new float[Eprun];
        eigenQS[i] = 0;
        for (ui j = 0; j < Eprun; j++)
        {
            eigenQ[i][j] = eigenVq1(i, j);
            if (eigenQ[i][j] > 0)
                eigenQS[i] += eigenQ[i][j];
        }
    }
    return PILOS(data_graph, query_graph, eigenQ, twohop, candidates, candidates_count, EWeight, eigenVD1, alpha, beta, edge_matrix);
}

/*Assuming FCS[i][j].ID is sorted we cand find any j given IDC=i and IDQ=FCS[i][j].ID with Binary search.
**Keep in mind sorted properties have to be kept when removing elements.
*/
inline VertexID findIndBS(vector<vector<CSV>> &FCS, VertexID IDC, VertexID IDQ)
{
    int lo = 0, hi = FCS[IDQ].size() - 1;
    int mid;
    // This below check covers all cases , so need to check
    while (hi - lo > 1)
    {
        int mid = (hi + lo) / 2;
        if (FCS[IDQ][mid].ID < IDC)
        {
            lo = mid + 1;
        }
        else
        {
            hi = mid;
        }
    }
    if (FCS[IDQ][lo].ID == IDC)
    {
        return lo;
    }
    else if (FCS[IDQ][hi].ID == IDC)
    {
        return hi;
    }
    cout << "error Prob" << endl;
    cout << IDC << "IDC,IDQ" << IDQ << endl;
    return -10000;
}

/*Degree check for every node so we avoid computations in the future.
**Assuming that to be here the degree check is valid
**For every query vertex the node has to have at least 1 edge in
**a neigborhood candidate space. We start by adding all the count
**for every edge and every time we remove 1 edge we decrease the count.
**
**
*/
void fillEN(vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph)
{
    for (int i = 0; i < qsiz; i++)
    {
        ui de = query_graph->getVertexDegree(i);

        for (int j = 0; j < FCS[i].size(); j++)
        {
            FCS[i][j].Nedge = new int[qsiz];
            memset(FCS[i][j].Nedge, 0, sizeof(int) * qsiz);
            for (int kk = 0; kk < FCS[i][j].edges.size(); kk++)
            {
                FCS[i][j].Nedge[FCS[i][j].edges[kk].first]++;
            }
            ui sd = 0;
            for (int kk = 0; kk < qsiz; kk++)
            {
                if (FCS[i][j].Nedge[kk] != 0)
                    sd++;
            }
            // No needed as to be here it passed the Degree Check
            if (de > sd)
                FCS[i][j].NedgeC = true;
        }
    }
}

int PILOS(Graph *data_graph, Graph *query_graph, float **&eigenVq1, int twohop, ui **&candidates, ui *&candidates_count, float **&EWeight, float **&eigenVD1, int alpha, int beta, Edges ***edge_matrix)
{ 
    int MemSize = 0;
    int qsiz = query_graph->getVerticesCount();
    int dsiz = data_graph->getVerticesCount();
    vector<vector<CSV>> FCS;
    FCS.reserve(qsiz);
    vector<ui> DegreeK; // Discovered nodes for 2 hop
    vector<vector<pair<ui, int>>> QueryNlabel;
    vector<map<ui, int>> NLabel;  // Number of Labels 1hop
    vector<map<ui, int>> NLabel2; // Number of Labels 2hop
    // Exctract 1hop label information for query graph
    ui *flag = new ui[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), 0);
    ui *updated_flag = new ui[data_graph->getVerticesCount()];
    std::fill(updated_flag, updated_flag + data_graph->getVerticesCount(), 0);

    ui lb = query_graph->getLabelsCount();
    ExtractNImap(NLabel, query_graph, qsiz);

    // Extract 2hop label information for query graph
    int **VS = NULL;
    VS = new int *[qsiz];
    for (int aa = 0; aa < qsiz; aa++)
    {
        VS[aa] = new int[qsiz];
        for (int bb = 0; bb < qsiz; bb++)
        {
            VS[aa][bb] = -1;
        }
    }
    ExtractUI2h(DegreeK, NLabel2, query_graph, qsiz, VS);

    Vertices(FCS, qsiz, dsiz, data_graph, query_graph, eigenVq1, NLabel, eigenVD1);
    int count = 0;
    int Tcount = 0;
    // Add Edges between nodes in candidate space

    EdgesCSBasicRL(FCS, qsiz, dsiz, data_graph, query_graph, flag, updated_flag);

    // Get candidate nodes neigborhood information for fast pruning¨
    // Initial Pruning on Candidate Space
    //if (getValue1() > MemSize)
    //    MemSize = getValue1();
    while (InitPrunTCSR(FCS, qsiz, query_graph))
        clearWrong(FCS);

    fillEN(FCS, qsiz, query_graph);
    int GDegree = query_graph->getGraphMaxDegree();
    // Neigborhood Pruning
    for (int i = 0; i < qsiz; i++)
    {
        for (int j = 0; j < FCS[i].size(); j++)
        {
            std::sort(FCS[i][j].edges.begin(), FCS[i][j].edges.end(), [](const auto &a, const auto &b)
                      { return a.first < b.first; });
        }
    }
    int cc = 0;
    while (RefinementNV(NLabel, FCS, qsiz, query_graph, data_graph, GDegree, flag, updated_flag))
        clearWrong(FCS);
    ui mc = 3;
    if (twohop == 2)
    {
        while (RFNV(NLabel, NLabel2, FCS, qsiz, query_graph, eigenVq1, DegreeK, twohop, alpha, beta, flag, updated_flag, VS) && mc < 5)
        {
            mc++;
            clearWrong(FCS);
            while (RefinementNV(NLabel, FCS, qsiz, query_graph, data_graph, GDegree, flag, updated_flag))
            {

                clearWrong(FCS);
            }
        }
    }
    clearWrong(FCS);
    allocateBufferFCS1(FCS, query_graph, candidates, candidates_count, EWeight);
    int totalCand = 0;

    for (int i = 0; i < qsiz; i++)
    {
        for (int j = 0; j < FCS[i].size(); j++)
        {
            candidates[i][j] = FCS[i][j].ID;
        }
        candidates_count[i] = FCS[i].size();
        totalCand = candidates_count[i] + totalCand;
    }
    //if (getValue1() > MemSize)
    //    MemSize = getValue1();
    //cout << "MemSize" << MemSize / 1000 << endl;
    return totalCand;
}

/*Neighborhood NLF in candidate space.
 */

bool RefinementNV(vector<map<ui, int>> NLabel, vector<vector<CSV>> &FCS, int qsiz, Graph *query_graph, Graph *data_graph, ui GDegree, ui *&flag, ui *&updated_flag)
{
    bool returnhere = true;
    ui i;
    ui j;
    ui IDC;
    ui pos;
    ui NI;
    ui query_vertex_num = query_graph->getVerticesCount();
    ui data_vertex_num = data_graph->getVerticesCount();
    ui query_graph_max_degree = query_graph->getGraphMaxDegree();
    ui data_graph_max_degree = data_graph->getGraphMaxDegree();

    returnhere = false;
    int *flag1 = new int[data_graph->getVerticesCount()];
    std::fill(flag1, flag1 + data_graph->getVerticesCount(), -1);
    for (i = 0; i < qsiz; i++)
    {
        for (j = 0; j < FCS[i].size(); j++)
        {

            if (FCS[i][j].IPchange == false || FCS[i][j].deleted == true)
                continue;
            // 1 degree rule check

            if (!FCS[i][j].NedgeC)
                for (int ia = 0; ia < query_graph->getVerticesCount(); ia++)
                {

                    if (FCS[i][j].Nedge[ia] == 1 && NLabel[i][query_graph->getVertexLabel(ia)] > 1)
                    {
                        for (int aa = 0; aa < FCS[i][j].edges.size(); aa++)
                        {
                            if (ia == FCS[i][j].edges[aa].first)
                            {
                                IDC = FCS[i][j].edges[aa].second;
                                pos = aa;
                                break;
                            }
                        }
                        for (int aa = 0; aa < FCS[i][j].edges.size(); aa++)
                        {
                            if (FCS[i][j].edges[aa].second == IDC && aa != pos && FCS[i][j].edges[aa].first != 1000000)
                            {

                                ui queryN = FCS[i][j].edges[aa].first;
                                NI = findIndBS(FCS, IDC, queryN);
                                for (int ao = 0; ao < FCS[queryN][NI].edges.size(); ao++)
                                {
                                    if (FCS[queryN][NI].edges[ao].first == i && FCS[queryN][NI].edges[ao].second == FCS[i][j].ID)
                                    {
                                        // remove neigboor
                                        FCS[queryN][NI].edges[ao].first = 1000000;
                                        FCS[queryN][NI].Nedge[i]--;
                                        if (FCS[queryN][NI].Nedge[i] == 0)
                                            FCS[queryN][NI].NedgeC = true;
                                        FCS[queryN][NI].change = true;
                                        FCS[queryN][NI].IPchange = true;

                                        // remove edge
                                        FCS[i][j].Nedge[queryN]--;
                                        if (FCS[i][j].Nedge[queryN] == 0)
                                            FCS[i][j].NedgeC = true;
                                        FCS[i][j].edges[aa].first = 1000000;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }

            // Degree Check
            if (FCS[i][j].NedgeC)
            {
                removeVertexAndEgjesFKNP(FCS, i, j);
                returnhere = true;
            } // Neighborhood Check
            else if (!OneHopEigenVG(FCS[i][j], NLabel[i], query_graph, flag, updated_flag))
            {
                removeVertexAndEgjesFKNP(FCS, i, j);
                returnhere = true;
            }
            else
            {
                FCS[i][j].IPchange = false;
            }
        }
    }
    return returnhere;
}

/*Main function for Eigenvalue Pruning.
**Small Changes Added OMax and Omax2 First is limit for eigenvalues and second limit to twohop.
** We store the edges as triplets nd create a sparse Eigen MatrixXD
** We calculate eigenValues only if it passed the second pruning rule and the size of Matrix is less than oMax
*/

bool RFNV(vector<map<ui, int>> NLabel, vector<map<ui, int>> NLabel2, vector<vector<CSV>> &FCS,
          int qsiz, Graph *query_graph, float **&eigenVq1, vector<ui> DM, int twohop, int alpha, int beta, ui *&flag, ui *&updated_flag, int **&VS)
{
    vector<T> tripletList;
    std::map<int, int> count_uniques;
    std::set<std::pair<int, int>> seen;
    std::vector<Triplet<double>> unique_triplets;
    ui *flagq = new ui[query_graph->getVerticesCount()];
    std::fill(flagq, flagq + query_graph->getVerticesCount(), 0);
    ui *updated_flagq = new ui[query_graph->getVerticesCount()];
    std::fill(updated_flagq, updated_flagq + query_graph->getVerticesCount(), 0);
    int IDDLC[3] = {0, 0, 0};
    bool returnhere = false;
    VertexID vertex = 0;
    vector<VertexID> temp2;
    vector<pair<VertexID, VertexID>> q_curr;
    int Eprun = qsiz - 3;
    Eprun = 30;
    VertexID vertexDegree = 0;
    VertexID vertexlabel = 0;
    ui SIDDSize = 0;
    bool continueE = false;
    bool con = true;
    ui oMax;
    oMax = alpha;
    ui oMax2;
    float **LM = new float *[oMax + 1];
    ui *SIDN = new ui[qsiz];
    ui lb = query_graph->getLabelsCount();
    for (int i = 0; i <= oMax; i++)
    {
        LM[i] = new float[oMax + 1];
        memset(LM[i], 0, oMax + 1 * oMax + 1 * sizeof(float));
    }
    VectorXd evalues(Eprun);
    ui i;
    ui j;
    float sumD = 0;
    map<ui, int> NLabelT;
    map<ui, int> NLabelT1;

    for (int dd = 0; dd < qsiz; dd++)
    {
        i = dd;
        ui NDL = query_graph->getVertexDegree(i);

        for (j = 0; j < FCS[i].size(); j++)
        {
            if (FCS[i][j].deleted == true || FCS[i][j].change == false)
                continue;
            if (!FCS[i][j].NedgeC)
            {
                tripletList.clear();
                q_curr.clear();
                IDDLC[0] = 0;
                NLabelT = NLabel2[i];
                NLabelT1 = NLabel[i];
                IDDLC[1] = NLabel2[i].size();
                IDDLC[2] = NLabel[i].size();
                for (int aa = 0; aa < oMax; aa++)
                {
                    for (int bb = 0; bb < oMax; bb++)
                    {
                        LM[aa][bb] = 0;
                    }
                }
                for (int aa = 0; aa < qsiz; aa++)
                {
                    flagq[aa] = 0;
                }
                flagq[i] = 1;
                SIDDSize = 0;
                OHEPM(FCS[i][j], flag, flagq, updated_flag, updated_flagq, IDDLC, LM, NLabelT, NLabelT1,
                      query_graph, q_curr, oMax);

                if (IDDLC[0] <= oMax)
                {

                    if (beta == 0)
                        SHEigen(q_curr, flag, flagq, updated_flag, updated_flagq, NLabelT, IDDLC, FCS, LM, query_graph, oMax, i, FCS[i][j]);
                    else
                    {
                        SHEigenb(q_curr, flag, flagq, updated_flag, updated_flagq, NLabelT, IDDLC, FCS, LM, query_graph, oMax, i, FCS[i][j], beta);
                    }
                    for (int aa = 0; aa < qsiz; aa++)
                        if (flagq[aa] == 1)
                            SIDDSize++;
                }
                if (IDDLC[0] <= oMax)
                    if ((IDDLC[0]) < DM[i] || IDDLC[1] > 0 || SIDDSize < DM[i])
                    { 
                        removeVertexAndEgjesFK(FCS, i, j);
                        returnhere = true;
                        for (int aa = 0; aa <= IDDLC[0]; aa++)
                        {
                            flag[updated_flag[aa]] = 0;
                        }
                        IDDLC[0] = oMax + 1;
                    }
                ui s2 = IDDLC[0];
                ui count2 = 0;

                if (s2 <= oMax && NDL != 1)

                {
                    int largest = 0;
                    con = true;
                    bool endFast = false;
                    ui count1 = 0;
                    for (int l = 0; l < s2; l++)
                    {
                        count1 = 0;
                        for (int m = 0; m < s2; m++)
                        {
                            if (LM[l][m] == -1)
                            {
                                count1++;

                                tripletList.emplace_back(T(l, m, -1));
                            }
                        }
                        tripletList.emplace_back(T(l, l, count1));
                    }
                    int aa = 0;

                    if (tripletList.size() == s2 * s2)
                    {
                        evalues.resize(Eprun);

                        for (ui ss = 0; ss < Eprun; ss++)
                        {
                            if (ss < s2)
                            {
                                evalues(ss) = s2;
                            }

                            else if (ss == IDDLC[0] - 1)
                                evalues(ss) = 0;
                            else
                                evalues(ss) = 0;
                        }
                    }
                    else
                    {
                        SparseMatrix<double> M(s2, s2);
                        M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                                          { return b; });
                        M.makeCompressed();

                        calcEigens1(M, Eprun, evalues, s2);
                    }
                    con = true;
                    sumD = 0;

                    for (int dd = 0; dd < Eprun; dd++)
                    {
                        if (eigenVq1[i][dd] <= -1)
                            break;
                        if (evalues[dd] < eigenVq1[i][dd])
                        {
                            if ((eigenVq1[i][dd] - evalues[dd]) > 0.0001)
                            {
                                con = false;
                                break;
                            }
                        }
                        // Eigen Ordering If we want to eigenvalues uncomment.
                        // Add the values up
                        else
                            sumD += evalues[dd];
                    }
                    if (!con)
                    {

                        removeVertexAndEgjesFK(FCS, i, j);
                        returnhere = true;
                        for (int aa = 0; aa <= IDDLC[0]; aa++)
                        {
                            flag[updated_flag[aa]] = 0;
                        }
                    }

                    else
                    {
                        FCS[i][j].change = false;
                        // Eigen Ordering
                        FCS[i][j].ED = sumD;
                        for (int aa = 0; aa <= IDDLC[0]; aa++)
                        {
                            flag[updated_flag[aa]] = 0;
                        }
                    }
                }
                else
                {
                    FCS[i][j].change = false;
                    for (int aa = 0; aa <= IDDLC[0]; aa++)
                    {
                        flag[updated_flag[aa]] = 0;
                    }
                }
            }
            else
            {
                removeVertexAndEgjesFK(FCS, i, j);
                returnhere = true;
            }
        }
    }
    return returnhere;
}
