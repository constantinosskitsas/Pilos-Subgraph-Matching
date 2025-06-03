#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/SymEigsShiftSolver.h>

#include <iostream>
#include <cmath>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Eigenvalues>
// #include "GrM.h"
#include "eigenHelper.h"
#include "thread"
#include <queue>
#include <mutex>
std::mutex mtxx;
typedef Eigen::Triplet<double> T;
using namespace Eigen;
using namespace std;
using namespace Spectra;
void calcEigens12(SparseMatrix<double> M, int k, VectorXd &evalues, int count)
{
    // int sizek = k*2;
    int dek = k;

    if (count == 0 || count == 1)
    {
        evalues.resize(dek);

        evalues(0) = 1;
        for (int i = 1; i < dek; i++)
            evalues(i) = 0;

        return;
    }

    if (k >= count)
        k = count - 1;

    // SparseGenRealShiftSolver<double> op(M);
    // SymEigsSolver<SparseGenRealShiftSolver<double>> eigs(op, k, count);
    //  SparseGenMatProd<double> op(Lpl);
    // GenEigsSolver< double, LARGEST_MAGN, SparseGenMatProd<double> > eigs(&op, 3, 6);
    // SparseGenMatProd<double> op(M);
    SparseGenMatProd<double> op(M);
    SymEigsSolver<SparseGenMatProd<double>> eigs(op, k, count);

    eigs.init();
    int nconv = eigs.compute(SortRule::LargestAlge);
    // int nconv = eigs.compute(SortRule::SmalestAlge);
    if (eigs.info() == CompInfo::Successful)
        evalues = eigs.eigenvalues();
    if (eigs.info() != CompInfo::Successful)
    {
        int info = static_cast<int>(eigs.info());
        while (true)
            cout << "problem No Eigs" << endl;
        cout << info;
        return;
        // cout<<count<<endl;
        // cout<<k<<endl;
        // printSM(M);
        evalues.resize(k);
        for (int ss = 0; ss < k - 1; ss++)
            evalues(ss) = count;
        evalues(k - 1) = 0;
    }

    if (evalues.size() < dek)
    {
        int sz = evalues.size();
        evalues.conservativeResize(dek);
        evalues(sz) = 0;
        sz++;
        for (int i = sz; i < dek; i++)
            // evalues(i) = -1;
            evalues(i) = 0;
    }
}
void calcEigensEigenLib(SparseMatrix<double> M, int k, VectorXd &evalues, int count)
{
    int sizek = k * 2;
    int dek = k;
    MatrixXd dMat;
    dMat = M;
    // SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(M);
    // SelfAdjointEigenSolver<MatrixXd<double>> eigensolver(dMat);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(dMat);
    // setNbThreads(4);
    if (eigensolver.info() != Success)
    {
        std::cerr << "Eigenvalue computation failed!" << std::endl;
        while (true)
            cout << "jjjj" << endl;
        return;
    }
    VectorXd eigenvalues = eigensolver.eigenvalues();

    if (eigenvalues.size() < k)
    {
        evalues = eigenvalues.tail(eigenvalues.size()).reverse();
        int sz = evalues.size();
        evalues.conservativeResize(dek);
        evalues(sz) = 0;
        sz++;
        for (int i = sz; i < dek; i++)
            // evalues(i) = -1;
            evalues(i) = 0;
    }
    else
    {
        evalues = eigenvalues.tail(k).reverse();
    }
}
bool isLaplacianMatrix(const Eigen::SparseMatrix<double> &matrix)
{
    if (matrix.rows() != matrix.cols())
    {
        cout << "not square" << endl;
        return false; // Matrix should be square
    }

    int size = matrix.rows();

    // Check if the matrix is symmetric
    if (!(matrix.transpose().isApprox(matrix)))
    {
        cout << "not symmetric" << endl;
        return false;
    }

    // Check if the matrix satisfies the Laplacian properties
    for (int i = 0; i < size; ++i)
    {
        double diagonal = matrix.coeff(i, i);
        double sumOfRow = 0.0;

        for (Eigen::SparseMatrix<double>::InnerIterator it(matrix, i); it; ++it)
        {
            if (it.row() != it.col())
            {
                sumOfRow += it.value();
            }
        }

        if (diagonal != abs(sumOfRow))
        {
            cout << "Wrong Sum" << endl;
            return false;
        }
    }

    return true;
}
void printSM(SparseMatrix<double> &M);
// bool checkM2(SparseMatrix<double> mat);
bool checkM2(SparseMatrix<double> mat, int rowO, int ColO);
bool checkM(SparseMatrix<double> mat)
{
    int count = 0;
    int row = 0;
    int col = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {

        count = 0;
        for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            count = count + it.value();
            if (it.value() < -1)
                cout << "prob" << endl;
            row = it.row(); // row index
            col = it.col(); // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
            if (checkM2(mat, row, col) == false)
            {
                cout << "Problem" << endl;
            }
        }
        if (count != 0)
        {
            //            "Problem";
            return true;
        }
    }
    return false;
}

bool checkM2(SparseMatrix<double> mat, int rowO, int ColO)
{
    int row = 0;
    int col = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
    {
        // count=0;
        for (SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
            // count=count+it.value();

            row = it.row(); // row index
            col = it.col(); // col index (here it is equal to k)
            // it.row();   // row index
            // it.col();   // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
            if (row == ColO && col == rowO)
                return true;
        }
    }
    return false;
}
VertexID checkANX(vector<VertexID> ID, VertexID CID)
{
    for (int i = 0; i < ID.size(); i++)
        if (ID[i] == CID)
            return i;
    return 1000000;
}

VertexID checkA(VertexID *ID, VertexID vertex, int count)
{
    for (int i = 0; i <= count; i++)
    {
        if (ID[i] == vertex)
            return i;
    }
    return 1000000;
}

void printSM(SparseMatrix<double> &M)
{
    int count = 0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
        { // cout<<"query ID"<<i<<endl;
          // cout<<"test"<<endl;
            count = count + it.value();
            if (it.value() == -1)
            {

                // cout << "";
                // cout << " (" << it.row() << ",";
                // cout << it.col() << ")  ";
                // cout << "e " << it.col() << " " << it.row() << " "
                //     << "0" << endl;
                cout << "v" << it.row() << " " << it.col();
                // cout<<"v "<< it.col()<<" "<<"0";
                cout << " " << it.value() << endl;
            }
            // row index
            // col index (here it is equal to k)
            // it.index(); // inner index, here it is equal to it.row()
        } // cout<<" "<<endl;
        if (count != 0)
        {
            cout << "problem" << endl;
        }
        count = 0;
        // cout<<endl;
    }
    cout << " " << endl;
}

void ExtractAdL(SparseMatrix<double> &M, Graph *data_graph, int degree, int depth, VertexID vertex)
{
    VertexID *neighbors_;
    VertexID *ID; // add size then resize
    VertexID *IDL;
    ID = new VertexID[data_graph->getVerticesCount()];
    IDL = new VertexID[data_graph->getVerticesCount()];
    memset(ID, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    memset(IDL, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    int count = 0;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

    ID[0] = vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    M.insert(0, 0) = (float)u_nbrs_count;
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        count++;
        M.insert(0, count) = -1.0;
        M.insert(count, 0) = -1.0;
        ID[count] = u_nbrs[j];
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            vx1 = checkA(ID, vertex, count);
            M.insert(vx1, vx1) = (float)u_nbrs_count;
            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];

                vx2 = checkA(ID, vertexpair, count);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID[count] = vertexpair;
                    M.insert(vx1, vx2) = -1.0;
                    M.insert(vx2, vx1) = -1.0;
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    M.coeffRef(vx1, vx2) = -1.0;
                    M.coeffRef(vx2, vx1) = -1.0;
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = checkA(ID, vertex, count);

        M.insert(vx1, vx1) = IDL[vx1];
    }
}

void calcEigens1(SparseMatrix<double> M, int k, VectorXd &evalues, int count)
{
    int sizek = k * 2;
    int dek = k;
    //
    // M.makeCompressed();
    // if(!M.isApprox(M.transpose()))
    // cout<<"really?"<<endl;

    SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(M);

    // EigenSolver<SparseMatrix<double>> eigensolver(M);
    if (eigensolver.info() != Success)
    {
        std::cerr << "Eigenvalue computation failed!" << std::endl;
        return;
    }
    VectorXd eigenvalues = eigensolver.eigenvalues();

    if (eigenvalues.size() < k)
    {
        ui ts = eigenvalues.size();
        evalues = eigenvalues.tail(eigenvalues.size()).reverse();
        // cout<<evalues<<endl;
        // cout<<"count "<<count<<" NE "<<evalues.size();
        int sz = evalues.size();
        evalues.conservativeResize(dek);
        if (ts == count)
            evalues(sz) = -1;
        else
            evalues(sz) = 0;
        sz++;
        for (int i = sz; i < dek; i++)
            evalues(i) = -1;
        // evalues(i) = 0;
        // cout<<evalues<<endl;
    }
    else
    {
        evalues = eigenvalues.tail(k).reverse();
    }
}

void CompactADLEIGSet(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    // vector<VertexID> IDL;

    VertexID *neighbors_;
    unordered_map<ui, ui> SID;
    VertexID *IDL;

    IDL = new VertexID[data_graph->getVerticesCount()];

    VertexID vx1;
    VertexID vx2 = 0;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    ui count = 1;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = data_graph->getVerticesCount() - 1;
    if (k > 10)
        k = Eprun;
    if (k < 10)
        k = 10;
    if (u_nbrs_count > 150)
    {

        evalues.resize(k);
        // evalues<<500,500,500;
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    SID.insert({vertex, 0});
    IDL[0] = 1;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    tripletList.push_back(T(0, 0, (float)u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            cout << "error :" << vertex << endl;

        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        SID.insert({u_nbrs[j], count});

        IDL[count] = 1;
        count++;
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            if (u_nbrs_count > 150)
            {

                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }

            auto it = SID.find(vertex);
            vx1 = it->second;

            tripletList.push_back(T(vx1, vx1, (float)u_nbrs_count));

            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;

                auto [it1, success] = SID.emplace(vertexpair, count);

                if (success)
                {
                    vx2 = count;
                    q_next.push(vertexpair);
                    IDL[count] = 1;
                    count++;
                }
                else
                {
                    vx2 = SID[vertexpair];
                    if (i == depth - 1)
                        IDL[vx2]++;
                }
                tripletList.push_back(T(vx1, vx2, -1));
                tripletList.push_back(T(vx2, vx1, -1));
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = SID[vertex];
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }

    SparseMatrix<double> M(count, count);
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });
    M.makeCompressed();
    if (!isLaplacianMatrix(M))
        cout << "wrong Laplacian Matrix-Malakas" << endl;
    /*if(checkM(M)){
                   cout<<"wrong Laplacian Matrix-Malakas"<<endl;
               }

   */

    calcEigens1(M, Eprun, evalues, count);
    tripletList.clear();
    ID.clear();
    SID.clear();
    // delete[] IDL;
    // IDL.clear();
    // delete[]  neighbors_;
    // delete[]
    // delete[]
    // delete[] u_nbrs;
}
void CompactADLEIGNV(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun, int oMax)
{
    vector<T> tripletList;
    // vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    VertexID vx1;
    VertexID vx2;
    unordered_map<int, int> ID;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 0;
    int qsiz = data_graph->getVerticesCount();
    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = 30;
    int **LM = new int *[qsiz];
    ui *SIDN = new ui[qsiz];
    for (int i = 0; i < oMax; i++)
    {
        LM[i] = new int[oMax];
        memset(LM[i], 0, oMax * oMax * sizeof(int));
    }
    ID.insert({vertex, 0});
    count++;
    // ID[0]=vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            cout << "error :" << vertex << endl;

        LM[0][count] = -1;
        LM[count][0] = -1;
        ID.insert({u_nbrs[j], count});
        q_curr.push(u_nbrs[j]);
        count++;
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            // vx1 = checkANX(ID, vertex);
            vx1 = ID[vertex];

            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;
                auto result = ID.insert({vertexpair, count});

                if (result.second)
                {

                    vx2 = count;
                    LM[vx1][vx2] = -1;
                    LM[vx2][vx1] = -1;
                    count++;
                }
                else
                {
                    vx2 = ID[vertexpair];
                    LM[vx1][vx2] = -1;
                    LM[vx2][vx1] = -1;
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }

    ui cunt = 0;
    for (int wish = 0; wish < qsiz; wish++)
    {
        cunt = 0;
        for (int best = 0; best < qsiz; best++)
        {
            if (LM[wish][best] == -1)
            {
                cunt++;
                tripletList.emplace_back(T(wish, best, -1));
            }
        }
        tripletList.emplace_back(T(wish, wish, cunt));
    }
    if (tripletList.size() == count * count)
    {
        evalues.resize(k);

        for (int ss = 0; ss < k; ss++)
        { // count-1 or count?
            if (ss < count)
                evalues(ss) = count;
            // evalues(ss) = count - 1;
            else if (ss == count)
                evalues(ss) = 0;
            else
                //
                // evalues(ss) = -1;
                evalues(ss) = 0;
        }
    }
    else
    {
        SparseMatrix<double> M(count, count);
        M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                          { return b; });
        // checkM(M);

        calcEigens1(M, k, evalues, count);
    }
    tripletList.clear();
    ID.clear();
    IDL.clear();
}

void CompactADLEIGImpro(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    unordered_map<ui, ui> SID;

    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 1;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = 16;
    if (u_nbrs_count > 150)
    {
        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    SID.insert({vertex, count});
    IDL.push_back(1);
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;

    tripletList.push_back(T(0, 0, u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            continue;

        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        SID.insert({u_nbrs[j], count});
        IDL.push_back(1);
        count++;
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        { // or just counter and clear
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            if (u_nbrs_count > 450 || count > 4500)
            {
                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                // evalues<<500,500,500,500;
                return;
            }

            vx1 = SID[vertex];
            tripletList.push_back(T(vx1, vx1, u_nbrs_count));

            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];

                auto result = SID.insert({vertexpair, count});
                if (result.second)
                {

                    vx2 = count;
                    q_next.push(vertexpair);
                    IDL.push_back(1);
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    count++;
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    vx2 = result.first->second;
                    if (vx2 == vx1)
                    {
                        // cout<<"hola"<<endl;
                        continue;
                    }
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }

        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = SID[vertex];
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }

    SparseMatrix<double> M(count, count);
    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });

    /*if(checkM(M)){
                   cout<<"wrong Laplacian Matrix-Malakas"<<endl;
               }

   */

    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    ID.clear();
    IDL.clear();
}

void MTcalc12(Graph *data_graph, int degree, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
{

    auto AdJAdl1 = [](Graph *data_graph, int degree, VertexID svertex, VertexID evertex, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
    {
        VectorXd evalues;
        int depth = 2;
        for (int i = svertex; i < evertex; i++)
        {

            if (LE)
                CompactADLEIGQS(data_graph, degree, evalues, i, depth, Eprun, oMax);
            // CompactADLEIGNV(data_graph, degree, evalues, i, depth, Eprun,oMax);

            else
                CompactADJEIG(data_graph, degree, evalues, i, depth);

            eigenVD.row(i) = evalues;
            evalues.setZero();
        }
    };
    VectorXd evalues;
    int Tnum = 1;
    int siz = data_graph->getVerticesCount();
    if (Tnum == 1)
    {

        for (int i = 0; i < siz; i++)
        {
            CompactADLEIG(data_graph, degree, evalues, i, 2, Eprun, oMax);
            eigenVD.row(i) = evalues;
            evalues.setZero();
        }
    }
    else
    {
        int div = (int)(siz / Tnum);
        thread th[Tnum];
        for (int d = 0; d < Tnum - 1; d++)
        {
            th[d] = thread(AdJAdl1, data_graph, degree, div * d, div * (d + 1), ref(eigenVD), LE, Eprun, oMax);
        }
        th[Tnum - 1] = thread(AdJAdl1, data_graph, degree, div * (Tnum - 1), siz, ref(eigenVD), LE, Eprun, oMax);
        for (int d = 0; d < Tnum; d++)
            th[d].join();
    }
}

void MTcalc12A(Graph *data_graph, int degree, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
{

    auto AdJAdl1 = [](int *pos1, Graph *data_graph, int degree, VertexID d, VertexID siz, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
    {
        VectorXd evalues;
        int depth = 2;
        int svertex;
        int evertex;
        int opd = d;
        ui num = 1000000;
        while ((opd * num) < siz)
        {
            svertex = opd * num;
            evertex = svertex + num;
            if (evertex > siz)
                evertex = siz;
            for (int i = svertex; i < evertex; i++)
            {
                CompactADLEIG(data_graph, degree, evalues, i, depth, Eprun, oMax);
                eigenVD.row(i) = evalues;
                evalues.setZero();
            }
            mtxx.lock();
            (*pos1)++;
            opd = *pos1;
            cout << opd << endl;
            mtxx.unlock();
        }
    };
    VectorXd evalues;
    int Tnum = 1;
    int pos = Tnum - 1;

    int *pos1 = &pos;
    int siz = data_graph->getVerticesCount();
    int div = (int)(siz / Tnum);
    thread th[Tnum];
    for (int d = 0; d < Tnum; d++)
    {
        th[d] = thread(AdJAdl1, pos1, data_graph, degree, d, siz, ref(eigenVD), LE, Eprun, oMax);
    }
    for (int d = 0; d < Tnum; d++)
        th[d].join();
}


void MTcalc13A(Graph *data_graph, int degree, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
{

    auto AdJAdl1 = [](int *pos1, Graph *data_graph, int degree, VertexID d, VertexID siz, MatrixXd &eigenVD, bool LE, int Eprun, int oMax)
    {
        VectorXd evalues;
        int depth = 2;
        int svertex;
        int evertex;
        int opd = d;
        ui num = 1000000;
        while ((opd * num) < siz)
        {
            svertex = opd * num;
            evertex = svertex + num;
            if (evertex > siz)
                evertex = siz;
            CompactADLEIG_BOOST_twoH(data_graph,svertex,evertex,Eprun,oMax,eigenVD);
                //eigenVD.row(i) = evalues;
               // evalues.setZero();
            mtxx.lock();
            (*pos1)++;
            opd = *pos1;
            //cout << opd << endl;
            mtxx.unlock();
        }
    };
    VectorXd evalues;
    int Tnum = 1;
    int pos = Tnum - 1;

    int *pos1 = &pos;
    int siz = data_graph->getVerticesCount();
    int div = (int)(siz / Tnum);
    thread th[Tnum];
    for (int d = 0; d < Tnum; d++)
    {
        th[d] = thread(AdJAdl1, pos1, data_graph, degree, d, siz, ref(eigenVD), LE, Eprun, oMax);
    }
    for (int d = 0; d < Tnum; d++)
        th[d].join();
}

void CompactADJEIG(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth)
{
    vector<T> tripletList;

    VertexID *neighbors_;
    VertexID *ID; // add size then resize
    VertexID *IDL;
    ID = new VertexID[data_graph->getVerticesCount()];
    IDL = new VertexID[data_graph->getVerticesCount()];
    memset(ID, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    memset(IDL, 0, sizeof(VertexID) * (data_graph->getVerticesCount()));
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    int count = 0;
    int k = 30;
    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

    ID[0] = vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        count++;
        tripletList.push_back(T(0, count, 1));
        tripletList.push_back(T(count, 0, 1));
        ID[count] = u_nbrs[j];
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);

            vx1 = checkA(ID, vertex, count);

            for (ui j = 0; j < u_nbrs_count; ++j)
            {
                vertexpair = u_nbrs[j];
                vx2 = checkA(ID, vertexpair, count);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID[count] = vertexpair;
                    tripletList.push_back(T(vx1, vx2, 1));
                    tripletList.push_back(T(vx2, vx1, 1));
                }
                else
                {
                    tripletList.push_back(T(vx1, vx2, 1));
                    tripletList.push_back(T(vx2, vx1, 1));
                }
            }
        }

        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }

    SparseMatrix<double> M(count, count);

    M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                      { return b; });
    if (checkM(M))
    {
        cout << "wrong Laplacian Matrix-Malakas" << endl;
    }

    calcEigens1(M, k, evalues, count);
    tripletList.clear();
    delete[] neighbors_;
    delete[] ID;
    delete[] IDL;
}
void CompactADLEIGQS(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun, int oMax)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 0;

    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = data_graph->getVerticesCount() - 1;
    if (k >= 10)
        k = Eprun;
    if (k <= 15)
        k = 10;

    k = Eprun;
    if (u_nbrs_count > oMax)
    {

        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    ID.push_back(vertex);
    IDL.push_back(1);
    // ID[0]=vertex;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    // cout<<u_nbrs_count<<" NC "<<endl;
    // cout<<vertex<<" vertex "<<endl;
    tripletList.push_back(T(0, 0, (float)u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
        if (u_nbrs[j] == vertex)
            cout << "error :" << vertex << endl;
        count++;
        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        ID.push_back(u_nbrs[j]);
        IDL.push_back(1);
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            if (u_nbrs_count > oMax)
            {

                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }
            vx1 = checkANX(ID, vertex);
            tripletList.push_back(T(vx1, vx1, (float)u_nbrs_count));

            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;

                vx2 = checkANX(ID, vertexpair);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID.push_back(vertexpair);
                    IDL.push_back(1);
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
                }
                else
                {
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = checkANX(ID, vertex);
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }
    count++;
    map<int, int> count_uniques;
    set<pair<int, int>> seen;
    vector<Triplet<double>> unique_triplets;
    for (auto t : tripletList)
    {
        if (seen.count({t.row(), t.col()}) == 0)
        {
            unique_triplets.push_back(Triplet<double>(t.row(), t.col(), t.value()));
            seen.insert({t.row(), t.col()});
            count_uniques[t.row()]++;
        }
    }
    tripletList = unique_triplets;
    if (tripletList.size() == count * count)
    {
        // cout << "MEga prob" << endl;
        evalues.resize(k);

        for (int ss = 0; ss < k; ss++)
        {
            if (ss < count)
                // evalues(ss)=count-1;
                evalues(ss) = count;
            else if (ss == count)
                evalues(ss) = 0;
            else // evalues(ss)=-1;
                evalues(ss) = -1;
        }
    }
    else
    {
        SparseMatrix<double> M(count, count);
        M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                          { return b; });
        // checkM(M);

        calcEigens1(M, k, evalues, count);
    }
    tripletList.clear();
    ID.clear();
    IDL.clear();
}



void CompactADLEIG(Graph *data_graph, int degree, VectorXd &evalues, VertexID vertex, int depth, int Eprun, int oMax)
{
    vector<T> tripletList;
    vector<VertexID> ID; // add size then resize
    vector<VertexID> IDD;
    vector<VertexID> IDL;
    VertexID *neighbors_;
    VertexID vx1;
    VertexID vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    vertexprint = vertex;
    int count = 0;
    ui u_nbrs_count;
    ui u_nbrs_count1;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
    int k = Eprun;
    if (u_nbrs_count > oMax)
    {
        evalues.resize(k);
        for (int i = 0; i < k; i++)
            evalues(i) = 500;
        return;
    }
    ID.push_back(vertex);
    IDL.push_back(1);
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    tripletList.push_back(T(0, 0, (float)u_nbrs_count));
    for (int j = 0; j < u_nbrs_count; ++j)
    {
       // if (u_nbrs[j] == vertex)
        //    cout << "error :" << vertex << endl;
        count++;
        tripletList.push_back(T(0, count, -1));
        tripletList.push_back(T(count, 0, -1));
        ID.push_back(u_nbrs[j]);
        IDL.push_back(1);
        q_curr.push(u_nbrs[j]);
    }

    for (int i = 1; i < depth; i++)
    {
        while (!q_curr.empty())
        {
            vertex = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
            if (u_nbrs_count > oMax)
            {
                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }
            vx1 = checkANX(ID, vertex);
            tripletList.push_back(T(vx1, vx1, (float)u_nbrs_count));
            for (ui j = 0; j < u_nbrs_count; ++j)
            {

                vertexpair = u_nbrs[j];
                if (vertexpair == vertex)
                    cout << "problem again!!" << endl;

                vx2 = checkANX(ID, vertexpair);
                if (vx2 == 1000000)
                {
                    count++;
                    vx2 = count;
                    q_next.push(vertexpair);
                    ID.push_back(vertexpair);
                    IDL.push_back(1);
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2] = 1;
                    }
            if (count > oMax)
            {
                evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500;
                return;
            }
                }
                else
                {
                    tripletList.push_back(T(vx1, vx2, -1));
                    tripletList.push_back(T(vx2, vx1, -1));
                    if (i == depth - 1)
                    {
                        IDL[vx2]++;
                    }
                }
            }
        }
        if (!q_next.empty())
            q_curr.swap(q_next);
        else
        {
            i = depth;
        }
    }
    while (!q_curr.empty())
    {
        vertex = q_curr.front();
        q_curr.pop();
        vx1 = checkANX(ID, vertex);
        tripletList.push_back(T(vx1, vx1, IDL[vx1]));
    }
    count++;
    map<int, int> count_uniques;
    set<pair<int, int>> seen;
    vector<Triplet<double>> unique_triplets;
    for (auto t : tripletList)
    {
        if (seen.count({t.row(), t.col()}) == 0)
        {
            unique_triplets.push_back(Triplet<double>(t.row(), t.col(), t.value()));
            seen.insert({t.row(), t.col()});
            count_uniques[t.row()]++;
        }
    }
    tripletList = unique_triplets;
    if (tripletList.size() == count * count)
    {
        // cout << "MEga prob" << endl;
        evalues.resize(k);

        for (int ss = 0; ss < k; ss++)
        {
            if (ss < count)
                // evalues(ss)=count-1;
                evalues(ss) = count;
            else if (ss == count)
                evalues(ss) = 0;
            else // evalues(ss)=-1;
                evalues(ss) = 0;
        }
    }
    else
    {
        SparseMatrix<double> M(count, count);
        M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                          { return b; });
        // checkM(M);

       calcEigens1(M, k, evalues, count);
        //calcEigens12(M,k,evalues,count);
    }
    tripletList.clear();
    ID.clear();
    IDL.clear();
}


void CompactADLEIG_BOOST_twoH(Graph *data_graph, int sv,int ev, int Eprun, int oMax,MatrixXd &eigenVD)
{   
    VertexID vertex;
    VectorXd evalues;
    int *flag = new int[data_graph->getVerticesCount()];
    std::fill(flag, flag + data_graph->getVerticesCount(), -1);
    int *updated_flag = new int[data_graph->getVerticesCount()];
    std::fill(updated_flag, updated_flag + data_graph->getVerticesCount(), -1);
    
    
    
    int **LM = new int *[oMax];
    for (int i = 0; i < oMax; i++)
    {
        LM[i] = new int[oMax];
        memset(LM[i], 0, oMax  * sizeof(int));
    }
    
    vector<T> tripletList;
    vector<VertexID> IDD;
    VertexID *neighbors_;
    int vx1;
    int vx2;
    VertexID vertexpair;
    VertexID vertexprint;
    queue<VertexID> q_curr;
    queue<VertexID> q_next;
    queue<VertexID> q_temp;
    vertexprint = vertex;
    int count = 0;
    ui u_nbrs_count;
    ui u_nbrs_count1;
    int k = Eprun;
    
for (int vertex=sv;vertex<ev;vertex++){
    count=0;
    flag[vertex]=0;
    updated_flag[0]=vertex;
    const VertexID *u_nbrs = data_graph->getVertexNeighbors(vertex, u_nbrs_count);
        if (u_nbrs_count >= oMax-1)
    {count++;
        goto TerEarly;
    }
    for (int j = 0; j < u_nbrs_count; j++)
    {
        count++;
        LM[0][count]=-1;
        LM[count][0]=-1;
        flag[u_nbrs[j]]=count;
        updated_flag[count]=u_nbrs[j];
        q_curr.push(u_nbrs[j]);
    }
        while (!q_curr.empty())
        {
           int vertex1 = q_curr.front();
            q_curr.pop();
            u_nbrs = data_graph->getVertexNeighbors(vertex1, u_nbrs_count);
            if (u_nbrs_count >= oMax)
            {   count++;
                goto TerEarly;
            }
            vx1 = flag[vertex1];
            for (ui j = 0; j < u_nbrs_count; j++)
            {
                vertexpair = u_nbrs[j];
                vx2 = flag[vertexpair];
                if (vx2 == -1)
                {   
                    count++;
                    vx2 = count;
                    if (count >= oMax){
                        goto TerEarly; }
                    flag[vertexpair]=count;
                    updated_flag[count]=vertexpair;
                    LM[vx1][vx2]=-1;
                    LM[vx2][vx1]=-1;

                }
                else
                {
                    LM[vx1][vx2]=-1;
                    LM[vx2][vx1]=-1;
                }
            }
        }
    if (count>oMax){
    TerEarly:
       std::queue<VertexID> empty;
        std::swap(q_curr, empty);
       evalues.resize(k);
                for (int tt = 0; tt < k; tt++)
                    evalues(tt) = 500; 
        int s2=count;
        for (int l = 0; l < s2; l++)
    {    
        flag[updated_flag[l]]=-1;
        updated_flag[l]=-1;
        int count1 = 0;
        
        for (int m = 0; m < s2; m++)
        {
            if (LM[l][m] == -1)
            {
                count1++;
                LM[l][m]=0;  
            }
        }
    }
    }else{
    count++;
    int s2=count;
    ui count1 = 0;
    for (int l = 0; l < s2; l++)
    {   //if (updated_flag[l]==-1)
        //    cout<<"err2"<<endl;
        flag[updated_flag[l]]=-1;
        updated_flag[l]=-1;
        count1 = 0;
        for (int m = 0; m < s2; m++)
        {
            if (LM[l][m] == -1)
            {
                count1++;
                LM[l][m]=0;    
                tripletList.emplace_back(T(l, m, -1));
            }
        }
        tripletList.emplace_back(T(l, l, count1));
    }

    

    if (tripletList.size() == count * count)
    {
        // cout << "MEga prob" << endl;
        evalues.resize(k);

        for (int ss = 0; ss < k; ss++)
        {
            if (ss < count)
                // evalues(ss)=count-1;
                evalues(ss) = count;
            else if (ss == count)
                evalues(ss) = 0;
            else // evalues(ss)=-1;
                evalues(ss) = 0;
        }
    }
    else
    {
        SparseMatrix<double> M(count, count);
        M.setFromTriplets(tripletList.begin(), tripletList.end(), [](double a, double b)
                          { return b; });
        // checkM(M);

        calcEigens1(M, k, evalues, count);
    }
    tripletList.clear();
    }
  eigenVD.row(vertex) = evalues;
evalues.setZero();  
}

}