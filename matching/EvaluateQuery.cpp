
#include "EvaluateQuery.h"
#include "utility/computesetintersection.h"
#include <vector>
#include <cstring>
#include <chrono>
#include <time.h>
#include "utility/pretty_print.h"
#include <random>
#define CLOCKTOSEC(elapsed_time) ((elapsed_time) / (double)1000000)

#if ENABLE_QFLITER == 1
BSRGraph ***EvaluateQuery::qfliter_bsr_graph_;
int *EvaluateQuery::temp_bsr_base1_ = nullptr;
int *EvaluateQuery::temp_bsr_state1_ = nullptr;
int *EvaluateQuery::temp_bsr_base2_ = nullptr;
int *EvaluateQuery::temp_bsr_state2_ = nullptr;
#endif

#ifdef SPECTRUM
bool EvaluateQuery::exit_;
#endif

#ifdef DISTRIBUTION
size_t *EvaluateQuery::distribution_count_;
#endif

void calculateCellFN(const Graph *query_graph, Edges ***edge_matrix,
                   ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j,int ROQ[])
{

    ui u_nbrs_count;
    bool add = true;
    const VertexID *u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
    ui PO=ROQ[i];
    for (int k = 0; k < candidates_count[i]; k++)
    {
        add = true;
        if (k == j)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
            continue;
        }
        if (candidatesHC[i][k] != 0)
            continue;
        for (int d = 0; d < u_nbrs_count; d++)
        { 
            ui CVQ = u_nbrs[d];
            if (ROQ[CVQ]<PO)
            continue;
            ui SP = edge_matrix[i][CVQ]->offset_[j];
            ui EP = edge_matrix[i][CVQ]->offset_[j + 1];
            ui SPC = edge_matrix[i][CVQ]->offset_[k];
            ui EPC = edge_matrix[i][CVQ]->offset_[k + 1];
            int Times = EP - SP;
            if ((EP - SP) == (EPC - SPC))
            {
                for (int ee = 0; ee < Times; ee++)
                {
                    if (edge_matrix[i][CVQ]->edge_[SP + ee] != edge_matrix[i][CVQ]->edge_[SPC + ee])
                    {
                        add = false;
                        break;
                    }
                }
            }
            else
            {
                add = false;
            }
            if (add == false)
                break;
        }
        if (add == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}

int parseLine(const char *line)
{
    // This assumes the format "VmSize: <value> kB"
    const char *p = line;
    // Move p to the first digit after "VmSize:"
    while (*p < '0' || *p > '9')
        p++;
    // Convert the number part of the line to an integer
    int value = atoi(p);
    return value;
}

int getValue1()
{
    FILE *file = fopen("/proc/self/status", "r");
    if (file == nullptr)
    {
        perror("fopen");
        return -1;
    }

    int result = -1;
    char line[256]; // Increase the buffer size to handle longer lines

    while (fgets(line, sizeof(line), file) != NULL)
    {
        // Print the line for debugging purposes
        // cout << "Read line: " << line;
        // Check if the line starts with "VmSize:"
        if (strncmp(line, "VmSize:", 7) == 0)
        {
            // cout << "Found VmSize line: " << line;
            result = parseLine(line);
            break;
        }
    }

    fclose(file);
    return result;
}
void calculateCell(const Graph *query_graph, Edges ***edge_matrix,
                   ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j)
{

    ui u_nbrs_count;
    bool add = true;
    const VertexID *u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
    for (int k = 0; k < candidates_count[i]; k++)
    {
        add = true;
        if (k == j)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
            continue;
        }
        if (candidatesHC[i][k] != 0)
            continue;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[i][CVQ]->offset_[j];
            ui EP = edge_matrix[i][CVQ]->offset_[j + 1];
            ui SPC = edge_matrix[i][CVQ]->offset_[k];
            ui EPC = edge_matrix[i][CVQ]->offset_[k + 1];
            int Times = EP - SP;
            if ((EP - SP) == (EPC - SPC))
            {
                for (int ee = 0; ee < Times; ee++)
                {
                    if (edge_matrix[i][CVQ]->edge_[SP + ee] != edge_matrix[i][CVQ]->edge_[SPC + ee])
                    {
                        add = false;
                        break;
                    }
                }
            }
            else
            {
                add = false;
            }
            if (add == false)
                break;
        }
        if (add == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}

void calculateCellL(const Graph *query_graph, Edges ***edge_matrix,
                    ui *candidates_count, size_t **&candidatesHC, unordered_map<size_t, std::vector<VertexID>> *&idToValues, int count2, int i, int j, unordered_map<size_t, std::vector<VertexID>> *&idToValues1, size_t **&candidatesHC1)
{

    ui u_nbrs_count;
    bool add = true;
    vector<ui> VN;
    const VertexID *u_nbrs = query_graph->getVertexNeighbors(i, u_nbrs_count);
    VN = idToValues1[i][candidatesHC1[i][j]]; // keeps location in candidates matrix
    int siz = VN.size();
    for (int d = 0; d < siz; d++)
    {
        add = true;
        int k = VN[d];
        if (k == j)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
            continue;
        }
        if (candidatesHC[i][k] != 0)
            continue;
        for (int d = 0; d < u_nbrs_count; d++)
        {
            ui CVQ = u_nbrs[d];
            ui SP = edge_matrix[i][CVQ]->offset_[j];
            ui EP = edge_matrix[i][CVQ]->offset_[j + 1];
            ui SPC = edge_matrix[i][CVQ]->offset_[k];
            ui EPC = edge_matrix[i][CVQ]->offset_[k + 1];
            int Times = EP - SP;
            if ((EP - SP) == (EPC - SPC))
            {
                for (int ee = 0; ee < Times; ee++)
                {
                    if (edge_matrix[i][CVQ]->edge_[SP + ee] != edge_matrix[i][CVQ]->edge_[SPC + ee])
                    {
                        add = false;
                        break;
                    }
                }
            }
            else
            {
                add = false;
            }
            if (add == false)
                break;
        }
        if (add == true)
        {
            candidatesHC[i][k] = count2;
            idToValues[i][count2].push_back(k);
        }
    }
}
void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui *pivot, ui **&bn, ui *&bn_count)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr] && nbr != pivot[i])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

size_t
EvaluateQuery::exploreGraph(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                            ui *candidates_count, ui *order, ui *pivot, size_t output_limit_num, size_t &call_count)
{
    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidateIndex(data_graph, cur_depth, embedding, idx_embedding, idx_count,
                                            valid_candidate_idx,
                                            edge_matrix, visited_vertices, bn, bn_count, order, pivot, candidates);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

// Release the buffer.
EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

void EvaluateQuery::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i)
    {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num)
        {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    visited_vertices = new bool[data_vertices_num];
    temp_buffer = new ui[max_candidates_num];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        valid_candidate_idx[i] = new ui[max_candidates_num];
    }

    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}

void EvaluateQuery::generateValidCandidateIndex(const Graph *data_graph, ui depth, ui *embedding, ui *idx_embedding,
                                                ui *idx_count, ui **valid_candidate_index, Edges ***edge_matrix,
                                                bool *visited_vertices, ui **bn, ui *bn_cnt, ui *order, ui *pivot,
                                                ui **candidates)
{
    VertexID u = order[depth];
    VertexID pivot_vertex = pivot[depth];
    ui idx_id = idx_embedding[pivot_vertex];
    Edges &edge = *edge_matrix[pivot_vertex][u];
    ui count = edge.offset_[idx_id + 1] - edge.offset_[idx_id];
    ui *candidate_idx = edge.edge_ + edge.offset_[idx_id];

    ui valid_candidate_index_count = 0;

    if (bn_cnt[depth] == 0)
    {
        for (ui i = 0; i < count; ++i)
        {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v])
                valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
        }
    }
    else
    {
        for (ui i = 0; i < count; ++i)
        {
            ui temp_idx = candidate_idx[i];
            VertexID temp_v = candidates[u][temp_idx];

            if (!visited_vertices[temp_v])
            {
                bool valid = true;

                for (ui j = 0; j < bn_cnt[depth]; ++j)
                {
                    VertexID u_bn = bn[depth][j];
                    VertexID u_bn_v = embedding[u_bn];

                    if (!data_graph->checkEdgeExistence(temp_v, u_bn_v))
                    {
                        valid = false;
                        break;
                    }
                }

                if (valid)
                    valid_candidate_index[depth][valid_candidate_index_count++] = temp_idx;
            }
        }
    }

    idx_count[depth] = valid_candidate_index_count;
}

void EvaluateQuery::releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn,
                                  ui *bn_count)
{
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

enumResult
EvaluateQuery::LFTJB(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                     ui *candidates_count,
                     ui *order, size_t output_limit_num, size_t &call_count, float **&EWeight)
{

#ifdef DISTRIBUTION
    distribution_count_ = new size_t[data_graph->getVerticesCount()];
    memset(distribution_count_, 0, data_graph->getVerticesCount() * sizeof(size_t));
    size_t *begin_count = new size_t[query_graph->getVerticesCount()];
    memset(begin_count, 0, query_graph->getVerticesCount() * sizeof(size_t));
#endif

    enumResult s;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }
    sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + candidates_count[start_vertex], [&EWeight, start_vertex](const VertexID &a, const VertexID &b)
         { return EWeight[start_vertex][a] > EWeight[start_vertex][b]; });

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

#ifdef SPECTRUM
    exit_ = false;
#endif

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

#ifdef DISTRIBUTION
            begin_count[cur_depth] = embedding_cnt;
            // printf("Cur Depth: %d, v: %u, begin: %zu\n", cur_depth, v, embedding_cnt);
#endif

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef DISTRIBUTION
                distribution_count_[v] += 1;
#endif

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                            bn_count, order, temp_buffer);
                ui countCD = idx_count[cur_depth];
                u = order[cur_depth];
                // if(cur_depth<(max_depth/4) &&(query_graph->getVertexDegree(u)>1))
                if (cur_depth < (max_depth) && (query_graph->getVertexDegree(u) > 1))
                    sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + countCD, [&EWeight, u](const VertexID &a, const VertexID &b)
                         { return EWeight[u][a] > EWeight[u][b]; });
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

#ifdef SPECTRUM
        if (exit_)
        {
            goto EXIT;
        }
#endif

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;

#ifdef DISTRIBUTION
            distribution_count_[embedding[u]] += embedding_cnt - begin_count[cur_depth];
            // printf("Cur Depth: %d, v: %u, begin: %zu, end: %zu\n", cur_depth, embedding[u], begin_count[cur_depth], embedding_cnt);
#endif
        }
    }

    // Release the buffer.

#ifdef DISTRIBUTION
    if (embedding_cnt >= output_limit_num)
    {
        for (int i = 0; i < max_depth - 1; ++i)
        {
            ui v = embedding[order[i]];
            distribution_count_[v] += embedding_cnt - begin_count[i];
        }
    }
    delete[] begin_count;
#endif

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
            //            delete qfliter_bsr_graph_[i][j];
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif
    int true_cand_sum = 0;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        true_cand_sum += s.candidate_true[i].size();
    }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::LFTJA(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                     ui *candidates_count,
                     ui *order, size_t output_limit_num, size_t &call_count, float **&EWeight)
{

#ifdef DISTRIBUTION
    distribution_count_ = new size_t[data_graph->getVerticesCount()];
    memset(distribution_count_, 0, data_graph->getVerticesCount() * sizeof(size_t));
    size_t *begin_count = new size_t[query_graph->getVerticesCount()];
    memset(begin_count, 0, query_graph->getVerticesCount() * sizeof(size_t));
#endif

    enumResult s;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    // ui MR=(max_depth/2);
    ui MR = max_depth - 2;
    // MR=(max_depth/2);
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }
    sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + candidates_count[start_vertex], [&EWeight, start_vertex](const VertexID &a, const VertexID &b)
         { return EWeight[start_vertex][a] < EWeight[start_vertex][b]; });

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

#ifdef SPECTRUM
    exit_ = false;
#endif

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

#ifdef DISTRIBUTION
            begin_count[cur_depth] = embedding_cnt;
            // printf("Cur Depth: %d, v: %u, begin: %zu\n", cur_depth, v, embedding_cnt);
#endif

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef DISTRIBUTION
                distribution_count_[v] += 1;
#endif

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                            bn_count, order, temp_buffer);
                ui countCD = idx_count[cur_depth];
                u = order[cur_depth];
                // if(cur_depth<(max_depth/4) &&(query_graph->getVertexDegree(u)>1))

                //   if(cur_depth<(max_depth) &&(query_graph->getVertexDegree(u)>1) &&cur_depth<MR)
                //            sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] +  countCD, [&EWeight, u](const VertexID& a, const VertexID& b) {
                //           if (EWeight[u][b]<0.01)
                //                   return true;
                //           else if(EWeight[u][a]<0.01)
                //               return false;
                //           else return EWeight[u][a] < EWeight[u][b];
                //           });
                // else {
                //    if(cur_depth<(max_depth) &&(query_graph->getVertexDegree(u)>1) &&cur_depth>=MR)
                //         sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] +  countCD, [&EWeight, u](const VertexID& a, const VertexID& b) {
                //        if (EWeight[u][b]<0.01)
                //                return true;
                //        else if(EWeight[u][a]<0.01)
                //            return false;
                //        else return EWeight[u][a] > EWeight[u][b];
                //        });
                // }
                // if(cur_depth<(max_depth) &&(query_graph->getVertexDegree(u)>1) &&cur_depth<MR)
                if (cur_depth < (max_depth)-1)
                    sort(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + countCD, [&EWeight, u](const VertexID &a, const VertexID &b)
                         { return EWeight[u][a] < EWeight[u][b]; });

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

#ifdef SPECTRUM
        if (exit_)
        {
            goto EXIT;
        }
#endif

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;

#ifdef DISTRIBUTION
            distribution_count_[embedding[u]] += embedding_cnt - begin_count[cur_depth];
            // printf("Cur Depth: %d, v: %u, begin: %zu, end: %zu\n", cur_depth, embedding[u], begin_count[cur_depth], embedding_cnt);
#endif
        }
    }

    // Release the buffer.

#ifdef DISTRIBUTION
    if (embedding_cnt >= output_limit_num)
    {
        for (int i = 0; i < max_depth - 1; ++i)
        {
            ui v = embedding[order[i]];
            distribution_count_[v] += embedding_cnt - begin_count[i];
        }
    }
    delete[] begin_count;
#endif

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
            //            delete qfliter_bsr_graph_[i][j];
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif
    int true_cand_sum = 0;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        true_cand_sum += s.candidate_true[i].size();
    }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::LFTJR(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                     ui *candidates_count,
                     ui *order, size_t output_limit_num, size_t &call_count)
{

#ifdef DISTRIBUTION
    distribution_count_ = new size_t[data_graph->getVerticesCount()];
    memset(distribution_count_, 0, data_graph->getVerticesCount() * sizeof(size_t));
    size_t *begin_count = new size_t[query_graph->getVerticesCount()];
    memset(begin_count, 0, query_graph->getVerticesCount() * sizeof(size_t));
#endif

    enumResult s;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    default_random_engine rng(std::random_device{}());
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }
    shuffle(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + idx_count[cur_depth] - 1, rng);

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

#ifdef SPECTRUM
    exit_ = false;
#endif

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

#ifdef DISTRIBUTION
            begin_count[cur_depth] = embedding_cnt;
            // printf("Cur Depth: %d, v: %u, begin: %zu\n", cur_depth, v, embedding_cnt);
#endif

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;

#ifdef DISTRIBUTION
                distribution_count_[v] += 1;
#endif

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                //                             bn_count, order, temp_buffer);
                if (idx_count[cur_depth] > 1)
                    shuffle(valid_candidate_idx[cur_depth], valid_candidate_idx[cur_depth] + idx_count[cur_depth] - 1, rng);
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

#ifdef SPECTRUM
        if (exit_)
        {
            goto EXIT;
        }
#endif

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;

#ifdef DISTRIBUTION
            distribution_count_[embedding[u]] += embedding_cnt - begin_count[cur_depth];
            // printf("Cur Depth: %d, v: %u, begin: %zu, end: %zu\n", cur_depth, embedding[u], begin_count[cur_depth], embedding_cnt);
#endif
        }
    }

    // Release the buffer.

#ifdef DISTRIBUTION
    if (embedding_cnt >= output_limit_num)
    {
        for (int i = 0; i < max_depth - 1; ++i)
        {
            ui v = embedding[order[i]];
            distribution_count_[v] += embedding_cnt - begin_count[i];
        }
    }
    delete[] begin_count;
#endif

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
            //            delete qfliter_bsr_graph_[i][j];
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif
    int true_cand_sum = 0;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        true_cand_sum += s.candidate_true[i].size();
    }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count,
                    ui *order, size_t output_limit_num, size_t &call_count)
{    int MemSize = 0;


    enumResult s;
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }
    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {


            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth];         // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
 
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1; // found a match
                // cout<<"PL"<<embedding_cnt<<endl;
                visited_vertices[v] = false;

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                                    // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                    //                            bn_count, order, temp_buffer);
                                    // get the candidates for next depth
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                bn_count, order, temp_buffer, candidates, candidates_count, query_graph);

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.

EXIT:
 //if (getValue1() > MemSize)
  //      MemSize = getValue1();
//cout << "MemSize" << MemSize / 1000 << endl;
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

#if ENABLE_QFLITER == 1
    delete[] temp_bsr_base1_;
    delete[] temp_bsr_base2_;
    delete[] temp_bsr_state1_;
    delete[] temp_bsr_state2_;

    for (ui i = 0; i < max_depth; ++i)
    {
        for (ui j = 0; j < max_depth; ++j)
        {
            //            delete qfliter_bsr_graph_[i][j];
        }
        delete[] qfliter_bsr_graph_[i];
    }
    delete[] qfliter_bsr_graph_;
#endif
    int true_cand_sum = 0;

    // for (int i=0; i<query_graph->getVerticesCount();i++){
    //     true_cand_sum += s.candidate_true[i].size();
    // }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::LFTJSS(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                      ui *candidates_count,
                      ui *order, size_t output_limit_num, size_t &call_count)
{

    enumResult s;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth];         // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
            // cout<<"u,v :"<<u<<","<<v<<endl;
            // cout<<"or id u,v :"<<u<<","<<valid_idx<<endl;
            if (visited_vertices[v])
            {
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1; // found a match

                cout << "PL" << embedding_cnt << endl;
                visited_vertices[v] = false;

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                                    //  generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                    //                              bn_count, order, temp_buffer);
                                    // get the candidates for next depth

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    int true_cand_sum = 0;

    // for (int i=0; i<query_graph->getVerticesCount();i++){
    //     true_cand_sum += s.candidate_true[i].size();
    // }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

enumResult
EvaluateQuery::LFTJVEQCN(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t **candidatesHC, unordered_map<size_t, std::vector<VertexID>> *idToValues2)
{
//auto start = std::chrono::high_resolution_clock::now();
//double ens=10000;
    int MemSize = 0;
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VP[query_graph->getVerticesCount()];
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;



    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates
    int count2 = 1;
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, start_vertex, 0);
    count2++;
    VP[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];
    VN[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    int k;

    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)
            {
                idx[cur_depth] += 1;
                continue;
            }
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth]; // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, u, valid_idx);
                count2++;
            }
            candidate = candidatesHC[u][valid_idx];
           //         if (getValue1() > MemSize)
           // MemSize = getValue1();
            VP[u].clear();
            VN[u].clear();
            // if (cur_depth<(max_depth-2)){
            if (cur_depth != max_depth - 1)
            {
                VP[u] = idToValues2[u][candidatesHC[u][valid_idx]]; // keeps location in candidates matrix
                VN[u] = idToValues2[u][candidatesHC[u][valid_idx]];
            }

            Match_BA[u] = false; // set embedings found to empty for start

            if (visited_vertices[v])
            {   
                utemp=reverse_embedding[v];
                
                if (cur_depth != max_depth - 1&&Match_BA[utemp] == false)
                {
                        //ui valid_idxN = valid_candidate_idx[k][idx[k] - 1];
                       // vtemp = candidates[utemp][valid_idxN];
                        VNTemp.clear();
                            if (VN[u].size() > 1 && VN[utemp].size() > 1)
                            { 
                                int ia = 0;
                                int io = 0;
                                while (ia < VN[u].size() && io < VN[utemp].size())
                                { 
                                    if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                                    {
                                        VNTemp.push_back(VN[utemp][io]);
                                        ia++;
                                        io++;
                                    }
                                    else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                        ia++;
                                    else
                                        io++;
                                }
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();                                                               
                }
                else if(cur_depth == max_depth - 1&&Match_BA[utemp] == false){
                int ia=0;
                int io=0;
                VNTemp.clear();
                while (ia<idx_count[cur_depth]&&io<VN[utemp].size()){
                    ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                    if (valid_idxN == 10000000)
                    {
                        ia += 1;
                        continue;
                    }
                    if (VN[utemp][io]==valid_idxN){
                        //VNTemp.push_back(VN[u][io]);
                        VNTemp.push_back(VN[utemp][io]);
                        ia++;
                        io++;
                    }
                    
                    else if(VN[utemp][io]>valid_idxN)
                    ia++;
                    else
                    io++;
                  // ia++;
                }
                VN[utemp]=VNTemp; 
                VNTemp.clear();
                }

        if (cur_depth != max_depth - 1 && !Match_BA[utemp]&&false)
        {
    if (VN[u].size() > 1 && VN[utemp].size() > 1)
    {
        size_t ia = 0;                // iterator in VN[u]
        size_t io = 0;                // iterator in VN[utemp] (read)
        size_t write = 0;             // iterator in VN[utemp] (write)

        while (ia < VN[u].size() && io < VN[utemp].size())
        {
            auto cu = candidates[u]     [ VN[u]     [ia] ];
            auto cv = candidates[utemp][ VN[utemp][io] ];

            if (cu == cv)                               // keep
            {
                VN[utemp][write++] = VN[utemp][io];
                ++ia;
                ++io;
            }
            else if (cu < cv) ++ia;                     // advance ia
            else              ++io;                     // advance io
        }
        VN[utemp].resize(write);                        // discard tail
    }
}

/* ------------------------------------------------------------------ */
/*  Case 2:  cur_depth == max_depth‑1  (intersect VN[utemp] with       */
/*           valid_candidate_idx slice [0 .. idx_count‑1) )            */
else if (cur_depth == max_depth - 1 && !Match_BA[utemp]&&false)
{
    size_t ia = 0;                       // iterator in valid_candidate_idx
    size_t io = 0;                       // iterator in VN[utemp] (read)
    size_t write = 0;                    // iterator in VN[utemp] (write)

    while (ia < idx_count[cur_depth] && io < VN[utemp].size())
    {
        ui valid_idxN = valid_candidate_idx[cur_depth][ia];

        if (valid_idxN == 10000000) { ++ia; continue; }  // skip holes

        if (VN[utemp][io] == valid_idxN)                // keep
        {
            VN[utemp][write++] = VN[utemp][io];
            ++ia;
            ++io;
        }
        else if (VN[utemp][io] > valid_idxN) ++ia;      // advance ia
        else                                          ++io;  // advance io
    }
    VN[utemp].resize(write);                            // discard tail
}
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                NotifyALL = true;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                NotifyALL = false;
                int ao = cur_depth;//-1;
                // while(ao>=0 &&Match_BA[ao]==false){
                while (ao >= 0)
                {
                    // VN[order[cur_depth]].clear();
                    Match_BA[order[ao]] = true;
                    ao--;
                }

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {   
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer, candidates, candidates_count, query_graph);

                // get the candidates for next depth
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        int ao = cur_depth;

        if (NotifyALL == true)
        { // notify all Match_ that there was a match.
            // after visiting all the nodes of last level so it wont be repeated
            NotifyALL = false;
            while (ao >= 0)
            {
                Match_BA[ao] = true;
                ao--;
            }
        } // no embeddings case
        if (cur_depth == max_depth - 1)
            continue;
        ui vq = order[cur_depth];
        if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
        {
            if (VN[vq].size() > 1)
            { // something to remove
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (valid_idxn == VN[vq][dd])
                        {
                            // if(candidates[vq][valid_idxn]==candidates[vq][VN[vq][dd]]){
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
        }
        else if (false)
        {
            if (VP[order[cur_depth]].size() > 1)
            {
                int sp = 0;
                ui sLo = 0;
                for (int dd = 0; dd < VN[order[cur_depth]].size(); dd++)
                {
                    for (int kk = idx[cur_depth]; kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (candidates[order[cur_depth]][valid_idxn] == candidates[order[cur_depth]][VN[order[cur_depth]][dd]])
                        {
                            sLo = valid_candidate_idx[cur_depth][idx[cur_depth] + sp];
                            valid_candidate_idx[cur_depth][idx[cur_depth] + sp] = valid_candidate_idx[cur_depth][kk];
                            valid_candidate_idx[cur_depth][kk] = sLo;
                            sp++;
                            break;
                        } // ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];//valid id
                    }
                }
            }
        }

        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

EXIT:

    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    int true_cand_sum = 0;
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}



enumResult
EvaluateQuery::LFTJVEQ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t **candidatesHC, unordered_map<size_t, std::vector<VertexID>> *idToValues2)
{

    //int MemSize = 0;
    //if (getValue1() > MemSize)
    //    MemSize = getValue1();
    enumResult s;

    // for (int i=0; i<query_graph->getVerticesCount();i++){
    //     s.candidate_true.push_back(set<ui>());
    // }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VP[query_graph->getVerticesCount()];
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;

    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates
    int count2 = 1;
    //if (start_vertex >= max_depth)
    //    cout << "gia afto" << endl;
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, start_vertex, 0);
    count2++;
    VP[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];
    VN[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    int k;

    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)
            {
                idx[cur_depth] += 1;
                continue;
            }
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth]; // u is current depth qid
            if (u >= max_depth)
                cout << "gia afto" << endl;
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0)
            {
                calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, u, valid_idx);
                count2++;
            }
            candidate = candidatesHC[u][valid_idx];
            VP[u].clear();
            VN[u].clear();
            // if (cur_depth<(max_depth-2)){
            if (cur_depth != max_depth - 1)
            {
                VP[u] = idToValues2[u][candidatesHC[u][valid_idx]]; // keeps location in candidates matrix
                VN[u] = idToValues2[u][candidatesHC[u][valid_idx]];
            }
            // if (VN[u].size()==0)
            // cout<<"error";
            //}

            Match_BA[u] = false; // set embedings found to empty for start

            if (visited_vertices[v])
            {
                if (cur_depth != max_depth - 1)
                {
                    for (k = 0; k < cur_depth; k++)
                    { // cout<<"Phase B"<<endl;
                        if (Match_BA[order[k]] == true)
                            continue;
                        utemp = order[k];
                        // ui valid_idxN=valid_candidate_idx[k][idx[k]-1];
                        ui valid_idxN = valid_candidate_idx[k][idx[k] - 1];

                        vtemp = candidates[utemp][valid_idxN];

                        VNTemp.clear();
                        if (vtemp == v)
                        {
                            // if (utemp!=reverse_embedding[v]) //alternative
                            //    cout<<"problima";
                            // VNTemp.clear();
                            if (VN[u].size() > 1 && VN[utemp].size() > 1)
                            { // VP Negative-Positive is larger
                                // if (VN[u].size()>1&&VN[utemp].size()>1){
                                // if (true){
                                int ia = 0;
                                int io = 0;
                                // cout<<VN[utemp].size()<<endl;
                                while (ia < VN[u].size() && io < VN[utemp].size())
                                { // sorted?

                                    if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                                    {
                                        // cout<<VN[u][ia]<<","<<VN[utemp][io]<<endl;
                                        // VN[u].push_back(VN[u][ia]) ;
                                        // VN[utemp] = VNTemp;
                                        // cout<<candidates[u][VN[u][ia]]<<candidates[utemp][VN[utemp][io]]<<endl;
                                        // VNTemp.push_back(VN[u][ia]);
                                        VNTemp.push_back(VN[utemp][io]);

                                        ia++;
                                        io++;
                                    }
                                    else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                        ia++;
                                    else
                                        io++;
                                }
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                            break;
                        }
                    }
                }
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                NotifyALL = true;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;

                NotifyALL = false;
                int ao = cur_depth;//-1;
                // while(ao>=0 &&Match_BA[ao]==false){
                while (ao >= 0)
                {
                    // VN[order[cur_depth]].clear();
                    Match_BA[order[ao]] = true;
                    ao--;
                }

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {   
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                // get the candidates for next depth
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        int ao = cur_depth;

        if (NotifyALL == true)
        { // notify all Match_ that there was a match.
            // after visiting all the nodes of last level so it wont be repeated
            NotifyALL = false;
            while (ao >= 0)
            {
                Match_BA[ao] = true;
                ao--;
            }
        } // no embeddings case
        if (cur_depth == max_depth - 1)
            continue;
        ui vq = order[cur_depth];
        if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
        {
            if (VN[vq].size() > 1)
            { // something to remove
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (valid_idxn == VN[vq][dd])
                        {
                            // if(candidates[vq][valid_idxn]==candidates[vq][VN[vq][dd]]){
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
        }
        else if (false)
        {
            if (VP[order[cur_depth]].size() > 1)
            {
                int sp = 0;
                ui sLo = 0;
                for (int dd = 0; dd < VN[order[cur_depth]].size(); dd++)
                {
                    for (int kk = idx[cur_depth]; kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (candidates[order[cur_depth]][valid_idxn] == candidates[order[cur_depth]][VN[order[cur_depth]][dd]])
                        {
                            sLo = valid_candidate_idx[cur_depth][idx[cur_depth] + sp];
                            valid_candidate_idx[cur_depth][idx[cur_depth] + sp] = valid_candidate_idx[cur_depth][kk];
                            valid_candidate_idx[cur_depth][kk] = sLo;
                            sp++;
                            break;
                        } // ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];//valid id
                    }
                }
            }
        }

        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

EXIT:
    //if (getValue1() > MemSize)
    //    MemSize = getValue1();
   // cout << "MemSize" << MemSize / 1000 << endl;
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    int true_cand_sum = 0;
    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}
//FOrward Neighors
enumResult
EvaluateQuery::LFTJVEQFN(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                       ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t **candidatesHC, unordered_map<size_t, std::vector<VertexID>> *idToValues2)
{
  
    enumResult s;
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);
    vector<ui> VP[query_graph->getVerticesCount()];
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;
    int RQ[query_graph->getVerticesCount()];
    for (int i=0;i<query_graph->getVerticesCount();i++){
        RQ[order[i]]=i;
    }
    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates
    int count2 = 1;
    calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, start_vertex, 0);
    count2++;
    VP[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];
    VN[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    int k;

    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;
    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)
            {
                idx[cur_depth] += 1;
                continue;
            }
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth]; // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0&&cur_depth != max_depth - 1)
            {
                calculateCellFN(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, u, valid_idx,RQ);
                //calculateCell(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, u, valid_idx);
                count2++;
            }
            candidate = candidatesHC[u][valid_idx];
            VP[u].clear();
            VN[u].clear();

            if (cur_depth != max_depth - 1)//change here
            {
                 // keeps location in candidates matrix
                VN[u] = idToValues2[u][candidatesHC[u][valid_idx]];
            if (cur_depth != 0) {
                size_t write = 0;
                //size_t ia = idx[cur_depth];
                size_t ia =0;
                for (size_t read = 0; read < VN[u].size() && ia < idx_count[cur_depth]; ) {
                    ui valid = valid_candidate_idx[cur_depth][ia];
                    if (valid == 10000000) { ++ia; continue; }
                    if (VN[u][read] == valid) VN[u][write++] = VN[u][read++], ++ia;
                    else if (VN[u][read] > valid) ++ia;
                    else ++read;
                }
                VN[u].resize(write);
                
            }
            VP[u] = VN[u];
        }
            
            if (cur_depth != max_depth - 1&&false)//change here
            {
                VP[u] = idToValues2[u][candidatesHC[u][valid_idx]]; // keeps location in candidates matrix
                VN[u] = idToValues2[u][candidatesHC[u][valid_idx]];
                if (cur_depth!=0){
                int ia=idx[cur_depth];
                int io=0;
                VNTemp.clear();
                while (ia<idx_count[cur_depth]&&io<VN[u].size()){
                    ui valid_idx = valid_candidate_idx[cur_depth][ia];
                    if (valid_idx == 10000000)
                    {
                        ia += 1;
                        continue;
                    }
                    if (VN[u][io]==valid_idx){
                        VNTemp.push_back(VN[u][io]);
                        ia++;
                        io++;
                    }
                    
                    else if(VN[u][io]>valid_idx)
                    ia++;
                    else
                    io++;
                }
                VN[u]=VNTemp;
            }
            }


            Match_BA[u] = false; // set embedings found to empty for start
            
            if (visited_vertices[v]&false)
            {
                if (cur_depth != max_depth - 1)
                {
                    for (k = 0; k < cur_depth; k++)
                    { // cout<<"Phase B"<<endl;
                        if (Match_BA[order[k]] == true)
                            continue;
                        utemp = order[k];
                        // ui valid_idxN=valid_candidate_idx[k][idx[k]-1];
                        ui valid_idxN = valid_candidate_idx[k][idx[k] - 1];

                        vtemp = candidates[utemp][valid_idxN];

                        VNTemp.clear();
                        if (vtemp == v)
                        {

                            if (VN[u].size() > 1 && VN[utemp].size() > 1)
                            {
                                int ia = 0;
                                int io = 0;
                                // cout<<VN[utemp].size()<<endl;
                                while (ia < VN[u].size() && io < VN[utemp].size())
                                { // sorted?

                                    if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                                    {
                                        VNTemp.push_back(VN[utemp][io]);
                                        ia++;
                                        io++;
                                    }
                                    else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                        ia++;
                                    else
                                        io++;
                                }
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                            break;
                        }
                    }
                }
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            
            
           if (visited_vertices[v])
            {   
                utemp=reverse_embedding[v];
                
                if (cur_depth != max_depth - 1&&Match_BA[utemp] == false&&false)
                {
                        //ui valid_idxN = valid_candidate_idx[k][idx[k] - 1];
                       // vtemp = candidates[utemp][valid_idxN];
                        VNTemp.clear();
                            if (VN[u].size() > 1 && VN[utemp].size() > 1)
                            { 
                                int ia = 0;
                                int io = 0;
                                while (ia < VN[u].size() && io < VN[utemp].size())
                                { 
                                    if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                                    {
                                        VNTemp.push_back(VN[utemp][io]);
                                        ia++;
                                        io++;
                                    }
                                    else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                        ia++;
                                    else
                                        io++;
                                }
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();                                                               
                }
                else if(cur_depth == max_depth - 1&&Match_BA[utemp] == false&&false){
                int ia=0;
                int io=0;
                VNTemp.clear();
                while (ia<idx_count[cur_depth]&&io<VN[utemp].size()){
                    ui valid_idxN = valid_candidate_idx[cur_depth][ia];
                    if (valid_idxN == 10000000)
                    {
                        ia += 1;
                        continue;
                    }
                    if (VN[utemp][io]==valid_idxN){
                        //VNTemp.push_back(VN[u][io]);
                        VNTemp.push_back(VN[utemp][io]);
                        ia++;
                        io++;
                    }
                    
                    else if(VN[utemp][io]>valid_idxN)
                    ia++;
                    else
                    io++;
                  // ia++;
                }
                VN[utemp]=VNTemp; 
                VNTemp.clear();
                }
                

                if (cur_depth != max_depth - 1 && !Match_BA[utemp])
                //if (cur_depth != max_depth - 1&&false)
                {size_t write = 0; 
            if (VN[u].size() > 0 && VN[utemp].size() > 0)
            {
                size_t ia = 0;                // iterator in VN[u]
                size_t io = 0;                // iterator in VN[utemp] (read)
                            // iterator in VN[utemp] (write)
        
                while (ia < VN[u].size() && io < VN[utemp].size())
                {
                    auto cu = candidates[u]     [ VN[u]     [ia] ];
                    auto cv = candidates[utemp][ VN[utemp][io] ];
        
                    if (cu == cv)                               // keep
                    {
                        VN[utemp][write++] = VN[utemp][io];
                        ++ia;
                        ++io;
                    }
                    else if (cu < cv) ++ia;                     // advance ia
                    else              ++io;                     // advance io
                }
                VN[utemp].resize(write);                          // discard tail
            }
        }
        
        /* ------------------------------------------------------------------ */
        /*  Case 2:  cur_depth == max_depth‑1  (intersect VN[utemp] with       */
        /*           valid_candidate_idx slice [0 .. idx_count‑1) )            */
        else if (cur_depth == max_depth - 1 && !Match_BA[utemp])
        //else if (cur_depth == max_depth - 1 &&false)
        {
            size_t ia = 0;                       // iterator in valid_candidate_idx
            size_t io = 0;                       // iterator in VN[utemp] (read)
            size_t write = 0;                    // iterator in VN[utemp] (write)
        
            while (ia < idx_count[cur_depth] && io < VN[utemp].size())
            {
                ui valid_idxN = valid_candidate_idx[cur_depth][ia];
        
                if (valid_idxN == 10000000) { ++ia; continue; }  // skip holes
        
                if (VN[utemp][io] == valid_idxN)                // keep
                {
                    VN[utemp][write++] = VN[utemp][io];
                    ++ia;
                    ++io;
                }
                else if (VN[utemp][io] > valid_idxN) ++ia;      // advance ia
                else                                          ++io;  // advance io
            }
            VN[utemp].resize(write);                            // discard tail
        }
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                NotifyALL = true;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                NotifyALL = false;
                int ao = cur_depth;//-1;
                // while(ao>=0 &&Match_BA[ao]==false){
                while (ao >= 0)
                {
                    // VN[order[cur_depth]].clear();
                    Match_BA[order[ao]] = true;
                    ao--;
                }

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {   
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer, candidates, candidates_count, query_graph);
                // get the candidates for next depth
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        int ao = cur_depth;

        if (NotifyALL == true)
        { // notify all Match_ that there was a match.
            // after visiting all the nodes of last level so it wont be repeated
            NotifyALL = false;
            while (ao >= 0)
            {
                Match_BA[ao] = true;
                ao--;
            }
        } // no embeddings case
        if (cur_depth == max_depth - 1)
            continue;
        ui vq = order[cur_depth];
        if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
        {
            if (VN[vq].size() > 0)
            { // something to remove
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (valid_idxn == VN[vq][dd])
                        {
                            // if(candidates[vq][valid_idxn]==candidates[vq][VN[vq][dd]]){
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
        }
        else if (false)
        {
            if (VP[order[cur_depth]].size() > 1)
            {
                int sp = 0;
                ui sLo = 0;
                for (int dd = 0; dd < VN[order[cur_depth]].size(); dd++)
                {
                    for (int kk = idx[cur_depth]; kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (candidates[order[cur_depth]][valid_idxn] == candidates[order[cur_depth]][VN[order[cur_depth]][dd]])
                        {
                            sLo = valid_candidate_idx[cur_depth][idx[cur_depth] + sp];
                            valid_candidate_idx[cur_depth][idx[cur_depth] + sp] = valid_candidate_idx[cur_depth][kk];
                            valid_candidate_idx[cur_depth][kk] = sLo;
                            sp++;
                            break;
                        } // ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];//valid id
                    }
                }
            }
        }

        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

EXIT:
    //if (getValue1() > MemSize)
    //    MemSize = getValue1();
   // cout << "MemSize" << MemSize / 1000 << endl;
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    int true_cand_sum = 0;

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}





enumResult
EvaluateQuery::LFTJVEQL(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                        ui *candidates_count, ui *order, size_t output_limit_num, size_t &call_count, size_t **candidatesHC, unordered_map<size_t, std::vector<VertexID>> *idToValues2, unordered_map<size_t, std::vector<VertexID>> *idToValues, size_t **candidatesHC1)
{

    enumResult s;

    // for (int i=0; i<query_graph->getVerticesCount();i++){
    //     s.candidate_true.push_back(set<ui>());
    // }

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);
    bool NotifyALL = true;
    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    vector<ui> VP[query_graph->getVerticesCount()];
    vector<ui> VN[query_graph->getVerticesCount()];
    vector<ui> VNTemp;

    bool *embd_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(embd_vertices, embd_vertices + data_graph->getVerticesCount(), false);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];
    bool Match_BA[max_depth] = {false};
    idx[cur_depth] = 0;                                    // set the starting position to 0 for the current depth
    idx_count[cur_depth] = candidates_count[start_vertex]; // set the size for current depth
    // all the nodes are possible candidates
    int count2 = 1;
    calculateCellL(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, start_vertex, 0, idToValues, candidatesHC1);
    count2++;
    VP[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];
    VN[start_vertex] = idToValues2[start_vertex][candidatesHC[start_vertex][0]];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidate_idx[cur_depth][i] = i;
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, bn, bn_count, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif
    int k;

    VertexID utemp; // u is current depth qid
    VertexID vtemp;
    int counter = 0;

    while (true)
    {
        // current position smaller than size
        while (idx[cur_depth] < idx_count[cur_depth])
        {

            // from all the valid IDS
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]]; // valid id
            if (valid_idx == 10000000)
            {
                idx[cur_depth] += 1;
                continue;
            }
            // the position is the first element from the idx matrix
            VertexID u = order[cur_depth];         // u is current depth qid
            VertexID v = candidates[u][valid_idx]; // v is FCS[qid][vid] valid_idx-> position in CS
            size_t candidate = candidatesHC[u][valid_idx];
            if (candidate == 0)
            {
                calculateCellL(query_graph, edge_matrix, candidates_count, candidatesHC, idToValues2, count2, u, valid_idx, idToValues, candidatesHC1);
                count2++;
            }
            candidate = candidatesHC[u][valid_idx];
            VP[u].clear();
            VN[u].clear();
            // if (cur_depth<(max_depth-2)){
            if (cur_depth != max_depth - 1)
            {
                VP[u] = idToValues2[u][candidatesHC[u][valid_idx]]; // keeps location in candidates matrix
                VN[u] = idToValues2[u][candidatesHC[u][valid_idx]];
            }
            // if (VN[u].size()==0)
            // cout<<"error";
            //}

            Match_BA[u] = false; // set embedings found to empty for start

            if (visited_vertices[v])
            {
                if (cur_depth != max_depth - 1)
                {
                    for (k = 0; k < cur_depth; k++)
                    { // cout<<"Phase B"<<endl;
                        if (Match_BA[order[k]] == true)
                            continue;
                        utemp = order[k];
                        // ui valid_idxN=valid_candidate_idx[k][idx[k]-1];
                        ui valid_idxN = valid_candidate_idx[k][idx[k] - 1];

                        vtemp = candidates[utemp][valid_idxN];

                        VNTemp.clear();
                        if (vtemp == v)
                        {
                            // if (utemp!=reverse_embedding[v]) //alternative
                            //    cout<<"problima";
                            // VNTemp.clear();
                            if (VN[u].size() > 1 && VN[utemp].size() > 1)
                            { // VP Negative-Positive is larger
                                // if (VN[u].size()>1&&VN[utemp].size()>1){
                                // if (true){
                                int ia = 0;
                                int io = 0;
                                // cout<<VN[utemp].size()<<endl;
                                while (ia < VN[u].size() && io < VN[utemp].size())
                                { // sorted?

                                    if (candidates[u][VN[u][ia]] == candidates[utemp][VN[utemp][io]])
                                    {
                                        // cout<<VN[u][ia]<<","<<VN[utemp][io]<<endl;
                                        // VN[u].push_back(VN[u][ia]) ;
                                        // VN[utemp] = VNTemp;
                                        // cout<<candidates[u][VN[u][ia]]<<candidates[utemp][VN[utemp][io]]<<endl;
                                        // VNTemp.push_back(VN[u][ia]);
                                        VNTemp.push_back(VN[utemp][io]);

                                        ia++;
                                        io++;
                                    }
                                    else if (candidates[u][VN[u][ia]] < candidates[utemp][VN[utemp][io]])
                                        ia++;
                                    else
                                        io++;
                                }
                            }
                            VN[utemp] = VNTemp;
                            VNTemp.clear();
                            break;
                        }
                    }
                }
                idx[cur_depth] += 1;
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }
            embd_vertices[v] = true;
            counter++;
            embedding[u] = v;             // set the embeding of depth qid to v
            idx_embedding[u] = valid_idx; // reverse the id of the embdeing to the position
            visited_vertices[v] = true;
            idx[cur_depth] += 1; // next element

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                NotifyALL = true;
                visited_vertices[v] = false;
                embd_vertices[v] = false;
                counter--;
                NotifyALL = false;
                int ao = cur_depth;
                // while(ao>=0 &&Match_BA[ao]==false){
                while (ao >= 0)
                {
                    // VN[order[cur_depth]].clear();
                    Match_BA[order[ao]] = true;
                    ao--;
                }

#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;     // next depth
                idx[cur_depth] = 0; // set the element to 0
                                    // generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                    //                            bn_count, order, temp_buffer);
                                    // get the candidates for next depth

#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        int ao = cur_depth;

        if (NotifyALL == true)
        { // notify all Match_ that there was a match.
            // after visiting all the nodes of last level so it wont be repeated
            NotifyALL = false;
            while (ao >= 0)
            {
                Match_BA[ao] = true;
                ao--;
            }
        } // no embeddings case
        if (cur_depth == max_depth - 1)
            continue;
        ui vq = order[cur_depth];
        if (Match_BA[vq] == false && idx[cur_depth] < idx_count[cur_depth])
        {
            if (VN[vq].size() > 1)
            { // something to remove
                for (int dd = 0; dd < VN[vq].size(); dd++)
                {
                    for (int kk = (idx[cur_depth]); kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (valid_idxn == VN[vq][dd])
                        {
                            // if(candidates[vq][valid_idxn]==candidates[vq][VN[vq][dd]]){
                            valid_candidate_idx[cur_depth][kk] = 10000000;
                            break;
                        }
                    }
                }
            }
        }
        else if (true)
        {
            if (VP[order[cur_depth]].size() > 1)
            {
                int sp = 0;
                ui sLo = 0;
                for (int dd = 0; dd < VN[order[cur_depth]].size(); dd++)
                {
                    for (int kk = idx[cur_depth]; kk < idx_count[cur_depth]; kk++)
                    {
                        ui valid_idxn = valid_candidate_idx[cur_depth][kk];
                        if (valid_idxn == 10000000)
                            continue;
                        // VertexID vn = candidates[order[cur_depth]][valid_idxn];
                        if (candidates[order[cur_depth]][valid_idxn] == candidates[order[cur_depth]][VN[order[cur_depth]][dd]])
                        {
                            sLo = valid_candidate_idx[cur_depth][idx[cur_depth] + sp];
                            valid_candidate_idx[cur_depth][idx[cur_depth] + sp] = valid_candidate_idx[cur_depth][kk];
                            valid_candidate_idx[cur_depth][kk] = sLo;
                            sp++;
                            break;
                        } // ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];//valid id
                    }
                }
            }
        }

        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
            embd_vertices[embedding[u]] = false;
            counter--;
        }
    }

EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    int true_cand_sum = 0;

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;
    return s;
}

void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer, ui **candidates, ui *candidates_count, const Graph *query_graph)
{
    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

#if ENABLE_QFLITER == 1
    BSRGraph &bsr_graph = *qfliter_bsr_graph_[previous_bn][u];
    BSRSet &bsr_set = bsr_graph.bsrs[previous_index_id];

    if (bsr_set.size_ != 0)
    {
        offline_bsr_trans_uint(bsr_set.base_, bsr_set.states_, bsr_set.size_,
                               (int *)valid_candidate_index[depth]);
        // represent bsr size
        valid_candidates_count = bsr_set.size_;
    }

    if (bn_cnt[depth] > 0)
    {
        if (temp_bsr_base1_ == nullptr)
        {
            temp_bsr_base1_ = new int[1024 * 1024];
        }
        if (temp_bsr_state1_ == nullptr)
        {
            temp_bsr_state1_ = new int[1024 * 1024];
        }
        if (temp_bsr_base2_ == nullptr)
        {
            temp_bsr_base2_ = new int[1024 * 1024];
        }
        if (temp_bsr_state2_ == nullptr)
        {
            temp_bsr_state2_ = new int[1024 * 1024];
        }
        int *res_base_ = temp_bsr_base1_;
        int *res_state_ = temp_bsr_state1_;
        int *input_base_ = temp_bsr_base2_;
        int *input_state_ = temp_bsr_state2_;

        memcpy(input_base_, bsr_set.base_, sizeof(int) * bsr_set.size_);
        memcpy(input_state_, bsr_set.states_, sizeof(int) * bsr_set.size_);

        for (ui i = 1; i < bn_cnt[depth]; ++i)
        {
            VertexID current_bn = bn[depth][i];
            ui current_index_id = idx_embedding[current_bn];
            BSRGraph &cur_bsr_graph = *qfliter_bsr_graph_[current_bn][u];
            BSRSet &cur_bsr_set = cur_bsr_graph.bsrs[current_index_id];

            if (valid_candidates_count == 0 || cur_bsr_set.size_ == 0)
            {
                valid_candidates_count = 0;
                break;
            }

            valid_candidates_count = intersect_qfilter_bsr_b4_v2(cur_bsr_set.base_, cur_bsr_set.states_,
                                                                 cur_bsr_set.size_,
                                                                 input_base_, input_state_, valid_candidates_count,
                                                                 res_base_, res_state_);

            swap(res_base_, input_base_);
            swap(res_state_, input_state_);
        }

        if (valid_candidates_count != 0)
        {
            valid_candidates_count = offline_bsr_trans_uint(input_base_, input_state_, valid_candidates_count,
                                                            (int *)valid_candidate_index[depth]);
        }
    }

    idx_count[depth] = valid_candidates_count;

    // Debugging.
#ifdef YCHE_DEBUG
    Edges &previous_edge = *edge_matrix[previous_bn][u];

    auto gt_valid_candidates_count =
        previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];
    ui *gt_valid_candidate_index = new ui[1024 * 1024];
    memcpy(gt_valid_candidate_index, previous_candidates, gt_valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i)
    {
        VertexID current_bn = bn[depth][i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
            current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  gt_valid_candidate_index, gt_valid_candidates_count, temp_buffer,
                                                  temp_count);

        std::swap(temp_buffer, gt_valid_candidate_index);
        gt_valid_candidates_count = temp_count;
    }
    assert(valid_candidates_count == gt_valid_candidates_count);

    cout << "Ret, Level:" << bn_cnt[depth] << ", BSR:"
         << pretty_print_array(valid_candidate_index[depth], valid_candidates_count)
         << "; GT: " << pretty_print_array(gt_valid_candidate_index, gt_valid_candidates_count) << "\n";

    for (auto i = 0; i < valid_candidates_count; i++)
    {
        assert(gt_valid_candidate_index[i] == valid_candidate_index[depth][i]);
    }
    delete[] gt_valid_candidate_index;
#endif
#else
    if (previous_bn > query_graph->getVerticesCount())
    {
        cout << "sexy" << endl;
        valid_candidate_index[depth] = new ui[candidates_count[u]];
        idx_count[depth] = candidates_count[u];
        for (ui i = 0; i < idx_count[depth]; ++i)
        {
            valid_candidate_index[depth][i] = i;
        }
    }
    else if (query_graph->checkEdgeExistence(previous_bn, u) == false)
    {
        cout << "sexy1" << endl;
        valid_candidate_index[depth] = new ui[candidates_count[u]];
        idx_count[depth] = candidates_count[u];
        for (ui i = 0; i < idx_count[depth]; ++i)
        {
            valid_candidate_index[depth][i] = i;
        }
    }
    else
    {
        Edges &previous_edge = *edge_matrix[previous_bn][u];
        valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
        ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

        memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

        ui temp_count;
        for (ui i = 1; i < bn_cnt[depth]; ++i)
        {
            VertexID current_bn = bn[depth][i];
            Edges &current_edge = *edge_matrix[current_bn][u];
            ui current_index_id = idx_embedding[current_bn];

            ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
            ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

            ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                                                      temp_buffer, temp_count);
            /*cout<<"current Candidates matrix"<<endl;
            for (int aa=0;aa<current_candidates_count;aa++){
                cout<<current_candidates[aa]<<",";
            }cout<<"valid Candidates matrix"<<endl;
            for (int aa=0;aa<valid_candidates_count;aa++){
                cout<<valid_candidate_index[depth][aa]<<",";
            }cout<<"temps"<<endl;
            for (int aa=0;aa<temp_count;aa++){
                cout<<temp_buffer[aa]<<",";
            }cout<<endl;*/
            std::swap(temp_buffer, valid_candidate_index[depth]);
            valid_candidates_count = temp_count;
        }

        idx_count[depth] = valid_candidates_count;
    }
#endif
}

enumResult EvaluateQuery::exploreGraphQLStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                              ui *candidates_count, ui *order,
                                              size_t output_limit_num, size_t &call_count)
{
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    enumResult s;
    set<ui> true_sets[query_graph->getVerticesCount()];
    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        s.candidate_true.push_back(set<ui>());
    }

    // Generate the bn.
    ui **bn;
    ui *bn_count;

    bn = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        bn[i] = new ui[max_depth];
    }

    bn_count = new ui[max_depth];
    std::fill(bn_count, bn_count + max_depth, 0);

    std::vector<bool> visited_query_vertices(max_depth, false);
    visited_query_vertices[start_vertex] = true;
    for (ui i = 1; i < max_depth; ++i)
    {
        VertexID cur_vertex = order[i];
        ui nbr_cnt;
        const VertexID *nbrs = query_graph->getVertexNeighbors(cur_vertex, nbr_cnt);

        for (ui j = 0; j < nbr_cnt; ++j)
        {
            VertexID nbr = nbrs[j];

            if (visited_query_vertices[nbr])
            {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_query_vertices[cur_vertex] = true;
    }

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    for (ui i = 0; i < max_depth; ++i)
    {
        VertexID cur_vertex = order[i];
        ui max_candidate_count = candidates_count[cur_vertex];
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;

                for (int i = 0; i < max_depth; i++)
                {
                    s.candidate_true[order[i]].insert(embedding[i]);
                }

                for (int i = 0; i < max_depth; i++)
                {
                    true_sets[order[i]].insert(embedding[i]);
                }

                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, candidates, candidates_count);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

// Release the buffer.
EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    int true_cand_sum = 0;

    for (int i = 0; i < query_graph->getVerticesCount(); i++)
    {
        true_cand_sum += true_sets[i].size();
    }

    s.embedding_cnt = embedding_cnt;
    s.candidate_true_count_sum = true_cand_sum;

    return s;
}

void EvaluateQuery::generateValidCandidates(const Graph *data_graph, ui depth, ui *embedding, ui *idx_count,
                                            ui **valid_candidate, bool *visited_vertices, ui **bn, ui *bn_cnt,
                                            ui *order, ui **candidates, ui *candidates_count)
{
    VertexID u = order[depth];

    idx_count[depth] = 0;

    for (ui i = 0; i < candidates_count[u]; ++i)
    {
        VertexID v = candidates[u][i];

        if (!visited_vertices[v])
        {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j)
            {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v))
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

size_t EvaluateQuery::exploreQuickSIStyle(const Graph *data_graph, const Graph *query_graph, ui **candidates,
                                          ui *candidates_count, ui *order,
                                          ui *pivot, size_t output_limit_num, size_t &call_count)
{
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    // Generate the bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, pivot, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    VertexID **valid_candidate;
    bool *visited_vertices;

    idx = new ui[max_depth];
    idx_count = new ui[max_depth];
    embedding = new ui[max_depth];
    visited_vertices = new bool[data_graph->getVerticesCount()];
    std::fill(visited_vertices, visited_vertices + data_graph->getVerticesCount(), false);
    valid_candidate = new ui *[max_depth];

    ui max_candidate_count = data_graph->getGraphMaxLabelFrequency();
    for (ui i = 0; i < max_depth; ++i)
    {
        valid_candidate[i] = new VertexID[max_candidate_count];
    }

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];
    std::copy(candidates[start_vertex], candidates[start_vertex] + candidates_count[start_vertex],
              valid_candidate[cur_depth]);

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidate[cur_depth][idx[cur_depth]];
            embedding[u] = v;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(query_graph, data_graph, cur_depth, embedding, idx_count, valid_candidate,
                                        visited_vertices, bn, bn_count, order, pivot);
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
            visited_vertices[embedding[order[cur_depth]]] = false;
    }

// Release the buffer.
EXIT:
    delete[] bn_count;
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i)
    {
        delete[] bn[i];
        delete[] valid_candidate[i];
    }

    delete[] bn;
    delete[] valid_candidate;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(const Graph *query_graph, const Graph *data_graph, ui depth, ui *embedding,
                                            ui *idx_count, ui **valid_candidate, bool *visited_vertices, ui **bn,
                                            ui *bn_cnt,
                                            ui *order, ui *pivot)
{
    VertexID u = order[depth];
    LabelID u_label = query_graph->getVertexLabel(u);
    ui u_degree = query_graph->getVertexDegree(u);

    idx_count[depth] = 0;

    VertexID p = embedding[pivot[depth]];
    ui nbr_cnt;
    const VertexID *nbrs = data_graph->getVertexNeighbors(p, nbr_cnt);

    for (ui i = 0; i < nbr_cnt; ++i)
    {
        VertexID v = nbrs[i];

        if (!visited_vertices[v] && u_label == data_graph->getVertexLabel(v) &&
            u_degree <= data_graph->getVertexDegree(v))
        {
            bool valid = true;

            for (ui j = 0; j < bn_cnt[depth]; ++j)
            {
                VertexID u_nbr = bn[depth][j];
                VertexID u_nbr_v = embedding[u_nbr];

                if (!data_graph->checkEdgeExistence(v, u_nbr_v))
                {
                    valid = false;
                    break;
                }
            }

            if (valid)
            {
                valid_candidate[depth][idx_count[depth]++] = v;
            }
        }
    }
}

size_t EvaluateQuery::exploreDPisoStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                        Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                        ui **weight_array, ui *order, size_t output_limit_num,
                                        size_t &call_count)
{
    int max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }

    // Generate backward neighbors.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    VertexID start_vertex = order[0];
    std::vector<dpiso_min_pq> vec_rank_queue;

    for (ui i = 0; i < candidates_count[start_vertex]; ++i)
    {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
        reverse_embedding[v] = start_vertex;
#endif
        vec_rank_queue.emplace_back(dpiso_min_pq(extendable_vertex_compare));
        updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer, weight_array,
                               tree, start_vertex, extendable,
                               vec_rank_queue, query_graph);

        VertexID u = vec_rank_queue.back().top().first.first;
        vec_rank_queue.back().pop();

#ifdef ENABLE_FAILING_SET
        if (idx_count[u] == 0)
        {
            vec_failing_set[cur_depth] = ancestors[u];
        }
        else
        {
            vec_failing_set[cur_depth].reset();
        }
#endif

        call_count += 1;
        cur_depth += 1;
        order[cur_depth] = u;
        idx[u] = 0;
        while (cur_depth > 0)
        {
            while (idx[u] < idx_count[u])
            {
                ui valid_idx = valid_candidate_idx[u][idx[u]];
                v = candidates[u][valid_idx];

                if (visited_vertices[v])
                {
                    idx[u] += 1;
#ifdef ENABLE_FAILING_SET
                    vec_failing_set[cur_depth] = ancestors[u];
                    vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                    continue;
                }
                embedding[u] = v;
                idx_embedding[u] = valid_idx;
                visited_vertices[v] = true;
                idx[u] += 1;

#ifdef ENABLE_FAILING_SET
                reverse_embedding[v] = u;
#endif

                if (cur_depth == max_depth - 1)
                {
                    embedding_cnt += 1;
                    visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                    reverse_embedding.erase(embedding[u]);
                    vec_failing_set[cur_depth].set();
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];

#endif
                    if (embedding_cnt >= output_limit_num)
                    {
                        goto EXIT;
                    }
                }
                else
                {
                    call_count += 1;
                    cur_depth += 1;
                    vec_rank_queue.emplace_back(vec_rank_queue.back());
                    updateExtendableVertex(idx_embedding, idx_count, valid_candidate_idx, edge_matrix, temp_buffer,
                                           weight_array, tree, u, extendable,
                                           vec_rank_queue, query_graph);

                    u = vec_rank_queue.back().top().first.first;
                    vec_rank_queue.back().pop();
                    idx[u] = 0;
                    order[cur_depth] = u;

#ifdef ENABLE_FAILING_SET
                    if (idx_count[u] == 0)
                    {
                        vec_failing_set[cur_depth - 1] = ancestors[u];
                    }
                    else
                    {
                        vec_failing_set[cur_depth - 1].reset();
                    }
#endif
                }
            }

            cur_depth -= 1;
            vec_rank_queue.pop_back();
            u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
            restoreExtendableVertex(tree, u, extendable);
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[u] = idx_count[u];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
        }
    }

// Release the buffer.
EXIT:
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

void EvaluateQuery::updateExtendableVertex(ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                           Edges ***edge_matrix, ui *&temp_buffer, ui **weight_array,
                                           TreeNode *tree, VertexID mapped_vertex, ui *extendable,
                                           std::vector<dpiso_min_pq> &vec_rank_queue, const Graph *query_graph)
{
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0)
        {
            // generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
            //                             tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j)
            {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            vec_rank_queue.back().emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }
}

void EvaluateQuery::restoreExtendableVertex(TreeNode *tree, VertexID unmapped_vertex, ui *extendable)
{
    TreeNode &node = tree[unmapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] += 1;
    }
}

void EvaluateQuery::generateValidCandidateIndex(VertexID u, ui *idx_embedding, ui *idx_count, ui *&valid_candidate_index,
                                                Edges ***edge_matrix, ui *bn, ui bn_cnt, ui *&temp_buffer)
{
    VertexID previous_bn = bn[0];
    Edges &previous_edge = *edge_matrix[previous_bn][u];
    ui previous_index_id = idx_embedding[previous_bn];

    ui previous_candidates_count =
        previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui *previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    ui valid_candidates_count = 0;
    for (ui i = 0; i < previous_candidates_count; ++i)
    {
        valid_candidate_index[valid_candidates_count++] = previous_candidates[i];
    }

    ui temp_count;
    for (ui i = 1; i < bn_cnt; ++i)
    {
        VertexID current_bn = bn[i];
        Edges &current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count =
            current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui *current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index,
                                                  valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index);
        valid_candidates_count = temp_count;
    }

    idx_count[u] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, TreeNode *tree, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];
        ancestors[u].set(u);
        for (ui j = 0; j < u_node.bn_count_; ++j)
        {
            VertexID u_bn = u_node.bn_[j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

size_t EvaluateQuery::exploreDPisoRecursiveStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree,
                                                 Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                                                 ui **weight_array, ui *order, size_t output_limit_num,
                                                 size_t &call_count)
{
    int max_depth = query_graph->getVerticesCount();

    ui *extendable = new ui[max_depth];
    for (ui i = 0; i < max_depth; ++i)
    {
        extendable[i] = tree[i].bn_count_;
    }

    // Generate backward neighbors.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    // Evaluate the query.
    size_t embedding_cnt = 0;

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, tree, order, ancestors);

    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
    VertexID start_vertex = order[0];

    for (ui i = 0; i < candidates_count[start_vertex]; ++i)
    {
        VertexID v = candidates[start_vertex][i];
        embedding[start_vertex] = v;
        idx_embedding[start_vertex] = i;
        visited_vertices[v] = true;
        reverse_embedding[v] = start_vertex;
        call_count += 1;

        exploreDPisoBacktrack(max_depth, 1, start_vertex, tree, idx_embedding, embedding, reverse_embedding,
                              visited_vertices, idx_count, valid_candidate_idx, edge_matrix,
                              ancestors, dpiso_min_pq(extendable_vertex_compare), weight_array, temp_buffer, extendable,
                              candidates, embedding_cnt,
                              call_count, nullptr);

        visited_vertices[v] = false;
        reverse_embedding.erase(v);
    }

    // Release the buffer.
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>
EvaluateQuery::exploreDPisoBacktrack(ui max_depth, ui depth, VertexID mapped_vertex, TreeNode *tree, ui *idx_embedding,
                                     ui *embedding, std::unordered_map<VertexID, VertexID> &reverse_embedding,
                                     bool *visited_vertices, ui *idx_count, ui **valid_candidate_index,
                                     Edges ***edge_matrix,
                                     std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                     dpiso_min_pq rank_queue, ui **weight_array, ui *&temp_buffer, ui *extendable,
                                     ui **candidates, size_t &embedding_count, size_t &call_count,
                                     const Graph *query_graph)
{
    // Compute extendable vertex.
    TreeNode &node = tree[mapped_vertex];
    for (ui i = 0; i < node.fn_count_; ++i)
    {
        VertexID u = node.fn_[i];
        extendable[u] -= 1;
        if (extendable[u] == 0)
        {
            // generateValidCandidateIndex(u, idx_embedding, idx_count, valid_candidate_index[u], edge_matrix, tree[u].bn_,
            //                             tree[u].bn_count_, temp_buffer);

            ui weight = 0;
            for (ui j = 0; j < idx_count[u]; ++j)
            {
                ui idx = valid_candidate_index[u][j];
                weight += weight_array[u][idx];
            }
            rank_queue.emplace(std::make_pair(std::make_pair(u, query_graph->getVertexDegree(u)), weight));
        }
    }

    VertexID u = rank_queue.top().first.first;
    rank_queue.pop();

    if (idx_count[u] == 0)
    {
        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return ancestors[u];
    }
    else
    {
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> current_fs;
        std::bitset<MAXIMUM_QUERY_GRAPH_SIZE> child_fs;

        for (ui i = 0; i < idx_count[u]; ++i)
        {
            ui valid_index = valid_candidate_index[u][i];
            VertexID v = candidates[u][valid_index];

            if (!visited_vertices[v])
            {
                embedding[u] = v;
                idx_embedding[u] = valid_index;
                visited_vertices[v] = true;
                reverse_embedding[v] = u;
                if (depth != max_depth - 1)
                {
                    call_count += 1;
                    child_fs = exploreDPisoBacktrack(max_depth, depth + 1, u, tree, idx_embedding, embedding,
                                                     reverse_embedding, visited_vertices, idx_count,
                                                     valid_candidate_index, edge_matrix,
                                                     ancestors, rank_queue, weight_array, temp_buffer, extendable,
                                                     candidates, embedding_count,
                                                     call_count, query_graph);
                }
                else
                {
                    embedding_count += 1;
                    child_fs.set();
                }
                visited_vertices[v] = false;
                reverse_embedding.erase(v);

                if (!child_fs.test(u))
                {
                    current_fs = child_fs;
                    break;
                }
            }
            else
            {
                child_fs.reset();
                child_fs |= ancestors[u];
                child_fs |= ancestors[reverse_embedding[v]];
            }

            current_fs |= child_fs;
        }

        restoreExtendableVertex(tree, mapped_vertex, extendable);
        return current_fs;
    }
}

size_t
EvaluateQuery::exploreCECIStyle(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count)
{

    int max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (int i = 0; i < max_depth; ++i)
    {
        if (candidates_count[i] > max_valid_candidates_count)
        {
            max_valid_candidates_count = candidates_count[i];
        }
    }
    // Allocate the memory buffer.
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (int i = 0; i < max_depth; ++i)
    {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    // Evaluate the query.
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i)
    {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

#ifdef ENABLE_FAILING_SET
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);
#endif

    while (true)
    {
        while (idx[cur_depth] < idx_count[cur_depth])
        {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v])
            {
#ifdef ENABLE_FAILING_SET
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

#ifdef ENABLE_FAILING_SET
            reverse_embedding[v] = u;
#endif
            if (cur_depth == max_depth - 1)
            {
                embedding_cnt += 1;
                visited_vertices[v] = false;
#ifdef ENABLE_FAILING_SET
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
#endif
                if (embedding_cnt >= output_limit_num)
                {
                    goto EXIT;
                }
            }
            else
            {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);
#ifdef ENABLE_FAILING_SET
                if (idx_count[cur_depth] == 0)
                {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                }
                else
                {
                    vec_failing_set[cur_depth - 1].reset();
                }
#endif
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else
        {
            VertexID u = order[cur_depth];
#ifdef ENABLE_FAILING_SET
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0)
            {
                if (!vec_failing_set[cur_depth].test(u))
                {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                }
                else
                {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
#endif
            visited_vertices[embedding[u]] = false;
        }
    }

// Release the buffer.
EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (int i = 0; i < max_depth; ++i)
    {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return embedding_cnt;
}

void EvaluateQuery::generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates)
{

    VertexID u = order[depth];
    TreeNode &u_node = tree[u];
    idx_count[depth] = 0;
    ui valid_candidates_count = 0;
    {
        VertexID u_p = tree[u].parent_;
        VertexID v_p = embedding[u_p];

        auto iter = TE_Candidates[u].find(v_p);
        if (iter == TE_Candidates[u].end() || iter->second.empty())
        {
            return;
        }

        valid_candidates_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < valid_candidates_count; ++i)
        {
            valid_candidates[depth][i] = v_p_nbrs[i];
        }
    }
    ui temp_count;
    for (ui i = 0; i < tree[u].bn_count_; ++i)
    {
        VertexID u_p = tree[u].bn_[i];
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty())
        {
            return;
        }

        ui current_candidates_count = iter->second.size();
        ui *current_candidates = iter->second.data();

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  valid_candidates[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidates[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, ui **bn, ui *bn_cnt, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < bn_cnt[i]; ++j)
        {
            VertexID u_bn = bn[i][j];
            ancestors[u] |= ancestors[u_bn];
        }
    }
}

void EvaluateQuery::computeAncestor(const Graph *query_graph, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors)
{
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    // Compute the ancestor in the top-down order.
    for (ui i = 0; i < query_vertices_num; ++i)
    {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < i; ++j)
        {
            VertexID u_bn = order[j];
            if (query_graph->checkEdgeExistence(u, u_bn))
            {
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
}
