#ifndef _GEN_QUERY_HPP_
#define _GEN_QUERY_HPP_

#include "config/config.h"
#include "config/type.h"
#include "graph/graph.h"
#include "graph/pattern.h"
#include <random>
#include <vector>
#include <unordered_map>
#include <assert.h>

namespace generate_query
{

inline ui getRand(ui s, ui t) 
{
    return s + (rand() * (long long)rand() + rand() * (long long)rand()) % (t - s);
}

// querySize is in terms of # of edges
inline void generateOneSizeRandomQueries(ui querySize, ui numQueries, std::vector<std::vector<Pattern>> &queriesUniverse, 
                                         const Graph *graph)
{
    std::random_device rd;
    std::default_random_engine engine(rd());
    srand((unsigned)time(NULL));

    for(ui turn = 1; turn <= numQueries; turn++) {
        std::vector<std::vector<bool>> quickRef(querySize + 1, std::vector<bool>(querySize + 1, false));
        std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);
        Pattern pattern(querySize + 1, 0);
        std::unordered_map<ui, ui> nodeMap;
        std::vector<ui> vertices;

        ui startV = dis(engine);
        ui nodeId = 1;
        vertices.emplace_back(startV);
        nodeMap[startV] = 0;

        bool isFoundPattern = false;
        for(ui i = 0; i < querySize; i++) {
            bool isFoundEdge = false;

            for(ui tryTimes = 0; tryTimes < 100; tryTimes++) {
                std::uniform_int_distribution<ui> dis1(0, vertices.size() - 1);
                ui s = vertices[dis1(engine)];
                // 80% percentage for inner edge
                bool isInnerEdge = rand() % 100 < 80 ? true : false;

                ui e = s;
                ui count_ = 0;
                const ui *neighbors_ = graph->getVertexNeighbors(s, count_);
                
                ui rd = (ui)rand() % count_;
                if(isInnerEdge) {
                    for(ui j = 0; j < count_; j++) {
                        ui candidate = neighbors_[rd];
                        if(nodeMap.find(candidate) != nodeMap.end() && !quickRef[nodeMap[s]][nodeMap[candidate]]) {
                            e = candidate;
                            break;
                        }
                    }
                }
                else 
                    e = neighbors_[rd];

                s = nodeMap[s];
                if (nodeMap.find(e) == nodeMap.end()) {
                    vertices.emplace_back(e);
                    nodeMap[e] = nodeId++;
                }
                e = nodeMap[e];

                if (e == s || quickRef[s][e])
                    continue;

                pattern.addEdge(s, e);
                quickRef[e][s] = quickRef[s][e] = 1;
                isFoundEdge = true;
                break;
            }

            if (!isFoundEdge) {
                isFoundPattern = false;
                break;
            }
        }

        if (isFoundPattern) {
            pattern.numVertex = nodeId;
            queriesUniverse[querySize].emplace_back(pattern);
        }
        else 
            -- turn;
    }

    if(queriesUniverse[querySize].size() != numQueries) {
        std::cerr << "error in generating queries universe, size: " << querySize << " numQueries: " << numQueries
                  << "results size: " << queriesUniverse[querySize].size() <<  std::endl;
        exit(1);
    }
}


inline void generateRandomQueriesUniverse(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<std::vector<Pattern>> &queriesUniverse, 
                                          const Graph *graph)
{
    queriesUniverse.resize(sizeHigh + 1);
    for(ui size = sizeLow; size <= sizeHigh; size++) 
        generateOneSizeRandomQueries(size, numQueries, queriesUniverse, graph);

    std::cout << "queries universe |U|: " << numQueries * (sizeHigh - sizeLow + 1) << std::endl;
}


/*
 * apis used for generating queries sets
 */
inline void generateRandomQueries(const std::vector<std::vector<Pattern>> &queriesUniverse, ui sizeLow, ui sizeHigh,
                                  ui numQueries, std::vector<Pattern> &queries)
{
    if(sizeLow < sizeHigh)
        return;

    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui eachNumQueries = numQueries / numDivSizes;
    for (ui i = sizeLow; i <= sizeHigh; ++i)
        for (ui j = 0; j < eachNumQueries; ++j)
            queries.emplace_back(queriesUniverse[i][j]);
}


inline void genQueriesExtend(ui desiredSize, std::vector<int> &vis, int visFlag, ui &nodeId, std::vector<ui> &vertices,
                             std::unordered_map<ui, ui> &nodeMap, std::vector<std::vector<bool>> &quickRef,
                             Pattern &pattern, const Graph *graph) 
{
    ui tryTimes = 0;
    ui s = 0, e = 0;

    ++ desiredSize;

    while (desiredSize > 1 && tryTimes < 10) {
        if(nodeId != vertices.size()) {
            std::cerr << "extend error: " << nodeId << " " << vertices.size() << std::endl;
            exit(1);
        }

        bool find = false;
        if (rand() % 100 < 80) {  // 80% inner edge
            for (ui i = 0; i < nodeId && !find; i++) {
                for (ui j = i + 1; j < nodeId && !find; j++) {
                    if (!quickRef[i][j] && graph->checkEdgeExistence(vertices[i], vertices[j])) {
                        s = vertices[i];
                        e = vertices[j];
                        if(i != nodeMap[s] || j != nodeMap[e]) {
                            std::cerr << "error in nodeMap: " << i << " " << nodeMap[s] << " "
                                      << j << " " << nodeMap[e] << std::endl;
                            exit(1);
                        }
                        find = true;
                    }
                }
            }
        } 
        else {
            ui rdidx = rand() % (vertices.size());
            s = vertices[rdidx];

            ui count_ = 0;
            const ui *neighbors_ = graph->getVertexNeighbors(s, count_);

            if (!count_) {
                e = neighbors_[rand() % count_];
                if (vis[e] != visFlag) {
                    // node
                    vertices.emplace_back(e);
                    nodeMap[e] = nodeId ++;
                    vis[e] = visFlag;
                    find = true;
                }
            }
        }
        if (find) {
            // edge
            quickRef[nodeMap[s]][nodeMap[e]] = true;
            quickRef[nodeMap[e]][nodeMap[s]] = true;
            pattern.addEdge(nodeMap[s], nodeMap[e]);
            -- desiredSize;
            tryTimes = 0;
        } 
        else 
            ++ tryTimes;
    }
}


inline void generateTreeQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    for (ui size = sizeLow; size <= sizeHigh; ++size) {
        for (ui j = 1; j <= eachNumQueries; ++j) {
            Pattern pattern(size + 1, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 1, std::vector<bool>(size + 1, false));

            ui nodeId = 0;
            ui x = dis(engine);
            vertices.emplace_back(x);
            nodeMap[x] = nodeId ++;
            vis[x] = ++ visFlag;

            ui curSize = 0;
            ui tryTimes = 0;
            while (curSize + 1 < size && tryTimes < 10) {
                x = vertices[rand() % vertices.size()];

                ui count_ = 0;
                const ui *neighbors_ = graph->getVertexNeighbors(x, count_);

                if (!count_) {
                    ++ tryTimes;
                    continue;
                }

                ui next = neighbors_[rand() % count_];
                if (vis[next] != visFlag) {
                    vertices.emplace_back(next);
                    nodeMap[next] = nodeId ++;
                    quickRef[nodeMap[x]][nodeMap[next]] = true;
                    quickRef[nodeMap[next]][nodeMap[x]] = true;
                    pattern.addEdge(nodeMap[x], nodeMap[next]);
                    vis[next] = visFlag;
                    ++ curSize;
                    tryTimes = 0;
                } 
                else
                    ++ tryTimes;
            }

            if (curSize >= size / 2.0) {
                genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                                 nodeMap, quickRef, pattern, graph);
                pattern.numVertex = nodeId;
                queries.emplace_back(pattern);
            } 
            else 
                --j;
        }
    }
    
    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating tree queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


inline void generatePathQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    for (ui size = sizeLow; size <= sizeHigh; size++) {
        for (ui j = 1; j <= eachNumQueries; ++j) {
            Pattern pattern(size + 1, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 1, std::vector<bool>(size + 1, false));

            ui nodeId = 0;
            ui x = dis(engine);
            vertices.emplace_back(x);
            nodeMap[x] = nodeId ++;
            vis[x] = ++ visFlag;

            ui count_ = 0;
            const ui *neighbors_ = graph->getVertexNeighbors(x, count_);

            ui curSize = 0;
            ui tryTimes = 0;
            while (curSize + 1 < size && curSize < 15 && tryTimes < 10 && count_ > 0) {
                std::uniform_int_distribution<ui> dis2(0, count_ - 1);
                ui rdidx = dis2(engine);
                ui next = neighbors_[rdidx];

                if (vis[next] != visFlag) {
                    vertices.emplace_back(next);
                    nodeMap[next] = nodeId ++;
                    vis[next] = visFlag;
                    pattern.addEdge(nodeMap[x], nodeMap[next]);
                    quickRef[nodeMap[x]][nodeMap[next]] = true;
                    quickRef[nodeMap[next]][nodeMap[x]] = true;
                    ++ curSize;
                    tryTimes = 0;
                    x = next;
                } 
                else
                    ++ tryTimes;
            }

            if (curSize >= size / 2.0) {
                genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                                 nodeMap, quickRef, pattern, graph);
                pattern.numVertex = nodeId;
                queries.emplace_back(pattern);
            } 
            else 
                -- j;
        }
    }

    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating path queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


inline void generateStarQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    int try2 = 0;
    for (ui size = sizeLow; size <= sizeHigh; size++) {
        for (ui j = 1; j <= eachNumQueries; j++) {
            Pattern pattern(size + 1, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 1, std::vector<bool>(size + 1, false));

            ui nodeId = 0;
            ui x = dis(engine);
            vertices.emplace_back(x);
            nodeMap[x] = nodeId ++;
            vis[x] = ++ visFlag;

            ui count_ = 0;
            const ui *neighbors_ = graph->getVertexNeighbors(x, count_);

            ui curSize = 0;
            ui tryTimes = 0;

            while (curSize + 1 < size - 1 && tryTimes < 10 && count_ > 0) {  // -1 ==> not pure star
                ui next = neighbors_[rand() % count_];
                if (vis[next] != visFlag) {
                    vertices.emplace_back(next);
                    nodeMap[next] = nodeId ++;
                    quickRef[nodeMap[x]][nodeMap[next]] = true;
                    quickRef[nodeMap[next]][nodeMap[x]] = true;
                    pattern.addEdge(nodeMap[x], nodeMap[next]);
                    vis[next] = visFlag;
                    ++ curSize;
                    tryTimes = 0;
                } 
                else
                    ++ tryTimes;
            }

            if (curSize >= size / 2.0 || try2 > 10) {
                genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                                 nodeMap, quickRef, pattern, graph);
                pattern.numVertex = nodeId;
                queries.emplace_back(pattern);
                try2 = 0;
            } 
            else {
                -- j;
                ++ try2;
            }
        }
    }

    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating star queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


inline void generateDoubleStarQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    int try2 = 0;
    for (ui size = sizeLow; size <= sizeHigh; ++size) {
        for (ui j = 1; j <= eachNumQueries; ++j) {
            Pattern pattern(size + 1, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 1, std::vector<bool>(size + 1, false));

            ui nodeId = 0;
            ui x = dis(engine);
            vertices.emplace_back(x);
            nodeMap[x] = nodeId ++;
            vis[x] = ++ visFlag;

            ui count_ = 0;
            const ui *neighbors_ = graph->getVertexNeighbors(x, count_);

            ui curSize = 0;
            ui tryTimes = 0;
            bool secondCore = false;

            while (curSize + 1 < size && tryTimes < 10 && count_ > 0) {  // -1 ==> not pure star
                ui next = neighbors_[rand() % count_];
                if (vis[next] != visFlag) {
                    vertices.emplace_back(next);
                    nodeMap[next] = nodeId ++;
                    quickRef[nodeMap[x]][nodeMap[next]] = true;
                    quickRef[nodeMap[next]][nodeMap[x]] = true;
                    pattern.addEdge(nodeMap[x], nodeMap[next]);
                    vis[next] = visFlag;
                    ++ curSize;
                    tryTimes = 0;
                    if (curSize >= size / 2.0 && !secondCore && rand() % 100 >= 60) {
                        secondCore = true;
                        x = next;
                    }
                } 
                else
                    ++ tryTimes;
            }

            if (curSize >= size / 2.0 || try2 > 10) {
                genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                                nodeMap, quickRef, pattern, graph);
                pattern.numVertex = nodeId;
                queries.emplace_back(pattern);
                try2 = 0;
            } 
            else {
                -- j;
                ++ try2;
            }
        }
    }
    
    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating double star queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


inline void generateFlowerQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    int try2 = 0;
    for (ui size = sizeLow; size <= sizeHigh; size++) {
        for (ui j = 1; j <= eachNumQueries; ++j) {
            Pattern pattern(size + 10, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 10, std::vector<bool>(size + 10, false));

            ui nodeId = 0;
            ui x = dis(engine);

            ui count_ = 0;
            const ui *neighbors_ = graph->getVertexNeighbors(x, count_);

            if (!count_) {
                -- j;
                continue;
            }

            vertices.emplace_back(x);
            nodeMap[x] = nodeId ++;
            vis[x] = ++ visFlag;

            ui curSize = 0;
            ui tryTimes = 0;
            int k = -1;

            while (k == -1 && tryTimes < 10) {
                ui next = neighbors_[rand() % count_];
                ui nextCount_ = 0;
                const ui *nextNeighbors_ = graph->getVertexNeighbors(next, nextCount_);
                ui common = 0;
                compute_set_intersection::computeCandidates(neighbors_, count_, 
                                                            nextNeighbors_, nextCount_, common);
                if(common < 3) {
                    ++ tryTimes;
                    continue;
                }

                vertices.emplace_back(next);
                nodeMap[next] = ++ nodeId;
                quickRef[nodeMap[x]][nodeMap[next]] = true;
                quickRef[nodeMap[next]][nodeMap[x]] = true;
                vis[next] = visFlag;

                ui kRange = std::min((size + 1) / 2, common);
                k = rand() % ((int)kRange - 1) + 2;
                curSize += k * 2;
            }

            if (k == -1) {
                -- j;
                continue;
            }

            genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                             nodeMap, quickRef, pattern, graph);

            for (int i = 0; i < k; ++i) {
                pattern.addEdge(0, nodeId);
                pattern.addEdge(1, nodeId);
                ++ nodeId;
            }

            if (curSize >= size / 2 || try2 > 10) {
                pattern.numVertex = nodeId;
                queries.emplace_back(pattern);
                try2 = 0;
            } 
            else {
                -- j;
                ++ try2;
            }
        }
    }

    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating flower queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


inline void collectAllPaths(const std::vector<std::vector<ui>> &bfsTree, std::vector<std::vector<ui>> &tempPath, 
                            ui start, ui end, std::vector<bool> &vis, std::vector<ui> &onePath, ui sizeHigh, 
                            ui sizeLow, int depth)
{
    if(depth > sizeHigh)
        return;

    vis[start] = true;
    if(start == end) {
        if(depth >= sizeLow)
            tempPath.emplace_back(onePath);
        return;
    }

    for(const auto &v : bfsTree[start]) {
        if(!vis[v]) {
            onePath.emplace_back(v);
            vis[v] = true;
            collectAllPaths(bfsTree, tempPath, v, end, vis, onePath, sizeHigh, sizeLow, depth + 1);
            onePath.pop_back();
        }
    }
}


inline void generateCyclesUniverse(std::vector<std::vector<ui>> &cyclesUniverse, const Graph *graph, ui sizeLow, ui sizeHigh)
{
    ui numVertex = graph->numVertex;
    std::vector<std::vector<ui>> bfsTree(numVertex, std::vector<ui> ());
    std::vector<std::vector<ui>> nonTreeEdges(numVertex, std::vector<ui> ());
    std::vector<bool> vis(numVertex, false);
    std::queue<ui> q;

    for(ui i = 0; i < numVertex; i++) {
        ui count_ = 0;
        const ui *neighbors_ = graph->getVertexNeighbors(i, count_);
        for(ui j = 0; j < count_; j++)
            if(neighbors_[j] > i)
                nonTreeEdges[i].emplace_back(neighbors_[j]);
        std::sort(nonTreeEdges[i].begin(), nonTreeEdges[i].end());
    }

    q.emplace(0);
    vis[0] = true;

    while(!q.empty()) {
        ui pivot = q.front();
        q.pop();

        ui count_ = 0;
        const ui *neighbors_ = graph->getVertexNeighbors(pivot, count_);

        for(ui i = 0; i < count_; i++) {
            ui next = neighbors_[i];
            if(!vis[next]) {
                q.emplace(next);
                vis[next] = true;
                bfsTree[pivot].emplace_back(next);
                bfsTree[next].emplace_back(pivot);
                if(pivot > next) {
                    auto iter2 = std::lower_bound(nonTreeEdges[next].begin(), nonTreeEdges[next].end(), pivot);
                    nonTreeEdges[next].erase(iter2);
                }
                else {
                    auto iter1 = std::lower_bound(nonTreeEdges[pivot].begin(), nonTreeEdges[pivot].end(), next);
                    nonTreeEdges[pivot].erase(iter1);
                }
            }
        }
    }

    for(auto &graList : bfsTree)
        std::sort(graList.begin(), graList.end());

    std::vector<std::vector<ui>> tempPath;
    std::vector<ui> onePath;
    for(ui i = 0; i < numVertex; i++) {
        for(const auto &v : nonTreeEdges[i]) {
            onePath.clear();
            std::vector<bool> vis2(numVertex, false);
            onePath.emplace_back(i);
            collectAllPaths(bfsTree, tempPath, i, v, vis2, onePath, sizeHigh, sizeLow, 0);
            
            for(const auto &path : tempPath)
                cyclesUniverse.emplace_back(std::move(path));
        }
    }
}


inline void generateCycleQueries(ui sizeLow, ui sizeHigh, ui numQueries, std::vector<Pattern> &queries, const Graph *graph)
{
    std::vector<std::vector<ui>> cyclesUniverse;
    generateCyclesUniverse(cyclesUniverse, graph, sizeLow, sizeHigh);

    ui numDivSizes = sizeHigh - sizeLow + 1;
    if(numQueries % numDivSizes != 0)
        numQueries = ceil(numQueries * 1.0 / numDivSizes - 1e-5) * numDivSizes;

    ui initialSize = queries.size();
    ui eachNumQueries = numQueries / numDivSizes;

    std::vector<int> vis(graph->numVertex, -1);
    int visFlag = 0;

    std::random_device rd;
    std::default_random_engine engine(rd());
    std::uniform_int_distribution<ui> dis(0, graph->numVertex - 1);

    for (ui size = sizeLow; size <= sizeHigh; ++size) {
        for (ui j = 1; j < eachNumQueries; j++) {
            Pattern pattern(size + 1, 0);  // not real size
            std::vector<ui> vertices;
            std::unordered_map<ui, ui> nodeMap;
            std::vector<std::vector<bool>> quickRef(size + 1, std::vector<bool> (size + 1, false));

            ui nodeId = 0;
            ui curSize = 0;
            ui rdidx = rand() % cyclesUniverse.size();
            ui length = cyclesUniverse[rdidx].size();

            if (length > size) {
                -- j;
                continue;
            }

            ++ visFlag;
            for(const auto &x : cyclesUniverse[rdidx]) {
                vertices.emplace_back(x);
                nodeMap[x] = nodeId ++;
                vis[x] = visFlag;
            }
            
            for (ui i = 0; i < length; i++) {
                ui j = (i + 1) % length;
                ui x = cyclesUniverse[rdidx][i];
                ui next = cyclesUniverse[rdidx][j];
                quickRef[nodeMap[x]][nodeMap[next]] = true;
                quickRef[nodeMap[next]][nodeMap[x]] = true;
                pattern.addEdge(nodeMap[x], nodeMap[next]);
                ++ curSize;
            }
            assert(nodeId == curSize && nodeId == length);

            ui x = cyclesUniverse[rdidx][rand() % cyclesUniverse[rdidx].size()];
            genQueriesExtend(size - curSize, vis, visFlag, nodeId, vertices,
                             nodeMap, quickRef, pattern, graph);
            pattern.numVertex = nodeId;
            queries.emplace_back(pattern);
        }
    }

    if(queries.size() != numQueries + initialSize) {
        std::cerr << "error in generating cycle queries, required numQueries: " << numQueries + initialSize
                  << "results size: " << queries.size() <<  std::endl;
        exit(1);
    }
}


/*
 * you shoud call
 * universeEachSize: # of queries for each size
 */
inline void generateQueries(ui sizeLow, ui sizeHigh, ui universeEachSize, const std::unordered_map<std::string, ui> &numQueries, 
                            std::vector<Pattern> &queries, const Graph *graph)
{
    // generate universe
    std::vector<std::vector<Pattern>> queriesUniverse;
    generateRandomQueriesUniverse(sizeLow, sizeHigh, universeEachSize, queriesUniverse, graph);

    // generate desired queries
    ui numRandom = numQueries.at("random");
    std::cout << "generating random queries ..." << std::endl << std::flush;
    generateRandomQueries(queriesUniverse, sizeLow, sizeHigh, numRandom, queries);
    ui numPath = numQueries.at("path");
    std::cout << "generating path queries ..." << std::endl << std::flush;
    generatePathQueries(sizeLow, sizeHigh, numPath, queries, graph);
    ui numTree = numQueries.at("tree");
    std::cout << "generating tree queries ..." << std::endl << std::flush;
    generateTreeQueries(sizeLow, sizeHigh, numTree, queries, graph);
    ui numStar = numQueries.at("star");
    std::cout << "generating star queries ..." << std::endl << std::flush;
    generateStarQueries(sizeLow, sizeHigh, numStar, queries, graph);
    ui numDStar = numQueries.at("doublestar");
    std::cout << "generating double star queries ..." << std::endl << std::flush;
    generateDoubleStarQueries(sizeLow, sizeHigh, numDStar, queries, graph);
    ui numFlower = numQueries.at("flower");
    std::cout << "generating flower queries ..." << std::endl << std::flush;
    generateFlowerQueries(sizeLow, sizeHigh, numFlower, queries, graph);
    ui numCycle = numQueries.at("cycle");
    std::cout << "generating cycle queries ..." << std::endl << std::flush;
    generateCycleQueries(sizeLow, sizeHigh, numCycle, queries, graph);

    // report stat
    ui totalSize = numRandom + numPath + numTree + numStar + numDStar + numFlower + numCycle;
    std::cout << "### total generated queries: " << totalSize << std::endl;
    std::cout << "###          random queries: " << numRandom << std::endl;
    std::cout << "###            path queries: " << numPath << std::endl;
    std::cout << "###            tree queries: " << numTree << std::endl;
    std::cout << "###            star queries: " << numStar << std::endl;
    std::cout << "###     double star queries: " << numDStar << std::endl;
    std::cout << "###          flower queries: " << numFlower << std::endl;
    std::cout << "###           cycle queries: " << numCycle << std::endl;
    std::cout << "generation done" << std::endl;
}

};

#endif // _GEN_QUERY_HPP_