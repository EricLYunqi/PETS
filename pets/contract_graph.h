#ifndef _CONTRACT_GRAPH_H_
#define _CONTRACT_GRAPH_H_

#include "config/config.h"
#include "config/type.h"
#include "graph/graph.h"
#include "util/set_intersection.hpp"
#include <chrono>
#include <random>
#include <unordered_map>
#include <assert.h>

typedef void (*ContractionFunc) (const Graph*, std::vector<std::vector<ui>>&, bool*);  
typedef std::vector<std::pair<ContractionFunc, SuperNodeType>> ContractionOrder;


class ContractGraph
{
public:
    static ui MAXCHORDSIZE;
    static ui MAXCYCLESIZE;
    static ui MAXPATHSIZE;
    static ui MAXSTARSIZE;
    static std::unordered_map<SuperNodeType, std::string> nodeMap;

public:
    ContractGraph() = default;
    ~ContractGraph() = default;
    ContractGraph(const ContractGraph &other) = delete;
    ContractGraph(ContractGraph &&other) = delete;

private:
    // contraction
    static void extractChord(const Graph *graph, std::vector<std::vector<ui>> &chords, 
                             bool *visited);
    static void extractCycle(const Graph *graph, std::vector<std::vector<ui>> &cycles, 
                             bool *visited);
    static void extractPath(const Graph *graph, std::vector<std::vector<ui>> &paths, 
                            bool *visited);
    static void extractStar(const Graph *graph, std::vector<std::vector<ui>> &stars, 
                            bool *visited);
    static void extractSingleton(const Graph *graph, std::vector<std::vector<ui>> &singleton, 
                                 bool *visited);

    // contraction utils
    // Pathon's algorithm
    static void generateFundmentalCycles(const Graph *graph, std::vector<std::vector<ui>> &fcycles, 
                                         bool *visited);
    // random walk
    static void generateOnePath(const Graph *graph, std::vector<ui> &onepath, bool *visited, 
                                ui root, ui curLen, std::vector<bool> &vis, bool &finished, 
                                std::default_random_engine &eng);

    // contraction verification
    static bool verifyCycle(const Graph *graph, const std::vector<ui> &onecycle);
    static bool verifyPath(const Graph *graph, const std::vector<ui> &onepath);

    // memory management
    static void allocateBuffers(const Graph *graph, bool *&visited);
    static void releaseBuffers(const Graph *graph, bool *&visited);

    // building utils
    // update super node basic infos
    // cpv: value for checking whether different super nodes share the same vertex
    static void updateSuperNode(SuperGraph *superGraph, const std::vector<std::vector<ui>> &nodesets, 
                                ui &superNodeId, ui cpv);
    // update super node indexes
    static void updateSuperNodeIndex(SuperGraph *superGraph, const std::vector<std::vector<ui>> &nodesets,
                                     const Graph *graph, SuperNodeType sntype);

    static void generateContractionPlan(GraphType graphType, ContractionOrder &order);

    static bool selfCheck(const SuperGraph *superGraph);

public:
    static void buildSuperGraph(const Graph *graph, GraphType graphType, SuperGraph *superGraph);
};

#endif // _CONTRACT_GRAPH_H_