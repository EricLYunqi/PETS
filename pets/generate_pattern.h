#ifndef _GENERATE_PATTERN_H_
#define _GENERATE_PATTERN_H_

#include "config/config.h"
#include "config/type.h"
#include "graph/graph.h"
#include "graph/pattern.h"
#include <numeric>
#include <random>
#include <unordered_set>


class GeneratePattern
{
public:
    GeneratePattern() = default;
    ~GeneratePattern() = default;
    GeneratePattern(const GeneratePattern &other) = delete;
    GeneratePattern(GeneratePattern &&other) = delete;

public:
    static void allocateBuffers(const SuperGraph *superGraph, std::vector<Pattern> *&patterns);
    static void releaseBuffers(std::vector<Pattern> *&patterns);

    /*
     * Algorithm 2, E(v_{Hi}, v_{Hj}) <- f_D(v_{Hi}, v_{Hj})
     * At this moment, we only decontract one super edge each time
     * Each v_{Hi} has pattern set RP_i
     */
    static void generatePatterns(const SuperGraph *superGraph, std::vector<Pattern> *patterns);

private:
    static ui getRandomIdx(ui lower, ui upper);
    // grow a subgraph random walk
    static void growSubgraph(const SuperGraph *superGraph, std::vector<Pattern> &patterns,
                             SubGraph &g, std::vector<ui> &depth, bool newV, 
                             std::vector<std::vector<ui>> &neighbors_, 
                             std::vector<ui> &pivots);
    // grow a subgraph by adding all possible edges
    static void growSubgraphAll(const SuperGraph *superGraph, std::vector<Pattern> &patterns,
                                SubGraph &g, std::vector<ui> &depth, bool newV, 
                                std::vector<std::vector<ui>> &neighbors_, 
                                std::vector<ui> &pivots);

    // hash a edge pair<ui, ui> to size_t
    static size_t hashEdge(const std::pair<ui, ui> &edge);
    // maintain DRP_i
    // max coverage problem: greedily select the g_{i, x} who has the maximum # of uncovered edges
    static void selectTopK(std::vector<Pattern> &patterns);
};


#endif // _GENERATE_PATTERN_H_