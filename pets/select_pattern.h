#ifndef _SELECT_PATTERN_H_
#define _SELECT_PATTERN_H_

#include "config/config.h"
#include "config/type.h"
#include "graph/pattern.h"
#include "util/net_simile.hpp"
#include <bitset>


class SelectPattern
{
public:
    SelectPattern() = default;
    ~SelectPattern() = default;
    SelectPattern(const SelectPattern &other) = delete;
    SelectPattern(SelectPattern &&other) = delete;

public:
    static void allocateBuffers(Pattern *&patterns);
    static void releaseBuffers(Pattern *&patterns);

    // naive approach: f(P) = 1/(3|P|) (f_cov(P) - f_sim(P) - f_cog(P) + |P| + |R(P)|)
    // pattern frequency is estimated by its super node's frequency
    // no isomorphism check between patterns
    static void calculateCoverage(const SuperGraph *superGraph, std::vector<Pattern> *DRP);
    static void calculateCognition(const SuperGraph *superGraph, std::vector<Pattern> *DRP);
    // selected: a subset of super node id selected in UpdatePattern
    // if "selected" is nonempty, we will use the DRPs of hitted nodes as candidates
    static void greedySelect(std::vector<Pattern> *DRP, Pattern *patterns, ui numVertex, 
                             ui numAttribute, const std::vector<ui> &selected = std::vector<ui>(), 
                             ui requiredSize = 0);

    // report pattern set's scores
    static void reportPatternSetScores(const Pattern *patterns, std::string &filename);

private:
    // get pattern pos in DRP from id
    static std::pair<ui, ui> getPatternPos(ui id);
    // https://en.wikipedia.org/wiki/Planar_graph#Other_planarity_criteria
    // if we add the dropped edges back, may need to rethink the criteria
    // also do we need to consider all cases of 3-cycle?
    static bool isPlanar(ui numV, ui nodeV, ui numE, SuperNodeType sntype);
    // count # of 1 bit
    static int countBits(unsigned long long n);
};


#endif // SELECT_PATTERN_H