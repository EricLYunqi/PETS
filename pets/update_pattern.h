#ifndef _UPDATE_PATTERN_H_
#define _UPDATE_PATTERN_H_

#include "config/config.h"
#include "config/type.h"
#include "graph/pattern.h"
#include "pets/select_pattern.h"
#include <parallel/algorithm>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

// do we need to consider attribute & structure at the same time in "isSubgraph"?

class PMIndex;

class UpdatePattern
{
public:
    using BoostGraph = boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>;

    template <typename Graph>
    struct vf2_bool_callback 
    {
        vf2_bool_callback(const Graph &graph1, const Graph &graph2)
            : graph1_(graph1), graph2_(graph2) {}

        template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
        bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
            return false; 
        }

    private:
        const Graph &graph1_;
        const Graph &graph2_;
    };

public:
    UpdatePattern() = default;
    ~UpdatePattern() = default;
    UpdatePattern(const UpdatePattern &other) = delete;
    UpdatePattern(UpdatePattern &&other) = delete;

public:
    static BoostGraph reformatPattern2Boost(const Pattern &pattern);
    static BoostGraph reformatSuperNode2Boost(SuperNodeType nodeType, ui numVertex);
    // check if lhs is a subgraph of rhs
    static bool isSubgraph(const Pattern &lhs, const Pattern &rhs, 
                           const BoostGraph *boostLhs = nullptr, 
                           const BoostGraph *boostRhs = nullptr);
    static bool isSubgraph(const BoostGraph &lhs, const BoostGraph &rhs);
    // attribute hit and structure hit seperately
    // nodeId[0]: attribute hitted
    // nodeId[1]: structure hitted
    // nodeId[2]: intersection
    // nodeId[3]: union
    static void collectHitedSuperNode(const SuperGraph *superGraph, const Pattern &partialQuery, 
                                      std::vector<ui> nodeId[4], const PMIndex &index);
    static bool updatePatternOnDisplay(const SuperGraph *superGraph, const std::vector<ui> nodeId[4],
                                       std::vector<Pattern> *DRP, const std::vector<ui> &nodeOnDisplay, 
                                       ui *pattern2Node, Pattern *displayedPattern, Pattern *&newPattern,
                                       UpdateStrategy strategy);
};

class PMIndex
{
public:
    // attr id -> [node id1, node id2]
    std::vector<std::vector<ui>> invertedAttributeList;
    // node length (numVertex, ith pos corresponds to length i) -> [node id1, node id2]
    std::vector<ui> *invertedLengthListChord{nullptr};
    std::vector<ui> *invertedLengthListCycle{nullptr};
    std::vector<ui> *invertedLengthListPath{nullptr};
    std::vector<ui> *invertedLengthListStar{nullptr};
    // Boost graph
    UpdatePattern::BoostGraph *nodeBoostGraphChord{nullptr};
    UpdatePattern::BoostGraph *nodeBoostGraphCycle{nullptr};
    UpdatePattern::BoostGraph *nodeBoostGraphPath{nullptr};
    UpdatePattern::BoostGraph *nodeBoostGraphStar{nullptr};

public:
    PMIndex() = default;
    PMIndex(const SuperGraph *superGraph, const InverseFunction *fcprime, 
            const Synopsis *synopses) {
        // init
        invertedAttributeList.resize(superGraph->numAttribute);
        ui reserveSize = superGraph->numVertex / superGraph->numAttribute;
        for(auto &bucket: invertedAttributeList)
            bucket.reserve(reserveSize);
        // give enough space no matter the super node type
        invertedLengthListChord = new std::vector<ui>[MAX_PATTERN_SIZE + 2];
        invertedLengthListCycle = new std::vector<ui>[MAX_PATTERN_SIZE + 2];
        invertedLengthListPath = new std::vector<ui>[MAX_PATTERN_SIZE + 2];
        invertedLengthListStar = new std::vector<ui>[MAX_PATTERN_SIZE + 2];
        nodeBoostGraphChord = new UpdatePattern::BoostGraph[MAX_PATTERN_SIZE + 2];
        nodeBoostGraphCycle = new UpdatePattern::BoostGraph[MAX_PATTERN_SIZE + 2];
        nodeBoostGraphPath = new UpdatePattern::BoostGraph[MAX_PATTERN_SIZE + 2];
        nodeBoostGraphStar = new UpdatePattern::BoostGraph[MAX_PATTERN_SIZE + 2];

        // build
        ui numVertex = superGraph->numVertex;
        for(ui i = 0; i < numVertex; i++) {
            ui numv = fcprime[i].numVertex;
            auto nodeType = synopses[i].type;

            for(const auto &a : fcprime[i].distinctAttributes)
                invertedAttributeList[a].emplace_back(i);
            switch(nodeType) {
                case SuperNodeType::CHORD:
                    invertedLengthListChord[numv].emplace_back(i);
                    nodeBoostGraphChord[numv] = UpdatePattern::reformatSuperNode2Boost(nodeType, numv);
                    break;
                case SuperNodeType::CYCLE:
                    invertedLengthListCycle[numv].emplace_back(i);
                    nodeBoostGraphCycle[numv] = UpdatePattern::reformatSuperNode2Boost(nodeType, numv);
                    break;
                case SuperNodeType::PATH:
                    invertedLengthListPath[numv].emplace_back(i);
                    nodeBoostGraphPath[numv] = UpdatePattern::reformatSuperNode2Boost(nodeType, numv);
                    break;
                case SuperNodeType::STAR:
                    invertedLengthListStar[numv].emplace_back(i);
                    nodeBoostGraphStar[numv] = UpdatePattern::reformatSuperNode2Boost(nodeType, numv);
                    break;
            }
        }
    }
    ~PMIndex() {
        delete[] invertedLengthListChord;
        delete[] invertedLengthListCycle;
        delete[] invertedLengthListPath;
        delete[] invertedLengthListStar;
        delete[] nodeBoostGraphChord;
        delete[] nodeBoostGraphCycle;
        delete[] nodeBoostGraphPath;
        delete[] nodeBoostGraphStar;
    }
    PMIndex(const PMIndex &other) = delete;
    PMIndex(PMIndex &&other) = delete;
};

#endif // _UPDATE_PATTERN_H_