#include "pets/update_pattern.h"

UpdatePattern::BoostGraph UpdatePattern::reformatPattern2Boost(const Pattern &pattern)
{
    UpdatePattern::BoostGraph g(pattern.getNumVertex());
    for(const auto &edge : pattern.edges)
        boost::add_edge(edge.first, edge.second, g);
    return g;
}


UpdatePattern::BoostGraph UpdatePattern::reformatSuperNode2Boost(SuperNodeType nodeType, ui numVertex)
{
    UpdatePattern::BoostGraph g(numVertex);

    switch(nodeType) {
        case SuperNodeType::CHORD:
            boost::add_edge(0, 1, g);
            for(ui i = 2; i < numVertex; i++) {
                boost::add_edge(0, i, g);
                boost::add_edge(1, i, g);
            }
            break;
        case SuperNodeType::CYCLE:
            for(ui i = 0; i < numVertex; i++) {
                ui next = (i + 1) % numVertex;
                boost::add_edge(i, next, g);
            }
            break;
        case SuperNodeType::PATH:
            for(ui i = 0; i < numVertex - 1; i++)
                boost::add_edge(i, i + 1, g);
            break;
        case SuperNodeType::STAR: 
            for(ui i = 1; i < numVertex; i++)
                boost::add_edge(0, i, g);
            break;
    }

    return g;
}


bool UpdatePattern::isSubgraph(const Pattern &lhs, const Pattern &rhs, 
                               const UpdatePattern::BoostGraph *boostLhs, 
                               const UpdatePattern::BoostGraph *boostRhs)
{
    if(lhs.getNumVertex() > rhs.getNumVertex() || lhs.getNumEdge() > rhs.getNumEdge())
        return false;
    
    UpdatePattern::BoostGraph gLhs = boostLhs == nullptr ? reformatPattern2Boost(lhs) : *boostLhs;
    UpdatePattern::BoostGraph gRhs = boostRhs == nullptr ? reformatPattern2Boost(rhs) : *boostRhs;

    vf2_bool_callback<UpdatePattern::BoostGraph> callback(gLhs, gRhs);
    return boost::vf2_subgraph_mono(gLhs, gRhs, callback);
}


bool UpdatePattern::isSubgraph(const UpdatePattern::BoostGraph &lhs, const UpdatePattern::BoostGraph &rhs)
{
    if(boost::num_vertices(lhs) > boost::num_vertices(rhs) || boost::num_edges(lhs) > boost::num_edges(rhs))
        return false;

    vf2_bool_callback<UpdatePattern::BoostGraph> callback(lhs, rhs);
    return boost::vf2_subgraph_mono(lhs, rhs, callback);
}


void UpdatePattern::collectHitedSuperNode(const SuperGraph *superGraph, const Pattern &partialQuery, 
                                          std::vector<ui> nodeId[4], const PMIndex &index)
{
    UpdatePattern::BoostGraph gQuery = reformatPattern2Boost(partialQuery);

    std::vector<ui> attributeHitted;
    for(const auto &attr : partialQuery.distinctAttributes)
        for(const auto &nid : index.invertedAttributeList[attr])
            attributeHitted.emplace_back(nid);
    
    std::vector<ui> structureHitted;
    for(ui i = 0; i < MAX_PATTERN_SIZE + 1; i++) {
        if(!index.invertedLengthListChord[i].empty()) {
            const auto &gNode = index.nodeBoostGraphChord[i];
            bool nodeIsSubgraphQuery = isSubgraph(gNode, gQuery);
            bool queryIsSubgraphNode = isSubgraph(gQuery, gNode);
            if(nodeIsSubgraphQuery || queryIsSubgraphNode)
                for(const auto &nid : index.invertedLengthListChord[i])
                    structureHitted.emplace_back(nid);
        }
    }

    __gnu_parallel::sort(attributeHitted.begin(), attributeHitted.end());
    __gnu_parallel::sort(structureHitted.begin(), structureHitted.end());
    
    std::vector<ui> intersectionHitted;
    __gnu_parallel::set_intersection(attributeHitted.begin(), attributeHitted.end(), 
                                     structureHitted.begin(), structureHitted.end(), 
                                     std::back_inserter(intersectionHitted));
    
    std::vector<ui> unionHitted;
    __gnu_parallel::set_union(attributeHitted.begin(), attributeHitted.end(), 
                              structureHitted.begin(), structureHitted.end(), 
                              std::back_inserter(unionHitted));

    std::cout << "###     attribute hitted: " << attributeHitted.size() << std::endl;
    std::cout << "###     structure hitted: " << structureHitted.size() << std::endl;
    std::cout << "###  intersection hitted: " << intersectionHitted.size() << std::endl;
    std::cout << "###         union hitted: " << unionHitted.size() << std::endl;


    nodeId[0] = attributeHitted;
    nodeId[1] = structureHitted;
    nodeId[2] = intersectionHitted;
    nodeId[3] = unionHitted;
}


bool UpdatePattern::updatePatternOnDisplay(const SuperGraph *superGraph, const std::vector<ui> nodeId[4],
                                           std::vector<Pattern> *DRP, const std::vector<ui> &nodeOnDisplay, 
                                           ui *pattern2Node, Pattern *displayedPattern, Pattern *&newPattern,
                                           UpdateStrategy strategy)
{
    std::vector<ui> hittedNode;
    switch(strategy) {
        case UpdateStrategy::ATTRIBUTE: hittedNode = nodeId[0]; break;
        case UpdateStrategy::STRUCTURE: hittedNode = nodeId[1]; break;
        case UpdateStrategy::INTERSECTION: hittedNode = nodeId[2]; break;
        case UpdateStrategy::UNION: hittedNode = nodeId[3]; break;
    }

    ui hitted = 0;
    for(ui i = 0; i < nodeOnDisplay.size(); i++)
        hitted = std::count(hittedNode.begin(), hittedNode.end(), nodeOnDisplay[i]) != 0 ? hitted + 1 : hitted;
    double hittedRate = hitted * 1.0 / nodeOnDisplay.size();

    if(hittedRate < HIT_THRESHOLD) {
        std::cout << "hit rate: " << hittedRate << " is greater than threshold: " << HIT_THRESHOLD << std::endl;
        newPattern = nullptr;
        return false;
    }

    std::vector<ui> keptPatternId;
    for(ui i = 0; i < PATTERN_SET_SIZE; i++)
        if(std::count(hittedNode.begin(), hittedNode.end(), pattern2Node[i]) != 0)
            keptPatternId.emplace_back(i);
    ui requiredSize = PATTERN_SET_SIZE - keptPatternId.size();

    newPattern = new Pattern[PATTERN_SET_SIZE];
    SelectPattern::greedySelect(DRP, newPattern, superGraph->numVertex, superGraph->numAttribute, 
                                hittedNode, requiredSize);
    for(ui i = requiredSize; i < PATTERN_SET_SIZE; i++)
        newPattern[i] = displayedPattern[keptPatternId[i - requiredSize]];

    return true;
}