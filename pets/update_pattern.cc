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

    
    vf2_bool_callback<UpdatePattern::BoostGraph> callback(lhs, rhs);
    return boost::vf2_subgraph_mono(lhs, rhs, callback);
}


void UpdatePattern::collectHitedSuperNode(const SuperGraph *superGraph, const Pattern &partialQuery, 
                                          std::vector<ui> &nodeId, const PMIndex &index)
{
    UpdatePattern::BoostGraph gQuery = reformatPattern2Boost(partialQuery);

    std::vector<ui> attributeHited;
    for(const auto &attr : partialQuery.distinctAttributes)
        for(const auto &nid : index.invertedAttributeList[attr])
            attributeHited.emplace_back(nid);
    
    std::vector<ui> structureHited;
    for(ui i = 0; i < MAX_PATTERN_SIZE + 1; i++) {
        if(!index.invertedLengthListChord[i].empty()) {
            const auto &gNode = index.nodeBoostGraphChord[i];
            
        }
    }
}