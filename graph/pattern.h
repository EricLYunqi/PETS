#ifndef _PATTERN_H_
#define _PATTERN_H_

#include "config/config.h"
#include "config/type.h"
#include "graph/graph.h"
#include <vector>
#include <unordered_set>

#define NOTEXIST 0xffffffff
using PatEdge = std::pair<ui, ui>;

// The subgraph grown from a super node
class SubGraph
{
public:
    ui nodeId{0};
    ui numEdge{0};
    std::vector<ui> outvertices; // (s -> t), we only store t
    // 1. (s -> t) where s \in v_{H}
    // 2. (s -> t) where s < t && s, t \notin v_{H}
    std::vector<std::pair<ui, ui>> outedges;
    bool hasTriangle{false};

    SubGraph() = default;
    SubGraph(ui _nid, ui _numE)
    : nodeId(_nid), numEdge(_numE) { }
    ~SubGraph() = default;

    bool ifVertexExist(ui v) const {
        return std::count(outvertices.begin(), outvertices.end(), v) == 0 
        ? false : true;
    }

    bool ifEdgeExist(std::pair<ui, ui> e) const {
        return std::count(outedges.begin(), outedges.end(), e) == 0
        ? false : true;
    }

    // avoid re-build subgraph for the same super node
    void reset4SameSuperNode() {
        numEdge = 0;
        outvertices.clear();
        outedges.clear();
        hasTriangle = false;
    }
};

class Pattern
{
public:
    static ui curId;

    ui id{0};     // pattern id
    ui nodeId{0}; // super node id
    ui numVertex{0};
    int finalRank{0};
    bool hasTriangle{false};

    std::vector<PatEdge> edges;
    std::vector<PatEdge> outedges; // used for top K selection
    std::vector<std::vector<PatEdge>> edgesList;  // s --> [(e, edge_id)]
    std::vector<short> degree;  // now only use for linear star search
    std::vector<std::vector<ui>> attributes;
    std::unordered_set<ui> distinctAttributes;

    double scoreCoverage{0.0};
    double scoreCognition{0.0};
    double scoreDiversity{0.0};
    double scoreAttribute{0.0};
    double finalScore{0.0};

public:
    Pattern() = default;
    Pattern(ui n, ui nid) {
        numVertex = n;
        edgesList.resize(n);
        attributes.resize(n);
        
        id = curId ++;
        nodeId = nid;
        scoreCoverage = scoreCognition = scoreDiversity = 0.0;
        finalScore = 0.0;
        
        degree.clear();
    }
    ~Pattern() = default;

public:
    void addAtrribute(ui id, ui count_, const ui *attributes_) {
        for(ui i = 0; i < count_; i++) {
            attributes[id].emplace_back(attributes_[i]);
            distinctAttributes.insert(attributes_[i]);
        }
    }

    void addVertex(ui count_, const std::vector<ui> &attributes_) {
        ++ numVertex;
        attributes.emplace_back();
        (*std::back_inserter(attributes)) = attributes_;
    }

    void addEdge(ui s, ui e) {
        edgesList[s].emplace_back(e, edges.size());
        edgesList[e].emplace_back(s, edges.size());
        if(s > e) std::swap(s, e);
        edges.emplace_back(s, e);
    }

    ui getEdgeId(ui s, ui e) const {
        for(auto p : edgesList[s])
            if (p.first == e) 
                return p.second;
        return NOTEXIST;
    }

    // These only returns out vertices & edges
    ui getNumVertex() const { 
        return numVertex; 
    }
    ui getNumEdge() const { 
        return edges.size(); 
    }

    // construct pattern from a super node
    void updateInfo(const SuperGraph *superGraph) {
        ui nodeV =  superGraph->fcprime[nodeId].numVertex;

        // vertices
        for(ui i = 0; i < nodeV; i++)
            attributes[i] = superGraph->fcprime[nodeId].attributes[i];

        // edges
        auto sntype = superGraph->synopses[nodeId].type;
        switch(sntype) {
            case SuperNodeType::CHORD : {
                addEdge(0, 1);
                for(ui i = 2; i < nodeV; i++) {
                    addEdge(0, i);
                    addEdge(1, i);
                }
                hasTriangle = true;
                break;
            }
            case SuperNodeType::CYCLE : {
                for(ui i = 0; i < nodeV - 1; i++)  
                    addEdge(i, i + 1);
                addEdge(nodeV - 1, 0);
                hasTriangle = nodeV == 3 ? true : false;
                break;
            }
            case SuperNodeType::PATH : {
                for(ui i = 0; i < nodeV - 1; i++)  
                    addEdge(i, i + 1);
                break;
            }
            case SuperNodeType::STAR : {
                for(ui i = 1; i < nodeV; i++)
                    addEdge(0, i);
                break;
            }
            case SuperNodeType::SINGLETON :
                break;
        }
    }

    // construct pattern from a subgraph
    void updateInfo(const SuperGraph *superGraph, const SubGraph &g) {
        // super node
        updateInfo(superGraph);

        // out part
        std::unordered_map<ui, ui> outvmap; // map out vertex to its pos in pattern

        ui nodeV = superGraph->fcprime[nodeId].numVertex;
        for(ui i = nodeV; i < numVertex; i++) {
            ui outv = g.outvertices[i - nodeV];
            ui outVNodeId = superGraph->vertexSuperNodeId[outv];
            ui outVNodePos = superGraph->vertexSuperNodePos[outv];
            attributes[i] = superGraph->fcprime[outVNodeId].attributes[outVNodePos];
            outvmap[outv] = i;
        }

        ui outEsize = g.outedges.size();
        for(ui i = 0; i < outEsize; i++) {
            ui st = g.outedges[i].first;
            ui end = g.outedges[i].second;
            ui stpos = superGraph->vertexSuperNodeId[st] == nodeId
            ? superGraph->vertexSuperNodePos[st] : outvmap[st];
            ui enpos = outvmap[end];
            addEdge(stpos, enpos);
            if(st > end) std::swap(st, end);
            outedges.emplace_back(st, end);
        }
    }
};


#endif // _PATTERN_H_