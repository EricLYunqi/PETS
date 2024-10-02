#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "config/config.h"
#include "config/type.h"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <iostream>
#include <unordered_map>


// We do not allocate memory during construction
class Base
{
public:
    ui numVertex{0};
    ui numEdge{0};
    ui numAttribute{0};

    // CSR
    ui *offsets{nullptr};
    ui *neighbors{nullptr};

public:
    Base() = default;
    ~Base() {
        delete[] offsets;
        delete[] neighbors;
    }
    Base(const Base &other) = delete;
    Base(Base &&other) = delete;
};


class Graph : public Base
{
public:
    ui maxDegree{0};

    std::pair<ui, ui> *edges{nullptr};
    ui *attrOffsets{nullptr};
    ui *attributes{nullptr};
    double *attrWeights{nullptr};

public:
    Graph() = default;
    ~Graph() {
        delete[] edges;
        delete[] attrOffsets;
        delete[] attributes;
        delete[] attrWeights;
    }
    Graph(const Graph &other) = delete;
    Graph(Graph &&other) = delete;

public:
    ui getVertexDegree(const VertexID id) const {
        return offsets[id + 1] - offsets[id];
    }

    ui getVertexNumAttr(const VertexID id) const {
        return attrOffsets[id + 1] - attrOffsets[id];
    }

    const ui* getVertexNeighbors(const VertexID id, ui &count) const {
        count = offsets[id + 1] - offsets[id];
        return neighbors + offsets[id];
    }

    const ui* getVertexAttributes(const VertexID id, ui &count) const {
        count = attrOffsets[id + 1] - attrOffsets[id];
        return attributes + attrOffsets[id];
    }

    bool checkEdgeExistence(VertexID u, VertexID v) const {
        if (getVertexDegree(u) < getVertexDegree(v)) 
            std::swap(u, v);

        ui count = 0;
        const VertexID* neighbors_ = getVertexNeighbors(v, count);

        auto niter = std::lower_bound(neighbors_, neighbors_ + count, u);

        if(niter == neighbors_ + count)
            return false;
        else if(*niter == u)
            return true;
            
        return false;
    }

public:
    // vertex id may not be continous
    void loadGraphFromFiles(const std::string &statFile, const std::string &edgeFile, 
                            const std::string &attrFile);
    bool selfCheck() const;
    void printMetaData(const std::string &dataname) const;
};


struct InverseFunction
{
    ui id{0}; // super node id
    ui numVertex{0};
    std::vector<ui> vertices;
    std::vector<std::vector<ui>> attributes;
    std::unordered_set<ui> distinctAttributes;

    InverseFunction() = default;
    InverseFunction(ui _id, ui _numVertex)
    : id(_id), numVertex(_numVertex), vertices(_numVertex), attributes(_numVertex) { }
    InverseFunction(ui _id, const Graph *graph, const std::vector<ui> &_vertices)
    : id(_id), numVertex(_vertices.size()), vertices(_vertices), attributes(_vertices.size()) {
        for(ui i = 0; i < numVertex; i++) {
            ui v = vertices[i];
            ui count_ = 0;
            const ui *attributes_ = graph->getVertexAttributes(v, count_);
            for(ui j = 0; j < count_; j++) {
                attributes[i].emplace_back(attributes_[j]);
                distinctAttributes.insert(attributes_[j]);
            }
        }
    }

    void updateInfo(ui _id, const Graph *graph, const std::vector<ui> &_vertices) {
        id = _id;
        numVertex = _vertices.size();
        vertices = _vertices;
        attributes.resize(numVertex);

        for(ui i = 0; i < numVertex; i++) {
            ui v = vertices[i];
            ui count_ = 0;
            const ui *attributes_ = graph->getVertexAttributes(v, count_);
            for(ui j = 0; j < count_; j++)
                attributes[i].emplace_back(attributes_[j]);
        }
    }
};

struct DecontractionFunction
{
    ui id{0}; // super edge id
    ui numEdge{0};
    std::vector<ui> offsets; // for the s super node, degree stroed by pos
    std::vector<std::pair<ui, ui>> edges; // original vertex id

    DecontractionFunction() = default;
    DecontractionFunction(ui _id, ui _numEdge)
    : id(_id), numEdge(_numEdge), edges(_numEdge) { }
    DecontractionFunction(ui _id, const std::vector<ui> &_offsets, const std::vector<std::pair<ui, ui>> &_edges)
    : id(_id), numEdge(_edges.size()), offsets(_offsets), edges(_edges) { }

    void updateInfo(ui _id, const std::vector<ui> &_offsets, const std::vector<std::pair<ui, ui>> &_edges) {
        id = _id;
        numEdge = _edges.size();
        offsets = _offsets;
        edges = _edges;
    }
};


struct Synopsis
{
    SuperNodeType type;
    std::vector<ui> info;
    std::pair<ui, ui> cog;

    Synopsis() = default;
    Synopsis(SuperNodeType _type, const std::vector<ui> &_vertices)
    : type(_type) {
        switch(_type) {
            case SuperNodeType::CHORD : {
                info.emplace_back(_vertices[0]);
                info.emplace_back(_vertices[1]);
                cog = std::make_pair(_vertices.size(), 2 * (_vertices.size() - 2) + 1);
                break;
            }
            case SuperNodeType::CYCLE : {
                info = _vertices;
                cog = std::make_pair(_vertices.size(), _vertices.size());
                break;
            }
            case SuperNodeType::PATH : {
                info = _vertices;
                cog = std::make_pair(_vertices.size(), _vertices.size() - 1);
                break;
            }
            case SuperNodeType::STAR : {
                info.emplace_back(_vertices[0]);
                cog = std::make_pair(_vertices.size(), _vertices.size() - 1);
                break;
            }
            case SuperNodeType::SINGLETON : {
                info.emplace_back(_vertices[0]);
                cog = std::make_pair(1, 0);
                break;
            }
        }
    }

    void updateInfo(SuperNodeType _type, const std::vector<ui> &_vertices) {
        type = _type;
        switch(_type) {
            case SuperNodeType::CHORD : {
                info.emplace_back(_vertices[0]);
                info.emplace_back(_vertices[1]);
                cog = std::make_pair(_vertices.size(), 2 * (_vertices.size() - 2) + 1);
                break;
            }
            case SuperNodeType::CYCLE : {
                info = _vertices;
                cog = std::make_pair(_vertices.size(), _vertices.size());
                break;
            }
            case SuperNodeType::PATH : {
                info = _vertices;
                cog = std::make_pair(_vertices.size(), _vertices.size() - 1);
                break;
            }
            case SuperNodeType::STAR : {
                info.emplace_back(_vertices[0]);
                cog = std::make_pair(_vertices.size(), _vertices.size() - 1);
                break;
            }
            case SuperNodeType::SINGLETON : {
                info.emplace_back(_vertices[0]);
                cog = std::make_pair(1, 0);
                break;
            }
        }
    }
};


/*
 * synopses: S, synopsis for each super node
 * fcprime: f'_C, inverse contraction function for each super node
 * fd: f_D, decontraction for each super edge
 * Use CSR format to store the super edge, "edges" will also store edges
 * Not only CSR & edges but f_D needs twice of numEdge space
 * A super edge (undirected, s - t) will have 2 ids depending on the s super node
 * (super edge id = offsets[node] + pos of t super node)
 * Two ids correspond to the same decontraction information, hence we dont need extra spaces for super edge id
 * In f_D, edges stored like CSR format but directly stored as vector<pair<ss, tt>>, the vector is sorted as pos[ss] in the super node
 */
class SuperGraph : public Base
{
public:
    ui numChord{0};
    ui numCycle{0};
    ui numPath{0};
    ui numStar{0};
    ui numSingleton{0};

    // node size -> # of nodes of that size
    std::unordered_map<ui, ui> nodeLengthMap[5];
    
    ui numOriginalVertex{0};
    ui numOriginalEdge{0};

    std::pair<ui, ui> *edges{nullptr};
    ui *vertexSuperNodeId{nullptr};
    ui *vertexSuperNodePos{nullptr};

    Synopsis *synopses{nullptr};
    InverseFunction *fcprime{nullptr};
    DecontractionFunction *fd{nullptr}; // 2 * numEdge

public:
    SuperGraph() = default;
    ~SuperGraph() {
        delete[] edges;
        delete[] vertexSuperNodeId;
        delete[] vertexSuperNodePos;
        delete[] synopses;
        delete[] fcprime;
        delete[] fd;
    }
    SuperGraph(const SuperGraph &other) = delete;
    SuperGraph(SuperGraph &&other) = delete;

public:
    // super node id: [minId, maxId)
    std::pair<ui, ui> getChordNodes() const {
        ui minId = 0;
        ui maxId = numChord;
        return std::make_pair(minId, maxId);
    }

    std::pair<ui, ui> getCycleNodes() const {
        ui minId = numChord;
        ui maxId = numChord + numCycle;
        return std::make_pair(minId, maxId);
    }

    std::pair<ui, ui> getPathNodes() const {
        ui minId = numChord + numCycle;
        ui maxId = numChord + numCycle + numPath;
        return std::make_pair(minId, maxId);
    }

    std::pair<ui, ui> getStarNodes() const {
        ui minId = numChord + numCycle + numPath;
        ui maxId = numChord + numCycle + numPath + numStar;
        return std::make_pair(minId, maxId);
    }

    std::pair<ui, ui> getSingletonNodes() const {
        ui minId = numChord + numCycle + numPath + numStar;
        ui maxId = numChord + numCycle + numPath + numStar + numSingleton;
        return std::make_pair(minId, maxId);
    }

    // functions name with "Node": for super node
    // functions name with "Vertex": for original vertex
    ui getNodeDegree(const VertexID id) const {
        return offsets[id + 1] - offsets[id];
    }

    // note that for undirected graph, u -> v & v -> u share the same infos
    // thus two different ids map to the same f_D index
    ui getSuperEdgeId(const VertexID id, const ui pos) const {
        return offsets[id] + pos;
    }

    const ui* getNodeNeighbors(const VertexID id, ui &count) const {
        count = offsets[id + 1] - offsets[id];
        return neighbors + offsets[id];
    }

    // get vertex's neighbors within one super edge
    void getVertexOutNeighbors(ui vid, ui pos, std::vector<ui> &neighbors_) const {
        ui nodeid = vertexSuperNodeId[vid];
        ui vpos = vertexSuperNodePos[vid];
        ui eid = getSuperEdgeId(nodeid, pos);
        ui stidx = fd[eid].offsets[vpos];
        ui enidx = fd[eid].offsets[vpos];
        neighbors_.reserve(enidx - stidx);
        for(ui i = stidx; i < enidx; i++)
            neighbors_.emplace_back(fd[eid].edges[i].second);
    }   

    // get vertex's neighbors within all super edges
    void getVertexOutNeighbors(ui vid, std::vector<ui> &neighbors_) const {
        ui nodeid = vertexSuperNodeId[vid];
        ui vpos = vertexSuperNodePos[vid];
        ui count_ = getNodeDegree(nodeid);

        for(ui i = 0; i < count_; i++) {
            ui eid = getSuperEdgeId(nodeid, i);
            ui stidx = fd[eid].offsets[vpos];
            ui enidx = fd[eid].offsets[vpos];
            for(ui j = stidx; j < enidx; j++)
                neighbors_.emplace_back(fd[eid].edges[j].second);
        }
    }

    // get vertex's all neighbors
    void getVertexAllNeighbors(ui vid, std::vector<ui> &neighbors_) const {
        ui nodeid = vertexSuperNodeId[vid];
        ui vpos = vertexSuperNodePos[vid];
        auto sntype = synopses[nodeid].type;
        ui nodenum = fcprime[nodeid].numVertex;
        const auto &snvertices = fcprime[nodeid].vertices;

        switch(sntype) {
            case SuperNodeType::CHORD : {
                if(vpos == 0) {
                    for(ui i = 1; i < nodenum; i++)
                        neighbors_.emplace_back(snvertices[i]);
                }
                else if(vpos == 1){
                    neighbors_.emplace_back(snvertices[0]);
                    for(ui i = 2; i < nodenum; i++)
                        neighbors_.emplace_back(snvertices[i]);
                }
                else {
                    neighbors_.emplace_back(snvertices[0]);
                    neighbors_.emplace_back(snvertices[1]);
                }
                break;
            }
            case SuperNodeType::CYCLE : {
                ui prev = vpos == 0 ? nodenum - 1 : vpos - 1;
                ui next = vpos == nodenum - 1 ? 0 : vpos + 1;
                neighbors_.emplace_back(snvertices[prev]);
                neighbors_.emplace_back(snvertices[next]);
                break;
            }
            case SuperNodeType::PATH : {
                if(vpos == 0) 
                    neighbors_.emplace_back(snvertices[1]);
                else if(vpos == nodenum - 1)
                    neighbors_.emplace_back(snvertices[nodenum - 2]);
                else {
                    neighbors_.emplace_back(snvertices[vpos - 1]);
                    neighbors_.emplace_back(snvertices[vpos + 1]);
                }
                break;
            }
            case SuperNodeType::STAR : {
                if(vpos != 0)
                    neighbors_.emplace_back(snvertices[0]);
                else
                    for(ui i = 1; i < nodenum; i++)
                        neighbors_.emplace_back(snvertices[i]);
                break;
            }
            case SuperNodeType::SINGLETON :
                break;
        }

        getVertexOutNeighbors(vid, neighbors_);
    }


    bool checkSuperEdgeExistence(VertexID u, VertexID v) const {
        if (getNodeDegree(u) < getNodeDegree(v)) 
            std::swap(u, v);

        ui count = 0;
        const VertexID* neighbors_ = getNodeNeighbors(v, count);
        auto niter = std::lower_bound(neighbors_, neighbors_ + count, u);

        if(niter == neighbors_ + count)
            return false;
        else if(*niter == u)
            return true;

        return false;
    }


    // general api
    void printGraphMetaData() const;
    void saveGraph(const std::string &filename) const; // super graph's structure
    void saveIndex(const std::string &filename) const; // index: synopses, f'_C, f_D
};

#endif // _GRAPH_H_