#include "pets/generate_pattern.h"


void GeneratePattern::allocateBuffers(const SuperGraph *superGraph, std::vector<Pattern> *&patterns)
{
    patterns = new std::vector<Pattern>[superGraph->numVertex];
}

void GeneratePattern::releaseBuffers(std::vector<Pattern> *&patterns)
{
    delete[] patterns;
    patterns = nullptr;
}


void GeneratePattern::generatePatterns(const SuperGraph *superGraph, std::vector<Pattern> *patterns)
{
    ui numVertex = superGraph->numVertex;

    std::vector<std::vector<ui>> tmpNeighbors_;
    std::vector<ui> tmpNeigh_;
    tmpNeighbors_.reserve(MAX_PATTERN_SIZE);

    for(ui sid = 0; sid < numVertex; sid++) {
        // skip singleton
        if(superGraph->synopses[sid].type == SuperNodeType::SINGLETON)
            break;

        // std::cout << sid << std::endl;
        ui nodeV = superGraph->synopses[sid].cog.first;
        ui nodeE = superGraph->synopses[sid].cog.second;

        if(nodeE >= MIN_PATTERN_SIZE) {
            patterns[sid].emplace_back(nodeV, sid);
            patterns[sid][patterns[sid].size() - 1].updateInfo(superGraph);
        }

        SubGraph g;
        g.nodeId = sid;
        g.numEdge = nodeE;

        // decontract
        ui count_ = superGraph->getNodeDegree(sid);
        for(ui j = 0; j < count_; j++) {
            ui eid = superGraph->getSuperEdgeId(sid, j);
            int edgeNum = (int)superGraph->fd[eid].numEdge;
            int maxEnum = std::min(edgeNum, MAX_PATTERN_SIZE - (int)nodeE);
            
            const auto &originalEdges = superGraph->fd[sid].edges;

            for(int enumNum = 1; enumNum <= maxEnum; enumNum++) {
                // std::cout << enumNum << " " << maxEnum << std::endl;
                std::vector<int> comb(enumNum);
                std::iota(comb.begin(), comb.end(), 0);

                while(comb[0] != edgeNum - enumNum + 1) {
                    // original subgraph
                    g.reset4SameSuperNode();
                    tmpNeighbors_.clear();

                    for(const auto &c : comb) {
                        ui outv = originalEdges[c].second;
                        if(originalEdges[c].first >= superGraph->numOriginalVertex ||
                           originalEdges[c].second >= superGraph->numOriginalVertex) {
                            std::cerr << "error in: " << c << " " << sid << " " << eid << std::endl;
                            exit(1);
                        }
                        if(!g.ifVertexExist(outv))
                            g.outvertices.emplace_back(outv);

                        g.outedges.emplace_back(originalEdges[c]);
                        ++ g.numEdge;
                    }

                    // recursively grow
                    bool newV = false;
                    std::vector<ui> depth(g.outvertices.size(), 1);
                    std::vector<ui> pivots(g.outvertices.size(), 0);
                    for(const auto &v : g.outvertices) {
                        // tmpNeigh_.clear();
                        superGraph->getVertexAllNeighbors(v, tmpNeigh_);
                        tmpNeighbors_.emplace_back(std::move(tmpNeigh_));
                    }
                    growSubgraph(superGraph, patterns[sid], g, depth, newV, tmpNeighbors_, pivots);

                    // enumerate next
                    int pivot = enumNum - 1;
                    while(pivot >= 0 && comb[pivot] == edgeNum - enumNum + pivot)
                        -- pivot;
                    if(pivot < 0)
                        break;
                    int temp = comb[pivot];
                    for(int st = pivot; st < enumNum; st++)
                        comb[st] = temp + 1 + st - pivot;
                }
            }
        }

        // maintain top K
        selectTopK(patterns[sid]);
    }

    // rearrange id
    ui newpid = 0;
    for(ui sid = 0; sid < numVertex; sid++) {
        if(superGraph->synopses[sid].type == SuperNodeType::SINGLETON)
            break;
        for(auto &p : patterns[sid])
            p.id = newpid ++; 
    }
}


ui GeneratePattern::getRandomIdx(ui lower, ui upper)
{
    std::random_device dev;
    std::mt19937 rng(dev());
    // distribution in range [lower, upper]
    std::uniform_int_distribution<std::mt19937::result_type> dist(lower,upper); 

    return dist(rng);
}


void GeneratePattern::growSubgraph(const SuperGraph *superGraph, std::vector<Pattern> &patterns,
                                   SubGraph &g, std::vector<ui> &depth, bool newV, 
                                   std::vector<std::vector<ui>> &neighbors_, 
                                   std::vector<ui> &pivots)
{
    if(g.numEdge > MAX_PATTERN_SIZE)
        return;
    else if(g.numEdge >= MIN_PATTERN_SIZE) {
        ui sid = g.nodeId;
        ui nodeV = superGraph->fcprime[sid].numVertex + g.outvertices.size();
        patterns.emplace_back(nodeV, sid);
        patterns[patterns.size() - 1].updateInfo(superGraph, g);

        if(g.numEdge == MAX_PATTERN_SIZE)
            return;
    }

    // grow
    ui watchPointSize = g.outvertices.size();
    ui idx = 0;
    do {
        idx = getRandomIdx(0, watchPointSize - 1);
    } while(depth[idx] > MAX_GROW_DEPTH);

    ui v = g.outvertices[idx];
    ui curPivot = pivots[idx] + 1;
    const auto &curNeighbors_ = neighbors_[idx];
    ui curSize = curNeighbors_.size();

    for(ui l = curPivot; l < curSize; l++) {
        ui to = curNeighbors_[l];
        auto e = v > to ? std::make_pair(to, v) : std::make_pair(v, to);
        if(superGraph->vertexSuperNodeId[to] == g.nodeId || g.ifEdgeExist(e))
            continue;

        newV = false;
        
        // ++ depth[idx];
        if(!g.ifVertexExist(to)) {
            g.outvertices.emplace_back(to);
            depth.emplace_back(depth[idx] + 1);
            newV = true;
            std::vector<ui> tmpNeigh_;
            superGraph->getVertexAllNeighbors(to, tmpNeigh_);
            neighbors_.emplace_back(std::move(tmpNeigh_));
            pivots.emplace_back(0);
        }
        g.outedges.emplace_back(e);
        ++ g.numEdge;
        pivots[idx] = l;

        growSubgraph(superGraph, patterns, g, depth, newV, neighbors_, pivots);

        g.outedges.pop_back();
        -- g.numEdge;
        if(newV) {
            g.outvertices.pop_back();
            neighbors_.pop_back();
            depth.pop_back();
            pivots.pop_back();
        }
        // -- depth[idx];
    }
}


void GeneratePattern::growSubgraphAll(const SuperGraph *superGraph, std::vector<Pattern> &patterns,
                                      SubGraph &g, std::vector<ui> &depth, bool newV, 
                                      std::vector<std::vector<ui>> &neighbors_, 
                                      std::vector<ui> &pivots)
{
    if(g.numEdge > MAX_PATTERN_SIZE)
        return;
    else if(g.numEdge >= MIN_PATTERN_SIZE) {
        ui sid = g.nodeId;
        ui nodeV = superGraph->fcprime[sid].numVertex + g.outvertices.size();
        patterns.emplace_back(nodeV, sid);
        patterns[patterns.size() - 1].updateInfo(superGraph, g);

        if(g.numEdge == MAX_PATTERN_SIZE)
            return;
    }

    // grow
    ui watchPointSize = g.outvertices.size();
    for(ui idx = 0; idx < watchPointSize; idx++) {
        ui v = g.outvertices[idx];
        ui curPivot = pivots[idx] + 1;
        const auto &curNeighbors_ = neighbors_[idx];
        ui curSize = curNeighbors_.size();

        for(ui l = curPivot; l < curSize; l++) {
            ui to = curNeighbors_[l];
            auto e = v > to ? std::make_pair(to, v) : std::make_pair(v, to);
            if(superGraph->vertexSuperNodeId[to] == g.nodeId || g.ifEdgeExist(e))
                continue;

            newV = false;
            
            // ++ depth[idx];
            if(!g.ifVertexExist(to)) {
                g.outvertices.emplace_back(to);
                depth.emplace_back(depth[idx] + 1);
                newV = true;
                std::vector<ui> tmpNeigh_;
                superGraph->getVertexAllNeighbors(to, tmpNeigh_);
                neighbors_.emplace_back(std::move(tmpNeigh_));
                pivots.emplace_back(0);
            }
            g.outedges.emplace_back(e);
            ++ g.numEdge;
            pivots[idx] = l;

            growSubgraph(superGraph, patterns, g, depth, newV, neighbors_, pivots);

            g.outedges.pop_back();
            -- g.numEdge;
            if(newV) {
                g.outvertices.pop_back();
                neighbors_.pop_back();
                depth.pop_back();
                pivots.pop_back();
            }
            // -- depth[idx];
        }
    }
}


size_t GeneratePattern::hashEdge(const std::pair<ui, ui> &edge)
{
    std::hash<ui> hasher;
    size_t seed = 0;
    // 0x9e3779b9 = 2654435769
    seed = hasher(edge.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed = hasher(edge.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}


void GeneratePattern::selectTopK(std::vector<Pattern> &patterns)
{
    ui pivot;
    ui psize = patterns.size();
    if(psize <= DIVERSIFIED_K)
        return;

    // only store out edges
    using Coverage = std::vector<size_t>;
    Coverage tempDRP, maxCov;
    std::unordered_set<size_t> DRP;

    for(pivot = 0; pivot < DIVERSIFIED_K; pivot++) {
        ui selidx = 0;
        int curcov = -1;

        for(ui pid = pivot; pid < psize; pid++) {
            // intersect
            const auto &curp = patterns[pid];
            tempDRP.clear();
            
            int coverage = 0;
            for(const auto &e : curp.outedges) {
                size_t hv = hashEdge(e);
                if(DRP.count(hv) == 0) {
                    tempDRP.emplace_back(hv);
                    ++ coverage;
                }
            }
            
            // update
            if(coverage > curcov) {
                selidx = pid;
                curcov = coverage;
                maxCov = tempDRP;
            }
        }

        // append DRP
        std::iter_swap(patterns.begin() + pivot, patterns.begin() + selidx);
        for(const auto &e : maxCov)
            DRP.emplace(e);
    }

    patterns.resize(DIVERSIFIED_K);
}