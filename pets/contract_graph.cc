#include "pets/contract_graph.h"


std::unordered_map<SuperNodeType, std::string> ContractGraph::nodeMap = {
    { SuperNodeType::CHORD, "Chord" }, 
    { SuperNodeType::CYCLE, "Cycle" }, 
    { SuperNodeType::PATH, "Path" }, 
    { SuperNodeType::STAR, "Star" }, 
    { SuperNodeType::SINGLETON, "Singleton" }
};
ui ContractGraph::MAXCHORDSIZE = (MAX_PATTERN_SIZE - 1) / 2 + 2;
ui ContractGraph::MAXCYCLESIZE = MAX_PATTERN_SIZE;
ui ContractGraph::MAXPATHSIZE = MAX_PATTERN_SIZE + 1;
ui ContractGraph::MAXSTARSIZE = MAX_PATTERN_SIZE + 1;


void ContractGraph::extractChord(const Graph *graph, std::vector<std::vector<ui>> &chords, 
                                 bool *visited)
{
    ui maxDegree = graph->maxDegree;
    ui *common = new ui[maxDegree];
    std::fill(common, common + maxDegree, 0);
    ui commonCount = 0;

    ui numVertex = graph->numVertex;
    for(ui i = 0; i < numVertex; i++) {
        if(visited[i])
            continue;

        ui count_ = 0;
        const ui *neighbors_ = graph->getVertexNeighbors(i, count_);

        for(ui j = 0; j < count_; j++) {
            ui to = neighbors_[j];
            if(visited[to])
                continue;

            ui toCount = 0;
            const ui *toNeighbors = graph->getVertexNeighbors(to, toCount);

            compute_set_intersection::computeCandidates(neighbors_, count_, toNeighbors, toCount, 
                                                        common, commonCount);

            std::vector<ui> tempChords;
            tempChords.reserve(ContractGraph::MAXCHORDSIZE);
            
            tempChords.emplace_back(i);
            tempChords.emplace_back(to);
            for(ui l = 0; l < commonCount; l++) {
                if(!visited[common[l]]) {
                    tempChords.emplace_back(common[l]);
                    // visited[common[l]] = true;
                }
                if(tempChords.size() >= ContractGraph::MAXCHORDSIZE)
                    break;
            }

            if(tempChords.size() < MIN_SUPER_NODE_SIZE)
                continue;

            for(const auto &v : tempChords)
                visited[v] = true;
            chords.emplace_back(std::move(tempChords));
            break;
        }
    }

    delete[] common;
}


void ContractGraph::extractCycle(const Graph *graph, std::vector<std::vector<ui>> &cycles, 
                                 bool *visited)
{
    // fundmental set
    std::vector<std::vector<ui>> fcycles;
    generateFundmentalCycles(graph, fcycles, visited);
    std::cout << "Cardinality of fundmental set: " << fcycles.size() << std::endl;
    
    // enumerate
    for(auto cycle : fcycles) {
        if(cycle.size() < MIN_SUPER_NODE_SIZE || cycle.size() > ContractGraph::MAXCYCLESIZE)
            continue;

        for(const auto &v : cycle)
            visited[v] = true;
        cycles.emplace_back(cycle);
    }
}


void ContractGraph::extractPath(const Graph *graph, std::vector<std::vector<ui>> &paths,
                                bool *visited)
{
    ui numVertex = graph->numVertex;
    std::vector<bool> vis;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine eng(seed);

    for(ui i = 0; i < numVertex; i++) {
        if(visited[i])
            continue;

        vis.resize(numVertex, false);
        std::vector<ui> onepath;
        onepath.reserve(ContractGraph::MAXPATHSIZE);
        bool finished = false;
        generateOnePath(graph, onepath, visited, i, 0, vis, finished, eng);

        if(!finished || !verifyPath(graph, onepath))
            continue;

        for(const auto &v : onepath)
            visited[v] = true;
        paths.emplace_back(std::move(onepath));
    }
}


void ContractGraph::extractStar(const Graph *graph, std::vector<std::vector<ui>> &stars, 
                                bool *visited)
{
    ui numVertex = graph->numVertex;
    for(ui i = 0; i < numVertex; i++) {
        if(visited[i])
            continue;

        ui count;
        const ui *neighbors = graph->getVertexNeighbors(i, count);
        std::vector<std::pair<ui, ui>> candidates;
        for(ui j = 0; j < count; j++) {
            ui to = neighbors[j];
            ui degree = graph->getVertexDegree(to);
            if(!visited[to])
                candidates.emplace_back(to, degree);
        }

        std::sort(candidates.begin(), candidates.end(), [](const std::pair<ui, ui> &lhs, const std::pair<ui, ui> &rhs) {
            return lhs.second < rhs.second;
        });

        std::vector<ui> tempStar;
        tempStar.reserve(ContractGraph::MAXSTARSIZE);
        tempStar.emplace_back(i);
        for(const auto &p : candidates) {
            ui candidate = p.first;
            bool issuccess = true;
            for(ui j = 1; j < tempStar.size(); j++) {
                ui leaf = tempStar[j];
                if(graph->checkEdgeExistence(candidate, leaf)) {
                    issuccess = false;
                    break;
                }
            }
            if(issuccess)
                tempStar.emplace_back(candidate);
            if(tempStar.size() >= ContractGraph::MAXSTARSIZE)
                break;
        }

        if(tempStar.size() <= MIN_SUPER_NODE_SIZE)
            continue;
        for(const auto &v : tempStar)
            visited[v] = true;

        stars.emplace_back(std::move(tempStar));
    }
}


void ContractGraph::extractSingleton(const Graph *graph, std::vector<std::vector<ui>> &singletons, 
                                     bool *visited)
{
    ui numVertex = graph->numVertex;
    for(ui i = 0; i < numVertex; i++) {
        if(visited[i])
            continue;
        std::vector<ui> tempSingleton = {i};
        visited[i] = true;
        singletons.emplace_back(std::move(tempSingleton));
    }
}


void ContractGraph::generateFundmentalCycles(const Graph *graph, std::vector<std::vector<ui>> &fcycles, 
                                             bool *visited)
{
    ui numVertex = graph->numVertex;

    std::vector<bool> X(numVertex, false);
    std::vector<int> T;
    std::vector<int> pre(numVertex, -1);

    T.emplace_back(0);
    pre[0] = 0;
    
    while(true) {
        int z = 0;
        bool find = false;

        for(auto riter = T.rbegin(); riter != T.rend(); riter++) {
            if(!X[*riter] && !visited[*riter]) {
                z = *riter;
                find = true;
                break;
            }
        }

        if(!find)
            break;

        ui count_ = 0;
        const ui *neighbors_ = graph->getVertexNeighbors(z, count_);

        for(ui i = 0; i < count_; i++) {
            int w = (int)neighbors_[i];
            if(visited[w])
                continue;

            auto titer = std::find(T.begin(), T.end(), w);
            if(titer != T.end()) {
                std::vector<ui> tempFCycle;
                bool success = true;

                int v = w;
                do {
                    assert(v != -1);
                    success &= visited[v];
                    tempFCycle.emplace_back((ui)v);
                    v = pre[v];
                } while(v != pre[v]);

                if(success)
                    fcycles.emplace_back(std::move(tempFCycle));
            }
            else {
                T.emplace_back(w);
                pre[w] = z;
            }
        }

        X[z] = true;
    }
}


void ContractGraph::generateOnePath(const Graph *graph, std::vector<ui> &onepath, bool *visited, 
                                    ui root, ui curLen, std::vector<bool> &vis, bool &finished, 
                                    std::default_random_engine &eng)
{
    vis[root] = true;
    if(finished)
        return;
    if(curLen >= ContractGraph::MAXPATHSIZE) {
        finished = true;
        return;
    }

    ui count_ = 0;
    const ui *neighbors_ = graph->getVertexNeighbors(root, count_);

    std::vector<ui> shuffled;
    shuffled.reserve(count_);
    for(ui i = 0; i < count_; i++) {
        ui to = neighbors_[i];
        if(!visited[to] && !vis[to])
            shuffled.emplace_back(to);
    }

    // std::random_device rd;
    // std::mt19937 g(rd());
    std::shuffle(shuffled.begin(), shuffled.end(), eng);

    if(shuffled.empty() && curLen < MIN_SUPER_NODE_SIZE)
        return;
    else if(shuffled.empty() && curLen >= MIN_SUPER_NODE_SIZE) {
        finished = true;
        return;
    }
    else if(curLen >= MIN_SUPER_NODE_SIZE) {
        if(shuffled.front() % 2) {
            finished = true;
            return;
        }
    }

    for(const auto &to : shuffled) {
        onepath.emplace_back(to);
        generateOnePath(graph, onepath, visited, to, curLen + 1, vis, finished, eng);
        if(finished)
            break;
        onepath.pop_back();
    }

    return;
}


bool ContractGraph::verifyCycle(const Graph *graph, const std::vector<ui> &onecycle)
{

}


bool ContractGraph::verifyPath(const Graph *graph, const std::vector<ui> &onepath)
{
    ui size = onepath.size();
    if(size < MIN_SUPER_NODE_SIZE || size > ContractGraph::MAXPATHSIZE)
        return false;

    bool isconnected = true;
    for(ui i = 1; i < size - 1; i++) {
        bool prev = graph->checkEdgeExistence(onepath[i - 1], onepath[i]);
        bool next = graph->checkEdgeExistence(onepath[i + 1], onepath[i]);
        isconnected &= (prev & next);
        if(!isconnected)
            return false;
    }

    for(ui i = 0; i < size - 2; i++)
        for(ui j = i + 2; j < size; j++)
            if(graph->checkEdgeExistence(onepath[i], onepath[j]))
                return false;

    return true;
}


void ContractGraph::allocateBuffers(const Graph *graph, bool *&visited)
{
    visited = new bool[graph->numVertex];
    std::fill(visited, visited + graph->numVertex, false);
}


void ContractGraph::releaseBuffers(const Graph *graph, bool *&visited)
{
    delete[] visited;
    visited = nullptr;
}


void ContractGraph::updateSuperNode(SuperGraph *superGraph, const std::vector<std::vector<ui>> &nodesets,
                                    ui &superNodeId, ui cpv)
{
    for(const auto &nset : nodesets) {

        // check point
        auto nset_copy = nset;
        std::sort(nset_copy.begin(), nset_copy.end());
        auto uiter = std::unique(nset_copy.begin(), nset_copy.end());
        if(uiter != nset_copy.end()) {
            std::cerr << "Check point on super node containing duplicate vertices fails: " << superNodeId;
            for(auto iiter = uiter; iiter != nset_copy.end(); iiter++)
                std::cerr << *iiter << " ";
            std::cerr << std::endl;
        }

        ui vsize = nset.size();
        for(ui idx = 0; idx < vsize; idx++) {
            ui v = nset[idx];

            // check point
            if(superGraph->vertexSuperNodeId[v] != cpv || superGraph->vertexSuperNodePos[v] != cpv) {
                std::cerr << "Check point on different super nodes sharing the same vertex fails: " << v 
                << " " << superGraph->vertexSuperNodeId[v] << " " << superGraph->vertexSuperNodePos[v]
                << " " << superNodeId << " " << idx << std::endl;
                for(const auto &v : nset) std::cerr << v << " ";
                std::cerr << std::endl;
                for(const auto &v : nodesets[superGraph->vertexSuperNodeId[v]]) std::cerr << v << " ";
                std::cerr << std::endl;
                exit(1);
            }

            superGraph->vertexSuperNodeId[v] = superNodeId;
            superGraph->vertexSuperNodePos[v] = idx;
        }   
        ++ superNodeId;
    }
}


void ContractGraph::updateSuperNodeIndex(SuperGraph *superGraph, const std::vector<std::vector<ui>> &nodesets,
                                         const Graph *graph, SuperNodeType sntype)
{
    int idx = static_cast<int>(sntype);
    auto &lengthMap = superGraph->nodeLengthMap[idx];
    for(const auto &nset : nodesets) {
        ui sid = superGraph->vertexSuperNodeId[nset[0]];
        ui nsetsize = nset.size();
        superGraph->fcprime[sid].updateInfo(sid, graph, nset);
        superGraph->synopses[sid].updateInfo(sntype, nset);

        if(lengthMap.find(nsetsize) == lengthMap.end())
            lengthMap[nsetsize] = 1;
        else
            ++ lengthMap[nsetsize];
    }
}


void ContractGraph::generateContractionPlan(GraphType graphType, ContractionOrder &order)
{
    // assign according to graph's type
    switch(graphType) {
        case GraphType::SOCIAL : {
            order.emplace_back(&ContractGraph::extractCycle, SuperNodeType::CYCLE);
            order.emplace_back(&ContractGraph::extractStar, SuperNodeType::STAR);
            order.emplace_back(&ContractGraph::extractChord, SuperNodeType::CHORD);
            order.emplace_back(&ContractGraph::extractPath, SuperNodeType::PATH);
            order.emplace_back(&ContractGraph::extractSingleton, SuperNodeType::SINGLETON);
            break;
        }
        case GraphType::WEB : {
            order.emplace_back(&ContractGraph::extractCycle, SuperNodeType::CYCLE);
            order.emplace_back(&ContractGraph::extractChord, SuperNodeType::CHORD);
            order.emplace_back(&ContractGraph::extractStar, SuperNodeType::STAR);
            order.emplace_back(&ContractGraph::extractPath, SuperNodeType::PATH);
            order.emplace_back(&ContractGraph::extractSingleton, SuperNodeType::SINGLETON);
            break;
        }
        case GraphType::TRAFFIC : {
            order.emplace_back(&ContractGraph::extractCycle, SuperNodeType::CYCLE);
            order.emplace_back(&ContractGraph::extractPath, SuperNodeType::PATH);
            order.emplace_back(&ContractGraph::extractStar, SuperNodeType::STAR);
            order.emplace_back(&ContractGraph::extractChord, SuperNodeType::CHORD);
            order.emplace_back(&ContractGraph::extractSingleton, SuperNodeType::SINGLETON);
            break;
        }
    }

    // report
    printf("Contraction order: ");
    for(ui i = 0; i < 5; i++)
        printf("%s\t", nodeMap[order[i].second].c_str());
    printf("\n");
}


bool ContractGraph::selfCheck(const SuperGraph *superGraph)
{
    // node size
    for(ui i = 0; i < superGraph->numVertex; i++) {
        ui MAXNODESIZE = 0;
        auto nodetype = superGraph->synopses[i].type;
        switch(nodetype) {
            case SuperNodeType::CHORD : MAXNODESIZE = ContractGraph::MAXCHORDSIZE; break;
            case SuperNodeType::CYCLE : MAXNODESIZE = ContractGraph::MAXCYCLESIZE; break;
            case SuperNodeType::PATH : MAXNODESIZE = ContractGraph::MAXPATHSIZE; break;
            case SuperNodeType::STAR : MAXNODESIZE = ContractGraph::MAXSTARSIZE; break;
            case SuperNodeType::SINGLETON : MAXNODESIZE = 1; break;
        }
        if(superGraph->fcprime[i].numVertex > MAXNODESIZE || 
           superGraph->synopses[i].cog.second > MAX_PATTERN_SIZE) {
            std::cerr << "Super node size exceed max pattern size: "
            << superGraph->fcprime[i].numVertex << " "
            << superGraph->synopses[i].cog.second << std::endl;
            return false;
        }
        if(superGraph->fcprime[i].numVertex < MIN_SUPER_NODE_SIZE && 
           superGraph->synopses[i].type != SuperNodeType::SINGLETON) {
            std::cerr << "Super node size is smaller than min node size: "
            << superGraph->fcprime[i].numVertex << std::endl;
            return false;
        }
    }

    return true;
}


void ContractGraph::buildSuperGraph(const Graph *graph, GraphType graphType, SuperGraph *superGraph)
{
    // buffers
    bool *visited = nullptr;
    allocateBuffers(graph, visited);

    std::vector<std::vector<ui>> chords;
    std::vector<std::vector<ui>> cycles;
    std::vector<std::vector<ui>> paths;
    std::vector<std::vector<ui>> stars;
    std::vector<std::vector<ui>> singletons;

    // order
    ContractionOrder order;
    generateContractionPlan(graphType, order);

    // contraction
    // super nodes
    for(ui i = 0; i < 5; i++) {
        SuperNodeType sntype = order[i].second;
        ContractionFunc CFP = order[i].first;

        switch(sntype) {
            case SuperNodeType::CHORD : {
                CFP(graph, chords, visited);
                break;
            }
            case SuperNodeType::CYCLE : {
                CFP(graph, cycles, visited);
                break;
            }
            case SuperNodeType::PATH : {
                CFP(graph, paths, visited);
                break;
            }
            case SuperNodeType::STAR : {
                CFP(graph, stars, visited);
                break;
            }
            case SuperNodeType::SINGLETON : {
                CFP(graph, singletons, visited);
                break;
            }
        }
    }

    // build graph
    // first build derived class "SuperGraph" then base class "Base"
    // super node id follows the alphabetic order of their names
    // chord, cycle, path, star with singleton as the last
    ui superNodeId = 0, tempId = 0;
    superGraph->numChord = chords.size();
    superGraph->numCycle = cycles.size();
    superGraph->numPath = paths.size();
    superGraph->numStar = stars.size();
    superGraph->numSingleton = singletons.size();
    tempId = superGraph->numChord + superGraph->numCycle + superGraph->numPath
            + superGraph->numStar + superGraph->numSingleton;

    superGraph->vertexSuperNodeId = new ui[graph->numVertex];
    superGraph->vertexSuperNodePos = new ui[graph->numVertex];
    std::fill(superGraph->vertexSuperNodeId, superGraph->vertexSuperNodeId + graph->numVertex, tempId + 1);
    std::fill(superGraph->vertexSuperNodePos, superGraph->vertexSuperNodePos + graph->numVertex, tempId + 1);

    ui cpv = tempId + 1;
    updateSuperNode(superGraph, chords, superNodeId, cpv);
    updateSuperNode(superGraph, cycles, superNodeId, cpv);
    updateSuperNode(superGraph, paths, superNodeId, cpv);
    updateSuperNode(superGraph, stars, superNodeId, cpv);
    updateSuperNode(superGraph, singletons, superNodeId, cpv);

    // check point
    int cp1 = std::count(superGraph->vertexSuperNodeId, superGraph->vertexSuperNodeId + graph->numVertex, tempId + 1);
    int cp2 = std::count(superGraph->vertexSuperNodePos, superGraph->vertexSuperNodePos + graph->numVertex, tempId + 1);
    if(cp1 != 0 || cp2 != 0) {
        std::cerr << "Check point on building super nodes fails: " << cp1 << " " << cp2 << std::endl;
        exit(1);
    }

    ui numOriginalVertex = graph->numVertex;
    using ExpandListEntry = std::pair<ui, std::vector<std::pair<ui, ui>>>;
    std::vector<std::vector<ui>> vdegree(superNodeId, std::vector<ui>(MAX_PATTERN_SIZE * 2, 0));
    std::vector<std::vector<ui>> adjacencyList(superNodeId);
    std::vector<std::vector<ExpandListEntry>> expandList(superNodeId);

    for(ui i = 0; i < numOriginalVertex; i++) {
        ui inode = superGraph->vertexSuperNodeId[i];
        ui ipos = superGraph->vertexSuperNodePos[i];
        ui count_ = 0;
        const ui *neighbors_ = graph->getVertexNeighbors(i, count_);

        for(ui j = 0; j < count_; j++) {
            ui to = neighbors_[j];
            ui tonode = superGraph->vertexSuperNodeId[to];
            ui topos = superGraph->vertexSuperNodePos[to];

            if(inode == tonode)
                continue;

            auto fiter = std::find(adjacencyList[inode].begin(), adjacencyList[inode].end(), tonode);
            if(fiter != adjacencyList[inode].end()) {
                // edge
                int idis = std::distance(adjacencyList[inode].begin(), fiter);
                expandList[inode][idis].second.emplace_back(ipos, topos);
                // fiter = std::find(adjacencyList[tonode].begin(), adjacencyList[tonode].end(), inode);
                // assert(fiter != adjacencyList[tonode].end());
                // int todis = std::distance(adjacencyList[tonode].begin(), fiter);
                // expandList[tonode][todis].second.emplace_back(topos, ipos);

                // degree
                ++ vdegree[inode][ipos];
                // ++ vdegree[tonode][topos];
                continue;
            }

            // edge
            adjacencyList[inode].emplace_back(tonode);
            expandList[inode].emplace_back(tonode, std::vector<std::pair<ui, ui>>());
            expandList[inode][expandList[inode].size() - 1].second.emplace_back(ipos, topos);
            // adjacencyList[tonode].emplace_back(inode);
            // expandList[tonode].emplace_back(inode, std::vector<std::pair<ui, ui>>());
            // expandList[tonode][expandList[tonode].size() - 1].second.emplace_back(topos, ipos);

            //degree
            ++ vdegree[inode][ipos];
            // ++ vdegree[tonode][topos];
        }
    }

    ui ecount = 0;
    for(auto &alist : adjacencyList) {
        ecount += alist.size();
        std::sort(alist.begin(), alist.end());
    }
    for(auto &elist : expandList) {
        std::sort(elist.begin(), elist.end(), [](const ExpandListEntry &entry1, const ExpandListEntry &entry2) {
            return entry1.first < entry2.first;
        });
        for(auto &oelist : elist) {
            std::sort(oelist.second.begin(), oelist.second.end());
            // check point
            auto uiter = std::unique(oelist.second.begin(), oelist.second.end());
            if(uiter != oelist.second.end()) {
                std::cerr << "check point on expand list containing duplicate entities fails" << std::endl;
                for(const auto &e : oelist.second) std::cerr << e.first << " " << e.second << std::endl;
                exit(1);
            }
        }
    }
    assert(ecount % 2 == 0);

    superGraph->edges = new std::pair<ui, ui>[ecount];
    superGraph->offsets = new ui[superNodeId + 1];
    superGraph->neighbors = new ui[ecount];

    // both edges & CSR consume twice space as numEdge
    // edges
    ui ecnt = 0;
    for(ui i = 0; i < superNodeId; i++) {
        for(const auto &to : adjacencyList[i]) {
            superGraph->edges[ecnt].first = i;
            superGraph->edges[ecnt ++].second = to;
        }
    }

    std::cout << "building csr matrix" << std::endl;
    // CSR
    superGraph->offsets[0] = 0;
    for(ui i = 0; i < superNodeId; i++)
        superGraph->offsets[i + 1] = superGraph->offsets[i] + adjacencyList[i].size();

    std::vector<ui> tempOffsets(superNodeId, 0);
    for(ui i = 0; i < ecount; i++) {
        ui s = superGraph->edges[i].first;
        ui t = superGraph->edges[i].second;
        ui offset = superGraph->offsets[s] + tempOffsets[s];
        superGraph->neighbors[offset] = t;
        ++ tempOffsets[s];
    }

    // base's variables
    superGraph->numVertex = superNodeId;
    superGraph->numEdge = ecount / 2;
    superGraph->numAttribute = graph->numAttribute;
    superGraph->numOriginalVertex = graph->numVertex;
    superGraph->numOriginalEdge = graph->numEdge;

    for(ui i = 0; i < superNodeId; i++)
        std::sort(superGraph->neighbors + superGraph->offsets[i], superGraph->neighbors + superGraph->offsets[i + 1]);

    std::cout << "building node index" << std::endl;
    // build index
    // when build synopsis for singleton
    // the constructor needs argument of a vector wrapper of singleton rather than itself
    superGraph->synopses = new Synopsis[superNodeId];
    superGraph->fcprime = new InverseFunction[superNodeId];
    superGraph->fd = new DecontractionFunction[ecount];

    // super node index
    updateSuperNodeIndex(superGraph, chords, graph, SuperNodeType::CHORD);
    updateSuperNodeIndex(superGraph, cycles, graph, SuperNodeType::CYCLE);
    updateSuperNodeIndex(superGraph, paths, graph, SuperNodeType::PATH);
    updateSuperNodeIndex(superGraph, stars, graph, SuperNodeType::STAR);
    updateSuperNodeIndex(superGraph, singletons, graph, SuperNodeType::SINGLETON);

    std::cout << "building edge index" << std::endl;
    // super edge index
    std::vector<std::vector<ui>> _offsets(superNodeId);
    for(ui i = 0; i < superNodeId; i++) {
        ui numVertex_ = superGraph->fcprime[i].numVertex;
        _offsets[i].resize(numVertex_ + 1, 0);
        for(ui j = 0; j < numVertex_; j++)
            _offsets[i][j + 1] = _offsets[i][j] + vdegree[i][j];
    }

    for(ui i = 0; i < superNodeId; i++) {
        ui count_ = 0;
        const ui *neighbors_ = superGraph->getNodeNeighbors(i, count_);
        for(ui j = 0; j < count_; j++) {
            ui to = neighbors_[j];
            ui eid = superGraph->getSuperEdgeId(i, j);
            for(auto &oedge : expandList[i][j].second) {
                if(oedge.first >= superGraph->fcprime[i].vertices.size() || 
                    oedge.second >= superGraph->fcprime[to].vertices.size()) {
                    std::cerr << oedge.first << " " << oedge.second << " "
                    << superGraph->fcprime[i].vertices.size() << " "
                    << superGraph->fcprime[to].vertices.size() << " "
                    << i << " " << to << std::endl;
                    exit(1);    
                }
                oedge.first = superGraph->fcprime[i].vertices[oedge.first];
                oedge.second = superGraph->fcprime[to].vertices[oedge.second];
            }
            superGraph->fd[eid].updateInfo(eid, _offsets[i], expandList[i][j].second);
        }
    }

    bool checkres = selfCheck(superGraph);
    if(!checkres) {
        std::cerr << "Self check fails on building super graph" << std::endl;
        exit(1);
    }
    
    // release
    releaseBuffers(graph, visited);
}