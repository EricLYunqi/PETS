#include "graph/graph.h"


void Graph::loadGraphFromFiles(const std::string &statFile, const std::string &edgeFile, 
                               const std::string &attrFile)
{
    FILE *statFp = fopen(statFile.c_str(), "r");
    if(statFp == nullptr) {
        std::cerr << "Cannot open: " << statFile << std::endl;
        exit(1);
    }

    auto fres = fscanf(statFp, "%u %u %u\n", &numVertex, &numEdge, &numAttribute);
    // printf("|V|: %d\t|E|: %d\t|L|: %d\n", (int)numVertex, (int)numEdge, (int)numAttribute);

    fclose(statFp);

    // allocate
    edges = new std::pair<ui, ui>[numEdge];

    offsets = new ui[numVertex + 1];
    neighbors = new ui[numEdge * 2];

    ui maxNumAttr = numAttribute * numVertex;
    attrOffsets = new ui[numVertex + 1];
    attributes = new ui[maxNumAttr];
    attrWeights = new double[maxNumAttr];

    offsets[0] = 0;
    attrOffsets[0] = 0;

    // edge
    FILE *edgeFp = fopen(edgeFile.c_str(), "r");
    if(edgeFp == nullptr) {
        std::cerr << "Cannot open: " << edgeFile << std::endl;
        exit(1);
    }

    std::vector<ui> degree(numVertex, 0);
    std::vector<ui> tempOffsets(numVertex, 0);

    for(ui i = 0; i < numEdge; i++) {
        ui begin;
        ui end;
        fres = fscanf(edgeFp, "%u %u\n", &begin, &end);

        edges[i].first = begin;
        edges[i].second = end;
        ++ degree[begin];
        ++ degree[end];
    }

    fclose(edgeFp);

    for(ui i = 0; i < numVertex; i++) {
        maxDegree = maxDegree >= degree[i] ? maxDegree : degree[i];
        offsets[i + 1] = offsets[i] + degree[i];
    }
    
    for(ui i = 0; i < numEdge; i++) {
        ui begin = edges[i].first;
        ui end = edges[i].second;

        ui offset = offsets[begin] + tempOffsets[begin];
        neighbors[offset] = end;

        offset = offsets[end] + tempOffsets[end];
        neighbors[offset] = begin;

        ++ tempOffsets[begin];
        ++ tempOffsets[end];
    }

    // sort
    for(ui i = 0; i < numVertex; i++)
        std::sort(neighbors + offsets[i], neighbors + offsets[i + 1]);
    std::sort(edges, edges + numEdge);

    // attribute
    FILE *attrFp = fopen(attrFile.c_str(), "r");
    if(attrFp == nullptr) {
        std::cerr << "Cannot open: " << attrFile.c_str() << std::endl;
        exit(1);
    }

    std::fill(tempOffsets.begin(), tempOffsets.end(), 0);

    for(ui i = 0; i < numVertex; i++) {
        ui num;
        fres = fscanf(attrFp, "%u", &num);

        attrOffsets[i + 1] = attrOffsets[i] + num;

        for(ui j = 0; j < num; j++) {
            ui attr;
            double weights;
            fres = fscanf(attrFp, "%u %lf", &attr, &weights);

            ui offset = attrOffsets[i] + tempOffsets[i];
            attributes[offset] = attr;
            attrWeights[offset] = weights;

            ++ tempOffsets[i];
        }
    }

    fclose(attrFp);

    bool res = selfCheck();
    if(!res) {
        std::cerr << "Self check fails" << std::endl;
        exit(1);
    }
}


bool Graph::selfCheck() const
{
    // self circle
    for(ui i = 0; i < numEdge; i++) {
        if(edges[i].first == edges[i].second) {
            std::cerr << "self circle in: " << edges[i].first << std::endl;
            return false;
        }
    }

    // duplicate edge
    auto eiter = std::unique(edges, edges + numEdge);
    if(eiter != edges + numEdge) {
        ui num = std::distance(eiter, edges + numEdge);
        std::cerr << "Duplicate edges: " << num << std::endl;
        return false;
    }

    // vertex's id
    std::vector<bool> tag(numVertex, false);
    for(ui i = 0; i < numEdge; i++) {
        tag[edges[i].first] = true;
        tag[edges[i].second] = true;
    }
    for(ui i = 0; i < numVertex; i++) {
        if(!tag[i]) {
            std::cerr << "Id " << i << " do not exist" << std::endl;
            return false;
        }
    }

    // connectivity
    std::queue<ui> Q;
    std::vector<int> vis(numVertex, 0);

    Q.push(0);
    vis[0] = 1;

    while(!Q.empty()) {
        ui u = Q.front();
        Q.pop();

        ui count;
        const ui *neighbors_ = getVertexNeighbors(u, count);

        for(ui i = 0; i < count; i++) {
            ui to = neighbors_[i];
            if(vis[to] == 0) {
                Q.push(to);
                vis[to] = 1;
            }
        }
    }

    ui totCluster = 1;
    for(ui i = 0; i < numVertex; i++) 
        if(vis[i] == 0)
            ++ totCluster;

    printf("Total clusters: %u\n", totCluster);

    return true;
}


void Graph::printMetaData(const std::string &dataname) const
{
    printf("working on: %s\n", dataname.c_str());
    printf("|V|: %u\t|E|: %u\t|L|: %u\n", numVertex, numEdge, numAttribute);
}


// super graph
void SuperGraph::printGraphMetaData() const 
{
    double ratio = (numVertex + numEdge) * 1.0 / (numOriginalVertex + numOriginalEdge);

    printf("---------- super graph stat ----------\n");
    printf("|V|: %u\t|E|: %u\t|L|: %u\n", numVertex, numEdge, numAttribute);
    printf("original: |V|: %u\t|E|: %u\n", numOriginalVertex, numOriginalEdge);
    printf("|chord|: %u\t|cycle| :%u\t|path|: %u\t|star| :%u\t|sinlgeton| :%u\n", 
           numChord, numCycle, numPath, numStar, numSingleton);
    printf("contraction ratio: %.4lf\n", ratio);
    printf("-------------- stat end --------------\n");
}