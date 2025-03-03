#include "pets/contract_graph.h"
#include "pets/generate_pattern.h"
#include "pets/select_pattern.h"
#include "pets/select_pattern.h"
#include "util/gen_query.hpp"
#include "util/pattern_io.hpp"
#include <sys/time.h>
// #define NDEBUG
#include <assert.h>

ui Pattern::curId = 0;

// timer
timeval ioBegin, ioEnd;
timeval conBegin, conEnd;
timeval genBegin, genEnd;
timeval selBegin, selEnd;

/*
 * global variable
 */
// main
Graph *graph = nullptr;
GraphType graphType;
SuperGraph *superGraph = nullptr;
std::vector<Pattern> *patterns = nullptr;
Pattern *topkPatterns = nullptr;
// queries
std::vector<Pattern> queries;
std::unordered_map<std::string, ui> numQueries = {
    { "random", 500 }, 
    { "path", 100 },
    { "tree", 100 },
    { "star", 50 },
    { "doublestar", 50 },
    { "cycle", 100 },
    { "flower", 100 } 
};
ui sizeLow = 5, sizeHigh = 30, universeEachSize = 100;


int main(int argc, char *argv[])
{
    const std::string dataname = argv[1];
    const std::string graphtype = argv[2];

    if(graphtype == "social") graphType = GraphType::SOCIAL;
    else if(graphtype == "web") graphType = GraphType::WEB;
    else if(graphtype == "traffic") graphType = GraphType::TRAFFIC;

    // load original graph
    const std::string statfile = "../datasets/" + dataname + "/stat.txt";
    const std::string edgefile = "../datasets/" + dataname + "/edgelist.txt";
    const std::string attrfile = "../datasets/" + dataname + "/attrs.txt";

    graph = new Graph;
    gettimeofday(&ioBegin, NULL);
    graph->loadGraphFromFiles(statfile, edgefile, attrfile);
    gettimeofday(&ioEnd, NULL);
    graph->printMetaData(dataname);

    // contraction
    // superGraph = new SuperGraph;
    // gettimeofday(&conBegin, NULL);
    // ContractGraph::buildSuperGraph(graph, graphType, superGraph);
    // gettimeofday(&conEnd, NULL);
    // superGraph->printGraphMetaData();

    // // generate
    // GeneratePattern::allocateBuffers(superGraph, patterns);
    // gettimeofday(&genBegin, NULL);
    // GeneratePattern::generatePatterns(superGraph, patterns);
    // gettimeofday(&genEnd, NULL);

    // // select
    // ui patternNumVertex = superGraph->numChord + superGraph->numCycle 
    //                     + superGraph->numPath + superGraph->numStar;
    // SelectPattern::allocateBuffers(topkPatterns);
    // gettimeofday(&selBegin, NULL);
    // SelectPattern::calculateCoverage(superGraph, patterns);
    // SelectPattern::calculateCognition(superGraph, patterns);
    // SelectPattern::greedySelect(patterns, topkPatterns, patternNumVertex, graph->numAttribute);
    // gettimeofday(&selEnd, NULL);

    // generate a set of queries
    // update pattern set according to them
    generate_query::generateQueries(sizeLow, sizeHigh, universeEachSize, numQueries, queries, graph);
    
    // save
    const std::string dirNamePattern = "output/" + dataname + "/patterns";
    const std::string dirNameQuery = "output/" + dataname + "/queries";
    std::vector<Pattern> outputQueries;
    outputQueries.emplace_back(queries[0]);
    outputQueries.emplace_back(queries[500]);
    outputQueries.emplace_back(queries[600]);
    outputQueries.emplace_back(queries[700]);
    outputQueries.emplace_back(queries[750]);
    outputQueries.emplace_back(queries[800]);
    outputQueries.emplace_back(queries[900]);
    pattern_io::savePatternsGraphvizFormat(dirNamePattern, topkPatterns);
    pattern_io::savePatternsGraphvizFormat(dirNameQuery, queries);

    // end
    double ioTime = ioEnd.tv_sec - ioBegin.tv_sec + (ioEnd.tv_usec - ioBegin.tv_usec) / 1e6;
    double conTime = conEnd.tv_sec - conBegin.tv_sec + (conEnd.tv_usec - conBegin.tv_usec) / 1e6;
    double genTime = genEnd.tv_sec - genBegin.tv_sec + (genEnd.tv_usec - genBegin.tv_usec) / 1e6;
    double selTime = selEnd.tv_sec - selBegin.tv_sec + (selEnd.tv_usec - selBegin.tv_usec) / 1e6;
    printf("###                     io time: %.4lf\n", ioTime);
    printf("###            contraction time: %.4lf\n", conTime);
    printf("###   pattern generatation time: %.4lf\n", genTime);
    printf("###      pattern selection time: %.4lf\n", selTime);

    // release
    GeneratePattern::releaseBuffers(patterns);
    SelectPattern::releaseBuffers(topkPatterns);

    return 0;
}