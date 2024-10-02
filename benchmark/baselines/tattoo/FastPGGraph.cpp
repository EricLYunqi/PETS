/**
 * Generate patterns from in truss subgraph
 */
#include "FastPGGraph.h"

#include <cstdlib>
#include <fstream>
#include <queue>

using namespace std;

namespace {
  const int COV_FUNC = 2;
  const int COG_FUNC = 6;
}

//#define DEBUG_MODE
int FastPGGraph::MAX_PATTERN_SIZE = 15;
int FastPGGraph::MIN_STAR_DEGREE = 4;  // stars minimum degree

int FastPGGraph::NUM_V;  // the number of vertex of the undirected graph
int FastPGGraph::NUM_E;  // the number of edge of the undirected graph
vector<Pattern> FastPGGraph::patterns;       // the candidate patterns set
vector<graph_t> FastPGGraph::boostPatterns;  // the candidate patterns set in
                                             // boost graph form
vector<PGEdge> FastPGGraph::edges;           //	the edges of the graph
vector<vector<vector<int>>> FastPGGraph::patIdMap = vector<vector<vector<int>>>(
    MAX_COMBINED_TRUSS_CLASS + 1,
    vector<vector<int>>(
        MAX_COMBINED_TRUSS_CLASS + 1,
        vector<int>(3, 0)));  // patterns map to their id,for calculating scores
unordered_map<int, int>
    FastPGGraph::outsideStarPatMap;  // when search stars, stars map to id.
map<pair<int, int>, int>
    FastPGGraph::outsideCombinedStarPatMap;  // when search combined stars,  map
                                             // to id.
unordered_map<int, int>
    FastPGGraph::outsideLinePatMap;  // when search lines,  map to id.
unordered_map<int, int>
    FastPGGraph::outsideCirclePatMap;  ////when search circle,  map to id.
unordered_map<size_t, int>
    FastPGGraph::outsideLinearCombinedStarPatMap;  // when search Linear
                                                   // combined stars,  map to
                                                   // id.
boost::mutex FastPGGraph::io_mutex;    // for lock when using muti threads.
vector<char> FastPGGraph::edgeType;    // the type of each edge
vector<int> FastPGGraph::binStart;     // for bin sort
vector<int> FastPGGraph::binSize;      ////for bin sort
vector<int> FastPGGraph::curLength;    // for bin sort
vector<int> FastPGGraph::sortedEdges;  // //for bin sort
vector<int> FastPGGraph::sortedPos;    // for bin sort
vector<bool>
    FastPGGraph::edgeAlive;  // for marking edges when doing truss marking.
vector<int> FastPGGraph::trussEdgeCnt;  // trussEdgeCnt[i]= the number of edges
                                        // of truss i
vector<vector<vector<SmallComponentInfo>>> FastPGGraph::smallComponents =
    vector<vector<vector<SmallComponentInfo>>>(
        MAX_INDIVIDUAL_NODE_COUNT + 1,
        vector<vector<SmallComponentInfo>>(MAX_SMALL_PATTERN_EDGE_COUNT + 1));
SummaryInfo
    FastPGGraph::summaryInfo;  // for recording summary infomation such as
                               // running time, edges ditrubution and so on.
vector<vector<NeighbVE>>
    FastPGGraph::neibours;            // get the neighbours of  each edges.
vector<int> FastPGGraph::GoutDegree;  // the degree of node in g_out
vector<vector<int>> FastPGGraph::NNH;
vector<int> FastPGGraph::support;
vector<char> FastPGGraph::trussClass;        // the trussness of edges
vector<int> FastPGGraph::finalPatternIdSet;  // the final choosed patterns id
queue<int> FastPGGraph::candidatePatternsQueue;  // for muti-threads, store the
                                                 // candidate patterns.
vector<int> bitCountTable16;   // bitCountTable16[i]= the number of 1 in the
                               // binary form of i.
vector<int> patternUseDistrb;  // for output information, patternUseDistrb[i]=
                               // the number of fianl patterns of size i.
vector<vector<Pattern>> FastPGGraph::queryPatterns;  // only first time to gen
vector<vector<int>> FastPGGraph::circleRecord;
vector<Pattern> FastPGGraph::query2Patterns;
vector<int> FastPGGraph::query2Setting = {
    500, 100, 100, 50,
    50,  100, 100, 1000};  // 0: random, 1: path, 2: tree, 3: single_star, 4:
                           // double_star, 5: cycle, 6: flower, 7: total
vector<string> FastPGGraph::query2SettingName = {
    "random",      "path",  "tree",   "single_star",
    "double_star", "cycle", "flower", "total"};
vector<SummaryInfo> FastPGGraph::newInfo;  // unused now
vector<vector<Pattern>>
    FastPGGraph::finalSet4Datastores;   //   when test com patterns
vector<int> FastPGGraph::tempfinalset;  // temp when
unordered_map<ull, int> ullUsed2patternId;
bool FastPGGraph::LOADED = false;
vector<int> FastPGGraph::patternUseClassCnt;
vector<vector<int>> FastPGGraph::Step4curQuery;
ofstream totout("data/tot_result.txt");
// double g_cogmin, g_cogmax;

void FastPGGraph::init() {
  srand(time(0));
  edges.clear();
  bitCountTable16.clear();
  bitCountTable16.emplace_back(0);
  for (ull i = 1; i < 1ULL << 16; ++i)
    bitCountTable16.emplace_back((i & 1) + bitCountTable16[i >> 1]);
  neibours.clear();
  newInfo.clear();
  trussClass.clear();
  patternUseClassCnt.clear();
  patternUseClassCnt.resize(7, 0);
  Step4curQuery.clear();
  Step4curQuery.resize(31);
  for (int i = 0; i < 31; i++) Step4curQuery[i].resize(100, 0);
  /*others */
  edgeType.clear();
  patterns.clear();
  patIdMap = vector<vector<vector<int>>>(
      MAX_COMBINED_TRUSS_CLASS + 1,
      vector<vector<int>>(
          MAX_COMBINED_TRUSS_CLASS + 1,
          vector<int>(3,
                      0)));  // patterns map to their id,for calculating scores
  outsideStarPatMap.clear();
  outsideCombinedStarPatMap.clear();
  outsideLinePatMap.clear();
  outsideCirclePatMap.clear();
  outsideLinearCombinedStarPatMap.clear();
#ifndef linux
  io_mutex.initialize();
#endif
  GoutDegree.clear();
  smallComponents = vector<vector<vector<SmallComponentInfo>>>(
      MAX_INDIVIDUAL_NODE_COUNT + 1,
      vector<vector<SmallComponentInfo>>(MAX_SMALL_PATTERN_EDGE_COUNT + 1));
  queryPatterns.clear();
  /*get final patterns*/
  summaryInfo.init();
  boostPatterns.clear();
  finalPatternIdSet.clear();
  patternUseDistrb.clear();
}

void FastPGGraph::init_part1() {
  vector<Pattern>().swap(patterns);
  vector<char>().swap(edgeType);
  vector<char>().swap(trussClass);
  patIdMap = vector<vector<vector<int>>>(
      MAX_COMBINED_TRUSS_CLASS + 1,
      vector<vector<int>>(
          MAX_COMBINED_TRUSS_CLASS + 1,
          vector<int>(3,
                      0)));  // patterns map to their id,for calculating scores
  unordered_map<int, int>().swap(outsideStarPatMap);
  map<pair<int, int>, int>().swap(outsideCombinedStarPatMap);
  unordered_map<int, int>().swap(outsideLinePatMap);
  unordered_map<int, int>().swap(outsideCirclePatMap);
#ifndef linux
  io_mutex.initialize();
#endif
  vector<int>().swap(GoutDegree);
  smallComponents = vector<vector<vector<SmallComponentInfo>>>(
      MAX_INDIVIDUAL_NODE_COUNT + 1,
      vector<vector<SmallComponentInfo>>(MAX_SMALL_PATTERN_EDGE_COUNT + 1));
  summaryInfo.init();
  vector<graph_t>().swap(boostPatterns);
#ifndef linux
  io_mutex.initialize();
#endif
  vector<int>().swap(finalPatternIdSet);
  patternUseDistrb.clear();
}

bool mycmp(const NeighbVE &a, const NeighbVE &b) {
  return a.neighbId < b.neighbId;
}

// pair hash to a size_t.
struct PairHash {
  size_t operator()(const std::pair<int, int> &v) const {
    std::hash<int> hasher;
    size_t seed = 0;
    seed = hasher(v.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed = hasher(v.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

// read the graph with num of node and edge
void FastPGGraph::loadGraphFromFileAndMap(string filename) {
  string edgeFileName = filename + "edgelist.txt";
  string statFileName = filename + "stat.txt";
  ifstream edgeFin(edgeFileName);
  ifstream statFin(statFileName);
  cout << "New style to input the graph..." << endl;
  int s, e, dummy;
  int edgeId = 0;
  trussEdgeCnt.resize(20, 0);
  statFin >> NUM_V >> NUM_E >> dummy;
  neibours.resize(NUM_V + 1);
  unordered_map<int, int> nodeMap;  // all nodes map to [0,numV-1](because the
                                    // node may not be continuously).
  unordered_set<pair<int, int>, PairHash> edgehash;  // against repeated edges
  int nodeId = 0;
  while (edgeFin >> s >> e) {
    auto itr = nodeMap.find(s);
    if (itr == nodeMap.end())
      itr = nodeMap.emplace(make_pair(s, nodeId++)).first;
    s = itr->second;
    itr = nodeMap.find(e);
    if (itr == nodeMap.end())
      itr = nodeMap.emplace(make_pair(e, nodeId++)).first;
    e = itr->second;

    if (s != e && edgehash.find(make_pair(s, e)) == edgehash.end()) {
      PGEdge edge = PGEdge(s, e);
      edges.emplace_back(edge);

      neibours[s].emplace_back(NeighbVE(e, edgeId));
      neibours[e].emplace_back(NeighbVE(s, edgeId));

      edgehash.emplace(make_pair(s, e));
      edgehash.emplace(make_pair(e, s));
      ++edgeId;
    }
  }
  NUM_E = edges.size();
  assert(nodeId == NUM_V);
  assert(edges.size() == NUM_E);
  // sort the neibours in order, for speed up when search.
  for (int i = 0; i < NUM_V; ++i) {
    sort(neibours[i].begin(), neibours[i].end(), mycmp);
    // neibours[i].shrink_to_fit();
  }
  // neibours.shrink_to_fit();
  // edges.shrink_to_fit();

  edgeFin.close();
  statFin.close();
}

/******************************  for output information
 * ***************************************************************************/
void FastPGGraph::outputGin4Gephi(string filename) {
  ofstream fout(filename);
  fout << "dl\nformat = edgelist1\nn=" << NUM_V << "\ndata:\n";
  for (int i = 0; i < NUM_E; ++i)
    if (edgeType[i] == EdgeType::RealTruss)
      fout << (edges[i].s + 1) << " " << (edges[i].e + 1) << "\n";
  fout.close();
}

void FastPGGraph::outputGout4Gephi(string filename) {
  ofstream fout(filename);
  fout << "dl\nformat = edgelist1\nn=" << NUM_V << "\ndata:\n";
  for (int i = 0; i < NUM_E; ++i)
    if (edgeType[i] != EdgeType::RealTruss)
      fout << (edges[i].s + 1) << " " << (edges[i].e + 1) << "\n";
  fout.close();
}

void FastPGGraph::outputGr4Gephi(string filename) {
  ofstream fout(filename);
  fout << "dl\nformat = edgelist1\nn=" << NUM_V << "\ndata:\n";
  for (int i = 0; i < NUM_E; ++i)
    if (edgeType[i] != EdgeType::RealTruss &&
        edgeType[i] != EdgeType::SimpleStar)
      fout << (edges[i].s + 1) << " " << (edges[i].e + 1) << "\n";
  fout.close();
}

/**
 * output outside edges
 * @param filename
 */
void FastPGGraph::outputOutsideEdgeInfo(string filename) {
  ofstream fout(filename);
  for (int i = 0; i < NUM_E; ++i)
    fout << edges[i].s << " " << edges[i].e << " " << trussClass[i] << "\n";
  fout.close();
}

/**
 * output remianed outside edges after removing simple patterns
 * @param filename
 */
void FastPGGraph::outputOutsideEdgeInfoRemoved(string filename) {
  ofstream fout(filename);
  for (int i = 0; i < NUM_E; ++i)
    if (trussClass[i] <= 2 && trussClass[i] > 0)
      fout << edges[i].s << " " << edges[i].e << "\n";
  fout.close();
}

/**
 * output remianed outside edges after removing star and combinedstar and linear
 * @param filename
 */
void FastPGGraph::outputOutsideEdgeInfoBesideStarsAndLines(string filename) {
  ofstream fout(filename);
  for (int i = 0; i < NUM_E; ++i)
    if (edgeType[i] == EdgeType::Default)
      fout << edges[i].s << " " << edges[i].e << "\n";
  fout.close();
}

/**
 * output patterns to file which can be read by the Davinci GUI
 * @param filename
 */
void FastPGGraph::outputPatterns(string filename) {
  ofstream fout(filename);
  for (int i = 0; i < patterns.size(); i++) {
    Pattern &pat = patterns[i];
    fout << "#id: " << pat.id << " |V|: " << pat.getNodesCnt()
         << " |E|: " << pat.getEdgeCnt() << " type: " << pat.getType()
         << " frequency: " << pat.getScore() << " cov:  " << pat.scoreCoverage
         << " cog: " << pat.scoreCognition << "\n";
  }
  fout.close();
}

void FastPGGraph::outputGraph4GUI(string filename) {
  ofstream fout(filename);
  fout << NUM_V << " " << NUM_E << endl;
  for (int i = 0; i < NUM_V; ++i)
    fout << "v " << i << " " << neibours[i].size() << endl;
  int c = 0;
  for (int i = 0; i < NUM_V; ++i)
    for (auto &e : neibours[i])
      if (i < e.neighbId) {
        fout << "e " << i << " " << e.neighbId << endl;
        ++c;
      }
  assert(c == NUM_E);
  fout.close();
}

/**
 * output patterns to file which can be read by the Davinci GUI
 * @param filename
 */
void FastPGGraph::outputFinalPatterns4GUI(string filename) {
  ofstream fout(filename);
  int score = 100;
  for (int i : finalPatternIdSet) {
    Pattern &pat = patterns[i];
    fout << "#" << i << " " << pat.getNodesCnt() << " " << pat.getEdgeCnt()
         << " " << pat.getType() << " " << score-- << endl;
    for (auto &edge : pat.edges) {
      fout << edge.s << " " << edge.e << endl;
    }
  }
  fout.close();
}

void createFolder(const string &folder) {
#ifdef linux
  { auto _ = system(("mkdir -p " + folder).c_str()); }
#elif defined(_WIN32) || defined(_WIN64)
  string nfolder = folder;
  for (int i = 0; i < nfolder.length(); ++i)
    if (nfolder[i] == '/') nfolder[i] = '\\';
  { auto _ = system(("mkdir " + nfolder + " 2> nul").c_str()); }
#endif
}

void FastPGGraph::outputPngPattern(Pattern &pattern, const string &folder) {
  createFolder(folder);
  ofstream fout(folder + "dotfile");
  fout << "graph \"result\" {\n";
  fout << "graph [fontsize = 30, ratio = 0.835, dpi = 100, size = \"15,15\" "
          "];\n";
  fout << "node [label = \"\\N\", shape = doublecircle, sides = 4, color = "
          "skyblue, style = filled ];\n";
  for (int i = 0; i < pattern.nodesCnt; ++i)
    fout << i << " [shape = doublecircle, label = \"" << i << "\"];\n";
  for (auto &e : pattern.edges) fout << e.s << " -- " << e.e << ";\n";
  fout << "}\n";
  fout.close();
  string command =
      "dot -Tpng  -Kneato -Gepsilon=0.0001 -Goverlap=false " + folder +
      "dotfile -o " + folder + to_string(pattern.finalRank) + "-" +
      to_string(pattern.type) + "-" + to_string(pattern.getEdgeCnt()) + "-" +
      to_string(pattern.id) + "-" + to_string(pattern.getScore()) + ".png";
  { auto _ = system(command.c_str()); }
}

void FastPGGraph::outputFinalPngPattern(Pattern &pattern, const string &folder,
                                        string filename) {
  if (pattern.getNodesCnt() > 50) return;
  createFolder(folder);
  ofstream fout(folder + "dotfile");
  fout << "graph \"result\" {\n";
  fout << "graph [fontsize = 30, ratio = 0.835, dpi = 100, size = \"15,15\" "
          "];\n";
  fout << "node [label = \"\\N\", shape = doublecircle, sides = 4, color = "
          "skyblue, style = filled ];\n";
  for (int i = 0; i < pattern.nodesCnt; ++i)
    fout << i << " [shape = doublecircle, label = \"" << i << "\"];\n";
  for (auto &e : pattern.edges) fout << e.s << " -- " << e.e << ";\n";
  fout << "}\n";
  fout.close();
  string command = "dot -Tpng  -Kneato -Gepsilon=0.0001 -Goverlap=false " +
                   folder + "dotfile -o " + folder + filename + ".png";
  { auto _ = system(command.c_str()); }
}

/**
 * png file name format: type-edgeCnt-id-score
 * @param filename
 */
void FastPGGraph::outputPngPatterns(string filename) {
#ifdef linux
  { auto _ = system(("mkdir " + filename).c_str()); }
#endif
  printf("outputPngPatterns......\n");
  for (auto &pattern : patterns)
    if (pattern.score > 1 && pattern.getEdgeCnt() <= 22) {
      outputPngPattern(pattern, filename);
    }
}

void FastPGGraph::outputSummaryInfo(string filename) {
  outputSummaryInfo(filename, "", false);
}

void FastPGGraph::outputSummaryInfo(string filename, string datasetName,
                                    bool append) {
  ofstream fout;
  if (append)
    fout.open(filename, ios_base::app);
  else
    fout.open(filename);
  fout << "===============================" << datasetName
       << "===============================\n";

  auto totalTime =
      summaryInfo.readGraphTime + summaryInfo.calTrussClass4EdgesTime +
      summaryInfo.calScore4CombinedTrussPatternsTime +
      summaryInfo.searchStarsTime +
      summaryInfo.searchLinesAndCirclesAndIndividualsTime +
      summaryInfo.calCovAndCogScoreTime + summaryInfo.finalPatternsGenTime;
  fout << "Nodes number: " << NUM_V << "\n";
  fout << "Edges number: " << edges.size() << "\n";
  fout << "G_in edges number: " << summaryInfo.insideEdgeCnt << "("
       << summaryInfo.insideEdgeCnt * 1.0 / edges.size() << "%)\n";
  fout << "G_out edges number: " << summaryInfo.outsideEdgeCnt << "("
       << summaryInfo.outsideEdgeCnt * 1.0 / edges.size() << "%)\n";
  fout << "candidate patterns number: " << patterns.size() << "\n";
  fout << "total time : " << totalTime << "s\n";
  fout << "\n";

  fout << "Time of reading graph: " << summaryInfo.readGraphTime << "s("
       << summaryInfo.readGraphTime / totalTime << "%)\n";
  fout << "Time of calculating truss number for edges: "
       << summaryInfo.calTrussClass4EdgesTime << "s("
       << summaryInfo.calTrussClass4EdgesTime / totalTime << "%)\n";
  fout << "Time of calculating score for truss: "
       << summaryInfo.calScore4CombinedTrussPatternsTime << "s("
       << summaryInfo.calScore4CombinedTrussPatternsTime / totalTime << "%)\n";
  fout << "Time of searching stars: " << summaryInfo.searchStarsTime << "s("
       << summaryInfo.searchStarsTime / totalTime << "%)\n";
  fout << "Time of searching lines, circles and individuals: "
       << summaryInfo.searchLinesAndCirclesAndIndividualsTime << "s("
       << summaryInfo.searchLinesAndCirclesAndIndividualsTime / totalTime
       << "%)\n";
  fout << "Time of cal Cov and Cog: " << summaryInfo.calCovAndCogScoreTime
       << "s(" << summaryInfo.calCovAndCogScoreTime / totalTime << "%)\n";
  fout << "Time of generating final patterns : "
       << summaryInfo.finalPatternsGenTime << "s("
       << summaryInfo.finalPatternsGenTime / totalTime << "%)\n";
  fout << "\n";

  fout << "G_out\t\t\tEdge No.\tComponent No.\n";
  fout << "Stars:\t\t\t" << summaryInfo.starsEdgeCnt << "\t"
       << summaryInfo.starsComponentCnt << "\n";
  fout << "Lines:\t\t\t" << summaryInfo.linesEdgeCnt << "\t"
       << summaryInfo.linesComponentCnt << "\n";
  fout << "Circles:\t\t" << summaryInfo.circlesEdgeCnt << "\t"
       << summaryInfo.circlesComponentCnt << "\n";
  fout << "Individuals:\t\t" << summaryInfo.individualsEdgeCnt << "\t"
       << summaryInfo.individualsComponentCnt << "\n";
  fout << "Uncovered:\t\t" << summaryInfo.uncoveredEdgeCnt << "\t"
       << summaryInfo.uncoveredComponentCnt << "\n";
  fout.close();
}

/*************************************generate some candidate
 * patterns*******************************************************************************/
/**
 * generate default patterns
 */
void FastPGGraph::genDefaultPatterns() {
  Pattern patOneEdge = Pattern(2);
  patOneEdge.addEdge(0, 1);
  patOneEdge.setType(PatType::OneEdge);
  patOneEdge.detilType = make_pair(1, -1);
  patterns.emplace_back(patOneEdge);

  Pattern patTwoEdge = Pattern(3);
  patTwoEdge.addEdge(0, 1);
  patTwoEdge.addEdge(1, 2);
  patTwoEdge.setType(PatType::TwoEdge);
  patTwoEdge.detilType = make_pair(2, -1);
  patterns.emplace_back(patTwoEdge);

  Pattern patTriangle = Pattern(3);
  patTriangle.addEdge(0, 1);
  patTriangle.addEdge(1, 2);
  patTriangle.addEdge(2, 0);
  patTriangle.setType(PatType::Triangle);
  patTriangle.detilType = make_pair(3, -1);
  patterns.emplace_back(patTriangle);

  Pattern patRectangle = Pattern(4);
  patRectangle.addEdge(0, 1);
  patRectangle.addEdge(1, 2);
  patRectangle.addEdge(2, 3);
  patRectangle.addEdge(3, 0);
  patRectangle.setType(PatType::Triangle);  // keep order, not OutsideCircle
  patRectangle.detilType = make_pair(4, -1);
  patterns.emplace_back(patRectangle);
}

/**
 * for each k-truss, generate the truss-generated patterns (one edge with k - 2
 * triangles around)
 * @param k
 */
void FastPGGraph::genTrussPattern(int k) {
  Pattern pattern = Pattern(k);
  pattern.addEdge(0, 1);
  for (int id = 2; id < k; id++) {
    pattern.addEdge(0, id);
    pattern.addEdge(1, id);
  }
  pattern.setType(PatType::Truss);
  pattern.setScore(trussEdgeCnt[k]);
  pattern.setTrussClass(k);
  pattern.detilType = make_pair(k, -1);
  patterns.emplace_back(pattern);
}

/**
 * generate all the truss-generated patterhs for 4 <= k <= 9
 */
void FastPGGraph::genTrussPatterns() {
  // currently hard code the max truss class to 9
  for (int k = 4; k <= 9; k++) {
    genTrussPattern(k);
  }
}
/**
 * for k1, k2, generate four kinds of combined truss
 * @param k1
 * @param k2
 * @param type, type == 1: k1-side-k2-main; type == 2:  k1-side-k2-side (2
 * kinds)
 */
void FastPGGraph::genCombinedPattern(int k1, int k2, int type) {
  if (k1 == 3 && k2 == 3) return;
  if (type == 1)  // k1-side-k2-main;
  {
    Pattern pattern = Pattern(k1 + k2 - 2);
    // the first pattern
    pattern.addEdge(0, 1);
    for (int i = 2; i < k1; i++) {
      pattern.addEdge(0, i);
      pattern.addEdge(1, i);
    }
    // second pattern based on (0,2)
    for (int i = k1; i <= k1 + k2 - 3; i++) {
      pattern.addEdge(0, i);
      pattern.addEdge(2, i);
    }
    // main with side
    pattern.setType(PatType::TrussCombined1);
    // pattern.setSurroundings();					//only
    // used when it is TrussExtend
    pattern.setScore(0);  // init it 0
    int id = patterns.size();
    patIdMap[k1][k2][0] = id;
    pattern.detilType = make_pair(k1, k2);
    patterns.emplace_back(pattern);
  } else {
    // side to side, main edge connected
    Pattern pattern = Pattern(k1 + k2 - 2);
    // the first pattern
    pattern.addEdge(0, 1);
    for (int i = 2; i < k1; i++) {
      pattern.addEdge(0, i);
      pattern.addEdge(1, i);
    }
    // the second pattern based on (0, k1) and shares an edge with (0, 2)
    pattern.addEdge(0, k1);
    pattern.addEdge(2, k1);
    for (int i = k1 + 1; i <= k1 + k2 - 3; i++) {
      pattern.addEdge(0, i);
      pattern.addEdge(k1, i);
    }
    pattern.setType(PatType::TrussCombined2);
    // pattern.setSurroundings();
    pattern.setScore(0);
    int id = patterns.size();
    patIdMap[k1][k2][1] = id;
    pattern.detilType = make_pair(k1, k2);
    patterns.emplace_back(pattern);

    // side to side, main edge not connected
    pattern = Pattern(k1 + k2 - 2);
    // the first pattern
    pattern.addEdge(0, 1);
    for (int i = 2; i < k1; i++) {
      pattern.addEdge(0, i);
      pattern.addEdge(1, i);
    }
    // the second pattern based on (2, k1) and shares an edge with (0, 2)
    pattern.addEdge(2, k1);
    pattern.addEdge(0, k1);
    for (int i = k1 + 1; i <= k1 + k2 - 3; i++) {
      pattern.addEdge(2, i);
      pattern.addEdge(k1, i);
    }
    pattern.setType(PatType::TrussCombined2);
    // pattern.setSurroundings();
    pattern.setScore(0);
    id = patterns.size();
    patIdMap[k1][k2][2] = id;
    pattern.detilType = make_pair(k1, k2);
    patterns.emplace_back(pattern);
  }
}

/**
 * generate all the combined truss sharing one edge
 */

void FastPGGraph::genCombinedTrussPatterns() {
  genFourMainThreeSidePattern();   // *one special case for combined pattern 4+3
  for (int k1 = 4; k1 <= 9; k1++)  // * as a matter of fact,when k1+k2>11 break
  {                                // * only gen （4,4/5/6/7）,(5,5/6)
    for (int k2 = k1; k2 <= 9; k2++) {
      int size1 = (k1 - 2) * 2 + 1;  // * ??
      int size2 = (k2 - 2) * 2 + 1;
      // the maximum size is set to 15
      if (size1 + size2 - 1 > MAX_PATTERN_SIZE) break;
      // small side + big main
      genCombinedPattern(k1, k2, 1);
      if (k1 != k2) {
        // big side + small main
        genCombinedPattern(k2, k1, 1);
      }
    }
  }
  for (int k1 = 4; k1 <= 9; k1++) {
    for (int k2 = k1; k2 <= 9; k2++) {
      int size1 = (k1 - 2) * 2 + 1;
      int size2 = (k2 - 2) * 2 + 1;
      if (size1 + size2 - 1 > MAX_PATTERN_SIZE) break;
      // side + side
      genCombinedPattern(k1, k2, 2);
    }
  }
}

/**
 * one special case for combined truss 4+3
 */
void FastPGGraph::genFourMainThreeSidePattern() { genCombinedPattern(4, 3, 1); }

/*
input (s,e) ,return the edge_id.
 because neibours is sorted, so ,we can using binary search
 if not find ,return -1.
*/
inline int FastPGGraph::getEdgebyNode(int s, int e) {
  if (neibours[s].size() > neibours[e].size()) {
    int t = s;
    s = e;
    e = t;
  }
  int l = 0, r = neibours[s].size(), mid;
  while (l <= r) {
    mid = (l + r) >> 1;
    if (neibours[s][mid].neighbId == e)
      return neibours[s][mid].edgeId;
    else if (neibours[s][mid].neighbId < e)
      l = mid + 1;
    else
      r = mid - 1;
  }
  assert(0);
  return -1;
}

/**
 * get the Node Neighborhood of an edge.  what is more ,the array are indexed .
 * @param edgeId
 * @return
      we using one-d array to save memery.
      for each edge, nnh[edge][k-3] means the exatly truss k NNH.
 */
void FastPGGraph::getNNH(int &edgeId) {
  int s = edges[edgeId].s;
  int e = edges[edgeId].e;
  int i = 0, j = 0;
  int lens = neibours[s].size();
  int lene = neibours[e].size();
  int minTruss;
  int numCnt = 0;
  vector<int> prefixSum(MAX_COMBINED_TRUSS_CLASS + 2);
  vector<int> cntP(MAX_COMBINED_TRUSS_CLASS + 1, 0);
  vector<int> nodesTruss;
  vector<int> oldNNH;
  while (j < lene && i < lens)  //
  {
    if (neibours[s][i].neighbId == neibours[e][j].neighbId) {
      minTruss = trussClass[neibours[s][i].edgeId];
      int tmp = trussClass[neibours[e][j].edgeId];
      if (tmp < minTruss) minTruss = tmp;
      nodesTruss.emplace_back(minTruss);
      oldNNH.emplace_back(neibours[s][i].neighbId);
      if (minTruss > 2)  // lessthan 2 is meaningless
      {
        ++cntP[minTruss];
        ++numCnt;
      }
      ++i;
      ++j;
    } else if (neibours[s][i].neighbId < neibours[e][j].neighbId)
      ++i;
    else
      ++j;
  }

  prefixSum[2] = 0;
  for (int i = 3; i < MAX_COMBINED_TRUSS_CLASS + 2; ++i) {
    prefixSum[i] = prefixSum[i - 1] + cntP[i - 1];
    cntP[i - 1] = 0;
  }
  vector<int> indexedNNH(MAX_COMBINED_TRUSS_CLASS + 2 + numCnt - 3);
  for (int i = 0; i < oldNNH.size(); ++i) {
    int w = oldNNH[i];
    minTruss = nodesTruss[i];
    if (minTruss > 2)  // lessthan 2 is meaningless
    {
      int local = MAX_COMBINED_TRUSS_CLASS + 1 + prefixSum[minTruss] +
                  cntP[minTruss] + 1;
      indexedNNH[local - 3] = w;
      ++cntP[minTruss];
    }
  }
  for (int i = 3; i < MAX_COMBINED_TRUSS_CLASS + 2; ++i) {
    indexedNNH[i - 3] = MAX_COMBINED_TRUSS_CLASS + 1 + prefixSum[i] + 1;
  }
  boost::mutex::scoped_lock lock(
      io_mutex);  // must lock:  many thread write it  at the same time
  NNH[edgeId] = indexedNNH;
}

/**
 * k1-side-k2-main pattern
 * @param k, k of the first truss
 * @param originalEdgeId, main edge of the first pattern
 * @param edgeId, the shared edgeId
 * @param nodeSet, the node neighborhood of originEdge
 * @param nodeSet, the node neighborhood of originEdge
 * @param Tid, the current thread id
 */
void FastPGGraph::addScoreForComPatternWithComEdge1(int k, int originalEdgeId,
                                                    int sharedEdgeId,
                                                    unordered_set<int> &nodeSet,
                                                    int &Tid) {
  int oriS = edges[originalEdgeId].s;
  int oriE = edges[originalEdgeId].e;

  int startK = 4;
  int endK = 4;
  if (k == 4)  // only need to consider 4 side + 3 main, otherwise it has been
               // calculated with 4 main k-1 side
    startK = 3;
  // equals to size1+size2-1<= MAX_PATTERN_SIZE
  //	endK = min(trussClass[sharedEdgeId], ((MAX_PATTERN_SIZE+7)>>1) - k);
  endK = trussClass[sharedEdgeId];
  int tmp = ((MAX_PATTERN_SIZE + 7) >> 1) - k;
  if (tmp < endK) endK = tmp;
  int k1 = k;
  int cnt1 = nodeSet.size();

  if (NNH[sharedEdgeId].size() == 0) getNNH(sharedEdgeId);
  vector<int> &mynnh = NNH[sharedEdgeId];

  int diff = 0;
  int common = 0;
  int scanTo = MAX_COMBINED_TRUSS_CLASS + 1;
  bool flag = false;
  for (int k2 = endK; k2 >= startK; --k2) {
    while (k2 < scanTo && !flag) {
      --scanTo;
      int gobein = mynnh[scanTo - 3];
      int goend = mynnh[scanTo + 1 - 3];  //[begin,end)
      for (int go = gobein; go < goend; ++go) {
        int id = mynnh[go - 3];
        if (id != oriS && id != oriE) {
          if (nodeSet.find(id) != nodeSet.end())
            common++;
          else
            diff++;
          if (diff >= k2 - 2) {
            flag = true;
            break;
          }
        }
      }
    }
    //(cnt1 + (common + diff) - common >= k1 - 2 + k2 - 2 && common + diff >= k2
    //- 2)
    if (flag || cnt1 + diff >= k1 + k2 - 4 && common + diff >= k2 - 2) {
      for (int kk = k2; kk >= startK;
           kk--)  // there is a pruning, if it satisfy, less than it must
                  // satisfy.
      {
        // 3 main/side with k + 1 side same to 4 main with k side, processes
        // separately to avoid double count
        // int patId = kk == 3 ? patIdMap[4][k1 - 1][0] : patIdMap[k1][kk][0];
        int patId = patIdMap[k1][kk][0];
        Pattern &pattern = patterns[patId];
        ++pattern.scorePart[Tid];
      }
      return;
    }
  }
}

/**   calculate score for combined truss
 * k1-side-k2-side combined truss .
 * @param k, k of the first pattern
 * @param originalEdgeId, main edge of the first pattern
 * @param edgeId, the shared side edge
 * @param nodeSet, node neighbourhood of the first edge
 * @param sharedS, the src node of the shared edge
 * @param sharedE, the dest node of the shared edge
 * @param Tid, the current thread id
 */
void FastPGGraph::addScoreForComPatternWithComEdge2(int k, int originalEdgeId,
                                                    int edgeId,
                                                    unordered_set<int> &nodeSet,
                                                    int &sharedS, int &sharedE,
                                                    int &Tid) {
  // connected main edge
  int startK = k;
  int endK = k;
  int size = (startK - 2) * 2 + 1;
  if (size * 2 - 1 > MAX_PATTERN_SIZE) return;
  if (NNH[edgeId].size() == 0) getNNH(edgeId);
  vector<int> &mynnh = NNH[edgeId];
  // //equals to size1+size2-1<= MAX_PATTERN_SIZE
  //	endK = min(trussClass[edgeId], ((MAX_PATTERN_SIZE + 7) >> 1) - k);
  endK = trussClass[edgeId];
  int tmp = ((MAX_PATTERN_SIZE + 7) >> 1) - k;
  if (tmp < endK) endK = tmp;

  bool flag1 = false;
  bool flag2 = false;
  for (int k2 = endK; k2 >= startK; --k2) {
    for (int j = k2; j <= MAX_COMBINED_TRUSS_CLASS; j++) {
      int gobein = mynnh[j - 3];
      int goend = mynnh[j + 1 - 3];  //[begin,end)
      for (int go = gobein; go < goend; ++go) {
        int id = mynnh[go - 3];
        int conEdgeId = getEdgebyNode(id, sharedS);
        int thirdEdgeId = getEdgebyNode(id, sharedE);
        if (conEdgeId == originalEdgeId || thirdEdgeId == originalEdgeId)
          continue;

        if (!flag1 && isOK(k, k2, originalEdgeId, conEdgeId, nodeSet, 1)) {
          for (int kk = k2; kk >= startK;
               kk--)  // if it satisfy, all smaller must be ok
          {
            int patId = k <= kk ? patIdMap[k][kk][1] : patIdMap[kk][k][1];
            Pattern &pattern = patterns[patId];
            ++pattern.scorePart[Tid];
          }
          flag1 = true;  // a prunning: if it satisfy, all smaller must be ok
        }

        if (!flag2 && isOK(k, k2, originalEdgeId, thirdEdgeId, nodeSet, 2)) {
          for (int kk = k2; kk >= startK; kk--) {
            int patId = k <= kk ? patIdMap[k][kk][2] : patIdMap[kk][k][2];
            Pattern &pattern = patterns[patId];
            ++pattern.scorePart[Tid];
          }
          flag2 = true;  // a prunning :if it satisfy, all smaller must be ok
        }
        if (flag1 && flag2) return;
      }
      if (flag1 && flag2) return;
    }
  }
}

/** calculate score for combined truss
 * share the side edge
 * @param k1
 * @param k2
 * @param oriEdgeId, main edge of the first pattern
 * @param edgeId, main edge of the second pattern
 * @param nodeSet, node neighborhood of the first main edge
 * @param  para  1:  two main edges share a node ; 2 : others
 * @return
 */
bool FastPGGraph::isOK(int k1, int k2, int oriEdgeId, int edgeId,
                       unordered_set<int> &nodeSet, int para) {
  int newU = edges[edgeId].s;
  int newV = edges[edgeId].e;
  int oriU = edges[oriEdgeId].s;
  int oriV = edges[oriEdgeId].e;

  int cnt1 = nodeSet.size();
  if (nodeSet.find(newU) != nodeSet.end()) cnt1--;
  if (nodeSet.find(newV) != nodeSet.end()) cnt1--;

  // remove common nodes, there are no k1 pattern
  if (para == 1 && cnt1 < k1 - 2) return false;
  if (para == 2 && cnt1 < k1 - 3) return false;

  if (trussClass[edgeId] < 3) return false;
  if (NNH[edgeId].size() == 0) getNNH(edgeId);
  vector<int> &mynnh = NNH[edgeId];

  int diff = 0;
  int common = 0;
  for (int i = k2; i <= MAX_COMBINED_TRUSS_CLASS; i++) {
    int gobein = mynnh[i - 3];
    int goend = mynnh[i + 1 - 3];  //[begin,end)
    for (int go = gobein; go < goend; ++go) {
      int w = mynnh[go - 3];
      if (w != oriU && w != oriV) {
        if (nodeSet.find(w) != nodeSet.end()) {
          if (w != newU && w != newV)
            common++;
          else
            diff++;
        } else
          diff++;
        if (para == 1) {
          //(common >= 1 && (cnt1 - 1) + (common + diff - 1) - (common - 1) >=
          // k1 - 3 + k2 - 3 && common + diff - 1 >= k2 - 3)
          if (common >= 1 && cnt1 + diff + 5 >= k1 + k2 &&
              common + diff + 2 >= k2)
            return true;
        } else {
          //(cnt1 + (common + diff) - (common) >= k1 - 3 + k2 - 3) && (diff +
          // common >= k2 - 3))
          if (cnt1 + diff >= k1 + k2 - 6 && diff + common >= k2 - 3)
            return true;
        }
      }
    }
  }
  return false;
}

/**
        calculate score for combined truss
 * in each thread, add score for the patterns the edges could contribute to
 * @param edgeIds: all the edges to be processed in this thread
 */
void FastPGGraph::mtAddScoreForComPattern(vector<int> &edgeIds, int &Tid) {
  for (int i = 0; i < edgeIds.size(); ++i) {
    int edgeId = edgeIds[i];
    /*if (!(i&0xffff))
            printf("threadId: %d : %d000 have done\n",Tid,i/10000);*/
    vector<int> &mynnh = NNH[edgeId];
    unordered_set<int> nodeSet;
    // the NNH  lg than trussClass[edgeId] must be satisfy
    for (int k = MAX_COMBINED_TRUSS_CLASS; k > trussClass[edgeId]; k--) {
      int gobein = mynnh[k - 3];
      int goend = mynnh[k + 1 - 3];  //[begin,end)
      for (int go = gobein; go < goend; ++go) {
        int w = mynnh[go - 3];
        nodeSet.emplace(w);
      }
    }

    int s = edges[edgeId].s;
    int e = edges[edgeId].e;
    for (int k = trussClass[edgeId]; k >= 4; k--) {
      int gobein = mynnh[k - 3];
      int goend = mynnh[k + 1 - 3];  //[begin,end)
      for (int go = gobein; go < goend; ++go) {
        int w = mynnh[go - 3];
        nodeSet.emplace(w);
      }

      for (auto comNode : nodeSet) {
        int edgeS = getEdgebyNode(s, comNode);
        int edgeE = getEdgebyNode(e, comNode);
        // p1-side-p2-main pattern
        if (trussClass[edgeS] >= 3)
          addScoreForComPatternWithComEdge1(k, edgeId, edgeS, nodeSet, Tid);
        if (trussClass[edgeE] >= 3)
          addScoreForComPatternWithComEdge1(k, edgeId, edgeE, nodeSet, Tid);

        // p1-side-p2-side pattern
        if (trussClass[edgeS] >= 3)
          addScoreForComPatternWithComEdge2(k, edgeId, edgeS, nodeSet, s,
                                            comNode, Tid);
        if (trussClass[edgeE] >= 3)
          addScoreForComPatternWithComEdge2(k, edgeId, edgeE, nodeSet, e,
                                            comNode, Tid);
      }
    }
  }
}

/**
    get the NNH first
 */

void FastPGGraph::preProceE(vector<int> &mybin, int &Tid) {
  for (int i = 0; i < mybin.size(); i++) {
    /*if (!(i & 0x3ffff))
            printf("Tid: %d, %d00000 edges have been pre-processed \n",Tid,
       i/100000);*/
    getNNH(mybin[i]);
  }
}

/*abandoned function*/
void FastPGGraph::preProceEdge() {
  assert(0);
  vector<int> bin[THREAD_NUM];
  int cur = 0;
  for (int i = 0; i < edges.size(); i++) {
    if (isSampledEdge(i) && edgeType[i] == EdgeType::RealTruss) {
      bin[cur++].emplace_back(i);
      if (cur >= THREAD_NUM) cur = 0;
    }
  }
  /*for (int i = 0; i < THREAD_NUM; i++)
          cout << "the bin size of " << i << " is " << bin[i].size() << endl;*/

  boost::thread *thread[THREAD_NUM];
  for (int i = 0; i < THREAD_NUM; i++)
    thread[i] = new boost::thread(boost::bind(preProceE, bin[i], i));

  for (int i = 0; i < THREAD_NUM; i++) {
    thread[i]->join();
    bin[i].clear();
    delete thread[i];
  }
  // vector<vector<int>> tmp;		//Rabbit died dog hum
  // unordered_NNH.swap(tmp);
}

/* the main function to calculate score for truss and combined truss*/

void FastPGGraph::calScore4CombinedTrussPatterns() {
  printf("calScore4CombinedTrussPatterns....\n");
  double startTime = time(nullptr);
  vector<int> bin[THREAD_NUM];
  NNH.resize(NUM_E);

  int cur = 0;
  for (int i = 0; i < edges.size(); i++)  // using muti-threads to get NNH
  {
    if (isSampledEdge(i) /* && trussClass[i] > 3*/)  // >3 means it is in truss
    {
      bin[cur].emplace_back(i);
      cur = (cur + 1) % THREAD_NUM;
    }
  }

  /*for (int i = 0; i < THREAD_NUM; i++)
  {
          cout << "the bin size of " << i << " is " << bin[i].size() << endl;
          //bin[i].shrink_to_fit();
  }*/

  boost::thread *thread[THREAD_NUM];
  for (int i = 0; i < THREAD_NUM; i++)
    thread[i] = new boost::thread(boost::bind(preProceE, bin[i], i));

  for (int i = 0; i < THREAD_NUM; i++) {
    thread[i]->join();
  }

  boost::thread *newThread[THREAD_NUM];  //// using muti-threads to cal score
                                         /// for combined truss
  for (int i = 0; i < THREAD_NUM; i++)
    newThread[i] =
        new boost::thread(boost::bind(mtAddScoreForComPattern, bin[i], i));

  for (int i = 0; i < THREAD_NUM; i++) newThread[i]->join();
  summaryInfo.calScore4CombinedTrussPatternsTime = time(nullptr) - startTime;

  // Rabbit died dog hum
  vector<vector<int>>().swap(NNH);
  for (int i = 0; i < patterns.size();
       i++)  // cal total scores for combinedtruss
  {
    for (int j = 0; j < THREAD_NUM; ++j)
      patterns[i].score += patterns[i].scorePart[j];
  }
}
/**
 * we are using sampling, now 10% , can adjust the percentage
 * @param edgeId
 * @return
 */
const double Pi = 0.1;
bool inline FastPGGraph::isSampledEdge(int edgeId)  // a naive sample way
{
  return true;  // if no sample
  // srand(time(NULL));
  if ((rand() % 100) * 0.01 < Pi)
    return true;
  else
    return false;
  /*if (edgeId % 10 == 0)
              return true;
  return false;*/
}

/*------------------------------------------------------------------------------------------------------------------------------------------
******************************************************************************************************************************************
                             233333333333    cal truss num for edges
2333333333333
******************************************************************************************************************************************
----------------------------------------------------------------------------------------------------------------------------------------*/

//
// cul the support for each edge.   not store the NNH to reduce the mem.
// neighbId of each edge is sorted, so it is fast

// cal the suppport for the given edgeid.
void FastPGGraph::getEdgeneighbor(int &edgeId) {
  int s = edges[edgeId].s;
  int e = edges[edgeId].e;
  int i = 0, j = 0;
  int lens = neibours[s].size();
  int lene = neibours[e].size();
  int mysup = 0;
  while (j < lene && i < lens) {
    if (neibours[s][i].neighbId == neibours[e][j].neighbId) {
      // mynb.emplace_back(neibours[s][i].neighbId);
      ++mysup;
      ++i;
      ++j;
    } else if (neibours[s][i].neighbId < neibours[e][j].neighbId)
      ++i;
    else
      ++j;
  }
  support[edgeId] = mysup;  // cal the support at the same time
  if (mysup == 0) trussClass[edgeId] = 2;
  // else
  //	unordered_NNH[edgeId] = mynb;
  // unordered_NNH[edgeId].shrink_to_fit();
  // boost::mutex::scoped_lock lock(io_mutex);	// must lock:  many thread write
  // it  at the same time
}

/**
 * sort the edges with their support using bin sort
 */
void FastPGGraph::arrangeSortedEdges() {
  cout << "Arrange sorted edges..." << endl;
  for (int i = 0; i < edges.size(); i++) {
    binSize[support[i]]++;
  }

  binStart[0] = 0;
  for (int i = 1; i < binStart.size(); i++) {
    binStart[i] = binStart[i - 1] + binSize[i - 1];
  }

  for (int i = 0; i < edges.size(); ++i) {
    int tsup = support[i];
    ++curLength[tsup];
    sortedPos[i] = binStart[tsup] + curLength[tsup] - 1;
    sortedEdges[sortedPos[i]] = i;
  }
}

/**
 * support(edgeId) -= 1, rearrange its position in the array
 * @param edgeId
 * @param k
 */
void inline FastPGGraph::decSupportForEdge(int edgeId) {
  // exchange order of cur_edge and binStart[ it's Support];
  int oldSupport = support[edgeId];
  int x = binStart[oldSupport];
  int pos = sortedPos[edgeId];
  int temp = sortedEdges[x];
  sortedEdges[x] = sortedEdges[pos];
  sortedEdges[pos] = temp;
  sortedPos[sortedEdges[x]] = x;
  sortedPos[sortedEdges[pos]] = pos;
  ++binStart[oldSupport];
  --support[edgeId];
}

/**
 * remove the edge with id edgeId, delete the support of the connected edges
 * which can form a triangle with it
 * @param edgeid
 *
 * @param k
 */
void FastPGGraph::decSupport(int edgeId, int k) {
  int s = edges[edgeId].s;
  int e = edges[edgeId].e;
  int i = 0, j = 0;
  int lens = neibours[s].size();
  int lene = neibours[e].size();
  int mysup = 0;
  while (j < lene && i < lens) {
    if (neibours[s][i].neighbId == neibours[e][j].neighbId) {
      int edgeId1 = neibours[s][i].edgeId;
      int edgeId2 = neibours[e][j].edgeId;
      // both of the two edges are alive, otherwise the triangle doesn't exists.
      if (edgeAlive[edgeId1] && edgeAlive[edgeId2]) {
        if (support[edgeId1] >= k - 2) decSupportForEdge(edgeId1);
        if (support[edgeId2] >= k - 2) decSupportForEdge(edgeId2);
      }
      ++i;
      ++j;
    } else if (neibours[s][i].neighbId < neibours[e][j].neighbId)
      ++i;
    else
      ++j;
  }

  /*vector<int>& mynb = unordered_NNH[edgeId];
  for(auto w :mynb)
  {
                  int edgeId1 = getEdgebyNode(edges[edgeId].s, w);
                  int edgeId2 = getEdgebyNode(edges[edgeId].e, w);
                  // both of the two edges are alive, otherwise the triangle
  doesn't exists. if (edgeAlive[edgeId1] && edgeAlive[edgeId2])
                  {
                          if (support[edgeId1] >= k - 2)
                                  decSupportForEdge(edgeId1);
                          if (support[edgeId2] >= k - 2)
                                  decSupportForEdge(edgeId2);
                  }
  }
  unordered_NNH[edgeId].swap(vector<int>());*/
}

/**
 * get K truss
 * @param k
 */
void FastPGGraph::getKTruss(int k) {
  // printf("-----------------------------------------\n");
  // cout << k-1 << "-Truss" << "\n";
  int cnt = 0;
  int sortedIdx;
  int goStart = binStart[k - 3];
  for (int i = goStart; i < binStart[k - 2]; ++i) {
    sortedIdx = sortedEdges[i];
    if (edgeAlive[sortedIdx]) {
      edgeAlive[sortedIdx] = false;
      // annotate the edge
      trussClass[sortedIdx] = k - 1;
      edgeType[sortedIdx] = EdgeType::RealTruss;  // mark it
      ++summaryInfo.insideEdgeCnt;
      decSupport(sortedIdx, k);
      ++cnt;
    }
  }
  // cout << "k:" << k-1  << " Cnt:" << cnt << "\n";
  // printf("-----------------------------------------\n");
}

/**
 * used for annotate each edge with trussness (truss number)
 * @param maxK, the largest k
 */
void FastPGGraph::annoEdgeWithTrussClass(int maxK) {
  GoutDegree.resize(NUM_V + 1, 0);
  for (int i = 4; i <= maxK; i++) getKTruss(i);

  printf("mark and cal truss and cal Gout_degree...\n");
  for (int i = 0; i < edges.size(); i++) {
    if (trussClass[i] == -1) {
      edgeType[i] = EdgeType::RealTruss;
      trussClass[i] = maxK;
      ++summaryInfo.insideEdgeCnt;
    }

    for (int j = 3; j <= trussClass[i];
         j++)  // cal truss.  less than 3 is meaningless
      ++trussEdgeCnt[j];

    if (edgeType[i] == EdgeType::Default)  // cal Goutdegree at the same time;
    {
      ++GoutDegree[edges[i].s];
      ++GoutDegree[edges[i].e];
      ++summaryInfo.outsideEdgeCnt;
    }
  }
}

/*get getEdgeneighbor using muti-threads */
void FastPGGraph::preProceE0(vector<int> &mybin, int &Tid) {
  for (int i = 0; i < mybin.size(); i++) {
    /*if (!(i & 0x1ffff))
            printf("Tid: %d, %d00000 edges have been pre-processed \n", Tid, i /
       100000);*/
    getEdgeneighbor(mybin[i]);
  }
}

/*cal trussNumber For all Edges */
void FastPGGraph::calTrussClass4Edges() {
  double startTime = time(nullptr);
  edgeType.resize(NUM_E, EdgeType::Default);
  support.resize(NUM_E, 0);
  trussClass.resize(NUM_E, -1);  //
  // unordered_NNH.resize(NUM_E);
  vector<int> bin[THREAD_NUM];
  int cur = 0;
  for (int i = 0; i < edges.size(); i++) {
    bin[cur++].emplace_back(i);
    if (cur >= THREAD_NUM) cur = 0;
  }

  /*	for (int i = 0; i < THREAD_NUM; i++)
          {
                  cout << "the bin size of " << i << " is " << bin[i].size() <<
     endl;
                  //bin[i].shrink_to_fit();
          }*/

  boost::thread *thread[THREAD_NUM];
  for (int i = 0; i < THREAD_NUM; i++)
    thread[i] = new boost::thread(boost::bind(preProceE0, bin[i], i));

  for (int i = 0; i < THREAD_NUM; i++) thread[i]->join();
  // now start to use
  sortedEdges.resize(NUM_E, 0);
  sortedPos.resize(NUM_E);
  edgeAlive.resize(NUM_E, true);
  binStart.resize(MAX_SUPPORT, 0);
  binSize.resize(MAX_SUPPORT, 0);
  curLength.resize(MAX_SUPPORT, 0);

  arrangeSortedEdges();
  annoEdgeWithTrussClass(MAX_COMBINED_TRUSS_CLASS);

  vector<int>().swap(binSize);
  vector<int>().swap(sortedEdges);
  vector<int>().swap(binStart);
  vector<int>().swap(curLength);
  vector<int>().swap(support);
  vector<int>().swap(sortedPos);
  vector<bool>().swap(edgeAlive);
  summaryInfo.calTrussClass4EdgesTime = time(nullptr) - startTime;
}

/*------------------------------------------------------------------------------------------------------------------------------------------
******************************************************************************************************************************************
233333333333   in G_out 2333333333333
******************************************************************************************************************************************
----------------------------------------------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------------------------------------------------
******************************************************************************************************************************************
233333333333    star patterns in G_out 2333333333333
******************************************************************************************************************************************
----------------------------------------------------------------------------------------------------------------------------------------*/

/*distribution sort in Gout by their Goutdegree   and  count the star
/*
   search stars and combined stars , and cal their scores

   abandand
*/
void FastPGGraph::searchStars() {
  double startTime = time(nullptr);
  printf("searching Stars patterns... \n");
  // for (auto & nodepair : nodesMap)
  for (int i = 0; i < NUM_V; ++i) {
    if (GoutDegree[i] >= MIN_STAR_DEGREE) {
      // travel my neibours
      // ******************** collect pattern ********************
      auto pat = outsideStarPatMap.find(GoutDegree[i]);
      if (pat == outsideStarPatMap.end()) {
        Pattern pattern = Pattern(GoutDegree[i] + 1);
        for (int j = 0; j < GoutDegree[i]; ++j) pattern.addEdge(0, j + 1);
        pattern.setType(PatType::OutsideStar);
        pattern.setScore(0);
        pat =
            outsideStarPatMap.emplace(make_pair(GoutDegree[i], patterns.size()))
                .first;
        pattern.detilType = make_pair(GoutDegree[i], -1);
        patterns.emplace_back(pattern);
      }
      ++patterns[pat->second].score;
      // ******************** collect pattern ********************
      ++summaryInfo.starsComponentCnt;
      // for (auto & mynb : node.neighbours)
      for (int j = 0; j < neibours[i].size(); ++j) {
        auto &type = edgeType[neibours[i][j].edgeId];
        if (type == EdgeType::Default) {
          type = EdgeType::SimpleStar;
          ++summaryInfo.starsEdgeCnt;
        } else if (type == EdgeType::SimpleStar) {
          type = EdgeType::CombinedStar;
          // ******************** collect pattern ********************
          int d1 = GoutDegree[i];
          int d2 = GoutDegree[neibours[i][j].neighbId];
          if (d1 > d2) swap(d1, d2);
          auto pat = outsideCombinedStarPatMap.find(make_pair(d1, d2));
          if (pat == outsideCombinedStarPatMap.end()) {
            Pattern pattern = Pattern(d1 + d2);
            pattern.addEdge(0, 1);
            for (int k = 0; k < d1 - 1; ++k) pattern.addEdge(0, 2 + k);

            for (int k = 0; k < d2 - 1; ++k) pattern.addEdge(1, d1 + 1 + k);

            pattern.setType(PatType::OutsideCombinedStar);
            pattern.setScore(0);
            pat = outsideCombinedStarPatMap
                      .emplace(make_pair(make_pair(d1, d2), patterns.size()))
                      .first;
            pattern.detilType = make_pair(d1, d2);
            pattern.degree.emplace_back(d1);
            pattern.degree.emplace_back(d2);
            patterns.emplace_back(pattern);
          }
          ++patterns[pat->second].score;
          // ******************** collect pattern ********************
        }
      }
    }
  }
  summaryInfo.searchStarsTime = time(nullptr) - startTime;
}

/*
search stars and combined stars , and cal their scores
linear stars vesion
*/
// used for hash vector to int
size_t vec2ull(const std::vector<short> &v) {
  std::hash<short> hasher;
  size_t seed = 0;
  for (auto i : v) {
    seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

size_t vec2ull_reverse(const std::vector<short> &v) {
  std::hash<short> hasher;
  size_t seed = 0;
  for (int i = v.size() - 1; i >= 0; --i) {
    seed ^= hasher(v[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}

void FastPGGraph::searchMultiLinearStars() {
  const int MAXLINEARSIZE = 6;  // level
  printf("searching Multi linear Stars patterns... \n");
  double startTime = time(nullptr);
  for (int i = 0; i < NUM_V; ++i) {
    if (GoutDegree[i] >= MIN_STAR_DEGREE) {
      // ******************** collect  single stars pattern ********************
      auto pat = outsideStarPatMap.find(GoutDegree[i]);
      if (pat == outsideStarPatMap.end()) {
        Pattern pattern = Pattern(GoutDegree[i] + 1);
        for (int j = 0; j < GoutDegree[i]; ++j) pattern.addEdge(0, j + 1);
        pattern.setType(PatType::OutsideStar);
        pattern.setScore(0);
        pat =
            outsideStarPatMap.emplace(make_pair(GoutDegree[i], patterns.size()))
                .first;
        pattern.detilType = make_pair(GoutDegree[i], -1);
        patterns.emplace_back(pattern);
      }
      ++patterns[pat->second].score;
      ++summaryInfo.starsComponentCnt;

      queue<StarState> q;
      StarState startS;
      startS.nodeId.push_back(i);
      startS.degree.push_back(GoutDegree[i]);
      startS.myv.emplace(i);
      startS.cursize = GoutDegree[i];
      q.push(startS);
      while (!q.empty()) {
        StarState cur = q.front();
        q.pop();
        int curv = cur.nodeId.back();
        for (int j = 0; j < neibours[i].size(); j++) {
          auto &type = edgeType[neibours[i][j].edgeId];
          if (type != EdgeType::Default && type != EdgeType::SimpleStar)
            continue;
          if (type == EdgeType::Default) {
            type = EdgeType::SimpleStar;
            ++summaryInfo.starsEdgeCnt;
          }

          int nextv = neibours[i][j].neighbId;
          if (cur.myv.find(nextv) == cur.myv.end() &&
              GoutDegree[nextv] >= MIN_STAR_DEGREE &&
              cur.nodeId.size() + 1 <= MAXLINEARSIZE &&
              cur.cursize + GoutDegree[nextv] - 1 <= MAX_PATTERN_SIZE) {
            StarState next = cur;
            next.myv.emplace(nextv);
            next.degree.push_back(GoutDegree[nextv]);
            next.nodeId.push_back(nextv);
            next.cursize += GoutDegree[nextv] - 1;
            q.push(next);

            /*     ******************** collect linear combined pattern
             * ********************      */
            size_t tmpId = vec2ull(next.degree);  //映射两次
            auto pat = outsideLinearCombinedStarPatMap.find(tmpId);
            if (pat == outsideLinearCombinedStarPatMap.end()) {
              int tnumv = 0;
              for (int k = 0; k < next.nodeId.size(); ++k)
                tnumv += next.degree[k];

              tnumv = tnumv - next.degree.size() + 2;  // 顶点数
              Pattern pattern = Pattern(tnumv);
              int tsums = next.degree.size();
              for (int k = 0; k < next.nodeId.size(); ++k) {
                if (k < next.nodeId.size() - 1)
                  pattern.addEdge(k, k + 1);  // main edge

                int tmp = next.degree[k] - 2;
                if (k == 0 || k == next.degree.size() - 1) ++tmp;

                for (int kk = 0; kk < tmp; ++kk)  // side edge
                  pattern.addEdge(k, tsums + kk);
                tsums += tmp;

                pattern.degree.emplace_back(
                    next.degree[k]);  // for record pattens degree in main node
              }

              pattern.setType(PatType::OutsideCombinedStar);
              pattern.setScore(0);
              outsideLinearCombinedStarPatMap.emplace(
                  make_pair(vec2ull_reverse(next.degree), patterns.size()));
              pat = outsideLinearCombinedStarPatMap
                        .emplace(make_pair(tmpId, patterns.size()))
                        .first;
              pattern.detilType =
                  make_pair(-1, -1);  // dose not make sence any more
              patterns.emplace_back(pattern);
            }
            ++patterns[pat->second].score;
          }
        }
      }
    }
  }
  /*for (auto p : patterns)
  {
          if (p.getType() == PatType::OutsideCombinedStar)
          {
                  outputPngPattern(p, "data/loc-brightkite/finalPatterns/");
          }
  }*/
  summaryInfo.searchStarsTime = time(nullptr) - startTime;
}

/**
 * search and mark lines, circles and individual components in G_out ,and cal
 * score for them
 */
void FastPGGraph::searchLinesAndCirclesAndIndividuals() {
  double startTime = time(nullptr);
  printf("searching Lines, Circles and Individuals patterns... \n");

  vector<bool> flag(NUM_V, false);
  queue<int> q;
  int nodeCnt, edgeCnt, d2Cnt, d1Cnt, degree;
  vector<int> edgeSet;
  circleRecord.clear();
  vector<int> circleRecordTmp;

  for (int i = 0; i < NUM_V; ++i) {
    if (GoutDegree[i] > 0 && !flag[i])  // means Gout
    {
      d2Cnt = d1Cnt = nodeCnt = edgeCnt = 0;  // reset counters
      edgeSet.clear();
      circleRecordTmp.clear();
      q.push(i);
      circleRecordTmp.emplace_back(i);
      flag[i] = true;
      ++nodeCnt;
      while (!q.empty())  // using bfs
      {
        int node = q.front();
        q.pop();
        degree = 0;
        for (auto ve : neibours[node]) {
          if (edgeType[ve.edgeId] == EdgeType::Default) {
            ++edgeCnt;
            ++degree;
            edgeSet.emplace_back(ve.edgeId);
            if (!flag[ve.neighbId]) {
              q.push(ve.neighbId);
              circleRecordTmp.emplace_back(ve.neighbId);
              flag[ve.neighbId] = true;
              ++nodeCnt;
            }
          }
        }
        if (degree == 2)
          ++d2Cnt;
        else if (degree == 1)
          ++d1Cnt;
      }
      assert((edgeCnt & 1) == 0);
      edgeCnt >>= 1;                       // every edge is counted twice
      if (nodeCnt == 4 && edgeCnt == 4) {  // 4-circle
        assert(nodeCnt == circleRecordTmp.size());
        circleRecord.emplace_back(circleRecordTmp);
      }
      if (d1Cnt == 2 && d2Cnt == nodeCnt - 2 && nodeCnt != 2 &&
          nodeCnt != 3)  // 1/2 edge pattern is default
      {                  // line
        // ******************** collect pattern ********************
        auto pat = outsideLinePatMap.find(edgeCnt);
        if (pat == outsideLinePatMap.end()) {
          Pattern pattern = Pattern(nodeCnt);
          for (int j = 0; j < edgeCnt; ++j) pattern.addEdge(j, j + 1);
          pattern.setType(PatType::OutsideLine);
          pattern.setScore(0);
          pat = outsideLinePatMap.emplace(make_pair(edgeCnt, patterns.size()))
                    .first;
          pattern.detilType = make_pair(nodeCnt - 1, -1);
          patterns.emplace_back(pattern);
        }
        ++patterns[pat->second].score;
        // ******************** collect pattern ********************
        ++summaryInfo.linesComponentCnt;
        summaryInfo.linesEdgeCnt += edgeCnt;
        for (auto edgeId : edgeSet) {
          edgeType[edgeId] = EdgeType::Line;
          //
        }
      } else if (d1Cnt == 0 && d2Cnt == nodeCnt && nodeCnt != 4 &&
                 nodeCnt != 3)  // trangle and  rectangle are default patterns
      {                         // circle
        // ******************** collect pattern ********************
        auto pat = outsideCirclePatMap.find(edgeCnt);
        if (pat == outsideCirclePatMap.end()) {
          Pattern pattern = Pattern(nodeCnt);
          for (int j = 0; j < nodeCnt - 1; ++j) pattern.addEdge(j, j + 1);
          pattern.addEdge(nodeCnt - 1, 0);
          pattern.setType(PatType::OutsideCircle);
          pattern.setScore(0);
          pat = outsideCirclePatMap.emplace(make_pair(edgeCnt, patterns.size()))
                    .first;
          pattern.detilType = make_pair(nodeCnt, -1);
          patterns.emplace_back(pattern);
        }
        ++patterns[pat->second].score;
        // ******************** collect pattern ********************
        ++summaryInfo.circlesComponentCnt;
        summaryInfo.circlesEdgeCnt += edgeCnt;
        for (auto edgeId : edgeSet) edgeType[edgeId] = EdgeType::Circle;

        assert(nodeCnt == circleRecordTmp.size());
        circleRecord.emplace_back(circleRecordTmp);
      } else if (nodeCnt > 3 && nodeCnt <= MAX_INDIVIDUAL_NODE_COUNT &&
                 edgeCnt != 4)  // individual
      {
        ++summaryInfo.individualsComponentCnt;
        summaryInfo.individualsEdgeCnt += edgeCnt;
        int c = 0;
        unordered_map<int, int> nodeMap;
        graph_t g(nodeCnt);
        for (auto edgeId : edgeSet) {
          edgeType[edgeId] = EdgeType::Individual;
          int s = edges[edgeId].s;
          auto itr = nodeMap.find(s);
          if (itr == nodeMap.end())
            itr = nodeMap.insert(make_pair(s, c++)).first;
          s = itr->second;
          int e = edges[edgeId].e;
          itr = nodeMap.find(e);
          if (itr == nodeMap.end())
            itr = nodeMap.emplace(make_pair(e, c++)).first;
          e = itr->second;
          if (s != e) add_edge(s, e, g);
        }
        assert(c == nodeCnt);
        bool found = false;
        std::vector<graph_t::vertex_descriptor> f(
            nodeCnt);  // judge using boost library
        for (auto &smallComponent : smallComponents[nodeCnt][edgeCnt])
          if (isomorphism(
                  g, smallComponent.graph,
                  boost::isomorphism_map(boost::make_iterator_property_map(
                      f.begin(), boost::get(boost::vertex_index, g))))) {
            found = true;
            ++smallComponent.score;
            ++patterns[smallComponent.patId].score;
            break;
          }
        if (!found) {
          Pattern pattern = Pattern(nodeCnt);
          for (auto ep = boost::edges(g); ep.first != ep.second; ++ep.first)
            pattern.addEdge(boost::source(*ep.first, g),
                            boost::target(*ep.first, g));
          smallComponents[nodeCnt][edgeCnt].emplace_back(
              SmallComponentInfo(g, 1, patterns.size()));
          pattern.setScore(1);
          pattern.setType(PatType::OutsideIndividual);
          pattern.detilType = make_pair(-1, -1);
          patterns.emplace_back(pattern);
        }
      } else if (nodeCnt > MAX_INDIVIDUAL_NODE_COUNT) {  // uncovered
        ++summaryInfo.uncoveredComponentCnt;
        summaryInfo.uncoveredEdgeCnt += edgeCnt;
      }
    }
  }
  summaryInfo.searchLinesAndCirclesAndIndividualsTime =
      time(nullptr) - startTime;
}

/*******************************************Pattern
 * Choose*****************************************************/

/* filterPatterns which is too big(our gui can not hold big patterns), or its
 * ferenquency is too low */

int FastPGGraph::filterPatterns(const Pattern &p) {
  // const int MAX_V2RECOGNISE = 15;

  if (p.score <= 2) return 1;
  if (p.getEdgeCnt() >= MAX_PATTERN_SIZE) return 2;
  return 0;
}

/*cal Coverage And Cognition  scores*/
void FastPGGraph::calCoverageAndCognition(const string &folder) {
  // filter large patterns which are useless.

  cout << "candidata patterns size: before filter:" << patterns.size() << "\n";
  double startTime = time(nullptr);
  vector<Pattern> filteredP;
  int cnt1 = 0, cnt2 = 0;
  for (const auto &p : patterns) {
    int isf = filterPatterns(p);
    if (isf == 0)
      filteredP.emplace_back(p);
    else if (isf == 1)
      cnt1++;
    else if (isf == 2)
      cnt2++;
  }
  // patterns = filteredP;
  patterns.swap(filteredP);
  cout << "removed by low frequence: " << cnt1 << "\n";
  cout << "removed by big size: " << cnt2 << "\n";
  cout << "candidata patterns size:" << patterns.size() << "\n";
  boostPatterns.resize(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)  // now used new id.
  {
    // patterns[j].id = j;
    boostPatterns[j] = pattern2boostGraph(patterns[j]);
  }

  switch (COV_FUNC) {
    case 1:
      calCoverage1(folder);
      break;
    case 2:
      calCoverage2(folder);  // final
      break;
    case 3:
      calCoverage3(folder);
      break;
    default:
      assert(0);
  }
  switch (COG_FUNC) {
    case 1:
      calCognition1(folder);
      break;
    case 2:
      calCognition2(folder);
      break;
    case 3:
      calCognition3(folder);
      break;
    case 4:
      calCognition4(folder);
      break;
    case 5:
      calCognition5(folder);
      break;
    case 6:
      calCognition6(folder);  // final
      break;
    default:
      assert(0);
  }

  // coverage of G_out patterns are too large, shrink them proportionally
  double prop_in = summaryInfo.insideEdgeCnt * 1.0 /
                   (summaryInfo.insideEdgeCnt + summaryInfo.outsideEdgeCnt);
  double prop_out = summaryInfo.outsideEdgeCnt * 1.0 /
                    (summaryInfo.insideEdgeCnt + summaryInfo.outsideEdgeCnt);
  cout << "prop_in=" << prop_in << " prop_out=" << prop_out << endl;
  for (auto &p : patterns)
    p.scoreCoverage *= (p.type <= PatType::TrussCombined2 ? prop_in : prop_out);

  summaryInfo.calCovAndCogScoreTime = time(nullptr) - startTime;
}

/*boost library use*/
graph_t FastPGGraph::pattern2boostGraph(const Pattern &p) {
  graph_t g(p.nodesCnt);
  for (auto edge : p.edges) add_edge(edge.s, edge.e, g);
  return g;
}

template <typename Graph>
struct mcgregor_vertex_callback {
  mcgregor_vertex_callback(const Graph &graph1, const Graph &graph2,
                           vector<int> &c)
      : graph1_(graph1), graph2_(graph2), common_vertices(c) {}

  template <typename CorrespondenceMapFirstToSecond,
            typename CorrespondenceMapSecondToFirst>
  bool operator()(
      CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
      CorrespondenceMapSecondToFirst correspondence_map_2_to_1,
      typename boost::graph_traits<Graph>::vertices_size_type subgraph_size) {
    if (common_vertices.empty()) {
      // Print out correspondences between vertices
      BGL_FORALL_VERTICES_T(vertex1, graph1_, Graph) {
        // Skip unmapped vertices
        if (get(correspondence_map_1_to_2, vertex1) !=
            boost::graph_traits<Graph>::null_vertex()) {
          //                        std::cout << vertex1 << " <-> " <<
          //                        get(correspondence_map_1_to_2, vertex1) <<
          //                        std::endl;
          common_vertices.emplace_back(vertex1);
          common_vertices.emplace_back(get(correspondence_map_1_to_2, vertex1));
        }
      }
      //            std::cout << "---" << std::endl;
    }
    return (false);
  }

 private:
  const Graph &graph1_;
  const Graph &graph2_;
  vector<int> &common_vertices;
};

template <typename Graph>
struct vf2_bool_callback {
  vf2_bool_callback(const Graph &graph1, const Graph &graph2)
      : graph1_(graph1), graph2_(graph2) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
    // Print (sub)graph isomorphism map
    //        BGL_FORALL_VERTICES_T(v, graph1_, Graph)std::cout << '(' <<
    //        get(boost::vertex_index_t(), graph1_, v) << ", "
    //                                                          <<
    //                                                          get(boost::vertex_index_t(),
    //                                                          graph2_, get(f,
    //                                                          v)) << ") ";
    //
    //        std::cout << std::endl;

    return false;  // end
  }

 private:
  const Graph &graph1_;
  const Graph &graph2_;
};

template <typename Graph>
struct vf2_vertex_callback {
  vf2_vertex_callback(const Graph &graph1, const Graph &graph2,
                      vector<vector<int>> &m)
      : graph1_(graph1), graph2_(graph2), matches(m) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1) const {
    if (matches.size() > 500) return false;

    // Print (sub)graph isomorphism map
    vector<int> m((int)boost::num_vertices(graph1_));
    BGL_FORALL_VERTICES_T(v, graph1_, Graph)
    m[get(boost::vertex_index_t(), graph1_, v)] =
        get(boost::vertex_index_t(), graph2_, get(f, v));
    matches.emplace_back(m);

    return true;  // not end, need all matching schemes
  }

 private:
  const Graph &graph1_;
  const Graph &graph2_;
  vector<vector<int>> &matches;
};

/**
 * maximum common subgraph of pattern p1 and pattern p2
 * @param p1 first pattern index
 * @param p2 second pattern index
 * @return {vertex count, edge count}
  abandaned
 */

vector<int> FastPGGraph::preJudge4MaxComSubGraph(Pattern &pattern1,
                                                 Pattern &pattern2) {
  assert(0);
  // special judge to speed up.
  int type1 = pattern1.type, type2 = pattern2.type;
  if (type1 == PatType::Truss && type2 == PatType::Truss ||
      type1 == PatType::OutsideStar && type2 == PatType::OutsideStar ||
      type1 == PatType::OutsideLine && type2 == PatType::OutsideLine ||
      type1 == PatType::OutsideCircle && type2 == PatType::OutsideCircle)
    return pattern1.nodesCnt <= pattern2.nodesCnt
               ? vector<int>({1, pattern1.nodesCnt, pattern1.getEdgeCnt()})
               : vector<int>({1, pattern2.nodesCnt, pattern2.getEdgeCnt()});

  if (type1 == PatType::Truss &&
      (type2 == PatType::TrussCombined1 || type2 == PatType::TrussCombined2)) {
    int tmpCom1 = min(pattern1.detilType.first, pattern2.detilType.first);
    int tmpCom2 = min(pattern1.detilType.first, pattern2.detilType.second);
    int maxCom = max(tmpCom1, tmpCom2);
    return vector<int>({1, maxCom, 2 * max(tmpCom1, tmpCom2) - 3});
  }

  if ((type1 == PatType::TrussCombined2 || type1 == PatType::TrussCombined1) &&
      (type2 == PatType::TrussCombined2 || type2 == PatType::TrussCombined1)) {
    int tmpCom1 = min(pattern1.detilType.first, pattern2.detilType.first) +
                  min(pattern1.detilType.second, pattern2.detilType.second);
    int tmpCom2 = min(pattern1.detilType.first, pattern2.detilType.second) +
                  min(pattern1.detilType.second, pattern2.detilType.first);
    int maxCom = max(tmpCom1, tmpCom2);
    return vector<int>({1, maxCom - 2,
                        2 * max(tmpCom1, tmpCom2) - 6 -
                            1});  // -1 means shared edge ,-2 mean share;
  }

  if (type1 == PatType::OutsideCombinedStar &&
      type2 == PatType::OutsideCombinedStar) {
    int tmpCom1 = min(pattern1.detilType.first, pattern2.detilType.first) +
                  min(pattern1.detilType.second, pattern2.detilType.second);
    int tmpCom2 = min(pattern1.detilType.first, pattern2.detilType.second) +
                  min(pattern1.detilType.second, pattern2.detilType.first);
    int maxCom = max(tmpCom1, tmpCom2);
    return vector<int>({1, maxCom, maxCom - 1});
  }

  if (type1 == PatType::OutsideStar && type2 == PatType::OutsideCombinedStar) {
    int tmpCom1 = min(pattern1.detilType.first, pattern2.detilType.first);
    int tmpCom2 = min(pattern1.detilType.first, pattern2.detilType.second);
    int maxCom = max(tmpCom1, tmpCom2);
    return vector<int>({1, maxCom + 1, maxCom});
  }

  if (type1 == PatType::OutsideStar && type2 == PatType::Truss) {
    int maxCom = min(pattern1.detilType.first, pattern2.detilType.first - 1);
    return vector<int>({1, maxCom + 1, maxCom});
  }

  if (type1 == PatType::OutsideStar && type2 == PatType::TrussCombined1) {
    int maxCom = min(pattern1.detilType.first,
                     pattern2.detilType.first + pattern2.detilType.second - 3);
    return vector<int>({1, maxCom + 1, maxCom});
  }

  if (type1 == PatType::OutsideCombinedStar &&
      type2 == PatType::TrussCombined1) {
    int tmpCom1 = min(pattern1.detilType.first, pattern2.detilType.second + 1) +
                  min(pattern1.detilType.second,
                      pattern2.detilType.first + pattern2.detilType.second - 3);
    int tmpCom2 =
        min(pattern1.detilType.first,
            pattern2.detilType.first + pattern2.detilType.second - 3) +
        min(pattern1.detilType.second, pattern2.detilType.second + 1);
    int maxCom = max(tmpCom1, tmpCom2);
    return vector<int>({1, maxCom, maxCom - 1});
  }

  /*stars/combined with others  */
  vector<int> tmpD(pattern2.nodesCnt, 0);
  for (auto edge : pattern2.edges) ++tmpD[edge.e], ++tmpD[edge.s];
  int maxD = 0, maxDoubleD = 0;
  if (type1 == PatType::OutsideStar) {
    for (int v : tmpD)
      if (v > maxD) maxD = v;
    int maxCom = min(maxD, pattern1.detilType.first);
    return vector<int>({1, maxCom + 1, maxCom});
  }

  if (type1 == PatType::OutsideCombinedStar) {
    int maxCom = 0;
    for (auto edge : pattern2.edges) {
      int tmpCom1 = min(pattern1.detilType.first, tmpD[edge.s]) +
                    min(pattern1.detilType.second, tmpD[edge.e]);
      int tmpCom2 = min(pattern1.detilType.first, tmpD[edge.e]) +
                    min(pattern1.detilType.second, tmpD[edge.s]);
      maxCom = max(maxCom, max(tmpCom1, tmpCom2));
    }
    return vector<int>({1, maxCom, maxCom - 1});
  }
  return vector<int>({0, 0, 0});
}

/**
 * maximum common subgraph
 * @param p1 id of pattern1
 * @param p2 id of pattern2
 * @return {common_node_cnt, common_edge_cnt}
  abandant
 */
vector<int> FastPGGraph::patternCommonSubgraph(int p1, int p2) {
  // some prejudge

  assert(0);
  vector<int> preJudgeResults =
      preJudge4MaxComSubGraph(patterns[p1], patterns[p2]);
  if (preJudgeResults[0])
    return vector<int>({preJudgeResults[1], preJudgeResults[2]});

  vector<int> preJudgeResultsReversed =
      preJudge4MaxComSubGraph(patterns[p2], patterns[p1]);
  if (preJudgeResultsReversed[0])
    return vector<int>(
        {preJudgeResultsReversed[1], preJudgeResultsReversed[2]});

  cout << "Their type are : " << patterns[p1].type
       << " size: " << patterns[p1].nodesCnt << "  --   " << patterns[p2].type
       << " size: " << patterns[p2].nodesCnt << "\n";
  graph_t &g1 = boostPatterns[p1];
  graph_t &g2 = boostPatterns[p2];
  vector<int> common_vertices;
  mcgregor_vertex_callback<graph_t> my_callback(g1, g2, common_vertices);
  mcgregor_common_subgraphs_maximum_unique(g1, g2, true, my_callback);
  vector<int> match(patterns[p1].getNodesCnt(), -1);
  int edgeCnt = 0;
  for (int i = 0; i < common_vertices.size(); i += 2)
    match[common_vertices[i]] = common_vertices[i + 1];
  for (auto edge : patterns[p1].edges)
    if (match[edge.s] != -1 && match[edge.e] != -1 &&
        patterns[p2].getEdgeId(match[edge.s], match[edge.e]) != -1)
      ++edgeCnt;
  return vector<int>({int(common_vertices.size() >> 1), edgeCnt});
}

bool isCombinedStarsSubGraph(const Pattern &p1, const Pattern &p2) {
  int l1 = p1.degree.size();
  int l2 = p2.degree.size();
  bool found;
  for (int j = 0; j <= l2 - l1; ++j) {
    found = true;
    for (int i = 0; i < l1; ++i)
      if (p1.degree[i] > p2.degree[j + i]) {
        found = false;
        break;
      }
    if (found) return true;
    found = true;
    for (int i = 0; i < l1; ++i)
      if (p1.degree[i] > p2.degree[l2 - j - i]) {
        found = false;
        break;
      }
    if (found) return true;
  }
  return false;
}

/**
 * judge if p1 is a subgraph of p2
 * @param p1 id of pattern1 (smaller)
 * @param p2 id of pattern2 (larger)
 * @return bool
 */
bool FastPGGraph::isSubGraph(int p1, int p2) {
  if (patterns[p1].getNodesCnt() > patterns[p2].getNodesCnt() ||
      patterns[p1].getEdgeCnt() > patterns[p2].getEdgeCnt())
    return false;
  if (patterns[p1].type == PatType::OutsideCombinedStar &&
      patterns[p2].type == PatType::OutsideCombinedStar)
    return isCombinedStarsSubGraph(patterns[p1], patterns[p2]);

  graph_t &g1 = boostPatterns[p1];
  graph_t &g2 = boostPatterns[p2];
  vf2_bool_callback<graph_t> callback(g1, g2);
  return vf2_subgraph_mono(g1, g2, callback);
}

/**
 * judge if p1 is a subgraph of p2
 * @param p1  (smaller)
 * @param p2  (larger)
 * @return bool
 */
bool FastPGGraph::isSubGraph(const Pattern &p1, const Pattern &p2) {
  if (p1.getNodesCnt() > p2.getNodesCnt() || p1.getEdgeCnt() > p2.getEdgeCnt())
    return false;
  if (p1.type == PatType::OutsideCombinedStar &&
      p2.type == PatType::OutsideCombinedStar)
    return isCombinedStarsSubGraph(p1, p2);
  graph_t g1 = pattern2boostGraph(p1);
  graph_t g2 = pattern2boostGraph(p2);
  vf2_bool_callback<graph_t> callback(g1, g2);
  return vf2_subgraph_mono(g1, g2, callback);
}

/**
 * all subgraph matching schemes on pattern2
 * @param id of pattern1 (smaller)
 * @param id of pattern2 (larger)
 * @return unordered_set<unsigned long long> unique edge sets on pattern2
 */
unordered_set<unsigned long long> FastPGGraph::getAllSubGraphEdgeSets(int p1,
                                                                      int p2) {
  // cout <<"doing "<< patterns[p1].type <<" size:
  // "<<patterns[p1].getEdgeCnt()<<"---" << patterns[p2].type << " size: " <<
  // patterns[p2].getEdgeCnt()<<"\n";
  return getAllSubGraphEdgeSets(patterns[p1], boostPatterns[p1], patterns[p2],
                                boostPatterns[p2]);
}

/**
 * all subgraph matching schemes on pattern2
 * @param id of pattern1 (smaller)
 * @param pattern2 (larger)
 * @return unordered_set<unsigned long long> unique edge sets on pattern2
 */
unordered_set<unsigned long long> FastPGGraph::getAllSubGraphEdgeSets(
    int p1, const Pattern &pattern2) {
  return getAllSubGraphEdgeSets(patterns[p1], boostPatterns[p1], pattern2,
                                pattern2boostGraph(pattern2));
}

/**
 * pick leftC items from a
 * @param a edge list [(e, edge_id)]
 * @param i current index of a
 * @param leftA left items number of a
 * @param leftC number of items to be choosed
 * @param exclude id of forbidden edge
 * @param cur current edge set
 * @param ans edge set
 */
void enumerateEdgeSets(const vector<pair<int, int>> &a, int i, int leftA,
                       int leftC, const int &exclude_edge, ull cur,
                       vector<ull> &ans) {
  if (leftC == 0) {
    ans.emplace_back(cur);
    return;
  }
  // choose a[i]
  if (a[i].second != exclude_edge)
    enumerateEdgeSets(a, i + 1, leftA - 1, leftC - 1, exclude_edge,
                      cur | (1ULL << a[i].second), ans);
  // not choose a[i]
  if (leftA > leftC)
    enumerateEdgeSets(a, i + 1, leftA - 1, leftC, exclude_edge, cur, ans);
}

// with 2 exclude edges
void enumerateEdgeSets(const vector<pair<int, int>> &a, int i, int leftA,
                       int leftC, const int &exclude_edge1,
                       const int &exclude_edge2, ull cur, vector<ull> &ans) {
  if (leftC == 0) {
    ans.emplace_back(cur);
    return;
  }
  // choose a[i]
  if (a[i].second != exclude_edge1 && a[i].second != exclude_edge2)
    enumerateEdgeSets(a, i + 1, leftA - 1, leftC - 1, exclude_edge1,
                      exclude_edge2, cur | (1ULL << a[i].second), ans);
  // not choose a[i]
  if (leftA > leftC)
    enumerateEdgeSets(a, i + 1, leftA - 1, leftC, exclude_edge1, exclude_edge2,
                      cur, ans);
}

/**
 * enumerate all possible edge sets of a chain(multi-combined star)
 * @param ith node
 * @param node number
 * @param current edge set
 * @param all possible edge sets of each node
 * @param ans edge set
 */
void enumerateEdgeSets(int d, int max_d, ull tmp,
                       const vector<vector<ull>> &route_edge_sets,
                       unordered_set<unsigned long long> &ans) {
  if (d > max_d) {
    ans.emplace(tmp);
    return;
  }
  int cnt = route_edge_sets[d].size();
  if (cnt > 5) {
    for (int i = 0; i < 5; ++i)
      enumerateEdgeSets(d + 1, max_d, tmp | route_edge_sets[d][rand() % cnt],
                        route_edge_sets, ans);
  } else
    for (auto &a : route_edge_sets[d])
      enumerateEdgeSets(d + 1, max_d, tmp | a, route_edge_sets, ans);
}

void traverseCombinedStar(const Pattern &p, int cur_node, int pre_edge,
                          ull nodes_set, ull edges_set, int cur_degree,
                          const vector<short> &degree,
                          vector<vector<ull>> &route_edge_sets,
                          unordered_set<unsigned long long> &edgeSets) {
  /**
   * a     -- b   -- c   -- d
   * front -- mid -- mid -- back
   */
  if (cur_degree == degree.size() - 1) {
    route_edge_sets[cur_degree].clear();
    // back node: exclude 1 edge
    enumerateEdgeSets(p.edgesList[cur_node], 0, p.edgesList[cur_node].size(),
                      degree[cur_degree] - 1, pre_edge, 0,
                      route_edge_sets[cur_degree]);
    enumerateEdgeSets(0, cur_degree, edges_set, route_edge_sets, edgeSets);
    return;
  }
  for (auto &j : p.edgesList[cur_node]) {  // select next node
    int next_node = j.first;
    if (((1ULL << next_node) & nodes_set) == 0 &&
        p.edgesList[next_node].size() >= degree[cur_degree + 1]) {
      route_edge_sets[cur_degree].clear();
      if (pre_edge == -1)
        // front node: exclude 1 edge
        enumerateEdgeSets(p.edgesList[cur_node], 0,
                          p.edgesList[cur_node].size(), degree[cur_degree] - 1,
                          j.second, 0, route_edge_sets[cur_degree]);
      else
        // mid node: exclude 2 edges
        enumerateEdgeSets(p.edgesList[cur_node], 0,
                          p.edgesList[cur_node].size(), degree[cur_degree] - 2,
                          j.second, pre_edge, 0, route_edge_sets[cur_degree]);
      traverseCombinedStar(p,
                           next_node,                        // cur_node
                           j.second,                         // pre_edge
                           (1ULL << next_node) | nodes_set,  // nodes_set
                           (1ULL << j.second) | edges_set,   // edges_set
                           cur_degree + 1,                   // cur_degree
                           degree, route_edge_sets, edgeSets);
    }
  }
}

unordered_set<unsigned long long> FastPGGraph::getAllSubGraphEdgeSets(
    const Pattern &pattern1, const graph_t &graph1, const Pattern &pattern2,
    const graph_t &graph2) {
  unordered_set<unsigned long long> edgeSets;
  if (pattern1.edges.size() > pattern2.edges.size()) return edgeSets;
  /*printf("p1.type=%d, p2.type=%d\n", pattern1.type, pattern2.type);
  printf("p1.n=%d, p1.e=%d\n", pattern1.nodesCnt, pattern1.edges.size());
  for (auto edge : pattern1.edges)
          printf("%d,%d\n", edge.s, edge.e);
  printf("p2.n=%d, p2.e=%d\n", pattern2.nodesCnt, pattern2.edges.size());
  for (auto edge : pattern2.edges)
          printf("%d,%d\n", edge.s, edge.e);*/

  if (pattern1.type == PatType::OutsideStar) {
    int d1 = pattern1.detilType.first;
    for (int i = 0; i < pattern2.getNodesCnt(); ++i)
      if (pattern2.edgesList[i].size() >= d1) {
        vector<ull> tmp;
        enumerateEdgeSets(pattern2.edgesList[i], 0,
                          pattern2.edgesList[i].size(), d1, -1, 0, tmp);
        for (auto edgeSet : tmp) edgeSets.emplace(edgeSet);
      }
    return edgeSets;
  } else if (pattern1.type == PatType::OutsideCombinedStar) {
    vector<vector<ull>> route_edge_sets(pattern1.degree.size());
    for (int i = 0; i < pattern2.getNodesCnt(); ++i)
      if (pattern2.edgesList[i].size() >= pattern1.degree[0])
        traverseCombinedStar(pattern2, i, -1, 1ULL << i, 0, 0, pattern1.degree,
                             route_edge_sets, edgeSets);
    return edgeSets;
  }

  vector<vector<int>> matches;
  vf2_vertex_callback<graph_t> callback(graph1, graph2, matches);
  vf2_subgraph_mono(graph1, graph2, callback);
  // printf("matches size=%d\n", matches.size());
  for (vector<int> &m : matches) {
    unsigned long long edgeSet = 0;  // edge set
    for (auto edge : pattern1.edges)
      edgeSet |= 1ULL << pattern2.getEdgeId(m[edge.s], m[edge.e]);
    edgeSets.emplace(edgeSet);
  }
  return edgeSets;
}

/* generate the final patterns,  size is the number of patterns that user
 * required */
void FastPGGraph::patternChooseByRedundancy(int size, string filename) {
  printf("patternChooseByRedundancy....\n");
  genDefaultPatterns();
  for (int i = patterns.size() - 4; i < patterns.size(); i++)  // calCognition3
    calCognition3Single(patterns[i]);
  boostPatterns.resize(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)
    boostPatterns[j] = pattern2boostGraph(patterns[j]);

  for (int i = patterns.size() - 4; i < patterns.size(); i++)
    finalPatternIdSet.emplace_back(i);  //最后放入的四个默认pattern
  vector<int> flag(patterns.size() + 1, 0);
  for (int idx : finalPatternIdSet) flag[idx] = 1;

  int curMinStep = 3;
  for (int i = 4; i < size; i++) {
    // cout << "choosing no." << i + 1 << " pattern..." << endl;
    double maxps = -1;
    int maxp = -1;

    // *************************** multi-thread implementation
    // ***************************
    vector<double> maxpsVec(THREAD_NUM, -1);
    vector<int> maxpVec(THREAD_NUM, -1);
    for (int j = 0; j < patterns.size(); ++j)
      if (!flag[j]) candidatePatternsQueue.emplace(j);
    boost::thread *thread[THREAD_NUM];
    for (int j = 0; j < THREAD_NUM; ++j)
      thread[j] = new boost::thread(calCandidatePatternsFinalScore, curMinStep,
                                    &maxpsVec[j], &maxpVec[j], j);

    // cout << "muti-thread task finished\n";
    for (int j = 0; j < THREAD_NUM; ++j) {
      thread[j]->join();
      delete thread[j];
      if (maxpsVec[j] > maxps) {
        maxps = maxpsVec[j];
        maxp = maxpVec[j];
      }
    }
    // *************************** multi-thread implementation
    // ***************************

    // cout <<"com e :"<< patterns[maxp].historyMaxComE << "\n";
    if (maxp == -1) {
      curMinStep--;
      --i;
      // printf("curMinStep decreased to %d\n", curMinStep);
      if (curMinStep <= 0) {
        break;
      }
      continue;
    }

    if (patterns[maxp].getEdgeCnt() > 4)  // penaty similar patterns
    {
      int mark = 0;
      for (int fid : finalPatternIdSet) {
        // vector<int>com = patternCommonSubgraph(maxp, fid);
        // int comEdgeSize = com[1];
        // if (patterns[maxp].getNodesCnt() == 4 && patterns[maxp].getEdgeCnt()
        // == 3 && patterns[maxp].type == 10) { 	cout <<
        // (patterns[maxp].type ==
        // PatType::OutsideLine) << endl;
        //}
        if (patterns[maxp].getEdgeCnt() <= MAX_PATTERN_SIZE / 2 + 1 ||
            patterns[maxp].type == PatType::OutsideStar ||
            patterns[maxp].type == PatType::OutsideCombinedStar ||
            patterns[maxp].type == PatType::OutsideLine ||
            patterns[maxp].type == PatType::OutsideCircle) {
          if (isSubGraph(fid, maxp) && patterns[fid].getEdgeCnt() >=
                                           patterns[maxp].getEdgeCnt() - 1 ||
              isSubGraph(maxp, fid) && patterns[maxp].getEdgeCnt() >=
                                           patterns[fid].getEdgeCnt() - 1) {
            mark = 1;
            break;
          }
        } else {
          if (isSubGraph(fid, maxp) && patterns[fid].getEdgeCnt() >=
                                           patterns[maxp].getEdgeCnt() - 2 ||
              isSubGraph(maxp, fid) && patterns[maxp].getEdgeCnt() >=
                                           patterns[fid].getEdgeCnt() - 2) {
            mark = 1;
            break;
          }
        }
      }
      if (mark) {
        flag[maxp] = 1;
        --i;
        // printf("too similar, choose again\n");
        continue;
      }
    }

    // printf("pattern %d is choosed\n", maxp);
    // cout << "choose " << patterns[maxp].getNodesCnt() << " " <<
    // patterns[maxp].getEdgeCnt() << endl; cout << "type " <<
    // patterns[maxp].type << endl; for (auto e : patterns[maxp].edges)
    // cout
    // <<
    // e.s << "-" << e.e << endl;
    flag[maxp] = 1;
    finalPatternIdSet.emplace_back(maxp);
  }

  /*int dis = size - finalPatternIdSet.size();
          for (int j = 0; j < patterns.size(); ++j)
                  if (!flag[j])
                  {
                          if (dis <= 0)
                                  break;
                          finalPatternIdSet.push_back(j);
                  }
  */

  /*when less than using finalscore sort*/
  flag.clear();
  flag.resize(patterns.size(), 0);
  for (int idx : finalPatternIdSet) flag[idx] = 1;
  int cursize = finalPatternIdSet.size();
  int id = 0;
  for (int i = cursize; i < size;) {
    if (flag[id] == 0) {
      i++;
      finalPatternIdSet.push_back(id);
      flag[id] = 1;
    }
    id++;
    if (id >= patterns.size()) break;
  }
  /*when less than using finalscore sort*/

  printf("finalPatternIdSet size: %lu\n", finalPatternIdSet.size());

  // if you want to output the final patterns for view,
  // for output final patterns for view!!!
  for (int i = 0; i < finalPatternIdSet.size(); i++) {
    patterns[finalPatternIdSet[i]].finalRank = i;
    //		outputPngPattern(patterns[finalPatternIdSet[i]], filename);
    outputFinalPngPattern(
        patterns[finalPatternIdSet[i]], filename,
        "#" + to_string(i + 1) + "-" + to_string(finalPatternIdSet[i]) + "=" +
            to_string(patterns[finalPatternIdSet[i]].type) + "=" +
            to_string(patterns[finalPatternIdSet[i]].scoreCoverage) + "-" +
            to_string(patterns[finalPatternIdSet[i]].scoreCognition));

    string typeClass = "others";
    switch (patterns[finalPatternIdSet[i]].type) {
      case 0:
      case 1:
      case 2:
        typeClass = "default";
        break;
      case 3:
      case 4:
      case 5:
      case 6:
        typeClass = "chord-type";
        break;
      case 8:
      case 9:
        typeClass = "star-type";
        break;
      case 10:
      case 11:
      case 12:
        typeClass = "remain-type";
        break;
    }
    outputFinalPngPattern(
        patterns[finalPatternIdSet[i]], filename + typeClass + "/",
        "#" + to_string(i + 1) + "-" + to_string(finalPatternIdSet[i]) + "=" +
            to_string(patterns[finalPatternIdSet[i]].type) + "=" +
            to_string(patterns[finalPatternIdSet[i]].scoreCoverage) + "-" +
            to_string(patterns[finalPatternIdSet[i]].scoreCognition) + "-" +
            to_string(patterns[finalPatternIdSet[i]].scoreDiversity));
  }
}

void FastPGGraph::patternChooseByGreedy(int size, string filename) {
  printf("patternChooseByGreedy....\n");
  boostPatterns.resize(patterns.size());
  vector<vector<double>> netSimileSignature(patterns.size());
  for (int j = 0; j < patterns.size(); ++j) {
    boostPatterns[j] = pattern2boostGraph(patterns[j]);
    netSimileSignature[j] = net_simile::Signature(patterns[j]);
    /*cout << "id=" << patterns[j].id << " size=" <<
    netSimileSignature[j].size() << endl; for (auto a : netSimileSignature[j])
            cout << a << " ";
    cout << endl;*/
  }

  vector<int> flag(patterns.size() + 1, 0);
  for (int idx : finalPatternIdSet) flag[idx] = 1;

  for (int i = 4; i < size; i++) {
    // cout << "choosing no." << i + 1 << " pattern..." << endl;
    double maxps = -1e8;
    int maxp = -1;

    // *************************** multi-thread implementation
    // ***************************
    vector<double> maxpsVec(THREAD_NUM, -1);
    vector<int> maxpVec(THREAD_NUM, -1);
    for (int j = 0; j < patterns.size(); ++j)
      if (!flag[j]) {
        auto cov = patterns[j].scoreCoverage;
        auto cog = patterns[j].scoreCognition;
        double sim = 0.0;
        for (int k = 0; k < finalPatternIdSet.size(); ++k) {
          double tmp = net_simile::CamberraDistance(netSimileSignature[j],
                                                    netSimileSignature[k]);
          if (tmp > sim) sim = tmp;
        }
        // cout << "sim=" << sim << endl;
        auto delta_s = cov - sim - cog;
        if (delta_s > maxps) {
          maxps = delta_s;
          maxp = j;
        }
      }
    if (maxp == -1) break;

    // printf("pattern %d is choosed\n", maxp);
    // cout << "choose " << patterns[maxp].getNodesCnt() << " " <<
    // patterns[maxp].getEdgeCnt() << endl; cout << "type " <<
    // patterns[maxp].type << endl; for (auto e : patterns[maxp].edges)
    // cout
    // <<
    // e.s << "-" << e.e << endl;
    flag[maxp] = 1;
    finalPatternIdSet.emplace_back(maxp);
  }

  /*int dis = size - finalPatternIdSet.size();
          for (int j = 0; j < patterns.size(); ++j)
                  if (!flag[j])
                  {
                          if (dis <= 0)
                                  break;
                          finalPatternIdSet.push_back(j);
                  }
  */

  /*when less than using finalscore sort*/
  flag.clear();
  flag.resize(patterns.size(), 0);
  for (int idx : finalPatternIdSet) flag[idx] = 1;
  int cursize = finalPatternIdSet.size();
  int id = 0;
  for (int i = cursize; i < size;) {
    if (flag[id] == 0) {
      i++;
      finalPatternIdSet.push_back(id);
      flag[id] = 1;
    }
    id++;
    if (id >= patterns.size()) break;
  }
  /*when less than using finalscore sort*/

  printf("finalPatternIdSet.size: %lu\n", finalPatternIdSet.size());
  for (int i = 0; i < finalPatternIdSet.size(); ++i) {
    patterns[finalPatternIdSet[i]].finalRank = i;
    double sim = 0.0;
    for (int j = 0; j < finalPatternIdSet.size(); ++j)
      if (i != j) {
        double tmp = net_simile::CamberraDistance(
            netSimileSignature[finalPatternIdSet[i]],
            netSimileSignature[finalPatternIdSet[j]]);
        if (tmp > sim) sim = tmp;
      }
    patterns[finalPatternIdSet[i]].scoreDiversity = 1.0 / sim;
  }
}

void FastPGGraph::selectFinalPatterns(string filename, double lambda, int size,
                                      int datasetid) {
  double startTime = time(nullptr);
  // for (auto &pattern : patterns)
  // 	pattern.finalScore = lambda * pattern.scoreCoverage + (1 - lambda) *
  // pattern.scoreCognition; for (auto &pattern : patterns)
  // pattern.finalScore = pattern.scoreCoverage * pattern.scoreCognition;
  for (auto &pattern : patterns)
    pattern.finalScore = pattern.scoreCoverage / pattern.scoreCognition;
  sort(patterns.begin(), patterns.end(), [](Pattern a, Pattern b) -> bool {
    return a.finalScore > b.finalScore;
  });
  // patternChooseByRedundancy(size, filename);
  genDefaultPatterns();  // add 4 patterns
  for (int i = patterns.size() - 4; i < patterns.size(); i++) {
    patterns[i].scoreCoverage = 1.0;
    calCognition6Single(patterns[i]);
  }
  finalPatternIdSet.clear();
  for (int i = patterns.size() - 4; i < patterns.size(); i++)
    finalPatternIdSet.emplace_back(i);  //最后放入的四个默认pattern
  patternChooseByGreedy(size, filename);
  summaryInfo.finalPatternsGenTime = time(nullptr) - startTime;

  // if you want to output the final patterns for view,
  // for output final patterns for view!!!
  if (!filename.empty()) {
    for (int i = 0; i < finalPatternIdSet.size(); i++) {
      //		outputPngPattern(patterns[finalPatternIdSet[i]],
      // filename);
      outputFinalPngPattern(
          patterns[finalPatternIdSet[i]], filename,
          "#" + to_string(i + 1) + "-" + to_string(finalPatternIdSet[i]) + "=" +
              to_string(patterns[finalPatternIdSet[i]].type) + "=" +
              to_string(patterns[finalPatternIdSet[i]].scoreCoverage) + "-" +
              to_string(patterns[finalPatternIdSet[i]].scoreCognition) + "-" +
              to_string(patterns[finalPatternIdSet[i]].scoreDiversity));

      string typeClass = "others";
      switch (patterns[finalPatternIdSet[i]].type) {
        case 0:
        case 1:
        case 2:
          typeClass = "default";
          break;
        case 3:
        case 4:
        case 5:
        case 6:
          typeClass = "chord-type";
          break;
        case 8:
        case 9:
          typeClass = "star-type";
          break;
        case 10:
        case 11:
        case 12:
          typeClass = "remain-type";
          break;
      }
      outputFinalPngPattern(
          patterns[finalPatternIdSet[i]], filename + typeClass + "/",
          "#" + to_string(i + 1) + "-" + to_string(finalPatternIdSet[i]) + "=" +
              to_string(patterns[finalPatternIdSet[i]].type) + "=" +
              to_string(patterns[finalPatternIdSet[i]].scoreCoverage) + "-" +
              to_string(patterns[finalPatternIdSet[i]].scoreCognition) + "-" +
              to_string(patterns[finalPatternIdSet[i]].scoreDiversity));
    }
  }

  for (int i = 0; i < finalPatternIdSet.size(); i++) {
    finalSet4Datastores[datasetid].emplace_back(patterns[finalPatternIdSet[i]]);
  }

  // outputPatterns(filename + "patterns.txt");
}

void FastPGGraph::selectFinalPatternsByCov(int size) {
  sort(patterns.begin(), patterns.end(), [](Pattern a, Pattern b) -> bool {
    return a.scoreCoverage > b.scoreCoverage;
  });
  genDefaultPatterns();
  boostPatterns.resize(patterns.size());
  finalPatternIdSet.clear();
  for (int i = patterns.size() - 4; i < patterns.size(); i++)
    finalPatternIdSet.emplace_back(i);  //最后放入的四个默认pattern
  for (int i = 0; i < size - 4; ++i) finalPatternIdSet.emplace_back(i);
}

/**
 * 1 groups
 * @author Zifeng Yuan
 */
void FastPGGraph::calCoverage1(const string &folder) {
  printf("evaluating patterns coverage scores... \n");

  // Two patterns A and B, if A is a subgraph of B, A's frequency is definitely
  // larger than B's frequency score == frequency for patterns of G_out
  vector<int> outside;
  for (int i = 0; i < patterns.size(); ++i)
    if (patterns[i].type > PatType::TrussCombined2) outside.emplace_back(i);

  int len = outside.size();
  vector<int> extraFrequency(len, 0);
  for (int i = 0; i < len; ++i)
    for (int j = i + 1; j < len; ++j)
      if (!(patterns[outside[i]].type == PatType::OutsideStar &&
            patterns[outside[j]].type == PatType::OutsideCombinedStar) &&
          !(patterns[outside[i]].type == PatType::OutsideCombinedStar &&
            patterns[outside[j]].type == PatType::OutsideStar)) {
        if (isSubGraph(outside[i], outside[j])) {
          // cout <<"v: "<<patterns[outside[i]].nodesCnt << "--- e:
          // "<<patterns[outside[i]].getEdgeCnt() << "-----vs----- "; cout << "
          // v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() << "\n";
          extraFrequency[i] += patterns[outside[j]].score;
        } else if (isSubGraph(outside[j], outside[i])) {
          // cout << "  v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() <<"-----vs----- " ; cout << "v: "
          // << patterns[outside[i]].nodesCnt << "--- e: " <<
          // patterns[outside[i]].getEdgeCnt() << "\n";
          extraFrequency[j] += patterns[outside[i]].score;
        }
      }

  for (int i = 0; i < len; ++i) patterns[outside[i]].score += extraFrequency[i];

  // score == frequency for patterns of G_out
  // coverage = frequency * size
  for (auto &pattern : patterns)
    // pattern.scoreCoverage = pattern.score;
    pattern.scoreCoverage = pattern.type <= PatType::TrussCombined2
                                ? pattern.score
                                : pattern.score * pattern.getEdgeCnt();

  // normalization
  vector<int> group(20);
  group[PatType::Truss] = group[PatType::TrussCombined1] =
      group[PatType::TrussCombined2] = 0;
  group[PatType::OutsideStar] = group[PatType::OutsideCombinedStar] = 0;
  group[PatType::OutsideLine] = group[PatType::OutsideCircle] =
      group[PatType::OutsideIndividual] = 0;

  vector<double> covmin(3, 1e18), covmax(3, 0);
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    if (pattern.scoreCoverage > covmax[g]) covmax[g] = pattern.scoreCoverage;
    if (pattern.scoreCoverage < covmin[g]) covmin[g] = pattern.scoreCoverage;
  }
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    pattern.scoreCoverage =
        (pattern.scoreCoverage - covmin[g] + 1) / (covmax[g] - covmin[g] + 1);
  }
}

void FastPGGraph::calCognition1(const string &folder) {
  double cogmin = 1e18, cogmax = 0;
  for (auto &pattern : patterns) {
    // cog cal
    double tmp = 0;
    vector<int> d(pattern.getNodesCnt() + 1, 0);
    int maxd = 0;
    for (PatEdge patEdge : pattern.edges) {
      ++d[patEdge.e];
      ++d[patEdge.s];
      if (d[patEdge.e] > maxd) maxd = d[patEdge.e];
      if (d[patEdge.s] > maxd) maxd = d[patEdge.s];
    }
    vector<int> numD(maxd + 1, 0);
    for (int j = 0; j <= pattern.getNodesCnt(); ++j) ++numD[d[j]];

    for (int j = 1; j <= maxd; ++j) {
      tmp += j * int(numD[j] > 0);
    }
    pattern.scoreCognition = -log(pattern.getNodesCnt() + tmp);
    if (pattern.scoreCognition > cogmax) cogmax = pattern.scoreCognition;
    if (pattern.scoreCognition < cogmin) cogmin = pattern.scoreCognition;
  }
  for (auto &pattern : patterns)
    pattern.scoreCognition =
        (pattern.scoreCognition - cogmin + 1) / (cogmax - cogmin + 1);
}

/**
 * 3 groups: (truss, combined_truss) (outside_star, outside_combined_star)
 * (others) or 4 groups: (truss) (combined_truss) (outside_star,
 * outside_combined_star) (others) or 5 groups: (truss) (combined_truss)
 * (outside_star) (outside_combined_star) (others)
 * @author Zifeng Yuan
 */
void FastPGGraph::calCoverage2(const string &folder) {
  printf("evaluating patterns coverage scores... \n");

  // Two patterns A and B, if A is a subgraph of B, A's frequency is definitely
  // larger than B's frequency score == frequency for patterns of G_out
  vector<int> outside;
  for (int i = 0; i < patterns.size(); ++i)
    if (patterns[i].type > PatType::TrussCombined2) outside.emplace_back(i);

  int len = outside.size();
  vector<int> extraFrequency(len, 0);
  for (int i = 0; i < len; ++i)
    for (int j = i + 1; j < len; ++j)
      if (!(patterns[outside[i]].type == PatType::OutsideStar &&
            patterns[outside[j]].type == PatType::OutsideCombinedStar) &&
          !(patterns[outside[i]].type == PatType::OutsideCombinedStar &&
            patterns[outside[j]].type == PatType::OutsideStar)) {
        if (isSubGraph(outside[i], outside[j])) {
          // cout <<"v: "<<patterns[outside[i]].nodesCnt << "--- e:
          // "<<patterns[outside[i]].getEdgeCnt() << "-----vs----- "; cout << "
          // v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() << "\n";
          extraFrequency[i] += patterns[outside[j]].score;
        } else if (isSubGraph(outside[j], outside[i])) {
          // cout << "  v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() <<"-----vs----- " ; cout << "v: "
          // << patterns[outside[i]].nodesCnt << "--- e: " <<
          // patterns[outside[i]].getEdgeCnt() << "\n";
          extraFrequency[j] += patterns[outside[i]].score;
        }
      }

  for (int i = 0; i < len; ++i) patterns[outside[i]].score += extraFrequency[i];

  // score == frequency for patterns of G_out
  // coverage = frequency * size
  for (auto &pattern : patterns)
    // pattern.scoreCoverage = pattern.score;
    pattern.scoreCoverage = pattern.type <= PatType::TrussCombined2
                                ? pattern.score
                                : pattern.score * pattern.getEdgeCnt();

  // normalization
  vector<int> group(20);
  group[PatType::Truss] = 0;
  group[PatType::TrussCombined1] = group[PatType::TrussCombined2] = 1;
  group[PatType::OutsideStar] = 2;
  group[PatType::OutsideCombinedStar] = 3;
  group[PatType::OutsideLine] = group[PatType::OutsideCircle] =
      group[PatType::OutsideIndividual] = 4;

  vector<double> covmin(5, 1e18), covmax(5, 0);
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    if (pattern.scoreCoverage > covmax[g]) covmax[g] = pattern.scoreCoverage;
    if (pattern.scoreCoverage < covmin[g]) covmin[g] = pattern.scoreCoverage;
  }
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    pattern.scoreCoverage =
        (pattern.scoreCoverage - covmin[g] + 1) / (covmax[g] - covmin[g] + 1);
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCoverage > patterns[b].scoreCoverage;
  });
  for (int i = 0; i < 50; ++i)
    outputFinalPngPattern(patterns[pr[i]], folder + "byCov2/",
                          "#" + to_string(i + 1) + "=" +
                              to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                              to_string(group[patterns[pr[i]].type]) + "=" +
                              to_string(patterns[pr[i]].scoreCoverage));
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::Truss) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov2/truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::TrussCombined1 ||
        patterns[pr[i]].type == PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov2/combined_truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type > PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov2/gout/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
    // ******************************* for debug only
    // *******************************
#endif
}

double cognitionF(double x) { return (exp(x) / (exp(x) + 1) - 0.5) * 2; }

void FastPGGraph::calCognition2Single(Pattern &pattern) {
  printf("evaluating single pattern cognition score... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density;
  int i = 0;
  density = (2.0 * pattern.getEdgeCnt()) /
            (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
  degree.resize(pattern.getNodesCnt(), 0);
  fill(degree.begin(), degree.end(), 0);
  for (auto &edge : pattern.edges) {
    ++degree[edge.s];
    ++degree[edge.e];
  }
  d1 = d2 = d3 = d4 = 0;
  for (auto d : degree) switch (d) {
      case 1:
        ++d1;
        break;
      case 2:
        ++d2;
        break;
      case 3:
        ++d3;
        break;
      default:
        ++d4;
    }

  double crossingPenalty = 0;
  if (pattern.getNodesCnt() >= 3) {
    if (pattern.getType() <= PatType::TrussCombined2 &&
        pattern.getEdgeCnt() > 3 * pattern.getNodesCnt() - 6)
      crossingPenalty = 0.8;
    if (pattern.getType() > PatType::TrussCombined2 &&
        pattern.getEdgeCnt() > 2 * pattern.getNodesCnt() - 4)
      crossingPenalty = 0.8;
  }

  if (pattern.getNodesCnt() <= 5 || pattern.getNodesCnt() <= 6 &&
                                        d3 + d4 <= 2 && d4 <= 1 &&
                                        density <= 0.40) {  // very simple
    pattern.scoreCognition = 1;
  } else if (pattern.getEdgeCnt() <= 9 ||
             pattern.type == PatType::OutsideStar ||
             pattern.type == PatType::OutsideCombinedStar ||
             pattern.type == PatType::OutsideLine ||
             pattern.type == PatType::OutsideCircle) {  // normal

    pattern.scoreCognition =
        1 - cognitionF((pattern.getNodesCnt() - 6) * 0.25 + d3 * 0.25 +
                       d4 * 0.5 + crossingPenalty);
  } else {
    pattern.scoreCognition =
        1 - cognitionF((2 * pattern.getNodesCnt() - 6) * 0.25 + d3 * 0.25 +
                       d4 * 0.5 + crossingPenalty);
  }
}

/**
 * @author Zifeng Yuan
 */
void FastPGGraph::calCognition2(const string &folder) {
  printf("evaluating patterns cognition scores... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density;
  int i = 0;
  for (auto &pattern : patterns) {
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    degree.resize(pattern.getNodesCnt(), 0);
    fill(degree.begin(), degree.end(), 0);
    for (auto &edge : pattern.edges) {
      ++degree[edge.s];
      ++degree[edge.e];
    }
    d1 = d2 = d3 = d4 = 0;
    for (auto d : degree) switch (d) {
        case 1:
          ++d1;
          break;
        case 2:
          ++d2;
          break;
        case 3:
          ++d3;
          break;
        default:
          ++d4;
      }

    double crossingPenalty = 0;
    if (pattern.getNodesCnt() >= 3) {
      if (pattern.getType() <= PatType::TrussCombined2 &&
          pattern.getEdgeCnt() > 3 * pattern.getNodesCnt() - 6)
        crossingPenalty = 0.8;
      if (pattern.getType() > PatType::TrussCombined2 &&
          pattern.getEdgeCnt() > 2 * pattern.getNodesCnt() - 4)
        crossingPenalty = 0.8;
    }

    if (pattern.getNodesCnt() <= 5 || pattern.getNodesCnt() <= 6 &&
                                          d3 + d4 <= 2 && d4 <= 1 &&
                                          density <= 0.40) {  // very simple
      pattern.scoreCognition = 1;
#ifdef DEBUG_MODE
      outputFinalPngPattern(
          pattern, folder + "byCog2/very-simple/",
          to_string(pattern.getNodesCnt()) + "=" + to_string(d1) + "-" +
              to_string(d2) + "-" + to_string(d3) + "-" + to_string(d4) + "=" +
              to_string(density) + "=" + to_string(pattern.scoreCognition));
#endif
    } else if (pattern.getEdgeCnt() <= 9 ||
               pattern.type == PatType::OutsideStar ||
               pattern.type == PatType::OutsideCombinedStar ||
               pattern.type == PatType::OutsideLine ||
               pattern.type == PatType::OutsideCircle) {  // normal

      pattern.scoreCognition =
          1 - cognitionF((pattern.getNodesCnt() - 6) * 0.25 + d3 * 0.25 +
                         d4 * 0.5 + crossingPenalty);
#ifdef DEBUG_MODE
      outputFinalPngPattern(
          pattern, folder + "byCog2/normal/",
          to_string(pattern.getNodesCnt()) + "=" + to_string(d1) + "-" +
              to_string(d2) + "-" + to_string(d3) + "-" + to_string(d4) + "=" +
              to_string(density) + "=" + to_string(pattern.scoreCognition));
#endif
    } else {
      pattern.scoreCognition =
          1 - cognitionF((2 * pattern.getNodesCnt() - 6) * 0.25 + d3 * 0.25 +
                         d4 * 0.5 + crossingPenalty);
#ifdef DEBUG_MODE
      outputFinalPngPattern(
          pattern, folder + "byCog2/sophisticated/",
          to_string(pattern.getNodesCnt()) + "=" + to_string(d1) + "-" +
              to_string(d2) + "-" + to_string(d3) + "-" + to_string(d4) + "=" +
              to_string(density) + "=" + to_string(pattern.scoreCognition));
#endif
    }

#ifdef DEBUG_MODE
    if (crossingPenalty > 0)
      outputFinalPngPattern(pattern, folder + "byCog2/crossing/",
                            to_string(pattern.getNodesCnt()) + "-" +
                                to_string(pattern.getEdgeCnt()) + "=" +
                                to_string(pattern.scoreCognition));
#endif
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCognition > patterns[b].scoreCognition;
  });
  for (int i = 0; i < patterns.size(); ++i) {
    outputFinalPngPattern(patterns[pr[i]], folder + "byCog2/",
                          "#" + to_string(i + 1) + "=" +
                              to_string(patterns[pr[i]].scoreCognition));
    outputFinalPngPattern(
        patterns[pr[i]],
        folder + "byCog2/" + to_string(patterns[pr[i]].type) + "/",
        "#" + to_string(i + 1) + "=" +
            to_string(patterns[pr[i]].scoreCognition));
  }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::Truss) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog2/truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::TrussCombined1 ||
        patterns[pr[i]].type == PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog2/combined_truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type > PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog2/gout/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
    // ******************************* for debug only
    // *******************************
#endif
}

void FastPGGraph::calCognition3Single(Pattern &pattern) {
  // printf("evaluating single pattern cognition score... \n");

  double density;
  density = (2.0 * pattern.getEdgeCnt()) /
            (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
  pattern.scoreCognition = density * pattern.getEdgeCnt();
  // pattern.scoreCognition = 1 - (pattern.scoreCognition - g_cogmin + 1) /
  // (g_cogmax - g_cogmin + 1);
}

/**
 * CATAPULT cognitive load
 * @author Zifeng Yuan
 */
void FastPGGraph::calCognition3(const string &folder) {
  printf("evaluating patterns cognition scores3... \n");

  double density;
  // double cogmin = 1e20, cogmax = -1e20;
  for (auto &pattern : patterns) {
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    pattern.scoreCognition = density * pattern.getEdgeCnt();
    // if (pattern.scoreCognition < cogmin)
    // 	cogmin = pattern.scoreCognition;
    // if (pattern.scoreCognition > cogmax)
    // 	cogmax = pattern.scoreCognition;
  }
  // g_cogmin = cogmin;
  // g_cogmax = cogmax;
  // for (auto &pattern : patterns)
  // 	pattern.scoreCognition = 1 - (pattern.scoreCognition - cogmin + 1) /
  // (cogmax - cogmin + 1);

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCognition > patterns[b].scoreCognition;
  });
  for (int i = 0; i < patterns.size(); ++i) {
    outputFinalPngPattern(patterns[pr[i]], folder + "byCog3/",
                          "#" + to_string(i + 1) + "=" +
                              to_string(patterns[pr[i]].scoreCognition));
    outputFinalPngPattern(
        patterns[pr[i]],
        folder + "byCog3/" + to_string(patterns[pr[i]].type) + "/",
        "#" + to_string(i + 1) + "=" +
            to_string(patterns[pr[i]].scoreCognition));
  }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::Truss) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog3/truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::TrussCombined1 ||
        patterns[pr[i]].type == PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog3/combined_truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type > PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog3/gout/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(patterns[pr[i]].scoreCognition));
    }
    // ******************************* for debug only
    // *******************************
#endif
}

/**
 * alternative CATAPULT cognitive load
 * @author Zifeng Yuan
 */
void FastPGGraph::calCognition4(const string &folder) {
  printf("evaluating patterns cognition scores... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density;
  int i = 0;
  for (auto &pattern : patterns) {
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    degree.resize(pattern.getNodesCnt(), 0);
    fill(degree.begin(), degree.end(), 0);
    for (auto &edge : pattern.edges) {
      ++degree[edge.s];
      ++degree[edge.e];
    }
    pattern.scoreCognition = 0;
    for (auto d : degree) pattern.scoreCognition += d - 1;
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCognition > patterns[b].scoreCognition;
  });
  for (int i = 0; i < patterns.size(); ++i) {
    outputFinalPngPattern(patterns[pr[i]], folder + "byCog4/",
                          "#" + to_string(i + 1) + "=" +
                              to_string(patterns[pr[i]].scoreCognition));
    outputFinalPngPattern(
        patterns[pr[i]],
        folder + "byCog4/" + to_string(patterns[pr[i]].type) + "/",
        "#" + to_string(i + 1) + "=" +
            to_string(patterns[pr[i]].scoreCognition));
  }
  // ******************************* for debug only
  // *******************************
#endif
}

void FastPGGraph::calCognition4Single(Pattern &pattern) {
  // printf("evaluating single pattern cognition score4... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density;
  int i = 0;

  density = (2.0 * pattern.getEdgeCnt()) /
            (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
  degree.resize(pattern.getNodesCnt(), 0);
  fill(degree.begin(), degree.end(), 0);
  for (auto &edge : pattern.edges) {
    ++degree[edge.s];
    ++degree[edge.e];
  }
  pattern.scoreCognition = 0;
  for (auto d : degree) pattern.scoreCognition += d - 1;
}

bool isPlanar(Pattern &pattern) {
  // ref: https://en.wikipedia.org/wiki/Planar_graph#Other_planarity_criteria
  auto v = pattern.getNodesCnt();
  auto e = pattern.getEdgeCnt();
  if (v < 3) return true;
  if (e > 3 * v - 6) return false;
  auto type = pattern.getType();
  // no cycles of length 3
  if (type != PatType::Triangle && type != PatType::Truss &&
      type != PatType::TrussCombined1 && type != PatType::TrussCombined2 &&
      e > 2 * v - 4)
    return false;
  /* wrong!!
  if (type == PatType::Truss && pattern.detilType.first - 2 > 2)
          return false;
  if (type == PatType::TrussCombined1 && pattern.detilType.first - 2 +
  pattern.detilType.second - 2 > 3) return false; if (type ==
  PatType::TrussCombined2 && pattern.detilType.first - 2 +
  pattern.detilType.second - 2 > 4) return false;
  */
  return true;
}

/**
 * 26/03/2020 update
 * @author Zifeng Yuan
 * larger ==> more complex
 */
void FastPGGraph::calCognition5(const string &folder) {
  printf("evaluating patterns cognition scores... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density, crossing;
  int i = 0;
  auto f1 = [](double x) -> double { return 1 - exp(-x); };
  auto f2 = [](double c1, double c2, double x) -> double {
    // f2(c2) = 0.5
    return 1 / (1 + exp(-c1 * (x - c2)));
  };
  for (auto &pattern : patterns) {
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    crossing = isPlanar(pattern)
                   ? 0
                   : pattern.getEdgeCnt() - 3 * pattern.getNodesCnt() + 6;
    /*1.0 * pattern.getEdgeCnt() * pattern.getEdgeCnt() * pattern.getEdgeCnt() /
     * pattern.getNodesCnt() / pattern.getNodesCnt() / 33.75 - 0.9 *
     * pattern.getNodesCnt();*/
    if (crossing < 0) crossing = 0;
    pattern.scoreCognition = (1.0 / 3) * f2(0.5, 10, pattern.getEdgeCnt()) +
                             (1.0 / 3) * f2(0.5, 0.6, density) +
                             (1.0 / 3) * f1(crossing);
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCognition < patterns[b].scoreCognition;
  });
  for (int i = 0; i < patterns.size(); ++i) {
    auto &pattern = patterns[pr[i]];
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    crossing = isPlanar(pattern)
                   ? 0
                   : pattern.getEdgeCnt() - 3 * pattern.getNodesCnt() + 6;
    if (crossing < 0) {
      cout << "crossing=" << crossing << " n=" << pattern.getNodesCnt()
           << " e=" << pattern.getEdgeCnt()
           << " id=" << to_string(i + 1) + "(" + to_string(pattern.id) + ")"
           << endl;
      crossing = 0;
    }
    outputFinalPngPattern(patterns[pr[i]], folder + "byCog5/",
                          "#" + to_string(i + 1) + "(" + to_string(pattern.id) +
                              ")" + "=" +
                              to_string(patterns[pr[i]].scoreCognition));
    outputFinalPngPattern(
        patterns[pr[i]],
        folder + "byCog5/" + to_string(patterns[pr[i]].type) + "/",
        "#" + to_string(i + 1) + "(" + to_string(pattern.id) + ")" + "=" +
            to_string(patterns[pr[i]].scoreCognition) + "=sz" +
            to_string(pattern.getEdgeCnt()) + "=d" + to_string(density) +
            "=cr" + to_string(crossing) + "=" +
            to_string(patterns[pr[i]].scoreCognition));
    if (pattern.type == PatType::Truss ||
        pattern.type == PatType::TrussCombined1 ||
        pattern.type == PatType::TrussCombined2)
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog5/",
                            "#" + to_string(i + 1) + "(" +
                                to_string(pattern.id) + ")" + "=" +
                                to_string(pattern.detilType.first) + "," +
                                to_string(pattern.detilType.second));
  }
  // ******************************* for debug only
  // *******************************
#endif
}

// larger ==> more complex
void FastPGGraph::calCognition5Single(Pattern &pattern) {
  double density, crossing;
  int i = 0;
  density = (2.0 * pattern.getEdgeCnt()) /
            (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
  crossing = isPlanar(pattern)
                 ? 0
                 : 1.0 * pattern.getEdgeCnt() * pattern.getEdgeCnt() *
                           pattern.getEdgeCnt() / pattern.getNodesCnt() /
                           pattern.getNodesCnt() / 33.75 -
                       0.9 * pattern.getNodesCnt();
  pattern.scoreCognition = (1.0 / 3) * (1 - exp(-pattern.getEdgeCnt())) +
                           (1.0 / 3) * (1 - exp(-density)) +
                           (1.0 / 3) * (1 - exp(-crossing));
}

/**
 * 14/04/2021 update: cog7 actually
 * @author Zifeng Yuan
 * larger ==> more complex
 */
void FastPGGraph::calCognition6(const string &folder) {
  printf("evaluating patterns cognition scores... \n");

  vector<int> degree;
  int d1, d2, d3, d4;
  double density, crossing;
  int i = 0;
  auto f1 = [](double x) -> double { return 1 - exp(-x); };
  auto f2 = [](double c1, double c2, double x) -> double {
    // f2(c2) = 0.5
    return 1 / (1 + exp(-c1 * (x - c2)));
  };
  for (auto &pattern : patterns) {
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    crossing = isPlanar(pattern)
                   ? 0
                   : pattern.getEdgeCnt() - 3 * pattern.getNodesCnt() + 6;
    /*1.0 * pattern.getEdgeCnt() * pattern.getEdgeCnt() * pattern.getEdgeCnt() /
     * pattern.getNodesCnt() / pattern.getNodesCnt() / 33.75 - 0.9 *
     * pattern.getNodesCnt();*/
    if (crossing < 0) crossing = 0;
    // pattern.scoreCognition = f1(pattern.getEdgeCnt() + density + crossing);
    // cog7
    pattern.scoreCognition =
        f2(0.5, 10, pattern.getEdgeCnt() + density + crossing);
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCognition < patterns[b].scoreCognition;
  });
  for (int i = 0; i < patterns.size(); ++i) {
    auto &pattern = patterns[pr[i]];
    density = (2.0 * pattern.getEdgeCnt()) /
              (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
    crossing = isPlanar(pattern)
                   ? 0
                   : pattern.getEdgeCnt() - 3 * pattern.getNodesCnt() + 6;
    if (crossing < 0) {
      cout << "crossing=" << crossing << " n=" << pattern.getNodesCnt()
           << " e=" << pattern.getEdgeCnt()
           << " id=" << to_string(i + 1) + "(" + to_string(pattern.id) + ")"
           << endl;
      crossing = 0;
    }
    outputFinalPngPattern(patterns[pr[i]], folder + "byCog6/",
                          "#" + to_string(i + 1) + "(" + to_string(pattern.id) +
                              ")" + "=" +
                              to_string(patterns[pr[i]].scoreCognition));
    outputFinalPngPattern(
        patterns[pr[i]],
        folder + "byCog6/" + to_string(patterns[pr[i]].type) + "/",
        "#" + to_string(i + 1) + "(" + to_string(pattern.id) + ")" + "=" +
            to_string(patterns[pr[i]].scoreCognition) + "=sz" +
            to_string(pattern.getEdgeCnt()) + "=d" + to_string(density) +
            "=cr" + to_string(crossing) + "=" +
            to_string(patterns[pr[i]].scoreCognition));
    if (pattern.type == PatType::Truss ||
        pattern.type == PatType::TrussCombined1 ||
        pattern.type == PatType::TrussCombined2)
      outputFinalPngPattern(patterns[pr[i]], folder + "byCog6/",
                            "#" + to_string(i + 1) + "(" +
                                to_string(pattern.id) + ")" + "=" +
                                to_string(pattern.detilType.first) + "," +
                                to_string(pattern.detilType.second));
  }
  // ******************************* for debug only
  // *******************************
#endif
}

// 14/04/2021 update: cog7 actually
// larger ==> more complex
void FastPGGraph::calCognition6Single(Pattern &pattern) {
  auto f1 = [](double x) -> double { return 1 - exp(-x); };
  auto f2 = [](double c1, double c2, double x) -> double {
    // f2(c2) = 0.5
    return 1 / (1 + exp(-c1 * (x - c2)));
  };
  double density = (2.0 * pattern.getEdgeCnt()) /
                   (pattern.getNodesCnt() * (pattern.getNodesCnt() - 1));
  double crossing = isPlanar(pattern)
                        ? 0
                        : pattern.getEdgeCnt() - 3 * pattern.getNodesCnt() + 6;
  /*1.0 * pattern.getEdgeCnt() * pattern.getEdgeCnt() * pattern.getEdgeCnt() /
   * pattern.getNodesCnt() / pattern.getNodesCnt() / 33.75 - 0.9 *
   * pattern.getNodesCnt();*/
  if (crossing < 0) crossing = 0;
  // pattern.scoreCognition = f1(pattern.getEdgeCnt() + density + crossing);
  // cog7
  pattern.scoreCognition =
      f2(0.5, 10, pattern.getEdgeCnt() + density + crossing);
}

/**
 * 2 groups  /just for test , finally we choose the way of  3 groups.
 * @author Zifeng Yuan
 */
void FastPGGraph::calCoverage3(const string &folder) {
  printf("evaluating patterns coverage scores... \n");

  // Two patterns A and B, if A is a subgraph of B, A's frequency is definitely
  // larger than B's frequency score == frequency for patterns of G_out
  vector<int> outside;
  for (int i = 0; i < patterns.size(); ++i)
    if (patterns[i].type > PatType::TrussCombined2) outside.emplace_back(i);
  int len = outside.size();
  vector<int> extraFrequency(len, 0);
  for (int i = 0; i < len; ++i)
    for (int j = i + 1; j < len; ++j)
      if (!(patterns[outside[i]].type == PatType::OutsideStar &&
            patterns[outside[j]].type == PatType::OutsideCombinedStar) &&
          !(patterns[outside[i]].type == PatType::OutsideCombinedStar &&
            patterns[outside[j]].type == PatType::OutsideStar)) {
        if (isSubGraph(outside[i], outside[j])) {
          // cout <<"v: "<<patterns[outside[i]].nodesCnt << "--- e:
          // "<<patterns[outside[i]].getEdgeCnt() << "-----vs----- "; cout << "
          // v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() << "\n";
          extraFrequency[i] += patterns[outside[j]].score;
        } else if (isSubGraph(outside[j], outside[i])) {
          // cout << "  v: " << patterns[outside[j]].nodesCnt << "--- e: " <<
          // patterns[outside[j]].getEdgeCnt() <<"-----vs----- " ; cout << "v: "
          // << patterns[outside[i]].nodesCnt << "--- e: " <<
          // patterns[outside[i]].getEdgeCnt() << "\n";
          extraFrequency[j] += patterns[outside[i]].score;
        }
      }

  for (int i = 0; i < len; ++i) patterns[outside[i]].score += extraFrequency[i];

  // score == frequency for patterns of G_out
  // coverage = frequency * size
  for (auto &pattern : patterns)
    // pattern.scoreCoverage = pattern.score;
    pattern.scoreCoverage = pattern.type <= PatType::TrussCombined2
                                ? pattern.score
                                : pattern.score * pattern.getEdgeCnt();

  // normalization
  vector<int> group(20);
  group[PatType::Truss] = group[PatType::TrussCombined1] =
      group[PatType::TrussCombined2] = 0;
  group[PatType::OutsideStar] = group[PatType::OutsideCombinedStar] = 1;
  group[PatType::OutsideLine] = group[PatType::OutsideCircle] =
      group[PatType::OutsideIndividual] = 1;

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  for (int i = 0; i < patterns.size(); ++i)
    outputFinalPngPattern(
        patterns[i],
        folder + "byCov3/group" + to_string(group[patterns[i].type]) + "/",
        to_string(patterns[i].scoreCoverage) + "=" +
            to_string(patterns[i].score));
    // ******************************* for debug only
    // *******************************
#endif

  vector<double> covmin(2, 1e18), covmax(2, 0);
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    if (pattern.scoreCoverage > covmax[g]) covmax[g] = pattern.scoreCoverage;
    if (pattern.scoreCoverage < covmin[g]) covmin[g] = pattern.scoreCoverage;
  }
  for (auto &pattern : patterns) {
    int g = group[pattern.type];
    pattern.scoreCoverage =
        (pattern.scoreCoverage - covmin[g] + 1) / (covmax[g] - covmin[g] + 1);
  }

#ifdef DEBUG_MODE
  // ******************************* for debug only
  // *******************************
  printf("exporting patterns by coverage score... \n");
  vector<int> pr;
  for (int i = 0; i < patterns.size(); ++i) pr.emplace_back(i);
  sort(pr.begin(), pr.end(), [&](const int &a, const int &b) -> bool {
    return patterns[a].scoreCoverage > patterns[b].scoreCoverage;
  });
  for (int i = 0; i < 50; ++i)
    outputFinalPngPattern(patterns[pr[i]], folder + "byCov3/",
                          "#" + to_string(i + 1) + "=" +
                              to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                              to_string(group[patterns[pr[i]].type]) + "=" +
                              to_string(patterns[pr[i]].scoreCoverage));
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::Truss) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov3/truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type == PatType::TrussCombined1 ||
        patterns[pr[i]].type == PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov3/combined_truss/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
  for (int i = 0, j = 0; i < pr.size() && j < 50; ++i)
    if (patterns[pr[i]].type > PatType::TrussCombined2) {
      ++j;
      outputFinalPngPattern(patterns[pr[i]], folder + "byCov3/gout/",
                            "#" + to_string(i + 1) + "=" +
                                to_string(patterns[pr[i]].getNodesCnt()) + "=" +
                                to_string(group[patterns[pr[i]].type]) + "=" +
                                to_string(patterns[pr[i]].scoreCoverage));
    }
    // ******************************* for debug only
    // *******************************
#endif
}

int bitCount(const ull &x) {
  return bitCountTable16[x & 0xffff] + bitCountTable16[(x >> 16) & 0xffff] +
         bitCountTable16[(x >> 32) & 0xffff] +
         bitCountTable16[(x >> 48) & 0xffff];
}

/*  we have 7 classes of canned patterns:
{ 0 simple k-truss patterns, 1 combined truss patterns,2 k-paths,
  3 k-cycles,4  k-stars,5  combined stars,6  other patterns}.
*/
int inline FastPGGraph::getFinalPatternClass(int patternId) {
  if (patterns[patternId].type <= 1 || patterns[patternId].type == 10)
    return 2;
  else if (patterns[patternId].type == 2 || patterns[patternId].type == 3)
    return 0;
  else if (patterns[patternId].type == 5 || patterns[patternId].type == 6)
    return 1;
  else if (patterns[patternId].type == 11)
    return 3;
  else if (patterns[patternId].type == 8)
    return 4;
  else if (patterns[patternId].type == 9)
    return 5;
  else
    return 6;
}

/*Approximate solution  using greedy strategy */
void FastPGGraph::setCoverSolve(vector<ull> &sets, const ull &aimset,
                                int &minstep, bool test) {
  ull status = 0;
  int cnt = 0;
  while (status != aimset && cnt < 100) {
    int maxcover = 0;
    ull maxcoverSet = 0;
    for (auto x : sets) {
      int addnum = bitCount(x & status ^ x);
      // printf("x=%d, addnum=%d\n", x, addnum);
      if (addnum > maxcover) {
        maxcoverSet = x;
        maxcover = addnum;
      }
    }
    cnt++;
    if (test) {
      int myclass = getFinalPatternClass(ullUsed2patternId[maxcoverSet]);
      ++patternUseClassCnt[myclass];
    }

    /* if(test) for test patternUseDistrb
            ++patternUseDistrb[bitCount(maxcoverSet)];*/
    status |= maxcoverSet;
  }
  minstep = cnt;

  /*modifed version
  ull status = 0;
  int cnt = 0;
  while (status != aimset)
  {
          int maxcover = 0;
          ull maxcoverSet = 0;
          for (auto x : sets)
          {
                  int addnum = bitCount(x & status ^ x);
                  if (addnum > maxcover)
                  {
                          maxcoverSet = x;
                          maxcover = addnum;
                  }
          }
          cnt=cnt+1+bitCount(maxcoverSet&status);
          status |= maxcoverSet;
  }
  minstep = cnt;*/
}

/* Get the best answer sulution  */
void FastPGGraph::setCoverSolve(vector<ull> &sets, const ull &aimset,
                                ull status, int curstep, int &minstep) {
  if (curstep >= minstep) return;
  if (status == aimset) {
    minstep = curstep;
    return;
  }

  sort(sets.begin(), sets.end(), [status](const ull &x, const ull &y) -> bool {
    return bitCount(x & status ^ x) > bitCount(y & status ^ y);
  });

  int maxSteps = bitCount(sets[0] & status ^ sets[0]);
  int minNeededStep = curstep + bitCount(status ^ aimset) / maxSteps;
  if (maxSteps == 1)  // a cut
  {
    if (minNeededStep < minstep) minstep = minNeededStep;
    return;
  }
  if (minNeededStep >= minstep) return;

  for (ull x : sets)
    setCoverSolve(sets, aimset, status | x, curstep + 1, minstep);

  /*Approximate solution  using greedy strategy */
  /*int maxcover=0;
  ull maxcoverSet=0;
  for (auto x : sets)
  {
          int addnum = bitCount(x & status ^ x);
          if (addnum > maxcover)
          {
                  maxcoverSet = x;
                  maxcover = addnum;
          }
  }
  setCoverSolve(sets, aimset, status | maxcoverSet, curstep + 1, minstep);*/
}

template <typename Query>
int FastPGGraph::minStep4filledQuery(const vector<int> &patternsId, Query query,
                                     int edgeCnt, bool test) {
  ullUsed2patternId.clear();
  vector<ull> edgeSets;
  for (auto id : patternsId) {
    // cout << "i am geting getAllSubGraphEdgeSets
    // "<<patterns[id].type<<"---"<<patterns[id].getEdgeCnt()<<" --- "<<
    // patterns[query].getEdgeCnt()<<"\n";
    unordered_set<ull> sets;
    // printf("id=%d\n", id);
    sets = getAllSubGraphEdgeSets(id, query);
    // printf("id=%d, sets size=%d\n", id, sets.size());
    // cout << "i got it! "<<sets.size()<<"\n";
    for (ull tset : sets) {
      edgeSets.emplace_back(tset);
      if (test) ullUsed2patternId[tset] = id;
    }
  }
  int ret = edgeCnt;
  // cout <<"set size:"<< edgeSets.size() << endl;
  // setCoverSolve(edgeSets, (1ULL << ret) - 1, 0, 0, ret);
  setCoverSolve(edgeSets, (1ULL << ret) - 1, ret, test);
  // cout << ret <<endl;
  return ret;
}

template int FastPGGraph::minStep4filledQuery<Pattern>(const vector<int> &patternsId, Pattern query,
                                     int edgeCnt, bool test);
template int FastPGGraph::minStep4filledQuery<int>(const vector<int> &patternsId, int query,
                                     int edgeCnt, bool test);

void FastPGGraph::calCandidatePatternsFinalScore(int curMinStep, double *maxps,
                                                 int *maxp, int threadId) {
  *maxps = -1;
  *maxp = -1;
  int c = 0;
  while (true) {
    int x;
    io_mutex.lock();
    // cout << candidatePatternsQueue.size() << endl;
    if (!candidatePatternsQueue.empty()) {
      x = candidatePatternsQueue.front();
      candidatePatternsQueue.pop();
    } else {
      // printf("thread %d: processed %d pattern(s) in total\n", threadId, c);
      io_mutex.unlock();
      return;
    }
    io_mutex.unlock();

    // if (c++ == 0)
    // printf("thread %d: processing no.%d pattern\n", threadId, x);

    /*vector<int>com = patternCommonSubgraph(j, newJoinPatternId);
    int comEdgeSize = com[1];
    int comNodeSize = com[0];*/

    // double newScore = double(patterns[j].getEdgeCnt() - comEdgeSize + 1) /
    // patterns[j].getEdgeCnt() * patterns[j].finalScore;
    int minStep =
        minStep4filledQuery(finalPatternIdSet, x, patterns[x].getEdgeCnt());

    // printf("thread %d: finish no.%d pattern\n", threadId, x);
    if (patterns[x].getEdgeCnt() >= MAX_PATTERN_SIZE / 2.0 + 1 &&
        minStep < curMinStep)  // big pattern but useless
    {
      continue;
    }
    /*if (newScore < patterns[j].dynamicScore)
        patterns[j].dynamicScore = newScore;*/

    double newScore =
        minStep * 1.0 / patterns[x].getEdgeCnt() * patterns[x].finalScore;
    /*for (int fid : finalPatternIdSet)
    {
        if (patterns[fid].type == patterns[j].type && patterns[j].type!=
    PatType::OutsideIndividual)
        {
            vector<int>com = patternCommonSubgraph(j, fid);
            int comEdgeSize = com[1];
            int comNodeSize = com[0];
            if (comEdgeSize >= min(patterns[fid].getEdgeCnt(),
    patterns[j].getEdgeCnt()) - 2)  //give the penaty
            {
                if (patterns[j].getEdgeCnt() <= 6 && comEdgeSize >=
    min(patterns[fid].getEdgeCnt(), patterns[j].getEdgeCnt()) - 1) newScore /=
    2; else newScore /= 3;
            }
        }
    }*/
    if (newScore > *maxps) {
      *maxps = newScore;
      *maxp = x;
    }
  }
}

//********************************** final test
//*******************************************************//

/*
random gen querys "times" times of size "size"
*/

void FastPGGraph::genQueryRandomly(string filename, int times, int size,
                                   bool gin_only) {
  random_device rd;
  default_random_engine engine(rd());
  srand((unsigned)time(NULL));
  int _times = times;

  while (times--) {
    vector<vector<int>> hash(size + 1, vector<int>(size + 1, 0));
    uniform_int_distribution<> dis(0, NUM_V - 1);
    int startV = dis(engine);
    Pattern pattern(size + 1);  // note : not real size. now used as a query
    unordered_map<int, int> nodeMap;
    int nodeId = 1;
    vector<int> Vset;
    Vset.emplace_back(startV);
    nodeMap[startV] = 0;
    bool found_pattern = true;  // for some small connect part.
    for (int i = 0; i < size; i++) {
      bool found_edge = false;
      for (int tri = 0; tri < 100; ++tri) {
        uniform_int_distribution<> dis1(0, Vset.size() - 1);
        int s = Vset[dis1(engine)];

        int e = s;
        // uniform_int_distribution<> dis2(0, neibours[s].size() - 1);
        int _size = neibours[s].size();
        if (rand() % 100 < 80) {  // 80% inner edge
          for (int z = 0; z < _size; ++z) {
            // int candidate_e = neibours[s][dis2(engine)].neighbId;
            int rd = rand() % _size;
            int candidate_e = neibours[s][rd].neighbId;
            if (/*edgeType[neibours[s][rd].edgeId] == EdgeType::RealTruss && */
                nodeMap.find(candidate_e) != nodeMap.end() &&
                hash[nodeMap[s]][nodeMap[candidate_e]] == 0) {
              if (gin_only) {
                // Gin
                if (trussClass[neibours[s][rd].edgeId] > 2) {
                  e = candidate_e;
                  break;
                }
              } else {
                e = candidate_e;
                break;
              }
            }
          }
        } else {
          e = neibours[s][rand() % _size].neighbId;
          // for (int z = 0; z < _size; ++z) {
          // 	int rd = rand() % _size;
          // 	if (edgeType[neibours[s][rd].edgeId] == EdgeType::RealTruss) {
          // 		e = neibours[s][rd].neighbId;
          // 		break;
          // 	}
          // }
        }

        s = nodeMap[s];
        if (nodeMap.find(e) == nodeMap.end()) {
          Vset.emplace_back(e);
          nodeMap[e] = nodeId++;
        }
        e = nodeMap[e];

        if (e == s || hash[s][e]) continue;
        pattern.addEdge(s, e);
        hash[e][s] = hash[s][e] = 1;
        found_edge = true;
        break;
      }
      if (!found_edge) {
        found_pattern = false;
        break;
      }
    }
    if (found_pattern) {
      pattern.nodesCnt = nodeId;
      queryPatterns[size].emplace_back(pattern);
      // outputPngPattern(pattern, filename);
    } else {
      ++times;
    }
  }
  assert(queryPatterns[size].size() == _times);
}

void FastPGGraph::genSubgraph(int size, string filename) {
  random_device rd;
  default_random_engine engine(rd());
  srand((unsigned)time(NULL));

  int t = 0;
  while (true) {
    vector<vector<int>> hash(size + 1, vector<int>(size + 1, 0));
    uniform_int_distribution<> dis(0, NUM_V - 1);
    int startV = dis(engine);
    Pattern pattern(size + 1);  // note : not real size. now used as a query
    unordered_map<int, int> nodeMap;
    int nodeId = 1;
    vector<int> Vset;
    Vset.emplace_back(startV);
    nodeMap[startV] = 0;
    bool found_pattern = true;  // for some small connect part.
    for (int i = 0; i < size; i++) {
      bool found_edge = false;
      for (int tri = 0; tri < 1000; ++tri) {
        uniform_int_distribution<> dis1(0, Vset.size() - 1);
        int s = Vset[dis1(engine)];

        int e = s;
        // uniform_int_distribution<> dis2(0, neibours[s].size() - 1);
        int _size = neibours[s].size();
        if (rand() % 100 < 80) {  // 80% inner edge
          for (int z = 0; z < _size; ++z) {
            // int candidate_e = neibours[s][dis2(engine)].neighbId;
            int rd = rand() % _size;
            int candidate_e = neibours[s][rd].neighbId;
            if (/*edgeType[neibours[s][rd].edgeId] == EdgeType::RealTruss && */
                nodeMap.find(candidate_e) != nodeMap.end() &&
                hash[nodeMap[s]][nodeMap[candidate_e]] == 0) {
              e = candidate_e;
              break;
            }
          }
        } else {
          e = neibours[s][rand() % _size].neighbId;
          // for (int z = 0; z < _size; ++z) {
          // 	int rd = rand() % _size;
          // 	if (edgeType[neibours[s][rd].edgeId] == EdgeType::RealTruss) {
          // 		e = neibours[s][rd].neighbId;
          // 		break;
          // 	}
          // }
        }

        s = nodeMap[s];
        if (nodeMap.find(e) == nodeMap.end()) {
          Vset.emplace_back(e);
          nodeMap[e] = nodeId++;
        }
        e = nodeMap[e];

        if (e == s || hash[s][e]) continue;
        pattern.addEdge(s, e);
        hash[e][s] = hash[s][e] = 1;
        found_edge = true;
        break;
      }
      if (!found_edge) {
        found_pattern = false;
        break;
      }
    }
    if (found_pattern) {
      ofstream fout(filename);
      for (const auto &e : pattern.edges)
        fout << e.s + 1 << " " << e.e + 1 << endl;
      fout.close();
      break;
    } else
      cout << "tried time:" << ++t << endl;
  }
}

/*
filename : for output results
timeseach : the number of each size

*/

double FastPGGraph::finalTest(string fileName, int testFromSize, int testToSize,
                              int timesEach) {
  // gen query  randomly  based on the data base.
  printf("final test from size %d to size %d with %d\n", testFromSize,
         testToSize, timesEach);
  // vector<vector<Pattern>> queryPatterns(testToSize+1);

  for (int size = testFromSize; size <= testToSize; ++size) {
    printf("size %d is geerating...\n", size);
    genQueryRandomly(fileName + "query/", timesEach, size);
  }

  patternUseDistrb.resize(testToSize + 1, 0);

  double avgTotal = 0;

  string dataName;
  int tmark = 1;
  for (char c : fileName) {
    if (c == '/') {
      if (tmark < 0) break;
      tmark = -tmark;
      continue;
    }
    if (tmark < 0) dataName += c;
  }

  ofstream out(fileName + dataName + "_results.txt");
  patternUseClassCnt.clear();
  patternUseClassCnt.resize(7, 0);

  for (int i = testFromSize; i <= testToSize; ++i) {
    printf("size %d is testing...\n", i);
    double avg4sizeI = 0;
    for (int j = 0; j < queryPatterns[i].size(); ++j) {
      int curstep = minStep4filledQuery(finalPatternIdSet, queryPatterns[i][j],
                                        queryPatterns[i][j].getEdgeCnt(), true);
      // out << curstep << "\n";
      avg4sizeI += curstep;
    }
    avg4sizeI /= timesEach;
    avgTotal += avg4sizeI;
    // out << "Querys with size " << i << "avg step : " << avg4sizeI << "\n";
    out << avg4sizeI / i << "\n";
  }
  avgTotal /= (testToSize - testFromSize + 1);
  // out << "All querys  avg step : " << avgTotal << "\n";

  for (int i = 1; i <= 15; ++i) {
    out << patternUseDistrb[i] << "\n";
  }
  return avgTotal;
}

void FastPGGraph::genQuery(int testFromSize, int testToSize, int timesEach,
                           bool gin_only) {
  printf("gen querys...\n");
  queryPatterns.resize(testToSize + 1);
  for (int size = testFromSize; size <= testToSize; ++size) {
    printf("size %d is generating...\n", size);
    genQueryRandomly("unused", timesEach, size, gin_only);
  }
  /*for(int i=4;i<31;i++)
          printf("%d\n",queryPatterns[i].size());*/
}

/*test avg 1-step/query_size  by given finalPatternIdSet */
double FastPGGraph::paratest(ofstream &out, int testFromSize, int testToSize,
                             int timesEach, vector<int> finalIdSet, int mark,
                             bool printIndividualRR) {
  printf("[func] paratest\n");
  assert(!queryPatterns.empty());
  cout << "\nfinal patterns set size: " << finalIdSet.size() << "\n";
  out << "\nfinal patterns set size: " << finalIdSet.size() << "\n";
  double avgTotal = 0;
  printf("testing...\n");
  double min_step_div_e = 1;
  double avg_step_div_e = 0;
  double min_step_div_en = 1;
  double avg_step_div_en = 0;
  // patternUseClassCnt.clear();
  // patternUseClassCnt.resize(7, 0);
  for (int i = testFromSize; i <= testToSize; ++i) {
    printf("%d,", i);
    double avg4sizeI = 0;
    // assert(queryPatterns[i].size() == timesEach);
    for (int j = 0; j < queryPatterns[i].size(); ++j) {
      int curstep;
      if (LOADED)
        curstep = Step4curQuery[i][j];
      else {
        curstep = minStep4filledQuery(finalIdSet, queryPatterns[i][j],
                                      queryPatterns[i][j].getEdgeCnt(), true);
        Step4curQuery[i][j] = curstep;
      }

      auto step_div_e = 1.0 * curstep / queryPatterns[i][j].getEdgeCnt();
      auto step_div_en =
          1.0 * curstep /
          (queryPatterns[i][j].getEdgeCnt() + queryPatterns[i][j].nodesCnt);
      if (step_div_e < min_step_div_e) min_step_div_e = step_div_e;
      avg_step_div_e += step_div_e;
      if (step_div_en < min_step_div_en) min_step_div_en = step_div_en;
      avg_step_div_en += step_div_en;
    }
    avg4sizeI /= timesEach;  // /100
    if (printIndividualRR) out << "size" << i << " avgRR=" << avg4sizeI << endl;
    avgTotal += avg4sizeI;
    // out << "Querys with size " << i << "avg step : " << avg4sizeI << "\n";
    // out << avg4sizeI / i << "\n";
  }
  avg_step_div_e /= (timesEach * 1.0 * (testToSize - testFromSize + 1));
  avg_step_div_en /= (timesEach * 1.0 * (testToSize - testFromSize + 1));

  double avgCov = 0, avgCog = 0;
  double avgDiv, minDiv, maxDiv;
  for (int id : finalIdSet) {
    avgCov += patterns[id].scoreCoverage;
    avgCog += patterns[id].scoreCognition;
  }
  avgCov /= finalIdSet.size();
  avgCog /= finalIdSet.size();
  collectDiversityScores(patterns, finalIdSet, &avgDiv, &minDiv, &maxDiv);

  cout << "\navgCov=" << avgCov << "\tavgCog=" << avgCog
       << "\tavgDiv=" << avgDiv << "\tminDiv=" << minDiv
       << "\tmaxDiv=" << maxDiv << "\n";
  out << "\navgCov=" << avgCov << "\tavgCog=" << avgCog << "\tavgDiv=" << avgDiv
      << "\tminDiv=" << minDiv << "\tmaxDiv=" << maxDiv << "\n";
  cout << "avg_step_div_e, min_step_div_e, avg_step_div_en, min_step_div_en:\n";
  out << "avg_step_div_e, min_step_div_e, avg_step_div_en, min_step_div_en:\n";
  cout << avg_step_div_e << " " << min_step_div_e << " " << avg_step_div_en
       << " " << min_step_div_en << "\n";
  out << avg_step_div_e << " " << min_step_div_e << " " << avg_step_div_en
      << " " << min_step_div_en << "\n";
  cout << "reversed= " << 1 - avg_step_div_e << " " << 1 - min_step_div_e << " "
       << 1 - avg_step_div_en << " " << 1 - min_step_div_en << "\n";
  out << "reversed= " << 1 - avg_step_div_e << " " << 1 - min_step_div_e << " "
      << 1 - avg_step_div_en << " " << 1 - min_step_div_en << "\n";

  /*double avgRR=avgTotal / (testToSize - testFromSize + 1);
  newInfo[mark].finalPatternsGenTime
  cout <<"\n"<< avgRR << " " << summaryInfo.finalPatternsGenTime << "\n\n"; //if
  you need output time. out << "AVGRR: " << avgRR <<"\n"; out << ret << "\n";
  out << "MAXrr: " << maxRR << "\n";
  double tottime=summaryInfo.calTrussClass4EdgesTime +
  summaryInfo.calScore4CombinedTrussPatternsTime
          + summaryInfo.searchStarsTime +
  summaryInfo.searchLinesAndCirclesAndIndividualsTime
          + summaryInfo.calCovAndCogScoreTime +
  summaryInfo.finalPatternsGenTime; out << "total time : " << tottime << "s\n";
  out << "total patternUse for each Class:\n";
  out << "avgCog,avgrr,maxrr,step_div_e,1-step_div_e:"<< "\n";

  double avgcog = 0;
  for (int id : finalIdSet)
          avgcog += patterns[id].scoreCognition;
  avgcog /= finalIdSet.size();
  out << avgcog << " " << avgRR << " " << maxRR << " " << step_div_e  << " " <<
  1-step_div_e << "\n";
  //totout << tottime << "\n";
  totout <<avgcog <<" "<< avgRR << " " << maxRR << " "<< step_div_e << " " << 1
  - step_div_e << "\n"; for (int i = 0; i < 7; i++)
  {
          out << patternUseClassCnt[i] << " ";
  }
  return avgRR;*/
  return 0.0;
}

void FastPGGraph::paratest2(ofstream &out, const vector<int> &finalIdSet) {
  printf("[func] paratest2\n");
  assert(!query2Patterns.empty());
  cout << "\nfinal patterns set size: " << finalIdSet.size() << "\n";
  out << "\nfinal patterns set size: " << finalIdSet.size() << "\n";
  /*{
          vector<int> size_cnt(50, 0);
          for (const auto &p : query2Patterns)
                  ++size_cnt[p.getEdgeCnt()];
          for (int i = 5; i < 30; ++i)
                  cout << "query_num[" << i << "]=" << size_cnt[i] << endl;
  }*/
  printf("testing...\n");
  vector<double> group_step_div_e(query2Setting.size(), 0);
  vector<double> group_step_div_en(query2Setting.size(), 0);
  double min_step_div_e = 1;
  double avg_step_div_e = 0;
  double min_step_div_en = 1;
  double avg_step_div_en = 0;
  int group_pre_cnt = 0;
  for (int i = 0, group_id = 0; i < query2Patterns.size(); ++i) {
    while (i - group_pre_cnt == query2Setting[group_id]) {
      group_pre_cnt = i;
      ++group_id;
    }
    if (i % 100 == 0) printf("i=%d, group_id=%d\n", i, group_id);
    const auto &p = query2Patterns[i];
    int cur_step = minStep4filledQuery(finalIdSet, p, p.getEdgeCnt(), true);
    // printf("p.edgecnt=%d, cur_step=%d\n", p.getEdgeCnt(), cur_step);
    assert(cur_step <= p.getEdgeCnt());
    // if (cur_step > p.getEdgeCnt()) cur_step = p.getEdgeCnt();
    double step_div_e = 1.0 * cur_step / p.getEdgeCnt();
    double step_div_en = 1.0 * cur_step / (p.getEdgeCnt() + p.getNodesCnt());
    if (step_div_e < min_step_div_e) min_step_div_e = step_div_e;
    avg_step_div_e += step_div_e;
    group_step_div_e[group_id] += step_div_e;
    if (step_div_en < min_step_div_en) min_step_div_en = step_div_en;
    avg_step_div_en += step_div_en;
    group_step_div_en[group_id] += step_div_en;
  }
  avg_step_div_e /= query2Patterns.size();
  avg_step_div_en /= query2Patterns.size();
  for (int i = 0; i < 7; ++i) {
    group_step_div_e[i] /= query2Setting[i];
    group_step_div_en[i] /= query2Setting[i];
    cout << "group_step_div_e/en[" << query2SettingName[i]
         << "]=" << group_step_div_e[i] << "," << group_step_div_en[i] << endl;
    out << "group_step_div_e/en[" << query2SettingName[i]
        << "]=" << group_step_div_e[i] << "," << group_step_div_en[i] << endl;
    if (i == 4) {
      cout << "avg of star: "
           << (group_step_div_e[3] + group_step_div_e[4]) * 0.5 << ","
           << (group_step_div_en[3] + group_step_div_en[4]) * 0.5 << endl;
      out << "avg of star: "
          << (group_step_div_e[3] + group_step_div_e[4]) * 0.5 << ","
          << (group_step_div_en[3] + group_step_div_en[4]) * 0.5 << endl;
    }
  }

  double avgCov = 0, avgCog = 0;
  double avgDiv, minDiv, maxDiv;
  for (int id : finalIdSet) {
    avgCov += patterns[id].scoreCoverage;
    avgCog += patterns[id].scoreCognition;
  }
  avgCov /= finalIdSet.size();
  avgCog /= finalIdSet.size();
  collectDiversityScores(patterns, finalIdSet, &avgDiv, &minDiv, &maxDiv);

  cout << "\navgCov=" << avgCov << "\tavgCog=" << avgCog
       << "\tavgDiv=" << avgDiv << "\tminDiv=" << minDiv
       << "\tmaxDiv=" << maxDiv << "\n";
  out << "\navgCov=" << avgCov << "\tavgCog=" << avgCog << "\tavgDiv=" << avgDiv
      << "\tminDiv=" << minDiv << "\tmaxDiv=" << maxDiv << "\n";
  cout << "avg_step_div_e, min_step_div_e, avg_step_div_en, min_step_div_en:\n";
  out << "avg_step_div_e, min_step_div_e, avg_step_div_en, min_step_div_en:\n";
  cout << avg_step_div_e << " " << min_step_div_e << " " << avg_step_div_en
       << " " << min_step_div_en << "\n";
  out << avg_step_div_e << " " << min_step_div_e << " " << avg_step_div_en
      << " " << min_step_div_en << "\n";
  cout << "reversed= " << 1 - avg_step_div_e << " " << 1 - min_step_div_e << " "
       << 1 - avg_step_div_en << " " << 1 - min_step_div_en << "\n";
  out << "reversed= " << 1 - avg_step_div_e << " " << 1 - min_step_div_e << " "
      << 1 - avg_step_div_en << " " << 1 - min_step_div_en << "\n";
}

/*
test by differ class of pattern.
1 only defualt pattens.
2 +truss
3 +combinedtruss
4 all

*/

void FastPGGraph::testByPatternClass(ofstream &out, int testFromSize,
                                     int testToSize, int timesEach) {
  vector<int> finalPatterIds;
  /*for defualt pattens*/
  // double avgcog = 0;
  for (int i = 0; i < 4; i++) {
    finalPatterIds.emplace_back(finalPatternIdSet[i]);
    calCognition3Single(patterns[finalPatternIdSet[i]]);
    // avgcog += patterns[finalPatternIdSet[i]].scoreCognition;
  }
  // totout << avgcog/4.0 << " ";
  paratest(out, testFromSize, testToSize, timesEach, finalPatterIds, true);
  /*for defualt+truss pattens*/
  int cnt = 4;
  for (int i = 4; i < finalPatternIdSet.size(); i++) {
    if (patterns[finalPatternIdSet[i]].type == PatType::Truss) {
      ++cnt;
      finalPatterIds.emplace_back(finalPatternIdSet[i]);
      // avgcog += patterns[finalPatternIdSet[i]].scoreCognition;
    }
  }

  // paratest(out, testFromSize, testToSize, timesEach, finalPatterIds, true);
  /*for defualt+truss+combined-truss pattens*/
  for (int i = 4; i < finalPatternIdSet.size(); i++) {
    if (patterns[finalPatternIdSet[i]].type == PatType::TrussCombined1 ||
        patterns[finalPatternIdSet[i]].type == PatType::TrussCombined2) {
      ++cnt;
      finalPatterIds.emplace_back(finalPatternIdSet[i]);
      // avgcog += patterns[finalPatternIdSet[i]].scoreCognition;
    }
  }
  // totout << avgcog/cnt << " ";
  paratest(out, testFromSize, testToSize, timesEach, finalPatterIds, true);

  /*all*/
  // avgcog = 0;
  // for (int i = 0; i < finalPatternIdSet.size(); i++)
  //{
  // if (patterns[finalPatternIdSet[i]].scoreCognition > maxcog)
  //	maxcog = patterns[finalPatternIdSet[i]].scoreCognition;
  // avgcog += patterns[finalPatternIdSet[i]].scoreCognition;
  //}
  // totout << avgcog/30 << " ";
  paratest(out, testFromSize, testToSize, timesEach, finalPatternIdSet, true);
  totout << "\n";
}

/* top-k common sub set test
finalSet4Datastores[i]:   dataset[i] finalpatterns.
*/
void FastPGGraph::getcomset4alldataset() {
  for (int i = 10; i <= 50; i += 10)  // top-10 20 30 40 50 test
  {
    cout << "top:" << i << endl;
    vector<Pattern> cur;
    for (int j = 0; j < i && j < finalSet4Datastores[0].size(); j++) {
      cur.emplace_back(finalSet4Datastores[0][j]);
    }

    for (int j = 1; j < 10; j++)  // travel 10 datasets
    {
      if (cur.size() == 4)  // only 4 default pattterns
        break;
      vector<Pattern> tcur;
      for (int k = 0; k < cur.size(); k++) {
        for (int q = 0; q < i && q < finalSet4Datastores[j].size(); q++) {
          if (isSubGraph(cur[k], finalSet4Datastores[j][q]) &&
              isSubGraph(finalSet4Datastores[j][q], cur[k]))
            tcur.emplace_back(cur[k]);
        }
      }
      cur = tcur;
    }
    cout << cur.size() << endl;
  }
}

/*use for test the result of final patterns only sorted by cov.*/
void FastPGGraph::predo4test() {
  sort(patterns.begin(), patterns.end(), [](Pattern a, Pattern b) -> bool {
    return a.scoreCoverage > b.scoreCoverage;
  });
  /*sort(patterns.begin(), patterns.end(), [](Pattern a, Pattern b)->bool {
          return a.score > b.score;
  });*/
  genDefaultPatterns();
  boostPatterns.resize(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)  // now used new id.
  {
    boostPatterns[j] = pattern2boostGraph(patterns[j]);
  }
}
/*use for test the result of final patterns only sorted by cov.*/
void FastPGGraph::choosefinalpatternsbycov(int topk) {
  finalPatternIdSet.clear();
  for (int i = patterns.size() - 4; i < patterns.size(); i++)
    finalPatternIdSet.emplace_back(i);  // put last 4 default patterns
  for (int i = 4; i < topk; i++) {
    finalPatternIdSet.push_back(i);
  }
}

/**
 * 3 results: number, diversity, cognition
 * 4 groups: default(ignored), truss, combinedtruss, G_out
 **/
void FastPGGraph::print3ResultsOf4Groups(const string &folder) {
  printf("printing 3 results of 4 groups... \n");
  vector<int> othersFinalSet;
  vector<int> num(4, 0);
  vector<double> min_cog_load(4, 100.0);
  vector<double> max_cog_load(4, -100.0);
  vector<int> min_diversity(4, 100);
  vector<int> max_diversity(4, -100);
  for (auto &pid : finalPatternIdSet) {
    int g = 0;
    // string prefix = "default";
    if (patterns[pid].type == PatType::Truss)
      // prefix = "truss";
      g = 1;
    else if (patterns[pid].type >= PatType::TrussCombined1 &&
             patterns[pid].type <= PatType::TrussCombined2)
      // prefix = "combinedtruss";
      g = 2;
    if (patterns[pid].type >= PatType::OutsideStar)
      // prefix = "Gout";
      g = 3;
    othersFinalSet.clear();
    for (auto opid : finalPatternIdSet)
      if (patterns[opid].getEdgeCnt() < patterns[pid].getEdgeCnt())
        othersFinalSet.emplace_back(opid);
    if (g == 0) {
      calCognition3Single(patterns[pid]);
    }
    double cog_load = patterns[pid].scoreCognition;
    int div = minStep4filledQuery(othersFinalSet, pid,
                                  patterns[pid].getEdgeCnt(), true);

    ++num[g];
    if (cog_load < min_cog_load[g]) min_cog_load[g] = cog_load;
    if (cog_load > max_cog_load[g]) max_cog_load[g] = cog_load;
    if (div < min_diversity[g]) min_diversity[g] = div;
    if (div > max_diversity[g]) max_diversity[g] = div;
    // cout << "pid=" << pid << " prefix=" << prefix << " minstep=" <<
    // minStep4filledQuery(othersFinalSet, pid, patterns[pid].getEdgeCnt(),
    // true) << endl; outputFinalPngPattern(patterns[pid], folder + prefix +
    // "_coverage/",
    // 	"#" + to_string(pid) + "=" + to_string(patterns[pid].getEdgeCnt()) + "="
    // + to_string(patterns[pid].scoreCoverage));
    // outputFinalPngPattern(patterns[pid], folder + prefix + "_cognition/",
    // 	"#" + to_string(pid) + "=" + to_string(patterns[pid].getEdgeCnt()) + "="
    // + to_string(minStep4filledQuery(othersFinalSet, pid,
    // patterns[pid].getEdgeCnt(), true))); outputFinalPngPattern(patterns[pid],
    // folder + prefix + "_diversity/",
    // 	"#" + to_string(pid) + "=" + to_string(patterns[pid].getEdgeCnt()) + "="
    // + to_string(patterns[pid].scoreCognition));
  }
  totout << "        "
            "\t\tnumber\t%\tmin_cog_load\tmax_cog_load\tmin_diversity\tmax_"
            "diversity\n";
  totout << "default:\t\t" << num[0] << "\t"
         << num[0] * 1.0 / finalPatternIdSet.size() << "\t" << min_cog_load[0]
         << "\t" << max_cog_load[0] << "\t" << min_diversity[0] << "\t"
         << max_diversity[0] << "\n";
  totout << "truss:\t\t" << num[1] << "\t"
         << num[1] * 1.0 / finalPatternIdSet.size() << "\t" << min_cog_load[1]
         << "\t" << max_cog_load[1] << "\t" << min_diversity[1] << "\t"
         << max_diversity[1] << "\n";
  totout << "combined_truss:\t" << num[2] << "\t"
         << num[2] * 1.0 / finalPatternIdSet.size() << "\t" << min_cog_load[2]
         << "\t" << max_cog_load[2] << "\t" << min_diversity[2] << "\t"
         << max_diversity[2] << "\n";
  totout << "Gout:\t\t" << num[3] << "\t"
         << num[3] * 1.0 / finalPatternIdSet.size() << "\t" << min_cog_load[3]
         << "\t" << max_cog_load[3] << "\t" << min_diversity[3] << "\t"
         << max_diversity[3] << "\n";
}

/*test vs with graphlets*/
void FastPGGraph::genGraphlet() {
  Pattern t2 = Pattern(2);
  t2.addEdge(0, 1);
  t2.setType(PatType::OneEdge);
  t2.detilType = make_pair(1, -1);
  patterns.emplace_back(t2);

  Pattern t3 = Pattern(3);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(3);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  /* 4 node */
  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(1, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 0);
  t3.addEdge(2, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(0, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(4);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(0, 2);
  t3.addEdge(1, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  /* 5 node */
  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(2, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(0, 2);
  t3.addEdge(0, 3);
  t3.addEdge(0, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 0);
  t3.addEdge(0, 3);
  t3.addEdge(1, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 4);
  t3.addEdge(2, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(0, 2);
  t3.addEdge(0, 3);
  t3.addEdge(0, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 4);
  t3.addEdge(0, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 4);
  t3.addEdge(2, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(0, 2);
  t3.addEdge(0, 3);
  t3.addEdge(0, 4);
  t3.addEdge(1, 2);
  t3.addEdge(3, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 4);
  t3.addEdge(1, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 4);
  t3.addEdge(4, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 4);
  t3.addEdge(4, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 0);
  t3.addEdge(4, 2);
  t3.addEdge(4, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(1, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 0);
  t3.addEdge(2, 3);
  t3.addEdge(4, 0);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(0, 2);
  t3.addEdge(1, 2);
  t3.addEdge(3, 1);
  t3.addEdge(1, 4);
  t3.addEdge(2, 3);
  t3.addEdge(4, 3);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(4, 0);
  t3.addEdge(4, 1);
  t3.addEdge(4, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 0);
  t3.addEdge(4, 0);
  t3.addEdge(4, 3);
  t3.addEdge(4, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(4, 0);
  t3.addEdge(4, 1);
  t3.addEdge(4, 3);
  t3.addEdge(4, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(1, 2);
  t3.addEdge(2, 3);
  t3.addEdge(3, 0);
  t3.addEdge(2, 0);
  t3.addEdge(1, 3);
  t3.addEdge(4, 0);
  t3.addEdge(4, 3);
  t3.addEdge(4, 2);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);

  t3 = Pattern(5);
  t3.addEdge(0, 1);
  t3.addEdge(0, 2);
  t3.addEdge(0, 3);
  t3.addEdge(0, 4);
  t3.addEdge(1, 2);
  t3.addEdge(1, 3);
  t3.addEdge(1, 4);
  t3.addEdge(2, 3);
  t3.addEdge(2, 4);
  t3.addEdge(3, 4);
  t3.setType(PatType::OneEdge);
  t3.detilType = make_pair(1, -1);
  patterns.emplace_back(t3);
}

/*test vs with graphlets*/
void FastPGGraph::graphletTest(int topk) {
  boostPatterns.resize(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)  // now used new id.
  {
    calCognition6Single(patterns[j]);
    boostPatterns[j] = pattern2boostGraph(patterns[j]);
  }
  finalPatternIdSet.clear();
  for (int i = 0; i < topk; i++) {
    finalPatternIdSet.push_back(i);
  }
  // re-assign coverage & cognition
  for (int i = 0; i < topk; ++i) {
    patterns[i].scoreCoverage = 1.0;
    calCognition6Single(patterns[i]);
  }
  // re-assign diversity
  vector<vector<double>> netSimileSignature(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)
    netSimileSignature[j] = net_simile::Signature(patterns[j]);
  for (int i = 0; i < topk; ++i) {
    double sim = 0.0;
    for (int j = 0; j < topk; ++j)
      if (i != j) {
        double tmp = net_simile::CamberraDistance(netSimileSignature[i],
                                                  netSimileSignature[j]);
        if (tmp > sim) sim = tmp;
      }
    patterns[i].scoreDiversity = 1.0 / sim;
  }
}

void FastPGGraph::peregrineTest(int size) {
  finalPatternIdSet.clear();
  // default patterns
  for (int i = patterns.size() - 4; i < patterns.size(); i++)
    finalPatternIdSet.emplace_back(i);
  patternChooseByGreedy(size, "");
  assert(finalPatternIdSet.size() == size);

  /*for (int i = 0; i < finalPatternIdSet.size(); i++) {
          //		outputPngPattern(patterns[finalPatternIdSet[i]],
  filename); outputFinalPngPattern(patterns[finalPatternIdSet[i]],
  "peregrine30/", "#" + to_string(i + 1)
                  + "-" + to_string(finalPatternIdSet[i])
                  + "=" + to_string(patterns[finalPatternIdSet[i]].type)
                  + "=" +
  to_string(patterns[finalPatternIdSet[i]].scoreCoverage)
                  + "-" +
  to_string(patterns[finalPatternIdSet[i]].scoreCognition)
                  + "-" +
  to_string(patterns[finalPatternIdSet[i]].scoreDiversity));
  }*/
}

/*test vs with random*/
void FastPGGraph::randomTest(int topk) {
  boostPatterns.resize(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)  // now used new id.
  {
    calCognition6Single(patterns[j]);
    boostPatterns[j] = pattern2boostGraph(patterns[j]);
  }
  finalPatternIdSet.clear();
  for (int i = 0; i < topk; i++) {
    finalPatternIdSet.push_back(i);
  }
}

/*for final Pattern Set Size*/
void FastPGGraph::prejob4finalPatternSetSize() {
  tempfinalset = finalPatternIdSet;
}
/*for final Pattern Set Size*/
void FastPGGraph::Test4finalPatternSetSize(ofstream &myout, int topk) {
  finalPatternIdSet.clear();
  double avgcog = 0;
  for (int i = 0; i < topk; i++) {
    avgcog += patterns[tempfinalset[i]].scoreCognition;
    finalPatternIdSet.push_back(tempfinalset[i]);
  }
  paratest(myout, 4, 30, 100, finalPatternIdSet, 0);
}

void FastPGGraph::saveFinalpatterns(string fname, int num) {
  ofstream sout(fname);
  assert(0 == num || finalPatternIdSet.size() == num);
  for (int i = 0; i < finalPatternIdSet.size(); i++) {
    Pattern &tmpp = patterns[finalPatternIdSet[i]];
    sout << tmpp.nodesCnt << " " << tmpp.getEdgeCnt() << "\n";
    for (int j = 0; j < tmpp.getEdgeCnt(); j++)
      sout << tmpp.edges[j].s << " " << tmpp.edges[j].e << "\n";
    sout << tmpp.scoreCognition << " " << tmpp.scoreCoverage << " "
         << tmpp.scoreDiversity /*finalScore*/ << " " << tmpp.finalRank << " "
         << tmpp.id << " " << tmpp.getType() << "\n";
    sout << tmpp.degree.size() << "\n";
    for (int j = 0; j < tmpp.degree.size(); j++) sout << tmpp.degree[j] << " ";
    sout << "\n\n";
  }
  sout.close();
}
void FastPGGraph::saveQueries(string filename, int testFromSize, int testToSize,
                              int timesEach) {
  ofstream sout(filename);
  for (int i = testFromSize; i <= testToSize; ++i) {
    assert(queryPatterns[i].size() == timesEach);
    for (int j = 0; j < timesEach; ++j) {
      Pattern &tmpp = queryPatterns[i][j];
      sout << tmpp.nodesCnt << " " << tmpp.getEdgeCnt() << "\n";
      for (int k = 0; k < tmpp.getEdgeCnt(); k++)
        sout << tmpp.edges[k].s << " " << tmpp.edges[k].e << "\n";
    }
    sout << "\n";
  }

  //	for (int i = testFromSize; i <= testToSize; ++i)
  //	{
  //		assert(queryPatterns[i].size() == timesEach);
  //		for (int j = 0; j < queryPatterns[i].size(); ++j)
  //		{
  //			sout << Step4curQuery[i][j] << " ";
  //		}
  //		sout << "\n";
  //	}
  //	sout << "\n";
  sout.close();
}

void FastPGGraph::loadFinalpatterns(string fname, int num) {
  assert(30 == num || 0 == num || 15 == num);
  printf("load final patterns...\n");
  ifstream lin(fname);
  patterns.clear();
  finalPatternIdSet.clear();
  boostPatterns.clear();
  int numv, nume, a, b, degreeSize;
  int i = 0;
  while (lin >> numv >> nume) {
    finalPatternIdSet.emplace_back(i);
    Pattern tmpp(numv);
    while (nume--) {
      lin >> a >> b;
      tmpp.addEdge(a, b);
    }
    lin >> tmpp.scoreCognition;
    lin >> tmpp.scoreCoverage;
    // lin >> tmpp.finalScore;
    lin >> tmpp.scoreDiversity;
    lin >> tmpp.finalRank;
    lin >> tmpp.id;
    lin >> tmpp.type;
    lin >> degreeSize;
    while (degreeSize--) {
      lin >> a;
      tmpp.degree.emplace_back(a);
    }
    patterns.emplace_back(tmpp);
    boostPatterns.emplace_back(pattern2boostGraph(patterns[i]));
    ++i;
  }
  assert(num == patterns.size() || 0 == num);
  // re-assign 4 default coverage
  for (int i = 0; i < 4; ++i) patterns[i].scoreCoverage = 1.0;
  // re-assign cognition
  for (int i = 0; i < patterns.size(); ++i) calCognition6Single(patterns[i]);
  // re-assign diversity
  vector<vector<double>> netSimileSignature(patterns.size());
  for (int j = 0; j < patterns.size(); ++j)
    netSimileSignature[j] = net_simile::Signature(patterns[j]);
  for (int i = 0; i < patterns.size(); ++i) {
    double sim = 0.0;
    for (int j = 0; j < patterns.size(); ++j)
      if (i != j) {
        double tmp = net_simile::CamberraDistance(netSimileSignature[i],
                                                  netSimileSignature[j]);
        if (tmp > sim) sim = tmp;
      }
    patterns[i].scoreDiversity = 1.0 / sim;
  }
}

void FastPGGraph::loadQueries(string filename, int testFromSize, int testToSize,
                              int timesEach) {
  printf("load querys and steps...\n");
  ifstream lin(filename);
  queryPatterns.clear();
  queryPatterns.resize(testToSize + 1);
  int numv, nume, a, b;
  for (int i = testFromSize; i <= testToSize; ++i) {
    for (int j = 0; j < timesEach; ++j) {
      lin >> numv >> nume;
      Pattern tmpp(numv);
      while (nume--) {
        lin >> a >> b;
        tmpp.addEdge(a, b);
      }
      queryPatterns[i].emplace_back(tmpp);
    }
    cout << queryPatterns[i].size() << "\n";
  }

  //	Step4curQuery.clear();
  //	Step4curQuery.resize(testToSize + 1);
  //	for (int i = testFromSize; i <= testToSize; ++i)
  //	{
  //		Step4curQuery[i].resize(timesEach);
  //		assert(queryPatterns[i].size() == timesEach);
  //		for (int j = 0; j < queryPatterns[i].size(); ++j)
  //		{
  //		   lin >> Step4curQuery[i][j];
  //		}
  //	}
  lin.close();
}

void FastPGGraph::saveQueries2(string filename, int num) {
  printf("[func] saveQueries2\n");
  ofstream sout(filename);
  assert(query2Patterns.size() == num);
  for (const auto &p : query2Patterns) {
    sout << p.nodesCnt << " " << p.getEdgeCnt() << "\n";
    for (int k = 0; k < p.getEdgeCnt(); k++)
      sout << p.edges[k].s << " " << p.edges[k].e << "\n";
  }
  sout.close();
}

void FastPGGraph::loadQueries2(string filename, int num) {
  printf("[func] loadQueries2\n");
  ifstream lin(filename);
  query2Patterns.clear();
  int numv, nume, a, b;
  for (int i = 0; i < num; ++i) {
    lin >> numv >> nume;
    Pattern p(numv);
    while (nume--) {
      lin >> a >> b;
      p.addEdge(a, b);
    }
    query2Patterns.emplace_back(p);
  }
  lin.close();
}

// [s, t)
int getRand(int s, int t) {
  // srand(time(NULL));
  return s +
         (rand() * (long long)rand() + rand() * (long long)rand()) % (t - s);
}

void FastPGGraph::genQueries2_genRandom(int num) {
  printf("[func] genQueries2_genRandom\n");
  assert(num % 25 == 0);
  int each_num = num / 25;
  for (int i = 5; i < 30; ++i)
    for (int j = 0; j < each_num; ++j)
      query2Patterns.emplace_back(queryPatterns[i][j]);
}

void FastPGGraph::genQueries2_genPath(int num) {
  printf("[func] genQueries2_genPath\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 1);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 1, vector<bool>(size + 1, false));
      int cur_node = 0;
      int x = getRand(0, NUM_V);
      printf("size=%d, j=%d, x=%d\n", size, j, x);
      nodeSet.emplace_back(x);
      nodeMap[x] = cur_node++;
      vis[x] = ++vis_flag;
      int cur_size = 0;
      int tri = 0;
      while (cur_size + 1 < size && cur_size < 15 && tri < 10) {
        if (neibours[x].size() == 0) break;
        int next = neibours[x][getRand(0, neibours[x].size())].neighbId;
        if (vis[next] != vis_flag) {
          nodeSet.emplace_back(next);
          nodeMap[next] = cur_node++;
          hasEdge[nodeMap[x]][nodeMap[next]] =
              hasEdge[nodeMap[next]][nodeMap[x]] = true;
          pattern.addEdge(nodeMap[x], nodeMap[next]);
          vis[next] = vis_flag;
          ++cur_size;
          tri = 0;
          x = next;
        } else
          ++tri;
      }

      if (cur_size >= size / 2.0) {
        genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                           nodeMap, hasEdge, pattern);
        pattern.nodesCnt = cur_node;
        query2Patterns.emplace_back(pattern);
      } else {
        --j;
      }
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_genTree(int num) {
  printf("[func] genQueries2_genTree\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 1);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 1, vector<bool>(size + 1, false));
      int cur_node = 0;
      int x = getRand(0, NUM_V);
      printf("size=%d, j=%d, x=%d\n", size, j, x);
      nodeSet.emplace_back(x);
      nodeMap[x] = cur_node++;
      vis[x] = ++vis_flag;
      int cur_size = 0;
      int tri = 0;
      while (cur_size + 1 < size && tri < 10) {
        x = nodeSet[getRand(0, nodeSet.size())];
        if (neibours[x].size() == 0) {
          ++tri;
          continue;
        }
        int next = neibours[x][getRand(0, neibours[x].size())].neighbId;
        if (vis[next] != vis_flag) {
          nodeSet.emplace_back(next);
          nodeMap[next] = cur_node++;
          hasEdge[nodeMap[x]][nodeMap[next]] =
              hasEdge[nodeMap[next]][nodeMap[x]] = true;
          pattern.addEdge(nodeMap[x], nodeMap[next]);
          vis[next] = vis_flag;
          ++cur_size;
          tri = 0;
        } else
          ++tri;
      }

      if (cur_size >= size / 2.0) {
        genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                           nodeMap, hasEdge, pattern);
        pattern.nodesCnt = cur_node;
        query2Patterns.emplace_back(pattern);
      } else {
        --j;
      }
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_genSingleStar(int num) {
  printf("[func] genQueries2_genSingleStar\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  int roadnet_tri = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 1);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 1, vector<bool>(size + 1, false));
      int cur_node = 0;
      int x = getRand(0, NUM_V);
      printf("size=%d, j=%d, x=%d\n", size, j, x);
      nodeSet.emplace_back(x);
      nodeMap[x] = cur_node++;
      vis[x] = ++vis_flag;
      int cur_size = 0;
      int tri = 0;
      while (cur_size + 1 < size - 1 && tri < 10) {  // -1 ==> not pure star
        if (neibours[x].size() == 0) {
          ++tri;
          continue;
        }
        int next = neibours[x][getRand(0, neibours[x].size())].neighbId;
        if (vis[next] != vis_flag) {
          nodeSet.emplace_back(next);
          nodeMap[next] = cur_node++;
          hasEdge[nodeMap[x]][nodeMap[next]] =
              hasEdge[nodeMap[next]][nodeMap[x]] = true;
          pattern.addEdge(nodeMap[x], nodeMap[next]);
          vis[next] = vis_flag;
          ++cur_size;
          tri = 0;
        } else
          ++tri;
      }

      if (cur_size >= size / 2.0 || roadnet_tri > 10) {
        genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                           nodeMap, hasEdge, pattern);
        pattern.nodesCnt = cur_node;
        query2Patterns.emplace_back(pattern);
        roadnet_tri = 0;
      } else {
        --j;
        ++roadnet_tri;
      }
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_genDoubleStar(int num) {
  printf("[func] genQueries2_genDoubleStar\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  int roadnet_tri = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 1);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 1, vector<bool>(size + 1, false));
      int cur_node = 0;
      int x = getRand(0, NUM_V);
      printf("size=%d, j=%d, x=%d\n", size, j, x);
      nodeSet.emplace_back(x);
      nodeMap[x] = cur_node++;
      vis[x] = ++vis_flag;
      int cur_size = 0;
      int tri = 0;
      bool second_core = false;
      while (cur_size + 1 < size && tri < 10) {  // -1 ==> not pure star
        if (neibours[x].size() == 0) {
          ++tri;
          continue;
        }
        int next = neibours[x][getRand(0, neibours[x].size())].neighbId;
        if (vis[next] != vis_flag) {
          nodeSet.emplace_back(next);
          nodeMap[next] = cur_node++;
          hasEdge[nodeMap[x]][nodeMap[next]] =
              hasEdge[nodeMap[next]][nodeMap[x]] = true;
          pattern.addEdge(nodeMap[x], nodeMap[next]);
          vis[next] = vis_flag;
          ++cur_size;
          tri = 0;
          if (cur_size >= size / 2.0 && !second_core && rand() % 100 >= 60) {
            second_core = true;
            x = next;
          }
        } else
          ++tri;
      }

      if (cur_size >= size / 2.0 || roadnet_tri > 10) {
        genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                           nodeMap, hasEdge, pattern);
        pattern.nodesCnt = cur_node;
        query2Patterns.emplace_back(pattern);
        roadnet_tri = 0;
      } else {
        --j;
        ++roadnet_tri;
      }
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_genCycle(int num) {
  printf("[func] genQueries2_genCycle\n");
  printf("circleRecord size=%lu\n", circleRecord.size());
  for (const auto &r : circleRecord) printf("%lu,", r.size());
  printf("\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 1);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 1, vector<bool>(size + 1, false));
      int cur_node = 0, cur_size = 0;
      int ind = getRand(0, circleRecord.size());
      int circle_length = circleRecord[ind].size();
      printf("size=%d, j=%d, circle_length=%d\n", size, j, circle_length);
      if (circle_length > size) {
        --j;
        continue;
      }

      ++vis_flag;
      for (int i = 0; i < circle_length; ++i) {
        int x = circleRecord[ind][i];
        nodeSet.emplace_back(x);
        nodeMap[x] = cur_node++;
        vis[x] = vis_flag;
      }
      for (int i = 0; i < circle_length; ++i) {
        int j = (i + 1) % circle_length;
        int x = circleRecord[ind][i], next = circleRecord[ind][j];
        hasEdge[nodeMap[x]][nodeMap[next]] =
            hasEdge[nodeMap[next]][nodeMap[x]] = true;
        pattern.addEdge(nodeMap[x], nodeMap[next]);
        ++cur_size;
      }
      assert(cur_node == cur_size && cur_node == circle_length);

      int x = circleRecord[ind][getRand(0, circleRecord[ind].size())];
      genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                         nodeMap, hasEdge, pattern);
      pattern.nodesCnt = cur_node;
      query2Patterns.emplace_back(pattern);
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_genFlower(int num) {
  printf("[func] genQueries2_genFlower\n");
  assert(num % 25 == 0);
  int before_num = query2Patterns.size();
  int each_num = num / 25;
  vector<int> vis(NUM_V, -1);
  int vis_flag = 0;
  int roadnet_tri = 0;
  for (int size = 5; size < 30; ++size)
    for (int j = 0; j < each_num; ++j) {
      Pattern pattern(size + 10);  // not real size
      vector<int> nodeSet;
      nodeSet.clear();
      unordered_map<int, int> nodeMap;
      nodeMap.clear();
      vector<vector<bool>> hasEdge(size + 10, vector<bool>(size + 10, false));
      int cur_node = 0;
      int x = getRand(0, NUM_V);
      printf("size=%d, j=%d, x=%d\n", size, j, x);
      if (neibours[x].size() == 0) {
        --j;
        continue;
      }
      nodeSet.emplace_back(x);
      nodeMap[x] = cur_node++;
      vis[x] = ++vis_flag;
      int cur_size = 0;
      int tri = 0;
      int k = -1;
      while (k == -1 && tri < 10) {
        auto neibour = neibours[x][getRand(0, neibours[x].size())];
        if (trussClass[neibour.edgeId] < 4) {
          ++tri;
          continue;
        }
        int next = neibour.neighbId;
        nodeSet.emplace_back(next);
        nodeMap[next] = cur_node++;
        assert(2 == cur_node && 0 == nodeMap[x] && 1 == nodeMap[next]);
        hasEdge[nodeMap[x]][nodeMap[next]] =
            hasEdge[nodeMap[next]][nodeMap[x]] = true;  // mocked
        // pattern.addEdge(nodeMap[x], nodeMap[next]);
        vis[next] = vis_flag;
        k = getRand(2, min((size + 1) / 2, trussClass[neibour.edgeId] -
                                               1));  // number of mocked edge
        cur_size += k * 2;                           // mocked
      }

      if (k == -1) {
        --j;
        continue;
      }
      genQueries2_extend(size - cur_size, vis, vis_flag, cur_node, nodeSet,
                         nodeMap, hasEdge, pattern, true);
      for (int i = 0; i < k; ++i) {
        pattern.addEdge(0, cur_node);
        pattern.addEdge(1, cur_node);
        ++cur_node;
      }
      if (cur_size >= size / 2 || roadnet_tri > 10) {
        pattern.nodesCnt = cur_node;
        query2Patterns.emplace_back(pattern);
        roadnet_tri = 0;
      } else {
        --j;
        ++roadnet_tri;
        continue;
      }
    }
  assert(before_num + num == query2Patterns.size());
}

void FastPGGraph::genQueries2_extend(int num, vector<int> &vis, int vis_flag,
                                     int &cur_node, vector<int> &nodeSet,
                                     unordered_map<int, int> &nodeMap,
                                     vector<vector<bool>> &hasEdge,
                                     Pattern &pattern, bool gout_only) {
  int tri = 0;
  int s, e;
  auto existsE = gout_only ? existsGoutEdge : existsEdge;
  while (num > 0 && tri < 10) {
    assert(cur_node == nodeSet.size());
    bool find = false;
    if (rand() % 100 < 80) {  // 80% inner edge
      for (int i = 0; i < cur_node && !find; ++i)
        for (int j = i + 1; j < cur_node && !find; ++j)
          if (!hasEdge[i][j] && existsE(nodeSet[i], nodeSet[j])) {
            s = nodeSet[i], e = nodeSet[j];
            assert(i == nodeMap[s] && j == nodeMap[e]);
            find = true;
          }
    } else {
      s = nodeSet[getRand(0, nodeSet.size())];
      if (!neibours[s].empty()) {
        e = neibours[s][getRand(0, neibours[s].size())].neighbId;
        if (vis[e] != vis_flag) {
          // node
          nodeSet.emplace_back(e);
          nodeMap[e] = cur_node++;
          vis[e] = vis_flag;
          find = true;
        }
      }
    }
    if (find) {
      // edge
      hasEdge[nodeMap[s]][nodeMap[e]] = hasEdge[nodeMap[e]][nodeMap[s]] = true;
      pattern.addEdge(nodeMap[s], nodeMap[e]);
      --num;
      tri = 0;
    } else {
      ++tri;
    }
  }
}

bool FastPGGraph::existsEdge(int s, int t) {
  for (const auto &pair : neibours[s])
    if (t == pair.neighbId) return true;
  return false;
}

bool FastPGGraph::existsGoutEdge(int s, int t) {
  for (const auto &pair : neibours[s])
    if (trussClass[pair.edgeId] < 3 && t == pair.neighbId) return true;
  return false;
}

void FastPGGraph::genQueries2(string query_filename, string filename) {
  printf("[func] genQueries2\n");
  loadQueries(query_filename, 4, 30, 100);
  assert(queryPatterns[4].size() == 100);
  query2Patterns.clear();
  genQueries2_genRandom(query2Setting[0]);
  genQueries2_genPath(query2Setting[1]);
  genQueries2_genTree(query2Setting[2]);
  genQueries2_genSingleStar(query2Setting[3]);
  genQueries2_genDoubleStar(query2Setting[4]);
  genQueries2_genCycle(query2Setting[5]);
  genQueries2_genFlower(query2Setting[6]);
  assert(query2Patterns.size() == query2Setting[7]);
  saveQueries2(filename, query2Setting[7]);
}

/**
 * search k-truss patterns
 * for test only.
 **/
void FastPGGraph::searchKTrussPatterns(int maxEdgeCnt) {
  printf("searching k-truss patterns for test...\n");
  vector<vector<vector<graph_t>>> patternSet = vector<vector<vector<graph_t>>>(
      MAX_INDIVIDUAL_NODE_COUNT + 1, vector<vector<graph_t>>(maxEdgeCnt + 1));
  // default patterns
  for (auto &p : patterns)
    patternSet[p.getNodesCnt()][p.getEdgeCnt()].emplace_back(
        pattern2boostGraph(p));

  vector<char> flag(NUM_V, 0);
  queue<int> q;
  int nodeCnt, edgeCnt;
  vector<int> edgeSet;
  for (int k = 4; k <= 9; ++k) {
    // k-clique
    if (k * (k - 1) / 2 <= maxEdgeCnt) {
      Pattern pat(k);
      graph_t g(k);
      for (int i = 0; i + 1 < k; ++i)
        for (int j = i + 1; j < k; ++j) {
          pat.addEdge(i, j);
          add_edge(i, j, g);
        }
      patterns.emplace_back(pat);
      patternSet[k][k * (k - 1) / 2].emplace_back(g);
    }
  }

  for (int k = 4; k <= 9; ++k) {
    for (int i = 0; i < NUM_V; ++i)
      if (flag[i] != k) {
        edgeSet.clear();
        flag[i] = k;
        nodeCnt = 1;
        edgeCnt = 0;
        q.push(i);
        while (!q.empty()) {
          int cur = q.front();
          q.pop();
          for (auto &ve : neibours[cur])
            if (trussClass[ve.edgeId] >= k) {
              ++edgeCnt;
              edgeSet.emplace_back(ve.edgeId);
              if (flag[ve.neighbId] != k) {
                flag[ve.neighbId] = k;
                ++nodeCnt;
                q.push(ve.neighbId);
              }
            }
        }
        assert((edgeCnt & 1) == 0);
        edgeCnt >>= 1;
        if (nodeCnt > 3 && nodeCnt <= MAX_INDIVIDUAL_NODE_COUNT &&
            edgeCnt <= maxEdgeCnt) {
          int c = 0;
          unordered_map<int, int> nodeMap;
          graph_t g(nodeCnt);
          for (auto edgeId : edgeSet) {
            int s = edges[edgeId].s;
            auto itr = nodeMap.find(s);
            if (itr == nodeMap.end())
              itr = nodeMap.insert(make_pair(s, c++)).first;
            s = itr->second;
            int e = edges[edgeId].e;
            itr = nodeMap.find(e);
            if (itr == nodeMap.end())
              itr = nodeMap.emplace(make_pair(e, c++)).first;
            e = itr->second;
            if (s != e) add_edge(s, e, g);
          }
          assert(c == nodeCnt);
          bool found = false;
          std::vector<graph_t::vertex_descriptor> f(nodeCnt);
          for (auto &bp : patternSet[nodeCnt][edgeCnt])
            if (isomorphism(
                    g, bp,
                    boost::isomorphism_map(boost::make_iterator_property_map(
                        f.begin(), boost::get(boost::vertex_index, g))))) {
              found = true;
              break;
            }
          if (!found) {
            Pattern pat(nodeCnt);
            for (auto ep = boost::edges(g); ep.first != ep.second; ++ep.first)
              pat.addEdge(boost::source(*ep.first, g),
                          boost::target(*ep.first, g));
            patterns.emplace_back(pat);
            patternSet[nodeCnt][edgeCnt].emplace_back(g);
          }
        }
      }
  }
}

void FastPGGraph::outputLocalPatterns4UI(const vector<Pattern> &patterns,
                                         const string &filename) {
  ofstream fout(filename);
  for (int i = 0; i < patterns.size(); ++i) {
    auto &p = patterns[i];
    int n = p.getNodesCnt(), m = p.getEdgeCnt();
    fout << "t # " << i << " " << n << "\n";
    for (int j = 0; j < n; ++j) fout << "v " << j << " 0\n";
    for (int j = 0; j < m; ++j)
      fout << "e " << p.edges[j].s << " " << p.edges[j].e << " 0\n";
    fout << "\n";
  }
  fout.close();
}

void FastPGGraph::loadPeregrine(string file) {
  patterns.clear();
  ifstream fin(file);
  // ofstream fout("clean_peregrine.txt");
  string line;
  vector<double> max_support(10, 0);
  while (getline(fin, line)) {
    int len = line.length();
    if (line[0] != '[') continue;

    int n = 0;
    vector<pair<int, int>> edges;
    int ea, eb;
    double support;
    bool edge_start = false;
    bool support_start = false;
    string tmp;
    int i;
    for (i = 0; i < len; ++i)
      if ('[' == line[i]) {
        edge_start = true;
        tmp = "";
      } else if ('-' == line[i] && edge_start) {
        ea = stoi(tmp);
        if (ea > n) n = ea;
        tmp = "";
      } else if (']' == line[i] && edge_start) {
        eb = stoi(tmp);
        if (eb > n) n = eb;
        edges.emplace_back(ea, eb);
        tmp = "";
        edge_start = false;
      } else if (':' == line[i]) {
        support_start = true;
        tmp = "";
      } else if ('0' <= line[i] && line[i] <= '9')
        tmp += line[i];
      else if (support_start && !tmp.empty())
        break;
    support = stod(tmp);
    if (n <= 4 || edges.size() <= 4) continue;
    if (support > max_support[n]) max_support[n] = support;
    // fout << line.substr(0, i - 1) << endl;

    Pattern p = Pattern(n);
    for (const auto &e : edges) p.addEdge(e.first - 1, e.second - 1);
    p.setType(PatType::OutsideIndividual);
    p.scoreCoverage = support;
    calCognition6Single(p);
    // if (p.nodesCnt > 4 || p.getEdgeCnt() > 4)
    patterns.emplace_back(p);
  }
  for (auto &p : patterns) {
    p.scoreCoverage /= max_support[p.nodesCnt];
    // printf("|V|=%d |E|=%d cov=%lf cog=%lf\n", p.nodesCnt, p.getEdgeCnt(),
    // p.scoreCoverage, p.scoreCognition);
  }
  genDefaultPatterns();
  for (int i = patterns.size() - 4; i < patterns.size(); i++) {
    patterns[i].scoreCoverage = 1.0;
    calCognition6Single(patterns[i]);
  }
  printf("size=%lu\n", patterns.size());
  fin.close();
  // fout.close();
}

void FastPGGraph::collectDiversityScores(const vector<Pattern> &patterns,
                                         const vector<int> &ids, double *avgDiv,
                                         double *minDiv, double *maxDiv) {
  vector<vector<double>> netSimileSignature(patterns.size());
  for (int j = 0; j < patterns.size(); ++j) {
    netSimileSignature[j] = net_simile::Signature(patterns[j]);
    /*printf("j=%d", j);
    for (int k = 0; k < netSimileSignature[j].size(); ++k)
            printf(" %lf", netSimileSignature[j][k]);
    printf("\n");*/
  }
  *avgDiv = 0.0;
  *minDiv = 1e10, *maxDiv = -1;
  for (int i = 0; i < ids.size(); ++i) {
    double sim = 0.0;
    for (int j = 0; j < ids.size(); ++j)
      if (i != j) {
        double tmp = net_simile::CamberraDistance(netSimileSignature[ids[i]],
                                                  netSimileSignature[ids[j]]);
        // printf("sim[%d][%d]=%lf\n", i, j, tmp);
        if (tmp > sim) sim = tmp;
      }
    double div = 1.0 / sim;
    *avgDiv += div;
    if (div < *minDiv) *minDiv = div;
    if (div > *maxDiv) *maxDiv = div;
  }
  *avgDiv /= ids.size();
}