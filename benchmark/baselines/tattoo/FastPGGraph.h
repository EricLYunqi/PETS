
#ifndef FASTPGGRAPH_H_
#define FASTPGGRAPH_H_

#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/mcgregor_common_subgraphs.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/thread/thread.hpp>
#include <cmath>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>

#include "base.h"
#include "utils/NetSimile.h"
#include "utils/json.hpp"

using namespace std;
using json = nlohmann::json;
#define ull unsigned long long
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS>
    graph_t;

struct SummaryInfo {
 public:
  int insideEdgeCnt, outsideEdgeCnt;
  int starsEdgeCnt, linesEdgeCnt, circlesEdgeCnt, individualsEdgeCnt,
      uncoveredEdgeCnt;
  int starsComponentCnt, linesComponentCnt, circlesComponentCnt,
      individualsComponentCnt, uncoveredComponentCnt;
  double readGraphTime, calTrussClass4EdgesTime,
      calScore4CombinedTrussPatternsTime, searchStarsTime,
      searchLinesAndCirclesAndIndividualsTime, calCovAndCogScoreTime,
      finalPatternsGenTime;

  SummaryInfo() { init(); }
  void init() {
    insideEdgeCnt = outsideEdgeCnt = 0;
    starsEdgeCnt = linesEdgeCnt = circlesEdgeCnt = individualsEdgeCnt =
        uncoveredEdgeCnt = 0;
    starsComponentCnt = linesComponentCnt = circlesComponentCnt =
        individualsComponentCnt = uncoveredComponentCnt = 0;
    readGraphTime = calTrussClass4EdgesTime =
        calScore4CombinedTrussPatternsTime = searchStarsTime =
            searchLinesAndCirclesAndIndividualsTime = calCovAndCogScoreTime =
                finalPatternsGenTime = 0;
  }
};

class SmallComponentInfo {
 public:
  graph_t graph;
  int score;
  int patId;

  SmallComponentInfo(const graph_t &graph, const int &score, const int &patId) {
    this->graph = graph;
    this->score = score;
    this->patId = patId;
  }
};

class FastPGGraph {
 public:
  // configurable params
  static int MAX_PATTERN_SIZE;
  static int MIN_STAR_DEGREE;

  static vector<int> trussEdgeCnt;  //* define the number of edges of truss[i].
                                    //(remember: truss[i+1] is also truss[i])
  static vector<PGEdge> edges;
  static vector<vector<int>> NNH;
  static vector<Pattern> patterns;
  static vector<graph_t> boostPatterns;
  static vector<vector<vector<int>>> patIdMap;
  static unordered_map<int, int> outsideStarPatMap;
  static map<pair<int, int>, int> outsideCombinedStarPatMap;
  static unordered_map<int, int> outsideLinePatMap;
  static unordered_map<int, int> outsideCirclePatMap;
  static vector<char> edgeType;
  static vector<int> binStart;
  static vector<int> binSize;
  static vector<int> curLength;
  static vector<int> sortedEdges;
  static vector<int> sortedPos;
  static vector<bool> edgeAlive;
  static vector<vector<vector<SmallComponentInfo>>> smallComponents;
  static vector<vector<NeighbVE>> neibours;
  static int NUM_V, NUM_E;
  static vector<int> GoutDegree;
  static vector<int> support;
  static vector<char> trussClass;
  static vector<int> finalPatternIdSet;
  static queue<int> candidatePatternsQueue;
  static vector<vector<Pattern>> queryPatterns;
  // query2
  static vector<vector<int>> circleRecord;
  static vector<Pattern> query2Patterns;
  static vector<int> query2Setting;
  static vector<string> query2SettingName;
  static vector<SummaryInfo> newInfo;
  static vector<vector<Pattern>> finalSet4Datastores;
  static vector<int> tempfinalset;
  static vector<int> patternUseClassCnt;
  static unordered_map<size_t, int> outsideLinearCombinedStarPatMap;
  static vector<vector<int>> Step4curQuery;
  static bool LOADED;

  FastPGGraph() {
    init();
  }
  
  FastPGGraph(string inputFilename) {
    double startTime = time(nullptr);
    init();
    // loadGraphFromEdgesFile(inputFilename);
    loadGraphFromFileAndMap(inputFilename);
    summaryInfo.readGraphTime = time(nullptr) - startTime;
  }
  static void genDefaultPatterns();
  static void outputGin4Gephi(string filename);
  static void outputGout4Gephi(string filename);
  static void outputGr4Gephi(string filename);
  static void outputOutsideEdgeInfo(string filename);
  static void outputOutsideEdgeInfoRemoved(string filename);
  static void outputOutsideEdgeInfoBesideStarsAndLines(string filename);
  static void outputSummaryInfo(string filename);
  static void outputSummaryInfo(string filename, string datasetName,
                                bool append);
  static void outputPatterns(string filename);
  static void outputGraph4GUI(string filename);
  static void outputFinalPatterns4GUI(string filename);
  static void outputPngPattern(Pattern &pattern, const string &folder);
  static void outputFinalPngPattern(Pattern &pattern, const string &folder,
                                    string filename);
  static void outputPngPatterns(string filename);
  static void genTrussPatterns();
  static void genCombinedTrussPatterns();
  static void preProceEdge();
  static void calScore4CombinedTrussPatterns();
  static void genFourMainThreeSidePattern();
  static int getEdgebyNode(int s, int e);

  static void searchStars();

  static void searchLinesAndCirclesAndIndividuals();

  static int filterPatterns(const Pattern &p);

  static void genQueryRandomly(string filename, int times, int size,
                               bool gin_only = false);
  static void genSubgraph(int size, string filename);

  static graph_t pattern2boostGraph(const Pattern &p);

  static vector<int> preJudge4MaxComSubGraph(Pattern &pattern1,
                                             Pattern &pattern2);
  static vector<int> patternCommonSubgraph(int p1, int p2);
  static bool isSubGraph(int p1, int p2);
  static bool isSubGraph(const Pattern &p1, const Pattern &p2);
  static unordered_set<unsigned long long> getAllSubGraphEdgeSets(int p1,
                                                                  int p2);
  static unordered_set<unsigned long long> getAllSubGraphEdgeSets(
      int p1, const Pattern &pattern2);
  static unordered_set<unsigned long long> getAllSubGraphEdgeSets(
      const Pattern &pattern1, const graph_t &graph1, const Pattern &pattern2,
      const graph_t &graph2);

  static void patternChooseByRedundancy(int size, string filename);
  static void patternChooseByGreedy(int size, string filename);

  static void selectFinalPatterns(string filename, double lambda, int size,
                                  int datasetid);
  static void selectFinalPatternsByCov(int size);

  static void calCoverageAndCognition(const string &folder);

  static void calCoverage1(const string &folder);
  static void calCoverage2(const string &folder);
  static void calCoverage3(const string &folder);
  static void calCognition1(const string &folder);
  static void calCognition2(const string &folder);
  static void calCognition2Single(Pattern &pattern);
  static void calCognition3(const string &folder);
  static void calCognition3Single(Pattern &pattern);
  static void calCognition4(const string &folder);
  static void calCognition4Single(Pattern &pattern);
  static void calCognition5(const string &folder);
  static void calCognition5Single(Pattern &pattern);
  static void calCognition6(const string &folder);
  static void calCognition6Single(Pattern &pattern);

  static void setCoverSolve(vector<ull> &sets, const ull &aimset, int &minstep,
                            bool test = false);

  static void setCoverSolve(vector<ull> &sets, const ull &aimset, ull status,
                            int curstep, int &minstep);

  static double finalTest(string fileName, int testFromSize, int testToSize,
                          int timesEach);

  static void genQuery(int testFromSize, int testToSize, int timesEach,
                       bool gin_only = false);

  static double paratest(ofstream &out, int testFromSize, int testToSize,
                         int timesEach, vector<int> finalIdSet, int mark,
                         bool printIndividualRR = false);
  static void paratest2(ofstream &out, const vector<int> &finalIdSet);
  static void testByPatternClass(ofstream &out, int testFromSize,
                                 int testToSize, int timesEach);
  static void getcomset4alldataset();

  static void predo4test();

  static void choosefinalpatternsbycov(int topk);

  static void print3ResultsOf4Groups(const string &folder);

  static void genGraphlet();

  static void graphletTest(int topk);
  static void peregrineTest(int size);
  static void randomTest(int topk);

  static void prejob4finalPatternSetSize();

  static void Test4finalPatternSetSize(ofstream &myout, int topk);

  static void getEdgeneighbor(int &edgeId);
  static void decSupportForEdge(int edgeId);
  static void getNNH(int &edgeId);
  static void arrangeSortedEdges();
  static void decSupport(int s, int k);
  static void calTrussClass4Edges();

  static void calCandidatePatternsFinalScore(int curMinStep, double *maxps,
                                             int *maxp, int threadId);
  static void init();
  static void init_part1();

  static void searchMultiLinearStars();

  static void saveFinalpatterns(string fname, int num);
  static void loadFinalpatterns(string fname, int num);
  static void saveQueries(string filename, int testFromSize, int testToSize,
                          int timesEach);
  static void loadQueries(string filename, int testFromSize, int testToSize,
                          int timesEach);
  // query2
  static void saveQueries2(string filename, int num);
  static void loadQueries2(string filename, int num);
  static void genQueries2(string query_filename, string filename);
  static void genQueries2_genRandom(int num);
  static void genQueries2_genPath(int num);
  static void genQueries2_genTree(int num);
  static void genQueries2_genSingleStar(int num);
  static void genQueries2_genDoubleStar(int num);
  static void genQueries2_genCycle(int num);
  static void genQueries2_genFlower(int num);
  static void genQueries2_extend(int num, vector<int> &vis, int vis_flag,
                                 int &cur_node, vector<int> &nodeSet,
                                 unordered_map<int, int> &nodeMap,
                                 vector<vector<bool>> &hasEdge,
                                 Pattern &pattern, bool gout_only = false);
  static bool existsEdge(int s, int t);
  static bool existsGoutEdge(int s, int t);

  static void searchKTrussPatterns(int maxEdgeCnt);
  static void outputLocalPatterns4UI(const vector<Pattern> &patterns,
                                     const string &filename);
  template <typename Query>
  static int minStep4filledQuery(const vector<int> &patternsId, Query query,
                                 int edgeCnt, bool test = false);

  static void loadPeregrine(string file);
  static void collectDiversityScores(const vector<Pattern> &patterns,
                                     const vector<int> &ids, double *avgDiv,
                                     double *minDiv, double *maxDiv);

  static SummaryInfo summaryInfo;

 private:
  static boost::mutex io_mutex;

  static int getFinalPatternClass(int patternId);
  static void loadGraphFromFileAndMap(string filename);

  static void genTrussPattern(int k);
  static void genCombinedPattern(int k1, int k2, int type);
  static void getKTruss(int k);
  static void annoEdgeWithTrussClass(int maxK);
  static void addScoreForComPatternWithComEdge1(int k, int originalEdgeId,
                                                int edgeId,
                                                unordered_set<int> &nodeSet,
                                                int &Tid);
  static bool isSampledEdge(int edgeId);
  static void preProceE0(vector<int> &mybin, int &Tid);

  static bool isOK(int k1, int k2, int oriEdgeId, int edgeId,
                   unordered_set<int> &nodeSet, int para);
  static void addScoreForComPatternWithComEdge2(int k, int originalEdgeId,
                                                int edgeId,
                                                unordered_set<int> &nodeSet,
                                                int &sharedS, int &sharedE,
                                                int &Tid);
  static void mtAddScoreForComPattern(vector<int> &edgeIds, int &Tid);
  static void preProceE(vector<int> &bin, int &Tid);
};

#endif  // FASTPGGRAPH_H_
