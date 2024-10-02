#ifndef BASE_H_
#define BASE_H_

#include <unordered_set>
#include <vector>

using namespace std;

constexpr int MAX_COMBINED_TRUSS_CLASS = 9;
constexpr int MAX_SUPPORT = 1e6;  // edge maximum support
constexpr int MAX_INDIVIDUAL_NODE_COUNT = 30;
constexpr int MAX_SMALL_PATTERN_EDGE_COUNT = 50;  // 450;

constexpr int THREAD_NUM = 6;  // thread number

class Edge {
 public:
  int s, e;
  int support;
  int id;
  bool alive;
  int pos;
  bool deleted;
  int trussClass;
  int originSupport;
  Edge(int edgeId, int i, int j) {
    id = edgeId;
    s = i;
    e = j;
    alive = true;
    deleted = false;
    trussClass = -1;
  }
};

enum PatType {
  OneEdge = 0,
  TwoEdge = 1,
  Triangle = 2,
  Truss = 3,           // truss-generated patterns
  TrussExtend = 4,     // now unused
  TrussCombined1 = 5,  // k1-side-k2-main
  TrussCombined2 = 6,  // k1-side-k2-side
  InRemain = 7,        // outside truss patterns
  OutsideStar = 8,
  OutsideCombinedStar = 9,
  OutsideLine = 10,
  OutsideCircle = 11,
  OutsideIndividual = 12
};

class StarState {
 public:
  vector<int> nodeId;
  vector<short> degree;
  int cursize;
  unordered_set<int> myv;
  StarState() {
    nodeId.clear();
    degree.clear();
    myv.clear();
    cursize = 0;
  }
};

class PatEdge {
 public:
  int s, e;
  PatEdge(int s, int e) {
    this->s = s;
    this->e = e;
  }
};

class Pattern {
 public:
  static int curId;
  int score = 0;
  int scorePart[10];
  int type;
  int id;
  vector<PatEdge> edges;
  int nodesCnt;
  int trussClass;
  double scoreCoverage, scoreCognition, scoreDiversity;
  double finalScore;
  pair<int, int> detilType;
  double dynamicScore;
  int historyMaxComE;
  int finalRank;
  vector<vector<pair<int, int>>> edgesList;  // s --> [(e, edge_id)]
  vector<short> degree;  // now only use for linear star search
  int outsideCoverage;

  Pattern(int n) {
    this->nodesCnt = n;
    edgesList.resize(n);
    for (int i = 0; i < 10; ++i) scorePart[i] = 0;
    this->id = curId++;
    scoreCoverage = scoreCognition = scoreDiversity = 0.0;
    finalScore = 0.0;
    dynamicScore = 0x3f3f3f3f;
    historyMaxComE = 0;
    // degree.resize(n, 0);
    degree.clear();
  }

  void setType(int type) { this->type = type; }
  void addEdge(int s, int e) {
    edgesList[s].emplace_back(make_pair(e, edges.size()));
    edgesList[e].emplace_back(make_pair(s, edges.size()));
    edges.emplace_back(s, e);
  }
  int getEdgeId(int s, int e) const {
    for (auto p : edgesList[s])
      if (p.first == e) return p.second;
    return -1;
  }
  void setScore(int score) { this->score = score; }
  void setTrussClass(int trussClass) { this->trussClass = trussClass; }
  int getNodesCnt() const { return nodesCnt; }
  int getEdgeCnt() const { return edges.size(); }
  int getScore() const { return score; }
  int getType() const { return type; }
  int getTrussClass() const { return trussClass; }
};

enum EdgeType {
  RealTruss = 0,  // in Gin
  Line = 1,
  Circle = 2,
  SimpleStar = 3,
  CombinedStar = 4,
  Individual = 5,
  Default = -1
};

struct NeighbVE {
  int neighbId;
  int edgeId;
  NeighbVE(int v, int e) {
    neighbId = v;
    edgeId = e;
  }
};

class PGEdge {
 public:
  int e, s;
  PGEdge(int i, int j) {
    s = i;
    e = j;
  }
};

#endif  // BASE_H_