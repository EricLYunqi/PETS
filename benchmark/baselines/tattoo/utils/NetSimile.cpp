#include "NetSimile.h"

#include <algorithm>
#include <cmath>
#include <iostream>
//#define DEBUG_MODE

// https://github.com/kristyspatel/Netsimile/blob/master/Netsimile.py
namespace net_simile {
// return 1 * 5
// code:
// http://basicodingfordummies.blogspot.com/2013/02/mean-median-variance-skewness-kurtosis.html
// code: https://www.johndcook.com/blog/skewness_kurtosis/
vector<double> AggregateFeatures(const vector<double> &a2) {
  int n = a2.size();
  auto a = a2;
  sort(a.begin(), a.end());
  auto sum = 0.0;
  for (auto v : a) sum += v;
  auto mean = sum / n;
  auto median = (n & 1) ? a[n / 2] : ((a[n / 2 - 1] + a[n / 2]) * 0.5);
  double variance = 0.0;
  for (auto v : a) variance += (v - mean) * (v - mean);
  variance /= n - 1;
#ifdef DEBUG_MODE
  cout << "Variance = " << variance << endl;
#endif
  auto std_deviation = sqrt(variance);
  double skewness = 0.0;
  if (std_deviation > 0.0) {
    for (auto v : a) skewness += (v - mean) * (v - mean) * (v - mean);
    skewness /= n * std_deviation * std_deviation * std_deviation;
  }
  double kurtosis = 0.0;
  if (std_deviation > 0.0) {
    for (auto v : a)
      kurtosis += (v - mean) * (v - mean) * (v - mean) * (v - mean);
    kurtosis /=
        n * std_deviation * std_deviation * std_deviation * std_deviation;
    kurtosis -= 3.0;
  }
  return {mean, median, std_deviation, skewness, kurtosis};
}

// return 1 * 35
vector<double> Signature(const Pattern &p) {
  int n = p.getNodesCnt();
  // (1) node degree
  vector<double> degree(n);
  for (int i = 0; i < n; ++i) degree[i] = p.edgesList[i].size();
  // (2) node coefficient
  vector<double> coefficient(n);
  for (int i = 0; i < n; ++i) {
    int triangle = 0, triple = 0;
    int edge_num = degree[i];
    for (int j = 0; j < edge_num; ++j)
      for (int k = j + 1; k < edge_num; ++k) {
        if (p.getEdgeId(p.edgesList[i][j].first, p.edgesList[i][k].first) != -1)
          ++triangle;
        ++triple;
      }
    coefficient[i] = triple == 0 ? 0.0 : triangle / (double)triple;
  }
#ifdef DEBUG_MODE
  cout << "coefficient=";
  for (auto a : coefficient) cout << a << " ";
  cout << endl;
#endif
  // (3) neighbour avg degree
  vector<double> neighbor_avg_degree(n);
  for (int i = 0; i < n; ++i) {
    neighbor_avg_degree[i] = 0.0;
    for (const auto &edge : p.edgesList[i])
      neighbor_avg_degree[i] += degree[edge.first];
    neighbor_avg_degree[i] /= degree[i];
  }
  // (4) neighbour avg coefficient
  vector<double> neighbor_avg_coefficient(n);
  for (int i = 0; i < n; ++i) {
    neighbor_avg_coefficient[i] = 0.0;
    for (const auto &edge : p.edgesList[i])
      neighbor_avg_coefficient[i] += coefficient[edge.first];
    neighbor_avg_coefficient[i] /= degree[i];
  }
  // (5) number of edges in node i’s egonet
  vector<double> m_egonet(n, 0);
  // (6) number of outgoing edges from ego(i)
  vector<double> out_m_egonet(n, 0);
  // (7) number of neighbors of ego(i)
  vector<double> n_egonet_neighbor(n, 0);
  for (int i = 0; i < n; ++i) {
    vector<bool> egonet(n, false);
    egonet[i] = true;
    int n_egonet = 1;
    for (const auto &edge : p.edgesList[i]) {
      egonet[edge.first] = true;
      ++n_egonet;
    }
    for (const auto &edge : p.edges)
      if (egonet[edge.s] && egonet[edge.e]) m_egonet[i] += 1;
    for (const auto &edge : p.edges)
      if (edge.s != i && edge.e != i && (egonet[edge.s] ^ egonet[edge.e]))
        out_m_egonet[i] += 1;
    vector<bool> egonet_neighbor(n, false);
    for (const auto &edge : p.edges)
      if (egonet[edge.s] && !egonet[edge.e])
        egonet_neighbor[edge.e] = true;
      else if (egonet[edge.e] && !egonet[edge.s])
        egonet_neighbor[edge.s] = true;
    for (int j = 0; j < n; ++j)
      if (j != i && egonet_neighbor[j]) n_egonet_neighbor[i] += 1;
  }

  vector<double> res;
  vector<double> f;
  f = AggregateFeatures(degree);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(coefficient);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(neighbor_avg_degree);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(neighbor_avg_coefficient);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(m_egonet);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(out_m_egonet);
  copy(f.begin(), f.end(), back_inserter(res));
  f = AggregateFeatures(n_egonet_neighbor);
  copy(f.begin(), f.end(), back_inserter(res));
  return res;
}
double CamberraDistance(const vector<double> &a, const vector<double> &b) {
  double res = 0;
  for (int i = 0; i < a.size(); ++i)
    if (std::abs(a[i] - b[i]) > 0)
      res += std::abs(a[i] - b[i]) / (std::abs(a[i]) + std::abs(b[i]));
  res /= a.size();
  return res;
}
}  // namespace net_simile
