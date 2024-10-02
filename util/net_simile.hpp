// https://github.com/kristyspatel/Netsimile/blob/master/Netsimile.py
// http://basicodingfordummies.blogspot.com/2013/02/mean-median-variance-skewness-kurtosis.html
// https://www.johndcook.com/blog/skewness_kurtosis/
#ifndef _NET_SIMILE_HPP_
#define _NET_SIMILE_HPP_

#include "config/config.h"
#include "config/type.h"
#include "graph/pattern.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace net_simile
{



};

class NetSimile
{
public:
    NetSimile() = default;
    ~NetSimile() = default;
    NetSimile(const NetSimile &other) = delete;
    NetSimile(NetSimile &&other) = delete;

public:
    // return 1 * 5
    static std::vector<double> aggregateFeatures(const std::vector<double> &a2) {
        int n = a2.size();
        auto a = a2;
        std::sort(a.begin(), a.end());
        double sum = 0.0;
        for (auto v : a)
            sum += v;
        double mean = sum / n * 1.0;
        double median = (n & 1) ? a[n / 2] * 1.0 : ((a[n / 2 - 1] + a[n / 2]) * 0.5);
        double variance = 0.0;
        for (auto v : a)
            variance += (v - mean) * (v - mean);
        variance /= n - 1;
#if DEBUG_MODE == 1
        std::cout << "Variance = " << variance << std::endl;
#endif
        double std_deviation = sqrt(variance);
        double skewness = 0.0;
        if (std_deviation > 0.0)
        {
            for (auto v : a)
                skewness += (v - mean) * (v - mean) * (v - mean);
            skewness /= n * std_deviation * std_deviation * std_deviation;
        }
        double kurtosis = 0.0;
        if (std_deviation > 0.0)
        {
            for (auto v : a)
                kurtosis += (v - mean) * (v - mean) * (v - mean) * (v - mean);
            kurtosis /=
                n * std_deviation * std_deviation * std_deviation * std_deviation;
            kurtosis -= 3.0;
        }

        return {mean, median, std_deviation, skewness, kurtosis};
    }

    // return 1 * 35
    static std::vector<double> signature(const Pattern &p) {
        auto n = (int)p.getNumVertex();

        // (1) node degree
        std::vector<double> degree(n);
        for (int i = 0; i < n; ++i)
            degree[i] = p.edgesList[i].size();

        // (2) node coefficient
        std::vector<double> coefficient(n);
        for (int i = 0; i < n; ++i) {
            int triangle = 0, triple = 0;
            int edge_num = degree[i];
            for (int j = 0; j < edge_num; ++j) {
                for (int k = j + 1; k < edge_num; ++k) {
                    if (p.getEdgeId(p.edgesList[i][j].first, p.edgesList[i][k].first) != NOTEXIST)
                        ++triangle;
                    ++triple;
                }
            }
            coefficient[i] = triple == 0 ? 0.0 : triangle / (double)triple;
        }
#if DEBUG_MODE == 1
        std::cout << "coefficient=";
        for (auto a : coefficient)
            std::cout << a << " ";
        std::cout << std::endl;
#endif

        // (3) neighbour avg degree
        std::vector<double> neighbor_avg_degree(n);
        for (int i = 0; i < n; ++i) {
            neighbor_avg_degree[i] = 0.0;
            for (const auto &edge : p.edgesList[i])
                neighbor_avg_degree[i] += degree[edge.first];
            neighbor_avg_degree[i] /= degree[i];
        }

        // (4) neighbour avg coefficient
        std::vector<double> neighbor_avg_coefficient(n);
        for (int i = 0; i < n; ++i) {
            neighbor_avg_coefficient[i] = 0.0;
            for (const auto &edge : p.edgesList[i])
                neighbor_avg_coefficient[i] += coefficient[edge.first];
            neighbor_avg_coefficient[i] /= degree[i];
        }

        // (5) number of edges in node i’s egonet
        std::vector<double> m_egonet(n, 0.0);
        // (6) number of outgoing edges from ego(i)
        std::vector<double> out_m_egonet(n, 0.0);
        // (7) number of neighbors of ego(i)
        std::vector<double> n_egonet_neighbor(n, 0.0);

        for (int i = 0; i < n; ++i) {
            std::vector<bool> egonet(n, false);
            egonet[i] = true;
            int n_egonet = 1;
            for (const auto &edge : p.edgesList[i]) {
                egonet[edge.first] = true;
                ++n_egonet;
            }
            for (const auto &edge : p.edges)
                if (egonet[edge.first] && egonet[edge.second])
                    m_egonet[i] += 1;
            for (const auto &edge : p.edges)
                if ((int)edge.first != i && (int)edge.second != i && (egonet[edge.first] ^ egonet[edge.second]))
                    out_m_egonet[i] += 1;
            std::vector<bool> egonet_neighbor(n, false);
            for (const auto &edge : p.edges)
                if (egonet[edge.first] && !egonet[edge.second])
                    egonet_neighbor[edge.second] = true;
                else if (egonet[edge.second] && !egonet[edge.first])
                    egonet_neighbor[edge.first] = true;
            for (int j = 0; j < n; ++j)
                if (j != i && egonet_neighbor[j])
                    n_egonet_neighbor[i] += 1;
        }

        std::vector<double> res;
        std::vector<double> f;
        // copy res to the end of f
        f = aggregateFeatures(degree);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(coefficient);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(neighbor_avg_degree);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(neighbor_avg_coefficient);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(m_egonet);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(out_m_egonet);
        std::copy(f.begin(), f.end(), std::back_inserter(res));
        f = aggregateFeatures(n_egonet_neighbor);
        std::copy(f.begin(), f.end(), std::back_inserter(res));

        return res;
    }

    // distance
    static double camberraDistance(const std::vector<double> &a, const std::vector<double> &b) {
        double res = 0.0;
        auto asize = (int)a.size();
        for (int i = 0; i < asize; ++i)
            if (std::abs(a[i] - b[i]) > 0)
                res += std::abs(a[i] - b[i]) / (std::abs(a[i]) + std::abs(b[i]));
        res /= asize;
        return res;
    }
};

#endif // NET_SIMILE_HPP