/*
 * This file can only used for testing set intersection algorithms
 * Do not include it in any other files
 */
#include "config/type.h"
#include "util/set_intersection.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <parallel/algorithm>
#include <sys/time.h>
#include <assert.h>


int main(int argc, char *argv[])
{
    int length = atoi(argv[1]);

    ui *bucket1 = new ui[length];
    ui *bucket2 = new ui[length];
    ui *bucket3 = new ui[length];
    ui *bucket4 = new ui[length];
    ui *vectorized = new ui[length];

    std::vector<ui> serial, parallel;

    srand((unsigned)time(NULL));
    
    ui prev1 = 0, prev2 = 0;
    for(ui i = 0; i < length; i++) {
        ui dice1 = ((long long)rand() * (long long)rand()) % 3 + 1;
        ui dice2 = ((long long)rand() * (long long)rand()) % 3 + 1;

        ui e1 = prev1 + dice1;
        ui e2 = prev2 + dice2;

        bucket1[i] = bucket3[i] = e1;
        bucket2[i] = bucket4[i] = e2;
        prev1 = e1;
        prev2 = e2;
    }

    timeval stdSerialBegin, stdParallelBegin, vectorizedBegin;
    timeval stdSerialEnd, stdParallelEnd, vectorizedEnd;

    gettimeofday(&stdSerialBegin, NULL);
    std::set_intersection(bucket1, bucket1 + length, bucket2, bucket2 + length, std::back_inserter(serial));
    gettimeofday(&stdSerialEnd, NULL);

    gettimeofday(&stdParallelBegin, NULL);
    __gnu_parallel::set_intersection(bucket1, bucket1 + length, bucket2, bucket2 + length, std::back_inserter(parallel));
    gettimeofday(&stdParallelEnd, NULL);

    gettimeofday(&vectorizedBegin, NULL);
    ui count_ = 0;
    compute_set_intersection::computeCandidates(bucket3, length, bucket4, length, vectorized, count_);
    gettimeofday(&vectorizedEnd, NULL);

    std::unordered_set<ui> verify;
    for(ui i = 0; i < length; i++)
        verify.insert(bucket3[i]);
    std::vector<ui> veri;
    for(ui i = 0; i < length; i++)
        if(verify.find(bucket4[i]) != verify.end())
            veri.emplace_back(bucket4[i]);
    std::sort(veri.begin(), veri.end());
    auto iter = std::unique(veri.begin(), veri.end());
    veri.resize(std::distance(veri.begin(), iter));

    std::cout << parallel.size() << " " << serial.size() << " " << count_ << " " << veri.size() << std::endl;
    for(ui i = 0; i < count_; i++) {
        if(serial[i] != vectorized[i]) {
            std::cerr << serial[i] << " " << vectorized[i] << std::endl;
            exit(1);
        }
    }

    double serialTime = stdSerialEnd.tv_sec - stdSerialBegin.tv_sec + (stdSerialEnd.tv_usec - stdSerialBegin.tv_usec) / 1e6;
    double parallelTime = stdParallelEnd.tv_sec - stdParallelBegin.tv_sec + (stdParallelEnd.tv_usec - stdParallelBegin.tv_usec) / 1e6;
    double vectorizedTime = vectorizedEnd.tv_sec - vectorizedBegin.tv_sec + (vectorizedEnd.tv_usec - vectorizedBegin.tv_usec) / 1e6;

    std::cout << "###     std serial time (s): " << serialTime << std::endl;
    std::cout << "###   std parallel time (s): " << parallelTime << std::endl;
    std::cout << "### std vectorized time (s): " << vectorizedTime << std::endl;

    delete[] bucket1;
    delete[] bucket2;
    delete[] bucket3;
    delete[] bucket4;
    delete[] vectorized;

    return 0;
}