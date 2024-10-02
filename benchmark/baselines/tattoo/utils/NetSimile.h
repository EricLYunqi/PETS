#ifndef NET_SIMILE_H_
#define NET_SIMILE_H_

#include <vector>

#include "../base.h"

using namespace std;

namespace net_simile {
// return 1 * 7
vector<double> AggregateFeatures(const vector<double>& a2);
// return 1 * 35
vector<double> Signature(const Pattern& p);
double CamberraDistance(const vector<double>& a, const vector<double>& b);
}  // namespace net_simile

#endif
