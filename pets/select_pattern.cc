#include "pets/select_pattern.h"

void SelectPattern::allocateBuffers(Pattern *&patterns)
{
    patterns = new Pattern[PATTERN_SET_SIZE];
}

void SelectPattern::releaseBuffers(Pattern *&patterns)
{
    delete[] patterns;
    patterns = nullptr;
}

void SelectPattern::calculateCoverage(const SuperGraph *superGraph, std::vector<Pattern> *DRP)
{
    ui numVertex = superGraph->numVertex;
    std::vector<double> covmin(5, 1e18);
    std::vector<double> covmax(5, 0);

    // |E| * freq(p)
    for (ui i = 0; i < numVertex; i++)
    {
        for (auto &p : DRP[i])
        {
            ui nodeid = p.nodeId;
            auto synopsis = superGraph->synopses[nodeid];
            auto sntype = static_cast<int>(synopsis.type);
            auto nmap = superGraph->nodeLengthMap[sntype];
            ui nodesizeV = synopsis.cog.first;
            ui nodesizeE = synopsis.cog.second;

            // here if we use superGraph->nodeLengthMap[sntype] the IDE throws an error
            // check it during compilation
            ui freq = nmap[nodesizeV];
            ui numE = nodesizeE + p.getNumEdge();
            p.scoreCoverage = (double)freq * numE * 1.0;

            covmin[sntype] = covmin[sntype] < p.scoreCoverage ? covmin[sntype] : p.scoreCoverage;
            covmax[sntype] = covmax[sntype] > p.scoreCoverage ? covmax[sntype] : p.scoreCoverage;
        }
    }

    // normalize
    for (ui i = 0; i < numVertex; i++)
    {
        for (auto &p : DRP[i])
        {
            auto sntype = static_cast<int>(superGraph->synopses[p.nodeId].type);
            p.scoreCoverage = (p.scoreCoverage - covmin[sntype] + 1) / (covmax[sntype] - covmin[sntype] + 1);
        }
    }
}


void SelectPattern::calculateCognition(const SuperGraph *superGraph, std::vector<Pattern> *DRP)
{
    ui numVertex = superGraph->numVertex;

    auto f1 = [](double x) -> double { 
        return 1 - exp(-x); 
    };
    auto f2 = [](double c1, double c2, double x) -> double {
        return 1 / (1 + exp(-c1 * (x - c2)));
    };

    for (ui i = 0; i < numVertex; i++) {
        for (auto &p : DRP[i]) {
            ui nodeid = p.nodeId;
            auto synopsis = superGraph->synopses[nodeid];
            ui nodesizeV = synopsis.cog.first;
            ui nodesizeE = synopsis.cog.second;

            ui numV = nodesizeV + p.getNumVertex();
            ui numE = nodesizeE + p.getNumEdge();
            double density = (2.0 * numE) / (numV * (numV - 1));
            double crossing = isPlanar(numV, nodesizeV, numE, synopsis.type) ? 0 : numE - 3 * numV + 6;

            /*1.0 * pattern.getEdgeCnt() * pattern.getEdgeCnt() * pattern.getEdgeCnt() /
             * pattern.getNodesCnt() / pattern.getNodesCnt() / 33.75 - 0.9 *
             * pattern.getNodesCnt();*/
            if (crossing < 0)
                crossing = 0;

            // cog7
            // pattern.scoreCognition = f1(pattern.getEdgeCnt() + density + crossing);
            p.scoreCognition = f2(0.5, 10, numE + density + crossing);
        }
    }
}


void SelectPattern::greedySelect(std::vector<Pattern> *DRP, Pattern *patterns, ui numVertex, 
                                 ui numAttribute)
{
    ui patternCount = 0;
    ui totpattern = numVertex * DIVERSIFIED_K;
    std::vector<std::vector<double>> signature(totpattern);

    // signature
    for(ui i = 0; i < numVertex; i++)
        for(const auto &p : DRP[i])
            signature[p.id] = NetSimile::signature(p);

    std::vector<int> flag(totpattern + 1, 0);
    // attr flag, sometimes # of attr may be 3*1e6, thus using bitset
    // unsigned long long attrflag = 0, tempaf = 0;
    std::bitset<3000000> af, tempaf;

    for (int i = 0; i < PATTERN_SET_SIZE; i++) {
        double maxps = -1e8;
        int maxnid = -1;
        int maxp = -1;

        for(int nid = 0; nid < (int)numVertex; nid++) {
            int pos = 0;
            for(auto &p : DRP[nid]) {
                if(flag[p.id])
                    continue;

                double cov = p.scoreCoverage;
                double cog = p.scoreCognition;
                double sim = 0.0;
                double attr = 0.0;
                tempaf.reset();

                // sim
                for (int k = 0; k < i; ++k) {
                    double tmp = NetSimile::camberraDistance(signature[p.id],
                                                             signature[patterns[k].id]);
                    if (tmp > sim)
                        sim = tmp;
                }
                // attr
                for(const auto &as : p.attributes)
                    for(const auto &a : as)
                        tempaf[a] = true;
                tempaf &= af;
                attr = (double)tempaf.count() / numAttribute * 1.0;
                p.scoreAttribute = attr;
                
                double delta_s = cov - sim - cog + attr;
                if (delta_s > maxps) {
                    maxps = delta_s;
                    maxnid = nid;
                    maxp = pos;
                }

                ++ pos;
            }
        }

        if (maxp == -1)
            break;

        // update
        const auto &selectp = DRP[maxnid][maxp];
        flag[selectp.id] = 1;
        patterns[i] = selectp;
        for(const auto &as : selectp.attributes)
            for(const auto &a : as)
                af[a] = true;
        ++ patternCount;
    }

    // when less than using finalscore sort
    flag.clear();
    flag.resize(totpattern, 0);
    for(ui i = 0; i < patternCount; i++)
        flag[patterns[i].id] = 1;

    int id = 0;
    for (int i = patternCount; i < PATTERN_SET_SIZE;) {
        if (flag[id] == 0) {
            ++ i;
            std::pair<ui, ui> pos = getPatternPos(id);
            patterns[i] = DRP[pos.first][pos.second];
            flag[id] = 1;
        }
        
        ++ id;
        if (id >= (int)totpattern)
            break;
    }

    for (int i = 0; i < PATTERN_SET_SIZE; ++i) {
        ui id = patterns[i].id;
        std::pair<ui, ui> pos = getPatternPos(id);

        patterns[i].finalRank = i;
        DRP[pos.first][pos.second].finalRank = i;

        double sim = 0.0;
        for (int j = i + 1; j < PATTERN_SET_SIZE; ++j) {
            double tmp = NetSimile::camberraDistance(signature[id], 
                                                     signature[patterns[j].id]);
            if (tmp > sim)
                sim = tmp;
        }

        patterns[i].scoreDiversity = 1.0 / sim;
        DRP[pos.first][pos.second].scoreDiversity = 1.0 / sim;
    }
}


void SelectPattern::reportPatternSetScores(const Pattern *patterns, std::string &filename)
{
    FILE *fp = fopen(filename.c_str(), "w");
    if(fp == nullptr) {
        std::cerr << "Cannot open: " << filename << std::endl;
        exit(1);
    }

    fprintf(fp, "Total patterns: %d\n", PATTERN_SET_SIZE);
    fprintf(fp, "_id\tcoverage\tcognition\tdiversity\tattribute\trank\n");

    for(ui i = 0; i < PATTERN_SET_SIZE; i++)
        fprintf(fp, "%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%d\n", patterns[i].id, 
                                                            patterns[i].scoreCoverage, 
                                                            patterns[i].scoreCognition, 
                                                            patterns[i].scoreDiversity, 
                                                            patterns[i].scoreAttribute, 
                                                            patterns[i].finalRank);

    fclose(fp);
}


// utils
std::pair<ui, ui> SelectPattern::getPatternPos(ui id)
{
    ui nodeid = id / DIVERSIFIED_K;
    ui nodepos = id % DIVERSIFIED_K;
    return std::make_pair(nodeid, nodepos);
}

bool SelectPattern::isPlanar(ui numV, ui nodeV, ui numE, SuperNodeType sntype)
{
    if (numV < 3)
        return true;
    if (numE > 3 * numV - 6)
        return false;
    // no cycles of length 3
    if (sntype != SuperNodeType::CHORD && (nodeV != 3 || sntype != SuperNodeType::CYCLE) && numE > 2 * numV - 4)
        return false;
    return true;
}

int SelectPattern::countBits(unsigned long long n)
{
    int count = 0;
    while (n) {
        count += n & 1;
        n >>= 1;
    }
    return count;
}