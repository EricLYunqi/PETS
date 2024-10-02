#ifndef _PATTERN_IO_HPP_
#define _PATTERN_IO_HPP_

#include "config/config.h"
#include "config/type.h"
#include "graph/pattern.h"
#include <fstream>

namespace pattern_io
{

/*
 * save patterns on disk
 * the normal format
 * |V|, |E|
 * (s1, t1)
 * (s2, t2)
 * ...
 * (s_{|E| - 1}, t_{|E| - 1})
 */
inline void savePatternNormalFormat(std::ofstream &fileStream, const Pattern &pattern, ui count)
{
    ui n = pattern.getNumVertex();
    ui m = pattern.getNumEdge();

    fileStream << "t # " << count << " " << n << "\n";
    for(ui j = 0; j < n; ++j) 
        fileStream << "v " << j << " 0\n";
    for(ui j = 0; j < m; ++j)
        fileStream << "e " << pattern.edges[j].first << " " << pattern.edges[j].second << " 0\n";
    fileStream << std::endl;
}

inline void savePatternsNormalFormat(const std::string &fileName, const std::vector<Pattern> &patterns)
{
    std::ofstream fileStream(fileName);
    for(ui i = 0; i < patterns.size(); i++)
        savePatternNormalFormat(fileStream, patterns[i], i);
    fileStream.close();
}

inline void savePatternsNormalFormat(const std::string &fileName, const Pattern *patterns)
{
    std::ofstream fileStream(fileName);
    for(ui i = 0; i < PATTERN_SET_SIZE; i++)
        savePatternNormalFormat(fileStream, patterns[i], i);
    fileStream.close();
}

/*
 * save patterns on disk
 * the graphviz format
 * and call command to draw patterns as png format
 */
inline bool savePatternGraphvizFormat(const std::string &dirName, const Pattern &pattern, ui count)
{
    if(pattern.getNumEdge() >= MAX_DISPLAY_PATTERN_SIZE) {
        std::cerr << "I can not display it: " << pattern.getNumEdge() << std::endl;
        return false;
    }

    const std::string dotFileName = dirName + "/dotfile";
    std::ofstream fileStream(dotFileName);

    fileStream << "graph \"result\" {\n";
    fileStream << "graph [fontsize = 30, ratio = 0.835, dpi = 100, size = \"15,15\" "
                  "];\n";
    fileStream << "node [label = \"\\N\", shape = doublecircle, sides = 4, color = "
                  "skyblue, style = filled ];\n";
    for(ui i = 0; i < pattern.getNumVertex(); ++i)
        fileStream << i << " [shape = doublecircle, label = \"" << i << "\"];\n";
    for(auto &e : pattern.edges) 
        fileStream << e.first << " -- " << e.second << ";\n";
    fileStream << "}\n";
    fileStream.close();

    const std::string pngFileName = dirName + "/pattern" + std::to_string(count) + ".png";
    const std::string drawCommand = "dot -Tpng  -Kneato -Gepsilon=0.0001 -Goverlap=false " +
                                    dotFileName + " -o " + pngFileName;

    return true;
}

inline void savePatternsGraphvizFormat(const std::string &dirName, const std::vector<Pattern> &patterns)
{
    ui success = 0;
    std::cout << "writing png to " << dirName << " ..." << std::endl;
    for(ui i = 0; i < patterns.size(); i++) {
        bool res = savePatternGraphvizFormat(dirName, patterns[i], i);
        success = res ? success + 1 : success;
    }
    std::cout << "successful written: " << success << std::endl;
}

inline void savePatternsGraphvizFormat(const std::string &dirName, const Pattern *patterns)
{
    ui success = 0;
    std::cout << "writing png to " << dirName << " ..." << std::endl;
    for(ui i = 0; i < PATTERN_SET_SIZE; i++) {
        bool res = savePatternGraphvizFormat(dirName, patterns[i], i);
        success = res ? success + 1 : success;
    }
    std::cout << "successful written: " << success << std::endl;
}

};

#endif // _PATTERN_IO_HPP_