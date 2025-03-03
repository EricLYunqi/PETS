#ifndef _TYPE_H_
#define _TYPE_H_

typedef unsigned int ui;
typedef unsigned int VertexID;

enum class GraphType
{
    SOCIAL = 0,
    WEB = 1,
    TRAFFIC = 2
};

enum class SuperNodeType
{
    CHORD = 0,
    CYCLE = 1,
    PATH = 2,
    STAR = 3,
    SINGLETON = 4
};

enum class UpdateStrategy
{
    ATTRIBUTE = 0,
    STRUCTURE = 1,
    INTERSECTION = 2,
    UNION = 3
};

#endif // _TYPE_H_