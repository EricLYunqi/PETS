#ifndef _CONFIG_H_
#define _CONFIG_H_

// If we need to change pattern size
// Remember to modify super node size
// as of the # of vertices
#define MIN_SUPER_NODE_SIZE 4
// deprecated
// #define MAX_SUPER_NODE_SIZE 15

// as of the # of edges
#define MIN_PATTERN_SIZE 4
#define MAX_PATTERN_SIZE 15

// max depth during dfs of a subgraph in Algorithm 2
#define MAX_GROW_DEPTH 3
// number of top k maintained for selection
#define DIVERSIFIED_K 2

// final pattern size: \gamma
#define PATTERN_SET_SIZE 30

/*
 * contract graph
 */
#define SUPER_NODE_VERIFICATION 1


/*
 * net simile
 */
#define DEBUG_MODE 0


// recomend to use these default settings
// since others are not tested
/**
 * Set intersection method.
 * 0: Hybrid method; 1: Merge based set intersections.
 */
#define HYBRID 0

/**
 * Accelerate set intersection with SIMD instructions.
 * 0: AVX2; 1: AVX512; 2: Basic;
 */
#define SI 0

/*
 * pattern io
 */
#define MAX_DISPLAY_PATTERN_SIZE 22

/*
 * Threshold for whether update
 */
#define HIT_THRESHOLD 0.5

#endif // _CONFIG_H_