1. Partial codes are adopted from TATTOO
2. Adjust the super node size within pattern size, thus could avoid judgment in Algorithm2 (generate pattern)
3. We drop some edges during contraction, like inner edges between circles, edges from internal vertices of a path, maybe add them in the furture.
   Current contraction scheme:
   - Cycle: drop internal edges
   - Chord: normal
   - Path: internal vertices can also has edges out
   - Star: no edges between leaves
4. Just like attribute community search, we could set an attribute constraint, i.e., given an input attribute set and select patterns having such attributes.
5. Judgement of planar graph is not sufficient & necessary
6. Remember to map node id if the node id is not continous
7. the Flower Queries Generation may not be correct
8. Why we need to consider attributes in selecting queries?

Shift + Option + F