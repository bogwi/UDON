# All purpose Graph map container in C++

## Description

This is largely a port of the adjacency map graph I did in ZIG earlier. 
All edges incident of a vertex are collected in a map,
using the adjacent vertex as a key. The graph can be initiated as one of four variants, 
subjects of graph type: `directed, undirected` and mode: `weighted, unweighted`. 
The only difference to the ZIG version is that this implementation uses pointers wherever possible 
reducing considerably memory overhead. (ZIG version must be updated 🙄). 
The code is written in one long file piled into a single class declaration.

The adjacency map structure uses hash maps as the sole storage container for outgoing and incoming edges:

`outGoing = map1`\
`inComing = outGoing* | map1`

where,

`map1: { key = origin_vertex; value = map2 }`\
`map2: { key = destination_vertex*; value = edge }`

If the graph is initiated as unordered, `inComing` simply points to `outGoing`. 
If directed, `inComing` becomes a separate instance.

Vertex is a simple two-fields class:
```c++
class Vertex { 
public:
    VTX_ID id; VTX_DAT data; 
    ... }
```
Edge stores only pointers, preventing double storage. Non-initiated weights are set to zero by default.
```c++
class Edge { 
public:
    const Vertex* origin; const Vertex* destination;
    edgeID id; long double weight=0; 
    ... }
```
The graph's template is next:
```c++
UDON<VTX_ID, VTX_DAT, edgeID, GraphType, GraphMode>
```
`VTX_ID` and `edgeID` assume `string | int`, but you can use any type you know how to hash.\
`VTX_DAT` - any type.

And a possible start may be:

```c++
auto graph = UDON<std::string, some_data*, std::string, undirected, weighted>();

auto A = graph.insertVertex("A", &data0, 0);
auto B = graph.insertVertex("B", &data1, 1);
auto AB = graph.insertEdge(A, B, "AB");
```
The biggest advantage of this type of graph is the virtually O(1) get query for any vertex or edge which is seen from the benchmark where the graph is iterated over with subsequent query back. The disadvantage is still some memory overhead as it is with hash maps.

## API
The graph has a number of methods to call. Worth mentioning are querying and extracting vertices by id only. There is no need to own vertices if you have their ids aldready. Having vertices ids only allows you to query any data thfor any vertex or edge the graph contains.\
The API is well populated with doc strings. LSPs like **clangd** will get info for each method. Please check `test.cc` for usage examples;
```
.clear
.mergeIntoSelf

.insertVertex
.insertVertexIfVertex
.updateVertexData
.updateVertexDataByPtr
.removeVertexByPtr
.removeVertexById
.containsVertex
.containsVertexByPtr
.containsVertexBytId
.extractVertexById
.extractVertexPtrByItd
.getVertexData
.vertexCount
.vertices -> VertexIterator
.verticesToVector
.adjacentVertices
.isOutgoingTo
.isOutogingToById
.isIncomingTo
.isIncomingToById

.insertEdge
.insertOrAssignEdge
.removeEdge
.removeEdgeByEndpointsPtrs
.removeEdgeByEndpointsIds
.containsEdge
.containsEdgeIfEdge
.extractEdgeIfEdge
.extractEdgeByEndpointsPtrs
.extractEdgeByEndpointsIds
.edgeCount
.degree
.incidentEdges
.mapEdgesToSet
.mapEndpointsToSet
.edges -> EdgeIterator

.traverseGraphAsTree
    breadth-first | preorder | postorder | inorder
.connectionTree
    path algorithms: depth-first search iterative | depth-first search recursive
    short path algorithm: breadth-first | Dijkstra
.topologicalSort
.primJarnikMST
.KruskalMST
```

## Performance
The accompanying `test.cc` and `bench.cc` have been tested with `clang++ -std=c++17` and `zig c++`.
`bench.cc` with flags `-Oz` or `-O3`. Most likely will not compile with gcc.

The benchmark is a simple synthetic test, connecting vertices in any order.
The resulting graph becomse very spidery. This is a speed test rather than an allocation test. 
The graph was started with 64-bit integer ids and 5-byte string data fields.

After compilation you can pass the number of vertices and edges per vertex to run. Example:
```
./a.out 10000 20
```
means 10k vertices with an average of 20 edges each. These are the default parameters.
The code and the bench have been tested for leaks with Valgrind.

#### Results
**apple M1 32GB RAM laptop.**

- All the tests are run on the whole graph by a single thread and the time is taken for each test.
- <u>Vertices and edges iteration and get query</u> tests give a measurable result in milliseconds when the graph count for vertices is > 1mil. 
Small graphs are in ns terittory. For example, when the graph has 1mil vertices and 20mil edges, 
iterating over **all the edges** with subsequent get query for each edge takes 379ms to complete the test.
- Short path algorithms start from an arbitrary origin vertex and compute a complete connection tree to all the other N-1 vertices in the graph.
- Dijsktra and Prim-Jarnik algorithms can be improved with faster priority queues.
- The benchmark can be considered a stress test for some laptop configurations when pushed hard, with over a million vertices and a factor of 10 for edges.  

```
BENCHING UDON GRAPH: connection random: try: ........10000:20
building undirected: .....................10000:199616:0.109s
Dijkstra alg: .........................................0.063s
Breadth-First Search alg: .............................0.039s
MST: Prim-Jarnik alg: .................................0.167s
MST: Kruskal alg: .....................................0.061s
Vertices iteration and get query: .........................0s
Remove all vertices: ..................................0.098s
dealocating: ..........................................0.001s

building directed: .......................10000:202142:0.055s
Dijkstra alg: .........................................0.006s
Breadth-First Search alg: .............................0.014s
Topological sort: .....................................0.017s
Vertices iteration and get query: .........................0s
Edges iteration and get query: ........................0.003s
Remove all edges: .....................................0.073s
Remove all vertices: ..................................0.007s
dealocating: ..........................................0.001s


BENCHING UDON GRAPH: connection random: try: .......100000:20
building undirected: ...................100000:1999623:1.208s
Dijkstra alg: .........................................1.061s
Breadth-First Search alg: .............................0.607s
MST: Prim-Jarnik alg: .................................2.251s
MST: Kruskal alg: .....................................0.911s
Vertices iteration and get query: .........................0s
Remove all vertices: ...................................1.43s
dealocating: ..........................................0.021s

building directed: .....................100000:2048296:1.212s
Dijkstra alg: .........................................0.105s
Breadth-First Search alg: ..............................0.11s
Topological sort: .....................................0.228s
Vertices iteration and get query: .....................0.001s
Edges iteration and get query: ........................0.036s
Remove all edges: .....................................1.129s
Remove all vertices: ..................................0.092s
dealocating: ..........................................0.021s


BENCHING UDON GRAPH: connection random: try: ......1000000:20
building undirected: ................1000000:19999583:18.734s
Dijkstra alg: .........................................16.07s
Breadth-First Search alg: .............................7.661s
MST: Prim-Jarnik alg: .................................28.51s
MST: Kruskal alg: ....................................15.845s
Vertices iteration and get query: .....................0.005s
Remove all vertices: .................................20.119s
dealocating: ..........................................0.168s

building directed: ..................1000000:20499186:20.132s
Dijkstra alg: .........................................1.278s
Breadth-First Search alg: .............................1.015s
Topological sort: .....................................3.397s
Vertices iteration and get query: .....................0.014s
Edges iteration and get query: ........................0.379s
Remove all edges: ....................................15.878s
Remove all vertices: ..................................1.163s
dealocating: ..........................................0.168s
```
