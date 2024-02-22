// MIT (c) bogwi@rakumail.jp
// This code uses tab=8 convention.
// Compiles with clang++ -std=c++17 or zig c++ compiler 

#include "udon.cc"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_set>

using std::cout;
using std::string;


void primJarnikKruskal () {
        cout << "UNDIRECTED GRAPH primJarnikKruskal() TEST .. " << std::flush;
        auto graph = UDON<string, string, string, undirected, weighted>();
        
        // Example of a planar graph:
        // https://en.wikipedia.org/wiki/Minimum_spanning_tree.

        auto A = graph.insertVertex("A", "any data type can be here");
        auto B = graph.insertVertex("B", "any data type can be here");
        auto C = graph.insertVertex("C", "any data type can be here");
        auto D = graph.insertVertex("D", "any data type can be here");
        auto E = graph.insertVertex("E", "any data type can be here");
        auto F = graph.insertVertex("F", "any data type can be here");
        auto G = graph.insertVertex("G", "any data type can be here");
        auto H = graph.insertVertex("H", "any data type can be here");
        auto I = graph.insertVertex("I", "any data type can be here");
        auto J = graph.insertVertex("J", "any data type can be here");

        graph.insertEdge(A, B, "AB", 3);
        graph.insertEdge(A, C, "AC", 4);
        graph.insertEdge(A, E, "AE", 10);
        graph.insertEdge(A, F, "AF", 18);

        graph.insertEdge(B, C, "BC", 1);
        graph.insertEdge(B, E, "BE", 9);
        graph.insertEdge(B, D, "BD", 5);

        graph.insertEdge(C, D, "CD", 4);

        graph.insertEdge(D, E, "DE", 7);
        graph.insertEdge(D, G, "DG", 9);
        graph.insertEdge(D, H, "DH", 9);

        graph.insertEdge(E, F, "EF", 8);
        graph.insertEdge(E, G, "EG", 8);
        graph.insertEdge(E, I, "EI", 9);

        graph.insertEdge(F, I, "FI", 9);
        graph.insertEdge(F, J, "FJ", 9);

        graph.insertEdge(G, H, "GH", 2);
        graph.insertEdge(G, I, "GI", 2);

        graph.insertEdge(H, I, "HI", 4);
        graph.insertEdge(H, J, "HJ", 6);

        graph.insertEdge(I, J, "IJ", 3);

        // The computed MST is a struct containing a variable .tree 
        // with endpoint pairs spanning the tree.
        auto mstPj = graph.primJarnikMST();
        auto mstK = graph.kruskalMST();
        assert(mstPj.cost == mstK.cost);

        // The graph above has only one MST, which is identical for both algorithms. 
        // We can assert that. However, on other graphs the MST may be different, but the cost is not. 

        assert(mstPj.contains(A, B) == mstK.contains(A, B));
        assert(mstPj.contains(B, C) == mstK.contains(B, C));
        assert(mstPj.contains(C, D) == mstK.contains(C, D));
        assert(mstPj.contains(D, E) == mstK.contains(D, E));
        assert(mstPj.contains(E, F) == mstK.contains(E, F));
        assert(mstPj.contains(E, G) == mstK.contains(E, G));
        assert(mstPj.contains(G, H) == mstK.contains(G, H));
        assert(mstPj.contains(G, I) == mstK.contains(G, I));
        assert(mstPj.contains(I, J) == mstK.contains(I, J));

        graph.clear();
        assert(graph.vertexCount() == 0); assert(graph.edgeCount() == 0);

        cout << "OK\n";
}

void topologicalSort() {
        struct TASK { string link; void followLink() { /* ... */ }};

        cout << "DIRECTED GRAPH topologicalSort() TEST .. " << std::flush;
        auto graph = UDON<string, TASK, string, directed, unweighted>();

        /** The model of the graph

        5->(11)                          7->(11, 8)                  3->(8, 10)
         
        11->(2, 9, 10)   8->(9)
        
        2                                9                           10

        */

        auto _5 = graph.insertVertex("5", { "https://example5.com"});
        auto _7 = graph.insertVertex("7", { "https://example7.com"});
        auto _3 = graph.insertVertex("3", { "https://example3.com"});
        auto _11 = graph.insertVertex("11", { "https://example11.com"});
        auto _8 = graph.insertVertex("8", { "https://example8.com"});
        auto _2 = graph.insertVertex("2", { "https://example2.com"});
        auto _9 = graph.insertVertex("9", { "https://example9.com"});
        auto _10 = graph.insertVertex("10", { "https://example10.com"});

        assert(graph.getVertexData("5")->link == "https://example5.com");
        assert(graph.getVertexData("7")->link == "https://example7.com");
        assert(graph.getVertexData("3")->link == "https://example3.com");
        assert(graph.getVertexData("11")->link == "https://example11.com");
        assert(graph.getVertexData("8")->link != "https://example5.com");
        assert(graph.getVertexData("2")->link != "https://example7.com");
        assert(graph.getVertexData("19") == None) ;
        assert(graph.getVertexData("100") == None);

        graph.insertEdge(_5, _11, "5-11");

        graph.insertEdge(_7, _11, "7-11");
        graph.insertEdge(_7, _8, "7-8");

        graph.insertEdge(_3, _8, "3-8");
        graph.insertEdge(_3, _10, "3-10");

        graph.insertEdge(_11, _2, "11-2");
        graph.insertEdge(_11, _9, "11-9");
        graph.insertEdge(_11, _10, "11-10");

        graph.insertEdge(_8, _9, "8-9");

        assert(graph.isOutgoingTo(_5, _11)); assert(graph.isIncomingTo(_11, _5));
        assert(graph.isOutgoingTo(_7, _11)); assert(graph.isIncomingTo(_11, _7));
        assert(graph.isOutgoingTo(_7, _8)); assert(graph.isIncomingTo(_8, _7));
        assert(graph.isOutgoingTo(_3, _8)); assert(graph.isIncomingTo(_8, _3));
        assert(graph.isOutgoingToById("3", "10")); assert(graph.isIncomingToById("10", "3"));
        assert(graph.isOutgoingToById("11", "2")); assert(graph.isIncomingToById("2", "11"));
        assert(graph.isOutgoingToById("11", "9")); assert(graph.isIncomingToById("9", "11"));
        assert(graph.isOutgoingToById("11", "10")); assert(graph.isIncomingToById("10", "11"));
        assert(graph.isOutgoingToById("8", "9")); assert(graph.isIncomingToById("9", "8"));
        assert(!graph.isIncomingTo(_9, _5)); assert(!graph.isOutgoingTo(_8, _11));
        assert(!graph.isIncomingToById("11", "11")); assert(!graph.isOutgoingToById("10", "17"));

        assert(graph.adjacentVertices(_11, outgoing).value().size() == 3);
        assert(graph.adjacentVertices(_11, incoming).value().size() == 2);
        assert(graph.adjacentVertices(_8, outgoing).value().size() == 1);
        assert(graph.adjacentVertices(_8, incoming).value().size() == 2);

        // All edges of a directed graph can be iterated over without double reporting.
        auto e_iter = graph.edges(); int c = 0;
        while (auto edge = e_iter.next()) {
                assert(graph.containsEdgeIfEdge(*edge.value())); c++;
        } 
        assert(c == graph.edgeCount());
        
        // reset() moves iterator to the beginning.
        e_iter.reset(); while (auto edge = e_iter.next()) {
                assert(graph.containsEdge(edge.value()->origin, edge.value()->destination));
        }

        // The graph sorts vertices topologically in arbitrary fashion, i.e the sorting scheme depends only 
        // on the precedence factor. No next vertex can have previous vertices as his dependencies. 
        // There are many topological sorts. This is the canonical one.
        // Read more at https://en.wikipedia.org/wiki/Topological_sorting.

        // Create the sort
        auto tS = graph.topologicalSort();

        // Assert the graph has no cycles
        assert(tS.acyclic);

        // Check dependecy coherence - whether the previous TASK does not depend on the next. 
                // temp container
        std::unordered_set<decltype(_5), std::hash<decltype(_5)>> previous;
        for (auto next : tS.topo) {
                // for every previous do the check
                for (auto buffered_task : previous) {
                        // there should be no outgoing edge from the next TASK to any of the previous
                        assert(!graph.isOutgoingTo(next, buffered_task));
                }
                // add the next to previous
                previous.insert({next});
        }

        // Intentially create a cycle
        graph.insertEdge(_9, _5, "9-5");
        tS = graph.topologicalSort();
        assert(!tS.acyclic);

        assert(graph.removeEdgeByEndpointsIds("5", "11")); assert(graph.removeEdgeByEndpointsIds("7", "11"));
        assert(graph.removeEdgeByEndpointsIds("7", "8")); assert(graph.removeEdgeByEndpointsIds("3", "8"));
        assert(graph.removeEdgeByEndpointsIds("3", "10")); assert(graph.removeEdgeByEndpointsIds("11", "2"));
        assert(graph.removeEdgeByEndpointsIds("11", "9")); assert(graph.removeEdgeByEndpointsIds("11", "10"));
        assert(graph.removeEdgeByEndpointsIds("8", "9")); assert(graph.removeEdgeByEndpointsIds("9", "5"));

        assert(graph.edgeCount() == 0); assert(graph.vertexCount() == 8);

        assert(graph.removeVertexById("5")); assert(!graph.containsVertexById("5"));
        assert(graph.removeVertexById("7")); assert(!graph.containsVertexById("7"));
        assert(graph.removeVertexById("3")); assert(!graph.containsVertexById("3"));
        assert(graph.removeVertexById("11")); assert(!graph.containsVertexById("11"));
        assert(graph.removeVertexById("8")); assert(!graph.containsVertexById("8"));
        assert(graph.removeVertexById("2")); assert(!graph.containsVertexById("2"));
        assert(graph.removeVertexById("9")); assert(!graph.containsVertexById("9"));
        assert(graph.removeVertexById("10")); assert(!graph.containsVertexById("10"));

        assert(graph.vertexCount() == 0);

        cout << "OK\n";

}

void shortestPath () {

        cout << "UNDIRECTED GRAPH shortestPath() TEST .. " << std::flush;
        string data = "@";
        auto graph = UDON<string, string, string, undirected, weighted>();

        auto A = graph.insertVertex("A", data); assert(graph.extractVertexPtrById("A").value() == A);
        auto B = graph.insertVertex("B", data); assert(graph.extractVertexById("B").value() == *B);
        auto C = graph.insertVertex("C", data); assert(graph.extractVertexPtrById("C").value() == C);
        auto D = graph.insertVertex("D", data); assert(graph.extractVertexById("D").value() == *D);
        auto E = graph.insertVertex("E", data); assert(graph.extractVertexPtrById("E").value() == E);
        auto F = graph.insertVertex("F", data); assert(graph.extractVertexById("F").value() == *F);
        auto G = graph.insertVertex("G", data); assert(graph.extractVertexPtrById("G").value() == G); 
        auto H = graph.insertVertex("H", data); assert(graph.extractVertexById("H").value() == *H);
        auto I = graph.insertVertex("I", data); assert(graph.extractVertexPtrById("I").value() == I);
        auto J = graph.insertVertex("J", data); assert(graph.extractVertexById("J").value() == *J);
        auto K = graph.insertVertex("K", data); assert(graph.extractVertexPtrById("K").value() == K);
        auto L = graph.insertVertex("L", data); assert(graph.extractVertexById("L").value() == *L);
        auto M = graph.insertVertex("M", data); assert(graph.extractVertexPtrById("M").value() == M);
        auto N = graph.insertVertex("N", data); assert(graph.extractVertexById("N").value() == *N);
        auto O = graph.insertVertex("O", data); assert(graph.extractVertexPtrById("O").value() == O);
        auto P = graph.insertVertex("P", data); assert(graph.extractVertexById("P").value() == *P);

        graph.insertEdge(A, B, "AB", 1); 
        graph.insertEdge(A, E, "AE", 1);  
                                                
        graph.insertEdge(A, F, "AF", 1);

        graph.insertEdge(B, C, "BC", 1);
        graph.insertEdge(B, F, "BF", 1);

        graph.insertEdge(C, D, "CD", 1);
        graph.insertEdge(C, G, "CG", 1);

        graph.insertEdge(D, G, "DG", 1);
        graph.insertEdge(D, H, "DH", 1);

        graph.insertEdge(E, F, "EF", 1);
        graph.insertEdge(E, I, "EI", 1);

        graph.insertEdge(F, I, "FI", 1);

        graph.insertEdge(G, J, "GJ", 1);
        graph.insertEdge(G, K, "GK", 1);
        graph.insertEdge(G, L, "GL", 1);

        graph.insertEdge(H, L, "HL", 1);

        graph.insertEdge(I, J, "IJ", 1);
        graph.insertEdge(I, M, "IM", 1);
        graph.insertEdge(I, N, "IN", 1);

        graph.insertEdge(J, K, "JK", 1);

        graph.insertEdge(K, N, "KN", 1);
        graph.insertEdge(K, O, "KO", 1);

        graph.insertEdge(L, P, "LP", 1);
        graph.insertEdge(M, N, "MN", 1);

        /**

        A ---- B ---- C ---- D
        |  \   |      |  /   |
        E ---- F      G      H
        |   /     /   |  \   |
        I ---- J ---- K      L
        |  \      /   |      |
        M ---- N      O      P

        */

        // Compute an open connection tree with root A to ALL the other vertices 
        // present in the graph using the Dijkstra( DIJ ) algorithm.
        auto origin = A;
        auto algorithm = DIJ;
        auto cT = graph.connectionTree(origin, algorithm);

        // Query the connection to the desired point
        assert(cT.getDistanceTo(B) == 1);
        assert(cT.getDistanceTo(E) == 1);
        assert(cT.getDistanceTo(F) == 1);
        assert(cT.getDistanceTo(C) == 2);
        assert(cT.getDistanceTo(I) == 2);
        assert(cT.getDistanceTo(D) == 3);
        assert(cT.getDistanceTo(J) == 3);
        assert(cT.getDistanceTo(M) == 3);
        assert(cT.getDistanceTo(G) == 3);
        assert(cT.getDistanceTo(N) == 3);
        assert(cT.getDistanceTo(K) == 4);
        assert(cT.getDistanceTo(H) == 4);
        assert(cT.getDistanceTo(L) == 4);
        assert(cT.getDistanceTo(O) == 5);
        assert(cT.getDistanceTo(P) == 5);

        // Since the weights are equal to 1, the Dijkstra algorithm will report the same routes as the
        // as breadth-first search (BFS). 
        origin = P;
        cT = graph.connectionTree(origin, algorithm);
        auto cT_BFS = graph.connectionTree(P, BFS);
        assert(cT.getDistanceTo(A) == cT_BFS.getDistanceTo(A));
        assert(cT.getDistanceTo(B) == cT_BFS.getDistanceTo(B));
        assert(cT.getDistanceTo(C) == cT_BFS.getDistanceTo(C));
        assert(cT.getDistanceTo(D) == cT_BFS.getDistanceTo(D));
        assert(cT.getDistanceTo(E) == cT_BFS.getDistanceTo(E));
        assert(cT.getDistanceTo(F) == cT_BFS.getDistanceTo(F));
        assert(cT.getDistanceTo(G) == cT_BFS.getDistanceTo(G));
        assert(cT.getDistanceTo(H) == cT_BFS.getDistanceTo(H));
        assert(cT.getDistanceTo(I) == cT_BFS.getDistanceTo(I));
        assert(cT.getDistanceTo(J) == cT_BFS.getDistanceTo(J));
        assert(cT.getDistanceTo(K) == cT_BFS.getDistanceTo(K));
        assert(cT.getDistanceTo(L) == cT_BFS.getDistanceTo(L));
        assert(cT.getDistanceTo(M) == cT_BFS.getDistanceTo(M));
        assert(cT.getDistanceTo(N) == cT_BFS.getDistanceTo(N));
        assert(cT.getDistanceTo(O) == cT_BFS.getDistanceTo(O));  

        // Overwriting the edge will of course affect the DIJ algorithm.
        // The connection tree must be recomputed. 
        origin = A;
        graph.insertOrAssignEdge(L, P, "LP", 100);
        cT = graph.connectionTree(origin, algorithm);
        assert(cT.getDistanceTo(P) == 104);

        graph.insertOrAssignEdge(L, P, "LP", 1);

        // Using the knockout set, some vertices can be skipped during the traverse. 
        // The additional parameter TreeTraverseControl must be set to [except].

        std::unordered_set<string> knockout;
        knockout.insert({"D", "C"});

        origin = B;
        cT = graph.connectionTree(origin, algorithm, &knockout, P, except);
        assert(cT.getDistanceTo(P) == 6);
        assert(cT.getDistanceTo(H) == 6);
        assert(cT.getDistanceTo(D) == 0);
        assert(cT.getDistanceTo(C) == 0);
        assert(cT.getDistanceTo(L) == 5);

        // If the target is provided, it will stop traversing as soon as it finds it.
        // However, depending on how far the destination is from the origin and the algorithm chosen,
        // additional routes will be calculated along the way.

        knockout.insert("F");
        origin = P;
        cT = graph.connectionTree(origin, algorithm, &knockout, B, except);
        assert(cT.getDistanceTo(B) == 7);
        assert(cT.getDistanceTo(A) == 6);

        // Setting the TreeTraverseControl parameter to [through], 
        // the walker will only move through
        // the vertices specified in the knockout set. 
        // Knockout set without the control parameter will be ignored.

        knockout.insert("B"); knockout.insert("H"); knockout.insert("L");
        cT = graph.connectionTree(origin, algorithm, &knockout, B, through);
        assert(cT.getDistanceTo(B) == 5);
        assert(cT.getDistanceTo(D) == 3);
        assert(cT.getDistanceTo(B) == 5);
        assert(cT.getDistanceTo(H) == 2);
        assert(cT.getDistanceTo(F) == 6);
        assert(cT.getDistanceTo(L) == 1);
        assert(cT.getDistanceTo(A) == 0);
        assert(cT.getDistanceTo(J) == 0);
        assert(cT.getDistanceTo(N) == 0);
        assert(cT.getDistanceTo(M) == 0);
        assert(cT.getDistanceTo(K) == 0);
        assert(cT.getDistanceTo(G) == 0);

        /** All of the above applies to the remaining two algorithms, DFSi and DFSr, which are
                depth-first iterative and recursive, respectively. 
                However, they do not guarantee a shortest path, 
                but rather explore the graph. Since the Udon graph does not preserve the order of insertion, 
                the DFS algorithms every time it is run, will give you a different connection tree, 
                exploring all the possible connections it could have. 
                The distances (lengths of the paths) can still be queried. */

        /** Paths are collected in reverse order in a vector container. 
                Can be iterated over or popped in O(1). The first node is the destination.
                Uncomment below to see the layout of the above graph. 
                Change the origin, the algorithm or the destination. */

        // struct printALG { string out(BestPathAlgorithm alg) { switch (alg) {
        //               case BFS: return "BFS"; case DIJ: return "DIJ"; 
        //               case DFSi:return "DFSi"; case DFSr:return "DFSr"; }}} PA;

        // origin = M;
        // algorithm = DFSr;
        // cT = graph.connectionTree(origin, algorithm, None, L);
        
        // cout << "print routes: origin: " << origin->id << " algorithm: " << PA.out(algorithm) << endl;;
        // auto iter = cT.discovered.begin();
        // auto iter_end = cT.discovered.end();
        // while (iter != iter_end) { 
        //       cout << "\t" << cT.getDistanceTo(iter->first) << ": ";
        //       for (auto step : *cT.getPathTo(iter->first)) {
        //               cout << " " << step->id;
        //       } cout << endl;
        //       iter++;
        // }

        
        /** connectionTree() IMPLICATIONS: 

                (1) The situation is different if the graph is directed, 
                but you want to compute backwards. 
                The walker can only walk along the direction of 
                the edges; such a connection tree will be empty.

                (2) There is also an experimental depth parameter. By default it is set to 0, zero,
                which means the full depth of the graph. Think of peeling the layers, how far into the graph it will go. */

        assert(graph.removeEdgeByEndpointsIds("A", "B")); assert(graph.removeEdgeByEndpointsPtrs(A, E));
        assert(graph.removeEdgeByEndpointsIds("A", "F")); assert(graph.removeEdgeByEndpointsPtrs(B, C));
        assert(graph.removeEdgeByEndpointsIds("B", "F")); assert(graph.removeEdgeByEndpointsPtrs(C, D));
        assert(graph.removeEdgeByEndpointsIds("C", "G")); assert(graph.removeEdgeByEndpointsPtrs(D, G));
        assert(graph.removeEdgeByEndpointsIds("D", "H")); assert(graph.removeEdgeByEndpointsPtrs(E, F));
        assert(graph.removeEdgeByEndpointsIds("E", "I")); assert(graph.removeEdgeByEndpointsPtrs(F, I));
        assert(graph.removeEdgeByEndpointsIds("G", "J")); assert(graph.removeEdgeByEndpointsPtrs(G, K));
        assert(graph.removeEdgeByEndpointsIds("G", "L")); assert(graph.removeEdgeByEndpointsPtrs(H, L));
        assert(graph.removeEdgeByEndpointsIds("I", "J")); assert(graph.removeEdgeByEndpointsPtrs(I, M));
        assert(graph.removeEdgeByEndpointsIds("I", "N")); assert(graph.removeEdgeByEndpointsPtrs(J, K));
        assert(graph.removeEdgeByEndpointsIds("K", "N")); assert(graph.removeEdgeByEndpointsPtrs(K, O));
        assert(graph.removeEdgeByEndpointsIds("L", "P")); assert(graph.removeEdgeByEndpointsPtrs(M, N));
        assert(graph.edgeCount() == 0);

        assert(graph.removeVertexByPtr(A)); assert(!graph.removeVertexById("A"));
        assert(graph.removeVertexById("B")); assert(!graph.removeVertexByPtr(B));
        assert(graph.removeVertexByPtr(C)); assert(!graph.removeVertexById("C"));
        assert(graph.removeVertexById("D")); assert(!graph.removeVertexByPtr(D));
        assert(graph.removeVertexByPtr(E)); assert(!graph.removeVertexById("E"));
        assert(graph.removeVertexById("F")); assert(!graph.removeVertexByPtr(F));
        assert(graph.removeVertexByPtr(G)); assert(!graph.removeVertexById("G"));
        assert(graph.removeVertexById("H")); assert(!graph.removeVertexByPtr(H));
        assert(graph.removeVertexByPtr(I)); assert(!graph.removeVertexById("I"));
        assert(graph.removeVertexById("J")); assert(!graph.removeVertexByPtr(J));
        assert(graph.removeVertexByPtr(K)); assert(!graph.removeVertexById("K"));
        assert(graph.removeVertexById("L")); assert(!graph.removeVertexByPtr(L));
        assert(graph.removeVertexByPtr(M)); assert(!graph.removeVertexById("M"));
        assert(graph.removeVertexById("N")); assert(!graph.removeVertexByPtr(N));
        assert(graph.removeVertexByPtr(O)); assert(!graph.removeVertexById("O"));
        assert(graph.removeVertexById("P")); assert(!graph.removeVertexByPtr(P));

        assert(graph.vertexCount() == graph.edgeCount());

        cout << "OK\n";
}

void general() {
        cout << "UNDIRECTED GRAPH general() TEST .. " << std::flush;

        auto graph = UDON<u64, char, int, undirected, weighted>();
        auto v1 = graph.insertVertex(111, 0);
        auto v2 = graph.insertVertex(222, 0);
        auto v3 = graph.insertVertex(333, 0);
        
        assert(v1->id == 111);
        assert(v2->id == 222);

        assert(graph.containsVertexByPtr(v3));
        assert(graph.vertexCount() == 3);

        auto e12 = graph.insertEdge(v1, v2, 12, 0.12);
        assert(e12.origin == v1);
        assert(e12.destination == v2);
        assert(e12.id == 12);
        assert(e12.weight == 0.12);
        assert(graph.edgeCount() == 1);

        assert(graph.containsEdge(v1, v2));
        assert(graph.containsEdgeIfEdge(e12));

        // MERGE
        decltype(graph) graph2;
        graph2.mergeIntoSelf(&graph);

        // Since graph2 was empty, had no same Vertices with graph, the latter will be emptied
        assert(graph2.vertexCount() == 3);
        assert(graph2.containsEdgeIfEdge(e12));
        assert(graph.vertexCount() == 0);

        auto beef = graph2.insertVertex(0xBEEF, 0);
        auto baaf = graph2.insertVertex(0xBAAF, 0);
        auto boof = graph2.insertEdge(beef, baaf, 0xBEEF + 0xBAAF, double(0xBEEF + 0xBAAF)/10);
        assert(graph2.vertexCount() == 3 + 2);
        assert(graph2.edgeCount() == 2);

        assert(graph2.containsVertexById(0xBEEF));
        assert(graph2.containsVertexById(0xBAAF));

        auto baba = graph2.insertVertex(0xBABA, 0);
        auto abba = graph2.insertVertex(0xABBA, 0);
        assert(graph.vertexCount() == 0);
        assert(graph2.vertexCount() == 3 + 2 + 2);

        assert(graph2.degree(beef, outgoing) == 1);
        assert(graph2.degree(baba, outgoing) = 0);

        auto hex = graph2.insertEdge(baba, abba, 0xAB, 0.8);
        auto e13 = graph2.insertEdge(baba, v1, 13, 0.13);
        auto e14 = graph2.insertEdge(baba, v2, 14, 0.14);
        auto e15 = graph2.insertEdge(baba, v3, 15, 0.15);     

        assert(graph2.extractEdgeByEndpointsPtrs(baba, abba).value() == hex);
        assert(graph2.extractEdgeByEndpointsPtrs(baba, v1).value() == e13);
        assert(graph2.extractEdgeByEndpointsPtrs(baba, v2).value() == e14);
        assert(graph2.extractEdgeByEndpointsPtrs(baba, v3).value() == e15);

        assert(graph2.extractEdgeByEndpointsIds(0xBEEF, 0xBAAF).value() == boof);
        assert(graph2.extractEdgeByEndpointsIds(0xBEEFEE, 0xBAAF) == None);
        assert(graph2.extractEdgeByEndpointsIds(0xBABA, 0xABBA).value() == hex);
        assert(graph2.extractEdgeByEndpointsIds(0xBABA, 0xABBAC) == None);

        assert(graph2.edgeCount() == 6);

        // Vertices can be iterated over using an iterator
        auto v_iter = graph2.vertices(); int c = 0;
        while ( auto vtx = v_iter.next()) {
                assert(graph2.containsVertexByPtr(vtx.value()));
                c++;
        } assert(c == graph2.vertexCount());
        
        // Use reset() to start over
        v_iter.reset(); while ( auto vtx = v_iter.next()) {
                assert(graph2.containsVertexById(vtx.value()->id));
        }

        // Query adjacent Vertices
        auto av = graph2.adjacentVertices(baba);
        assert(av->size() == 4);
        for (auto vtx : *av ) {
                bool as = (vtx == abba) | (vtx == v1) | (vtx == v2) | (vtx == v3);
                assert(as);
        }

        // Query Incident Edges
        auto ie = graph2.incidentEdges(baba, outgoing);
        assert(ie->size() == 4);
        for (auto edge : *ie) {
                bool as = (*edge == hex) | (*edge == e13) | (*edge == e14) | (*edge == e15);
                assert(as);
        
        }

        // All egdes can be gathered into a container using two methods
        // with respect to their ids
        auto edges = graph2.mapEdgesToSet();
        assert(edges->size() == graph2.edgeCount());
        for (auto Item : *edges) {
                assert(graph2.containsEdgeIfEdge(Item.second));
                assert(graph2.containsEdge(Item.second.origin, Item.second.destination));
        }

        // with respect to their endpoints
        auto edges2 = graph2.mapEndpointsToSet();
        assert(edges2->size() == graph2.edgeCount());
        for (auto endpoint : *edges2) {
                assert(graph2.containsEdge(endpoint.origin, endpoint.destination));
                assert(graph2.containsEdgeIfEdge(
                        graph2.extractEdgeByEndpointsPtrs(endpoint.origin, endpoint.destination).value())
                );
        }

        assert(graph2.removeVertexByPtr(baba)); assert(!graph2.removeVertexByPtr(baba));
        assert(!graph2.removeVertexById(0xBABA)); assert(!graph2.containsVertexByPtr(baba));
        assert(!graph2.containsVertexById(0xBABA)); assert(!graph2.containsEdge(baba, abba));
        assert(!graph2.containsEdgeIfEdge(e13)); assert(!graph2.containsEdgeIfEdge(e14));
        assert(!graph2.containsEdgeIfEdge(e15)); assert(graph2.edgeCount() == 2);
        assert(graph2.containsVertexByPtr(v1)); assert(graph2.containsVertexById(222));
        assert(graph2.containsVertexByPtr(v3)); assert(graph2.removeEdgeByEndpointsPtrs(beef, baaf));
        assert(!graph2.containsEdgeIfEdge(boof)); 
        assert(graph2.containsVertexById(0xBEEF) == graph2.containsVertexById(0xBAAF));

        
        // Create new Graph with string type ids
        char data = '0';

        auto graphSTR = UDON<string, char*, u64, undirected, weighted>();
        auto LAX = graphSTR.insertVertex("LAX", &data);
        auto TEX = graphSTR.insertVertex("TEX", &data);
        auto LT = graphSTR.insertEdge(LAX, TEX, 110011, 0.110011);

        assert(graphSTR.getVertexData("LAX") == graphSTR.getVertexData("TEX"));

        graphSTR.insertEdge(LAX, TEX, 110011, 0.110011);
        graphSTR.insertEdge(LAX, TEX, 110011, 0.110011);
        graphSTR.insertEdge(LAX, TEX, 110011, 0.110011);
        assert(graphSTR.edgeCount() == 1);

        assert(graphSTR.containsVertexById("LAX"));
        assert(graphSTR.containsVertexById("TEX"));

        assert(graphSTR.extractVertexById("LAX") == *LAX);
        assert(graphSTR.extractVertexPtrById("LAX") == LAX);
        assert(graphSTR.extractVertexPtrById("MUX") == None);

        assert(graphSTR.isOutgoingTo(LAX, TEX));
        assert(graphSTR.isIncomingTo(LAX, TEX));


        assert(graphSTR.extractEdgeByEndpointsIds("LAX", "TEX").value() == LT);

        // updating vertex data field
        char new_data = '1';
        graphSTR.updateVertexData("LAX", &new_data);
        graphSTR.updateVertexDataByPtr(TEX, &new_data);
        assert(graphSTR.extractVertexById("LAX")->data == &new_data);
        assert(graphSTR.extractVertexPtrById("TEX").value()->data == &new_data);
        assert(graphSTR.getVertexData("LAX") == graphSTR.getVertexData("TEX"));
        
        assert(graphSTR.containsEdge(LAX, TEX));
        assert(graphSTR.removeEdge(LT));
        assert(!graphSTR.removeEdge(LT));
        assert(graphSTR.containsVertexByPtr(LAX));
        assert(graphSTR.containsVertexByPtr(TEX));
        assert(graphSTR.removeVertexByPtr(LAX));
        assert(!graphSTR.removeVertexByPtr(LAX));
        assert(graphSTR.removeVertexByPtr(TEX));
        assert(!graphSTR.removeVertexByPtr(TEX));
        assert(graphSTR.vertexCount() == 0);

        cout << "OK\n";
}

int main() {
        general();
        shortestPath();
        topologicalSort();
        primJarnikKruskal();
        

        EXIT_SUCCESS;
}
