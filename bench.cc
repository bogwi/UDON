// MIT (c) bogwi@rakumail.jp
// This code uses tab=8 convention.
// Compiles with clang++ -std=c++17 or zig c++ compiler

#include "udon.cc"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include <vector>

using namespace std::chrono;
using CLOCK = steady_clock;
using std::cout;
using std::string;

using _ = void;
using u64 = unsigned long long;

std::random_device rd;  
std::mt19937 gen(rd());

template<class Graph>
_ run(Graph* graph, int width) {
        // Connection Tree
        
        // Obtain arbitrary vertex from the graph     
        auto iter = graph->vertices();
        auto origin = iter.next().value();

        // Calculate Connection tree using Dijkstra algorithm
        string title1 = "Dijkstra alg: ";
        cout << title1 << std::flush;

        auto begin = CLOCK::now();
        auto cT = graph->connectionTree(origin, DIJ);
        auto end = CLOCK::now();
        auto elapsed_time = double(duration_cast<milliseconds>(end - begin).count())/1000;
        cout << std::setfill('.') 
                << std::setw(width - title1.size())
                << elapsed_time << "s\n";

        // Calculate Connection tree using Breadth-first algorithm 
        string title2 = "Breadth-First Search alg: "; 
        cout << title2 << std::flush;

        begin = CLOCK::now();
        cT = graph->connectionTree(origin, BFS);
        end = CLOCK::now();
        elapsed_time = double(duration_cast<milliseconds>(end - begin).count())/1000;
        cout << std::setfill('.') 
                << std::setw(width - title2.size()) 
                << elapsed_time << "s\n";
}

auto undirectedRandomGraph (u64 vtx_n, int edges_n, int width) {

        /** Create randomly connected undirected graph 
                with edges_n outgoing edges on average for each vertex.
                This is a speed test rather than an allocation test.
                The Insertion speed depeneds strongly on the size of the vertex.
                The data field can be particurlary large. */

        string title = "building undirected: ";
        cout << title << std::flush;
        auto graph = UDON<u64, string, string, undirected, weighted>();

        auto begin = CLOCK::now();
        for(u64 i=0; i < vtx_n; i++) {
                graph.insertVertex(i, "data");
        } auto end = CLOCK::now();
        auto insertion_time = end - begin;

        std::uniform_int_distribution<u64> distribN(0, vtx_n - 1);
        std::uniform_int_distribution<u64> edgeWeight(1, edges_n);

        auto vertices = graph.verticesToVector();
        begin = CLOCK::now();
        for (auto origin : vertices) {
                for (int i=0; i < edges_n; i++) {
                         graph.insertEdge(origin, vertices[distribN(gen)], "id", edgeWeight(gen));
                }
        }
        end = CLOCK::now();
        vertices.clear();
 
        auto elapsed_time = duration_cast<milliseconds>(end - begin + insertion_time).count();
        auto printTime = double(elapsed_time)/1000;
        auto edge_count = graph.edgeCount();

        cout << std::setfill('.') 
                << std::setw(width 
                             - title.size()
                             - std::to_string(edge_count).size()
                             - std::to_string(printTime).size()
                             + 1)
                << graph.vertexCount() << ':' 
                << edge_count << ':';
                printf("%.3fs\n", printTime);           

        run(&graph, width);

        // Minimum Spanning Tree - Prim-Jarnik algorithm
        string title3 = "MST: Prim-Jarnik alg: "; 
        cout << title3 << std::flush;

        begin = CLOCK::now();
        auto pjMST = graph.primJarnikMST();
        end = CLOCK::now();
        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title3.size())
                << double(elapsed_time)/1000 << "s\n";

        // Minimum Spanning Tree - Kruskal algorithm
        string title4 = "MST: Kruskal alg: ";
        cout << title4 << std::flush;

        begin = CLOCK::now();
        auto kMST = graph.kruskalMST();
        end = CLOCK::now();
        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title4.size())
                << double(elapsed_time)/1000 << "s\n";


        // Iteration and query over all vertices
        string title6 = "Vertices iteration and get query: ";
        cout << title6 << std::flush;

        begin = CLOCK::now();
        auto v_iter = graph.vertices(); 
        while (auto vtx = v_iter.next()) {
                auto _= graph.extractVertexById(vtx.value()->id); 
        }
        end = CLOCK::now();

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title6.size())
                << double(elapsed_time)/1000 << "s\n";
         
        string title7 = "Remove all vertices: ";
        cout << title7 << std::flush;

        v_iter.reset();
        begin = CLOCK::now();
        while (auto vtx = v_iter.next()) {
                graph.removeVertexByPtr(vtx.value()); 
        }
        end = CLOCK::now();
        if (graph.vertexCount() != 0) { 
                cerr << "GRAPH:BENCH:ERROR -> .removeVertexByPtr() failure\n"; exit(EXIT_FAILURE); }

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title7.size())
                << double(elapsed_time)/1000 << "s\n";


        cout << "dealocating: " << std::flush;
        return CLOCK::now();
}

auto directedRandomGraph (u64 vtx_n, int edges_n, int width) {

        /** Create randomly connected directed acyclic graph 
                with edges_n outgoing edges on average for each vertex.
                This is a speed test rather than an allocation test.
                The Insertion speed depeneds strongly on the size of the vertex.
                The data field can be particurlary large. */

        string title = "building directed: ";
        cout << title << std::flush;
        auto graph = UDON<u64, string, string, directed, weighted>();
                
        auto begin = CLOCK::now();
        for(u64 i=0; i < vtx_n; i++) {
                graph.insertVertex(i, "data");
        } auto end = CLOCK::now();
        auto insertion_time = end - begin;

        std::uniform_int_distribution<u64> distribN(0, vtx_n - 1);
        std::uniform_int_distribution<u64> edges(1, 2 * edges_n);

        auto vertices = graph.verticesToVector();
        begin = CLOCK::now();
        for (int i=0; i < vertices.size(); i++) {
                auto origin = vertices[i];
                auto edge_num = edges(gen);
                
                for (int _=0; _ < edge_num; _++) {
                        std::uniform_int_distribution<u64> distribN(i, vtx_n - 1);
                        auto index = distribN(gen);
                        if (index == i) continue; // prevent connecting to itself

                        auto dest = vertices[index];
                        graph.insertEdge(origin, dest, "id", 1);
                }
        } end = CLOCK::now();
        vertices.clear();

        auto elapsed_time = duration_cast<milliseconds>(end - begin + insertion_time).count();
        auto printTime = double(elapsed_time)/1000;

        auto edge_count = graph.edgeCount();

        cout << std::setfill('.') 
                << std::setw(width 
                             - title.size()
                             - std::to_string(edge_count).size()
                             - std::to_string(printTime).size()
                             + 1)
                << graph.vertexCount() << ':' 
                << edge_count << ':';
                printf("%.3fs\n", printTime); 
                
        run(&graph, width);

        // Topological Sort
        string title2 = "Topological sort: ";
        cout << title2 << std::flush;

        begin = CLOCK::now();
        auto tS = graph.topologicalSort();
        end = CLOCK::now();

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width- title2.size())
                << double(elapsed_time)/1000 << "s\n";

        
        // Iteration and query over all vertices
        string title4 = "Vertices iteration and get query: ";
        cout << title4 << std::flush;

        auto v_iter = graph.vertices();
        begin = CLOCK::now();
        while (auto vtx = v_iter.next()) {
                auto _= graph.extractVertexById(vtx.value()->id); 
        } 
        end = CLOCK::now();

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title4.size())
                << double(elapsed_time)/1000 << "s\n";

        // Iteration over all edges in the graph with query
        string title6 = "Edges iteration and get query: ";
        cout << title6 << std::flush;

        auto e_iter = graph.edges();
        begin = CLOCK::now();
        while (auto edge = e_iter.next()) { 
                auto _= graph.extractEdgeByEndpointsPtrs(edge.value()->origin, edge.value()->destination);
        }
        end = CLOCK::now();
                
        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title6.size())
                << double(elapsed_time)/1000 << "s\n";

        string title7 = "Remove all edges: ";
        cout << title7 << std::flush;

        e_iter.reset();
        begin = CLOCK::now();
        while (auto edge = e_iter.next()) {
                graph.removeEdgeByEndpointsPtrs(edge.value()->origin, edge.value()->destination);
        }
        end = CLOCK::now();
        if (graph.edgeCount() != 0) { 
                cerr << "GRAPH:BENCH:ERROR -> .removeEdgeByEndpointsPtrs() failure\n"; exit(EXIT_FAILURE); }

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title7.size())
                << double(elapsed_time)/1000 << "s\n";

        string title8 = "Remove all vertices: ";
        cout << title8 << std::flush;

        v_iter.reset();
        begin = CLOCK::now();
        while (auto vtx = v_iter.next()) {
                graph.removeVertexByPtr(vtx.value()); 
        }
        end = CLOCK::now();
        if (graph.vertexCount() != 0) { 
                cerr << "GRAPH:BENCH:ERROR -> .removeVertexByPtr() failure\n"; exit(EXIT_FAILURE); }

        elapsed_time = duration_cast<milliseconds>(end - begin).count();
        cout << std::setfill('.')
                << std::setw(width - title8.size())
                << double(elapsed_time)/1000 << "s\n";


        cout << "dealocating: " << std::flush;
        return CLOCK::now();
}


int main(int argc, char** argv) {

        string help = "\n\tUDON graph bench usage example: ./a.out [A] [B]\n\
        -----------------------------------------------\n\
        [A] - integer number of vertices\n\
        [B][optional] - integer number of edges per vertex\n\
        Default A = 5000, B = 20\n\n\
        The benchmark connection scheme is probabilistic,\n\
        heavily depending on A and B, where the B parameter\n\
        is expected. If called with A=1000 B=1000,\n\
        the result will be approximately 1000 : 400, or\n\
        B always <= 40%̀ of the total number of vertices.\n\n";

        u64 vtx_n = 1e4;
        int edges_n = 20;
        bool has_v = false;

        if (argc > 3) { cout << help; return 0; }

        if (argc > 1) {
                std::vector<std::string> args(argv + 1, argv + argc);

                for (const auto& arg : args) {
                        if (arg == "-h") { cout << help; return 0; }
                        
                        for (char c : arg) {
                                if (c < '0' || c > '9') {
                                        { cout << help; return 0; }
                                }
                        }
                        u64 temp = 0;
                        for (int i = 0; i < arg.size(); i++) {
                                temp = temp * 10 + (arg[i] - '0');
                        }
                        if (!has_v) {
                                vtx_n = temp; has_v = true;
                        } else if (temp > vtx_n) {
                                edges_n = vtx_n;
                        } else edges_n = (temp < 2) ? 2 : int(temp);
                }
        }
        int width = 60;
        
        string title = "BENCHING UDON GRAPH: connection random: try: "; 
        cout << "\n" << title << std::setfill('.') 
                << std::setw(width
                             - title.size()
                             - std::to_string(edges_n).size()) 
                << vtx_n << ':' << edges_n << endl;

        auto dealloc_s = undirectedRandomGraph(vtx_n, edges_n, width);
        auto dealloc_e = CLOCK::now();
        auto udT = double(duration_cast<milliseconds>(dealloc_e - dealloc_s).count()) / 1000;
        cout << std::setfill('.') 
                << std::setw(width - 13) 
                << udT << 's' << endl << endl; 

        dealloc_s = directedRandomGraph(vtx_n, edges_n, width);
        dealloc_e = CLOCK::now();
        auto ddT = double(duration_cast<milliseconds>(dealloc_e - dealloc_s).count()) / 1000;
        cout << std::setfill('.') 
                << std::setw(width - 13) 
                << udT << 's' << endl << endl; 
        
        EXIT_SUCCESS;

}

