// MIT (c) bogwi@rakumail.jp
// This code uses tab=8 convention.
// Compiles with clang++ -std=c++17 or zig c++ compiler


#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <list>
#include <optional>
#include <ostream>
#include <queue>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::optional;
using std::cerr;
using std::endl;
using u64 = unsigned long long;

#define None std::nullopt // We use None for clarity of the code.
#define INF std::numeric_limits<long double>::max()

template<class VTX_ID>
std::hash<VTX_ID> hasher;


// Hash vertex by its id identifier.
template<class VTX_ID, class Vertex>
struct Hash {
        u64 operator()(Vertex v) const {
                return hasher<VTX_ID>(v.id);
        }
};

// Hash a pair or vertices.
template<class VTX_ID, class Endpoints>
struct Hash2 {
        u64 operator()(Endpoints pair) const {
                return hasher<VTX_ID>(pair.origin->id) 
                        ^ hasher<VTX_ID>(pair.destination->id);
        }
};


// Control graph's type and mode during initiation.
enum GraphType { undirected, directed };
enum GraphMode { unweighted, weighted };


// Used by `incidentEdges()` and `adjacentVertices()` methods.
enum EdgeClass { outgoing, incoming };


// Algorithms used when traversing the graph as a tree.
enum TreeTraverseAlgorithm { bfs, pre, post, ino };
enum TreeTraverseControl { all, except, through };

// Algorithms used to calculate the path from one vertex to another.
// BFS is breadth first. It produces the shortest path by counting 
// the number of vertices you have to go through.
// `DIJ` is the Dijkstra algorithm. It produces the shortest path 
// taking into account the weights of the edges.
// `DFSi` is depth first search.
// `DFSi` is recursive depth-first search.
// The latter two are more likely to explore the possible connections between vertices
// in the graph, rather than finding the shortest path. However, in some circumstances they
// can work as best-path algorithms.
enum PathAlgorithm { BFS, DFSi, DFSr, DIJ };


// The structure of the graph is next:
// Origin_vertex : { Destination_vertex_pointer0 : { Edge0 },
                 // Destination_vertex_pointer1 : { Edge1 },
                 // Destination_vertex_pointer2 : { Edge2 },
// ... }


// ADJACENCY MAP GRAPH VARIANT. 
// `VTX_ID` -> `integer` or `string` else provide your own hasher template.
// `VTX_DAT` -> any type.
// `edgeID` -> any type supporting the equals `==` operator.  
// GraphType -> `undirected` | `directed`.
// GraphMode -> `unweighted` | `weighted`.
template<typename VTX_ID,
        typename VTX_DAT,
        typename edgeID,
        enum GraphType graphType,
        enum GraphMode graphMode>
class UDON {
        public:
                // VTX_ID' -> by default string or integer, will be hashed.
                // If you want the ID parameter to be something else, change the hash function.
                // `VTX_DAT` -> Any type or pointer.
                // template<class VTX_ID, class VTX_DAT>
                class Vertex {
                        public:
                                VTX_ID id;
                                
                                VTX_DAT data;
                                
                                Vertex (VTX_ID id, VTX_DAT data ) : id(id), data(data) {};

                                bool operator==(const Vertex& other ) const {
                                        return id == other.id;
                                }
                                void operator=(const Vertex& other) {
                                        id = other.id;
                                        data = other.data;
                                }
                };

                // Class that manages the endpoints of the edge.
                class Endpoints { 
                        public:
                                const Vertex* origin; 
                                const Vertex* destination; 

                                Endpoints (const Vertex* origin, const Vertex* dest) : 
                                        origin(origin), destination(dest){};

                                bool operator==(const Endpoints& other) const {
                                        return (origin->id == other.origin->id) 
                                        && (destination->id == other.destination->id);
                                }
                                void operator=(const Endpoints& other) {
                                        origin = other.origin;
                                        destination = other.destination;
                                }

                };

                // The edge class does not store vertices, only their pointers.
                // The comparison is based on the vertices defined in the vertex class.
                // Of course edges can have any number of fields, add as many as you need.
                class Edge {
                        public:
                                const Vertex* origin; const Vertex* destination;
                                edgeID id; long double weight=0;

                                Edge () {};

                                Edge (
                                        const Vertex* origin, 
                                        const Vertex* destination, 
                                        edgeID id, 
                                        long double weight) : 
                                        origin(origin), 
                                        destination(destination),
                                        id(id), 
                                        weight(weight) {};

                                bool operator==(const Edge& other) const {
                                        return (origin == other.origin) 
                                                && (destination == other.destination)
                                                && (id == other.id);
                                }

                                void operator=(const Edge& other) {
                                        origin = other.origin;
                                        destination = other.destination;
                                        id = other.id;
                                        weight = other.weight;
                                }
                                
                                Endpoints endpoints() {
                                        Endpoints endpoints_ { origin, destination };
                                        return endpoints_;
                                }

                                const Vertex* opposite(const Vertex* v) {
                                        if (v == origin)
                                                return destination;
                                        return origin;
                                }
                };
        

        using HashMap2 = std::unordered_map<const Vertex*, Edge, std::hash<const Vertex*>>;
        using HashMap1 = std::unordered_map<Vertex, HashMap2, Hash<VTX_ID, Vertex>>;

        private:
                HashMap1* outGoing = new HashMap1;
                HashMap1* inComing; 

                // Release all allocated memory.
                // The function is linear, directly proportional to the size of the graph.
                void deinit() {
                        delete outGoing;
                        if (graphType == directed) {
                                delete inComing;
                        }
                }
        public:
                UDON () {

                        // If the graph is undirected, inComing will point to outGoing.
                        // If the graph is directed, inComing is a separate container.

                        if (graphType == undirected) {
                                inComing = outGoing;
                        } else inComing = new HashMap1;
                };

                ~UDON () {
                        deinit();
                }

                // Clear the graph of all the content.
                void clear() {
                        for ( auto iter = outGoing->begin(); iter != outGoing->end(); ++iter) {
                                iter->second.clear();
                        }
                        outGoing->clear();
                        if (graphType) {
                                for ( auto iter = inComing->begin(); iter != inComing->end(); ++iter) {
                                        iter->second.clear();
                                }
                                inComing->clear();
                        }
                }

                // Merge the other graph into itself.
                // The former will lose all its unique content. 
                void mergeIntoSelf(UDON* other) {
                        outGoing->merge(*other->outGoing);
                        if (graphType == directed) {
                                inComing->merge(*other->inComing);
                        } else inComing = outGoing;
                        
                }

                // Insert a vertex with given id and data, and return a `const` pointer to it.
                // Does not clobber an existing one with the same id.
                const Vertex* insertVertex(VTX_ID id, VTX_DAT data) {
                        Vertex vtx = Vertex(id, data);
                        HashMap2 m2;
                        outGoing->insert({vtx, m2});

                        if (graphType == directed) {
                                HashMap2 m2_;
                                inComing->insert({vtx, m2_});
                        }
                        return &outGoing->find(vtx)->first;
                }

                // Insert a given vertex created from the template. 
                // Does not clobber an existing one with the same id. 
                void insertVertexIfVertex(Vertex vtx) {
                        HashMap2 m2;
                        outGoing->insert({vtx, m2});
                        if (graphType == directed) {
                                HashMap2 m2_;
                                this->inComing->insert({vtx, m2_});
                        }
                }

                // Update the vertex associated with the given id with new data 
                // and return a pointer to the updated version.
                const Vertex* updateVertexData(VTX_ID id, VTX_DAT data) {
                        int codeLine = __LINE__ + 1;
                        if (containsVertexById(id)) {
                                Vertex query = Vertex(id, data);
                                auto nh = outGoing->extract(query);
                                nh.key() = query;
                                outGoing->insert(std::move(nh));
                                return &outGoing->find(query)->first;
                        }
                        cerr << "GRAPH:ERROR -> Missing Vertex with given id: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                }

                // Update the vertex associated with the given pointer with new data 
                // and return a pointer to the updated version.
                const Vertex* updateVertexDataByPtr(const Vertex* vertex, VTX_DAT data) {
                        int codeLine = __LINE__ + 1;
                        if (containsVertexByPtr(vertex)) {
                                Vertex query = Vertex(vertex->id, data);
                                auto nh = outGoing->extract(*vertex);
                                nh.key() = query;
                                outGoing->insert(std::move(nh));
                                return &outGoing->find(query)->first;
                        }
                        cerr << "GRAPH:ERROR -> Missing the given Vertex: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                }

                // Return data field belonging to the vertex associated with the given id,
                // or return None if no such vertex in the graph,      
                optional<VTX_DAT> getVertexData(VTX_ID id) {
                        Vertex query = Vertex(id, {});
                        auto search = outGoing->find(query);
                        if (search != outGoing->end())
                                return search->first.data; 
                        return None;
                }
                
                // Insert edge with given endpoints, id and optional weight. 
                // The method assumes that the graph ALREADY CONTAINS origin and destination vertices. 
                // Otherwise returns a missing vertex error. Will not clobber an existing edge.
                // To insert a vertex, use the `insertVertex()` method.
                Edge insertEdge(
                        const Vertex* origin, 
                        const Vertex* dest,
                        edgeID edge_id,
                        long double weight=0) {

                        Edge edge (origin, dest, edge_id, weight);

                        int codeLine = __LINE__ + 1;
                        if (auto search = outGoing->find(*origin); search != outGoing->end())
                                search->second.insert({dest, edge});
                        else {
                                cerr << "GRAPH:ERROR -> Missing origin Vertex: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                        }
                        codeLine = __LINE__ + 1; 
                        if (auto search = inComing->find(*dest); search != inComing->end())
                                search->second.insert({origin, edge});
                        else {
                                cerr << "GRAPH:ERROR -> Missing destination Vertex: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                        }
                        return edge; 
                }

                // Insert or assign an edge with given endpoints, id and optional weight.
                // This method will clobber an existing edge with new id and weight. 
                // Reports any missing vertices.
                Edge insertOrAssignEdge(
                        const Vertex* origin, 
                        const Vertex* dest,
                        edgeID edge_id,
                        long double weight=0) {
                        
                        Edge edge (origin, dest, edge_id, weight);

                        int codeLine = __LINE__ + 1;
                        if (auto search = outGoing->find(*origin); search != outGoing->end())
                                search->second.insert_or_assign(dest, edge);
                        else {
                                cerr << "GRAPH:ERROR -> Missing origin Vertex: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                        }
                        codeLine = __LINE__ + 1; 
                        if (auto search = inComing->find(*dest); search != inComing->end())
                                search->second.insert_or_assign(origin, edge);
                        else {
                                cerr << "GRAPH:ERROR -> Missing destination Vertex: code line: "
                                << codeLine << "\n" ;
                                exit(EXIT_FAILURE);
                        }
                        return edge; 
                }

                // Remove the given edge from the graph. 
                // Does not remove endpoints. Returns `false` if one of the endpoints is missing.
                // It is always good to check the return.
                bool removeEdge(Edge edge) {
                        bool origin, dest = false;
                        if (auto search = outGoing->find(*edge.origin); search != outGoing->end()){
                                origin = search->second.erase(edge.destination);
                        } else return false;
                        if (auto search = inComing->find(*edge.destination); search != inComing->end()){
                                dest = search->second.erase(edge.origin);
                        } else return false;
                        return (origin && dest) ? true : false;
                }

                // Remove the edge associated with given endpoints pointers from the graph. 
                // Does not remove endpoints. Returns `false` if one of the endpoints is missing.
                // It is always good to check the return.
                bool removeEdgeByEndpointsPtrs(const Vertex* origin, const Vertex* destination) {
                        bool origin_, dest_ = false;
                        if (auto search = outGoing->find(*origin); search != outGoing->end()){
                                origin_ = search->second.erase(destination);
                        } else return false;

                        if (auto search = inComing->find(*destination); search != inComing->end()){
                                dest_ = search->second.erase(origin);
                        } else return false;
                        return (origin_ && dest_) ? true : false;
                }

                // Remove the edge associated with given endpoints ids from the graph. 
                // Does not remove endpoints. Returns `false` if one of the endpoints is missing.
                // It is always good to check the return.
                bool removeEdgeByEndpointsIds(VTX_ID origin, VTX_ID destination) {
                        auto org = extractVertexPtrById(origin);
                        auto dest = extractVertexPtrById(destination);

                        if (org.has_value() && dest.has_value())
                                return removeEdgeByEndpointsPtrs(org.value(), dest.value());
                        return false;
                }

                // Remove the given vertex from the graph.
                // Returns `false` only if the vertex is missing but not an error. 
                // It is always good to check the return.
                bool removeVertexByPtr(const Vertex* vtx) {
                        Vertex vtx_ = *vtx;

                        if (auto adjacent = outGoing->find(vtx_); adjacent != outGoing->end()) {
                                HashMap2 adj2nd = adjacent->second;
                                for (auto Item : adj2nd) {
                                        auto adj = *Item.first;
                                        auto dest = inComing->find(adj);
                                        dest->second.erase(vtx);
                                        // if (auto dest = inComing->find(adj); dest != inComing->end()) {
                                        //         dest->second.erase(vtx);
                                        // }
                                }
                        } else return false;

                        outGoing->erase(vtx_);
                        
                        // We do not check if inComing has the given vertex.
                        // If it passes the first [if block] and the graph is directed, then OK.   
                        if (graphType == directed) {
                                auto adjacent = inComing->find(vtx_);
                                HashMap2 adj2nd = adjacent->second;
                                for (auto Item : adj2nd) {
                                        auto adj = *Item.first;
                                        auto dest = outGoing->find(adj);
                                        dest->second.erase(vtx);
                                        // if (auto dest = outGoing->find(adj); dest != outGoing->end()) {
                                        //         dest->second.erase(vtx);
                                        // }
                                }
                                inComing->erase(vtx_);
                        }
                        return true;   
                }

                // Remove the vertex associated with the given id from the graph.
                // Returns `false` only if the vertex is missing but not an error. 
                // It is always good to check the return.
                bool removeVertexById(VTX_ID id) {
                        auto query = extractVertexPtrById(id);

                        if (query.has_value()) {
                                return removeVertexByPtr(query.value());
                        }
                        return false;
                }

                // Return `true` if the graph contains the given vertex.
                bool containsVertex(Vertex vtx) {
                        auto search = outGoing->find(vtx);
                        return (search != outGoing->end()) ? true : false;
                }

                // Return `true` if the graph contains a vertex associated with the given const pointer.
                bool containsVertexByPtr(const Vertex* vtx) {
                        auto search = outGoing->find(*vtx);
                        return (search != outGoing->end()) ? true : false;
                }

                // Return `true` if the graph contains a vertex associated with the given id.
                bool containsVertexById(VTX_ID id) {
                        Vertex query = Vertex(id, {});
                        return containsVertex(query);
                }

                // Return a vertex associated with the given `id` or return `None` 
                // if there is no such vertex in the graph.
                optional<Vertex> extractVertexById(VTX_ID id) {
                        Vertex query = Vertex(id, {});
                        if (auto result = outGoing->find(query); result != outGoing->end())
                                return result->first;
                        return None;
                }

                // Return a pointer to the vertex associated with the given `id` or return `None` 
                // if there is no such vertex in the graph.
                optional<const Vertex*> extractVertexPtrById(VTX_ID id) {
                        Vertex query = Vertex(id, {});
                        if (auto result = outGoing->find(query); result != outGoing->end())
                                return &result->first;
                        return None;
                }

                // Return the number of all vertices present in the graph.
                u64 vertexCount() {
                        return outGoing->size();
                }

                using VERTICES = std::vector<const Vertex*>;

                struct VertexIterator {
                        UDON* contex;
                        typename HashMap1::iterator begin;
                        // return the next vertex
                        optional<const Vertex*> next() {
                                while (begin != contex->outGoing->end()) {
                                        return &begin++->first;
                                } return None;
                        }
                        // reset iterator to the beginning
                        void reset() { 
                                begin = contex->outGoing->begin();
                        }
                };
                   
                // Return `VertexIterator` to iterate over the vertices present in the graph.
                // Usage: 
                // `auto iter = graph.vertices()`; 
                // `while(auto vtx == iter.next()) { ... do stuff with vxt }`; 
                // Use `reset()` to get it back to the starting point. 
                // This method does not allocate anything.
                VertexIterator vertices() {
                        auto begin = outGoing->begin();
                        return { this, begin };
                }

                // Return `<vector>` filled with pointers to all the vertices present in the graph.
                VERTICES verticesToVector() {
                        VERTICES vertices;
                        for ( auto iter = outGoing->begin(); iter != outGoing->end(); ++iter) {
                                vertices.push_back(&iter->first);
                        }
                        return vertices;
                }

                // Return `true` if the graph contains an edge in between given endpoints.
                bool containsEdge(const Vertex* origin, const Vertex* dest) {
                        
                        if (auto search = outGoing->find(*origin); search != outGoing->end()) {
                                auto search2 = search->second.find(dest);
                                        if (search2 != search->second.end()) {
                                        return true;
                                } else {
                                        return false;
                                }
                        } else {
                                return false;
                        }
                }

                // Return `true` if the graph contains the given edge.
                bool containsEdgeIfEdge(Edge edge) {
                        return containsEdge(edge.origin, edge.destination);
                }

                // Return an edge with given endpoints' pointers or `None` if there is no such edge
                // in the graph or one of the endpoints is missing. This function can be used
                // instead of `containsEdge()` granted you need the return.
                optional<Edge> extractEdgeByEndpointsPtrs(const Vertex* origin, const Vertex* dest) {

                        if (auto search = outGoing->find(*origin); search != outGoing->end()) {
                                if (auto search2 = search->second.find(dest); 
                                        search2 != search->second.end()) {
                                                return search2->second;
                                } else return None;
                        } else return None;
                }

                // Return an edge with given endpoints' ids or `None` if there is no such 
                // in the graph or one of the endpoints is missing.
                optional<Edge> extractEdgeByEndpointsIds(VTX_ID origin_id, VTX_ID dest_id) {
                        Vertex origin = Vertex(origin_id, {});
                        optional<const Vertex*> dest_ptr = extractVertexPtrById(dest_id);
                        if (dest_ptr.has_value()) {
                                if (auto search = outGoing->find(origin); 
                                        search != outGoing->end()) {
                                                return search->second.find(dest_ptr.value())->second;
                                } else return None;
                        } return None;
                }

                // Return the number of all edges in the graph.
                u64 edgeCount() {
                        u64 numOfEdges = 0;
                        for (auto Item : *outGoing) {
                                numOfEdges += Item.second.size();
                        }
                        return (graphType == directed) ? numOfEdges : numOfEdges / 2;
                }

                // Return the number of all edges incident to the given vertex. 
                // The default parameter for edges is `outgoing`.
                optional<u64> degree(const Vertex* vtx, EdgeClass edge_class=outgoing) {
                        if (edge_class == incoming) {
                                if (auto search = inComing->find(*vtx); search != inComing->end()) {
                                        return search->second.size();
                                } else return None;
                        }
                        if (auto search = outGoing->find(*vtx); search != outGoing->end()) {
                                        return search->second.size();
                                } else return None;
                }
                
                using EDGES = std::vector<Edge*>;
                // Return `<vector>` with pointers to each incident edge
                // associated with the given vertex in the graph, 
                // The default parameter for `EdgeClass` is `outgoing`. 
                // Returns `None` if the given vertex is missing.
                optional<EDGES> incidentEdges(const Vertex* vtx, EdgeClass edge_class=outgoing) {
                        EDGES edges;
                        if (edge_class == incoming) {
                                if (auto incident = inComing->find(*vtx); 
                                incident != inComing->end()) {
                                        for (auto iter = incident->second.begin(); 
                                        iter != incident->second.end(); ++iter) {
                                                edges.push_back(&iter->second);
                                        } return edges;
                                } else return None;
                        }
                        if (auto incident = outGoing->find(*vtx); 
                                incident != outGoing->end()) {
                                        for (auto iter = incident->second.begin(); 
                                        iter != incident->second.end(); ++iter) {
                                                edges.push_back(&iter->second);
                                        } return edges;
                                } else return None;
                }

                // Return `<vector>` with pointers to each adjacent vertex
                // associated with the given vertex in the graph, 
                // The default parameter for `EdgeClass` is `outgoing`. 
                // Returns `None` if the given vertex is missing.
                optional<VERTICES> adjacentVertices(const Vertex* vtx, EdgeClass edge_class=outgoing) {
                        VERTICES vertices;
                        if (edge_class == incoming) {
                                if (auto adjacent = inComing->find(*vtx); adjacent != inComing->end()) {
                                        for ( auto Item: adjacent->second) {
                                                vertices.push_back(Item.first);
                                        }
                                        return vertices;
                                } else return None;
                        }
                        if (auto adjacent = outGoing->find(*vtx); adjacent != outGoing->end()) {
                                        for (auto Item: adjacent->second) {
                                                vertices.push_back(Item.first);
                                        }
                                return vertices;
                        } else return None;
                }

                // Return `true` if the `tested` vertex is outgoing to the `origin` vertex. 
                bool isOutgoingTo(const Vertex* origin, const Vertex* tested) {
                        if (auto out = outGoing->find(*origin); out != outGoing->end()) {
                                if (auto answer = out->second.find(tested);
                                                answer != out->second.end())
                                        return true;
                                return false;
                        } return false;
                }

                // Return `true` if the `tested` vertex is incoming to the `origin` vertex.
                bool isIncomingTo(const Vertex* origin, const Vertex* tested) {
                        if (auto in = inComing->find(*origin); in != inComing->end()) {
                                if (auto answer = in->second.find(tested);
                                                answer != in->second.end())
                                        return true;
                                return false;
                        } return false;
                }

                // Return `true` if the `tested` vertex is outgoing to the `origin` vertex,
                // using vertex's ids.
                bool isOutgoingToById(VTX_ID origin_id, VTX_ID tested_id) {
                        Vertex origin(origin_id, {});
                        auto tested = extractVertexPtrById(tested_id);
                        if (tested.has_value()) {
                                if (auto out = outGoing->find(origin); out != outGoing->end()) {
                                        if (auto answer = out->second.find(tested.value());
                                                        answer != out->second.end())
                                                return true;
                                        return false;
                                } return false;
                        } return false;
                }

                // Return `true` if the `tested` vertex is incoming to the `origin` vertex,
                // using vertex's ids.
                bool isIncomingToById(VTX_ID origin_id, VTX_ID tested_id) {
                        Vertex origin(origin_id, {});
                        auto tested = extractVertexPtrById(tested_id);
                        if (tested.has_value()) {
                                if (auto out = inComing->find(origin); out != inComing->end()) {
                                        if (auto answer = out->second.find(tested.value());
                                                        answer != out->second.end())
                                                return true;
                                        return false;
                                } return false;
                        } return false;
                }

                using EDGES_MAP = std::unordered_map<edgeID, Edge>;
                // Map all edges present in the graph with respect to their ids. 
                // Valid ids are integers and string literals. 
                // If the graph is `undirected`, only unique edges will be mapped.
                // Returns None if the graph has no edges.
                // This method is the only exception where the edge's id must be hashable.
                // Otherwise add a custom hasher to `EDGES_MAP` container or
                // check `mapEndpointsToSet()` method.      
                optional<EDGES_MAP> mapEdgesToSet() { 
                        static_assert (
                                std::is_integral_v<edgeID> == true or
                                typeid(edgeID) == typeid(std::string),
                                "Edge ID is ineligible for hashing. Provide custom hasher for EDGES_MAP container."
                        );

                        EDGES_MAP set;
                        // for every origin in the graph
                        for (auto origin : *outGoing) {
                                // for every destination in the graph
                                for (auto dest : origin.second) {
                                        // insert { edge.id, edge }
                                        set.insert( {dest.second.id, dest.second} );
                                }
                        }
                        if (set.size() > 0) return set;
                        return None;
                }

                using ENDPOINTS_SET = std::unordered_set<Endpoints, Hash2<VTX_ID, Endpoints>>;
                // Map all edges present in the graph with respect to their endpoints.
                // If the graph is `undirected`, only unique edges will be mapped.
                // Returns None if the graph has no edges.
                optional<ENDPOINTS_SET> mapEndpointsToSet() {
                        ENDPOINTS_SET set;
                        // for every origin in the graph
                        for (auto origin : *outGoing) {
                                // for every destination in the graph
                                for (auto dest : origin.second) {
                                        // insert { edge.endpoints() }
                                        set.insert( {dest.second.endpoints()} );
                                }
                        }
                        if (set.size() > 0) return set;
                        return None;     
                }
                
                struct EdgeIterator {
                        UDON* contex;
                        typename HashMap1::iterator begin;
                        typename HashMap2::iterator begin2;
                        // return the next edge
                        optional<Edge*> next() {
                                while (begin != contex->outGoing->end()) {
                                        while (begin2 != begin->second.end()) {
                                                return &begin2++->second;
                                        } 
                                        begin++;
                                        if (begin != contex->outGoing->end()) {
                                                begin2 = begin->second.begin(); 
                                        }
                                } return None;
                        }
                        // reset iterator to the beginning
                        void reset() {
                                begin = contex->outGoing->begin();
                                begin2 = begin->second.begin();
                        }
                };
                // Return `EdgeIterator` to iterate over the edges present in the graph.
                // Usage: 
                // `auto iter = graph.edges()`; 
                // `while(auto edge == iter.next()) { ... do stuff with edge }`; 
                // Use `reset()` to get it back to the starting point. 
                // Does not allocate anything.
                // *IMPORTANT* 
                // To prevent reporting each edge twise,
                // this method works only with directed graphs.
                EdgeIterator edges() {
                        static_assert(graphType == directed, "The graph is not directed.");
                        auto begin = outGoing->begin();
                        auto begin2 = begin->second.begin();
                        return { this, begin, begin2 };
                }


                /// CLASSIC TRAVERSING ALGORITHMS  ///

                using Path = std::vector<const Vertex*>;
                using Tq = std::list<const Vertex*>;
                
                // Traverse the graph from the given origin and return a path `<vector>` with pointers to vertices(nodes).
                // Asserts that the graph is directed. Cannot check if the graph contains cycles. 
                // Inorder traversing only takes first 2 child nodes, 
                // so general graphs will have skips with inorder algorithm.
                // Supports 4 classic algorithms where the codes `bf, pre, post, ino` are passed as parameters. 
                // Example: `auto path = graph.traverseGraphAsTree(A, bfs)`.
                // 
                // *Implications*. 
                // - This graph does not support insertion order and has no knowledge of child line-up.
                // - It will most likely not report child nodes in the exact way they were inserted.
                // - The report is rather right-handed.
                Path traverseGraphAsTree(
                        const Vertex* origin, 
                        TreeTraverseAlgorithm algorithm,
                        optional<const Vertex*> target=None) {

                        static_assert(graphType == directed,
                                      "This method is used only with directed graphs.");
                        
                        int codeLine = __LINE__ + 1;
                        
                        codeLine = __LINE__ + 1;
                        if (!containsVertexByPtr(origin)) {
                                cerr << "GRAPH:ERROR -> Missing origin Vertex: code line: "
                                << codeLine << "\n" ; 
                                exit(EXIT_FAILURE);
                        }
                        if (target.has_value()) {
                                codeLine = __LINE__ + 1;
                                if (!containsVertexByPtr(target.value())) {
                                        cerr << "GRAPH:ERROR -> Missing target Vertex: code line: "
                                        << codeLine << "\n" ; 
                                        exit(EXIT_FAILURE);
                                }
                        }

                        Path path;

                        TreeTraverse result = { &path, this, target };
                        switch (algorithm) {
                                case bfs: result.breadthFirst(origin); break;
                                case pre: result.preorder(origin); break;
                                case post: result.postOrder(origin); break;
                                case ino: result.inorder(origin); break;
                        } return path;
                }
        private:
                struct TreeTraverse {
                        Path* path;
                        UDON* cnt;
                        optional<const Vertex*> target;
                        bool stop = false;

                        void breadthFirst(const Vertex* origin) {
                                Tq queue;
                                queue.push_back(origin);

                                if (!target.has_value()) {
                                        while (queue.size() > 0) {
                                                auto elder = queue.front(); queue.pop_front();
                                                path->push_back(elder);
                                                optional<VERTICES> children = cnt->adjacentVertices(elder, outgoing);
                                                if (children.has_value()) {
                                                        for (auto vtx : children.value()) {
                                                                queue.push_back(vtx);
                                                        }
                                                }
                                        }
                                } else {
                                        while (queue.size() > 0) {
                                                auto elder = queue.front(); queue.pop_front();
                                                path->push_back(elder);
                                                if (elder->id == target.value()->id)
                                                        break;
                                                optional<VERTICES> children = cnt->adjacentVertices(elder, outgoing);
                                                if (children.has_value()) {
                                                        for (auto vtx : children.value()) {
                                                                queue.push_back(vtx);
                                                        }
                                                }
                                        }
                                }
                        }
                        void preorder(const Vertex* origin) {
                                path->push_back(origin);
                                optional<VERTICES> children = cnt->adjacentVertices(origin, outgoing);
                                if (children.has_value()) {
                                        for (auto vtx : children.value()) {
                                                preorder(vtx);
                                        }
                                }
                        }
                        void preorderIfTarget(const Vertex* origin) {
                                if ((origin->id != target.value()->id) && !stop) {
                                        path->push_back(origin);
                                } else if (origin->id == target.value()->id) {
                                        path->push_back(origin);
                                        stop = true;
                                        return;
                                } else return;

                                optional<VERTICES> children = cnt->adjacentVertices(origin, outgoing);
                                if (children.has_value()) {
                                        for (auto vtx : children.value()) {
                                                preorderIfTarget(vtx);
                                        }
                                }
                        }
                        void postOrder(const Vertex* origin) {
                                optional<VERTICES> children = cnt->adjacentVertices(origin, outgoing);
                                if (children.has_value()) {
                                        for (auto vtx : children.value()) {
                                                postOrder(vtx);
                                        }
                                        path->push_back(origin);
                                }
                        }
                        void postOrderIfTarget(const Vertex* origin) {
                                optional<VERTICES> children = cnt->adjacentVertices(origin, outgoing);
                                if (children.has_value()) {
                                        for (auto vtx : children.value()) {
                                                postOrder(vtx);
                                        }
                                        if ((origin->id != target.value()->id) && !stop) {
                                                path->push_back(origin);
                                        } else if (origin->id == target.value()->id) {
                                                path->push_back(origin);
                                                stop = true;
                                                return;
                                        } else return;
                                }
                        }
                        void inorder(const Vertex* origin) {
                                optional<VERTICES> kids = cnt->adjacentVertices(origin, outgoing);
                                if (kids.has_value()) {
                                        if (kids.value().size() > 0)
                                                inorder(kids.value()[0]);
                                        path->push_back(origin);
                                        if (kids.value().size() > 1)
                                                inorder(kids.value()[1]);
                                }
                        }
                        void inorderIfTarget(const Vertex* origin) {
                                optional<VERTICES> kids = cnt->adjacentVertices(origin, outgoing);
                                if (kids.has_value()) {
                                        if (kids.value().size() > 0)
                                                inorderIfTarget(kids.value()[0]);

                                        if ((origin->id != target.value()->id) && !stop) {
                                                path->push_back(origin);
                                        } else if (origin->id == target.value()->id) {
                                                path->push_back(origin);
                                                stop = true;
                                                return;
                                        } else return;

                                        if (kids.value().size() > 1)
                                                inorderIfTarget(kids.value()[1]);
                                }
                        }
                };
 
                // AUX types and containers.

                struct WAE { 
                        public:
                                // TODO: store edge pointers not edges.
                                long double weight;
                                Edge edge; 
                                WAE(long double weight, Edge edge) : 
                                        weight(weight), edge(edge) {};

                                void operator=(const WAE& other) {
                                        weight = other.weight;
                                        edge = other.edge;
                                }
                };
                using HashMap3 = std::unordered_map<const Vertex*, WAE, std::hash<const Vertex *>>;
                struct Element { 
                        public:
                                long double weight; 
                                const Vertex* vertex;

                                Element (long double weight, const Vertex* vertex) : 
                                        weight(weight), vertex(vertex) {}

                                bool operator<(const Element& other) const {
                                        return weight < other.weight;
                                }
                                bool operator>(const Element& other) const {
                                        return weight > other.weight;
                                }
                        
                };
                using minQ = std::priority_queue<Element, std::vector<Element>, std::greater<Element>>;

        public:
                struct Connections {
                        private:
                                const Vertex* last_lookup;
                                bool dij = false;
                                
                        public:
                                const Vertex* origin;
                                Path path{};
                                HashMap3 discovered;

                                Connections () {}
                                Connections (const Vertex* origin,
                                Path path,
                                HashMap3 discovered,
                                const Vertex* last_lookup, 
                                bool dij) : 
                                origin(origin), path(path), discovered(discovered), 
                                last_lookup(last_lookup), dij(dij) {}

                                // Return the distance from the origin vertex 
                                // at which connectionTree() was evaluated
                                // to the given destination vertex.
                                // The `DIJ` algorithm returns  the value with respect to edges' weight.
                                // The `BIJ` algorithm returns the length of the shortest path.
                                long double getDistanceTo(const Vertex* dest) {
                                        if (auto Item = discovered.find(dest); Item != discovered.end()) {
                                                checkPath(dest);
                                                if (dij) 
                                                        return Item->second.weight;
                                                return path.size();
                                        }
                                        return 0;
                                }

                                // Return `true` if the destination vertex has been found.
                                bool checkFound(const Vertex* dest) {
                                        checkPath(dest);
                                        return path.begin() != path.end();
                                }

                                // Return a path to the given destination vertex. 
                                // The `Path` container is filled backwards, so that
                                // it can be popped if individual vertices are to be consumed.
                                Path* getPathTo(const Vertex* dest) { 
                                        checkPath(dest);
                                        return &path;
                                }
                        private:
                                // Check if a path for the given destination exist 
                                // to avoid unnecessary re-calculation.
                                void checkPath(const Vertex* dest) {
                                        if (last_lookup != dest | path.size() == 0) {
                                                path.clear();
                                                constructPath(dest);
                                                last_lookup = dest;
                                                return;
                                        }
                                        return;
                                }
                                // Construct the actual path with respect to the destination vertex 
                                // and the algorithm used by the `connectionTree()` method.
                                void constructPath(const Vertex* dest) {
                                        if ( auto Item = discovered.find(dest); Item != discovered.end()) {
                                                path.push_back(dest);
                                                const Vertex* walk = dest;
                                                
                                                WAE* wae = &Item->second;
                                                // ? Dijkstra algorithm case : Breadth-first 
                                                if (dij) {
                                                        if (wae->weight != 0 and wae->weight != INF) {
                                                                while (walk != origin)  {
                                                                        if (auto Item = discovered.find(walk); 
                                                                        Item != discovered.end()) {
                                                                        Edge edge = Item->second.edge;
                                                                        const Vertex* parent = edge.opposite(walk);

                                                                        path.push_back(parent);
                                                                        walk = parent;
                                                                        } else { break;}
                                                                }
                                                        } else { Item->second.weight = 0; }
                                                } else {
                                                        while (walk != origin)  {
                                                                if (auto Item = discovered.find(walk); 
                                                                Item != discovered.end()) {
                                                                Edge edge = Item->second.edge;
                                                                const Vertex* parent = edge.opposite(walk);

                                                                path.push_back(parent);
                                                                walk = parent;
                                                                } else { break;}
                                                        }
                                                }
                                                path.pop_back();
                                        }
                                }
                };

                using VERTEX_ID_SET = std::unordered_set<VTX_ID>;
                // Compute a connection tree from origin to destination vertices
                // with respect to the given `PathAlgorithm`. See `test.cc` for usage.
                Connections connectionTree(
                        const Vertex* origin,
                        PathAlgorithm algorithm,
                        optional<VERTEX_ID_SET*> knockout=None,
                        optional<const Vertex*> target=None,
                        TreeTraverseControl control=all, 
                        int depth=0) {

                        Connections connections;
                        switch (algorithm) {
                                case BFS:
                                        connections = BFS_DFSi_ALL(
                                                origin, knockout,
                                                target, depth,
                                                control, algorithm); break;
                                case DFSi:
                                        connections = BFS_DFSi_ALL(
                                                origin, knockout,
                                                target, depth,
                                                control, algorithm); break;
                                case DIJ:
                                        connections = dijkstra_ALL(
                                                origin, knockout,
                                                target, depth,
                                                control); break;
                                case DFSr:
                                        connections = DFSr_ALL(
                                                origin, knockout,
                                                target, depth,
                                                control); break;
                        }
                        return connections;
                }

        private:
                Connections BFS_DFSi_ALL(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        optional<const Vertex*> target,
                        int depth,
                        TreeTraverseControl control,
                        PathAlgorithm algorithm
                ) { switch (control) { 
                                case all:
                                        return bfsDfs(origin, knockout, false, target, depth, algorithm);
                                        break;
                                case except:
                                        return bfsDfs(origin, knockout, false, target, depth, algorithm);
                                        break;
                                case through:
                                        return bfsDfs(origin, knockout, true, target, depth, algorithm);
                                        break;
                        }
                }
                Connections bfsDfs(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        bool reflect,
                        optional<const Vertex*> target,
                        int depth,
                        PathAlgorithm algorithm
                ) {
                        HashMap3 discovered;
                        Path path;
                        Tq found;
                        found.push_back(origin);

                        if (target.has_value()) {
                                if (depth == 0) {
                                        while (found.size() > 0) {
                                                const Vertex* origin__;
                                                switch (algorithm) {
                                                        case DFSi:
                                                                origin__ = found.back();
                                                                found.pop_back(); break;
                                                        default:
                                                                origin__ = found.front();
                                                                found.pop_front(); break;
                                                }
                                                bfsDfsLoop(origin__, knockout,  &discovered, reflect, &found);
                                                if (target.value() == origin__) {break;}
                                        }
                                } else {
                                        int depth_ = depth;
                                        while (found.size() > 0 and depth_ > 0) {
                                                const Vertex* origin__;
                                                switch (algorithm) {
                                                        case DFSi:
                                                                origin__ = found.back();
                                                                found.pop_back(); break;
                                                        default:
                                                                origin__ = found.front();
                                                                found.pop_front(); break;
                                                }
                                                bfsDfsLoop(origin__, knockout,  &discovered, reflect, &found);
                                                depth_ -= 1;
                                                if (target.value() == origin__) break;
                                        }
                                }
                        } else {
                                if (depth == 0) {
                                        while (found.size() > 0) {
                                                const Vertex* origin__;
                                                switch (algorithm) {
                                                        case DFSi:
                                                                origin__ = found.back();
                                                                found.pop_back(); break;
                                                        default:
                                                                origin__ = found.front();
                                                                found.pop_front(); break;
                                                }
                                                bfsDfsLoop(origin__, knockout,  &discovered, reflect, &found);
                                        }
                                } else {
                                        int depth_ = depth;
                                        while (found.size() > 0 and depth_ > 0) {
                                                const Vertex* origin__;
                                                switch (algorithm) {
                                                        case DFSi:
                                                                origin__ = found.back();
                                                                found.pop_back(); break;
                                                        default:
                                                                origin__ = found.front();
                                                                found.pop_front(); break;
                                                }
                                                bfsDfsLoop(origin__, knockout,  &discovered, reflect, &found);
                                                depth_ -= 1;
                                        }
                                }
                        }
                        discovered.erase(origin);
                        
                        Connections result( origin, path, discovered, origin, false ); 
                        return result;
                }
                void bfsDfsLoop(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        HashMap3* discovered,
                        bool reflect,
                        Tq* found
                ) {
                        EDGES edges = incidentEdges(origin, outgoing).value();
                        
                        for (u64 i=0; i < edges.size(); i++) {
                                const Vertex* vtx = edges[i]->opposite(origin);
                                if (knockout.has_value() 
                                        and !reflect 
                                        and (knockout.value()->find(vtx->id) != knockout.value()->end()))
                                        {continue;}
                                if (knockout.has_value() 
                                        and reflect 
                                        and !(knockout.value()->find(vtx->id) != knockout.value()->end()))
                                        {continue;}

                                if (discovered->find(vtx) == discovered->end()) {
                                        WAE wae(0, *edges[i]);
                                        discovered->insert({vtx, wae});
                                        found->push_back(vtx);
                                }
                        }
                }

                Connections dijkstra_ALL(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        optional<const Vertex*> target,
                        int depth,
                        TreeTraverseControl control
                ) { 
                        switch (control) { 
                                case all:
                                        return dijkstra(origin, knockout, false, target, depth);
                                case except:
                                        return dijkstra(origin, knockout, false, target, depth);
                                case through:
                                        return dijkstra(origin, knockout, true, target, depth);
                        }
                }
                Connections dijkstra(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        bool reflect,
                        optional<const Vertex*> target,
                        int depth
                ) {      
                        static_assert(graphMode == weighted,
                                      "Dijkstra algorithm works only with weighted graphs");

                        // possible run time protection
                        if (graphMode == unweighted) {
                                cerr << "GRAPH:ERROR -> Dijkstra Algorithms but unweighted graph\n";
                                exit(EXIT_FAILURE);
                        }

                        HashMap3 discovered;
                        Path path;
                        minQ minHeap;
                        long double inf = std::numeric_limits<long double >::max();

                        for (auto Item = outGoing->begin(); Item != outGoing->end(); ++Item) {
                                VTX_ID vtx_id = Item->first.id;
                                if (knockout.has_value() 
                                        and !reflect 
                                        and (knockout.value()->find(vtx_id) != knockout.value()->end()))
                                        {continue;}
                                if (knockout.has_value() 
                                        and reflect 
                                        and !(knockout.value()->find(vtx_id) != knockout.value()->end()))
                                        {continue;}

                                Edge edge;
                                WAE wae(inf, edge);
                                discovered.insert({ &Item->first, wae });
                        }
                        
                        Edge origin_edge;
                        WAE wae(0, origin_edge);
                        discovered.insert_or_assign( origin, wae );

                        Element element( 0, origin );
                        minHeap.push(element);

                        if (target.has_value()) {
                                if (depth == 0) {
                                        while (!minHeap.empty()) {
                                                auto fringe = minHeap.top().vertex;
                                                minHeap.pop();
                                                dijkstraLoop(fringe, knockout, &discovered, reflect, &minHeap);
                                                if (target.value() == fringe) break;
                                        }
                                } else {
                                        int depth_ = depth;
                                        while (!minHeap.empty() and depth_ > 0) {
                                                auto fringe = minHeap.top().vertex;
                                                minHeap.pop();
                                                dijkstraLoop(fringe, knockout, &discovered, reflect, &minHeap);
                                                if (target.value() == fringe) break;
                                                depth_ -= 1;
                                        }
                                }
                        } else { 
                                if (depth == 0) {
                                        while (!minHeap.empty()) {
                                                const Vertex* fringe = minHeap.top().vertex;
                                                minHeap.pop();
                                                dijkstraLoop(fringe, knockout, &discovered, reflect, &minHeap);
                                        }
                                } else {
                                        int depth_ = depth;
                                        while (!minHeap.empty() and depth_ > 0) {
                                                auto fringe = minHeap.top().vertex;
                                                minHeap.pop();
                                                dijkstraLoop(fringe, knockout, &discovered, reflect, &minHeap);
                                                depth_ -= 1;
                                        }
                                }
                        } 

                        discovered.erase(origin);
                        Connections result( origin, path, discovered, origin, true ); 
                        return result;

                };

                int codeLine = __LINE__ + 2;
                void dijkstraLoop(
                        const Vertex* fringe,
                        optional<VERTEX_ID_SET*> knockout,
                        HashMap3* discovered,
                        bool reflect,
                        minQ* minHeap) {
                                if (auto adjacent = outGoing->find(*fringe); adjacent != outGoing->end()) {
                                        for ( auto Item: adjacent->second) {
                                                const Vertex* vtx = Item.first;
                                                if (knockout.has_value() 
                                                        and !reflect
                                                        and (knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}
                                                if (knockout.has_value() 
                                                        and reflect
                                                        and !(knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}

                                                Edge edge = Item.second;
                                                auto distance = discovered->find(fringe)->second.weight 
                                                        + edge.weight;

                                                WAE* wae = &discovered->find(vtx)->second;
                                                if (wae->weight > distance) {
                                                        wae->weight = distance;
                                                        
                                                        wae->edge = edge;

                                                        Element element(distance, vtx);
                                                        minHeap->push(element);
                                                }
                                        }
                                } else {
                                        cerr << "GRAPH:ERROR -> Dijkstra Algorithm missing vertex: code line :" 
                                                << codeLine << endl;
                                        exit(EXIT_FAILURE);
                                }
                        }

                Connections DFSr_ALL(
                        const Vertex* origin,
                        optional<VERTEX_ID_SET*> knockout,
                        optional<const Vertex*> target,
                        int depth,
                        TreeTraverseControl control
                ) { 
                        Path path;
                        HashMap3 discovered;
                        switch (control) {
                                case all:
                                        dfsR(origin, None, &discovered, false, target, depth);
                                case except:
                                        dfsR(origin, knockout, &discovered, false, target, depth);
                                case through:
                                        dfsR(origin, knockout, &discovered, true, target, depth);
                        }

                        discovered.erase(origin);
                        Connections result( origin, path, discovered, origin, false ); 
                        return result;
                        }

                        void dfsR(
                                const Vertex* origin,
                                optional<VERTEX_ID_SET*> knockout,
                                HashMap3* discovered,
                                bool reflect,
                                optional<const Vertex*> target,
                                int depth
                        ) {
                                int depth_ = depth;
                                auto run = dfsRALL(
                                        this,
                                        knockout,
                                        discovered,
                                        reflect,
                                        target,
                                        &depth
                                );
                                if (target.has_value()) {
                                        if (depth == 0) {
                                                return run.recur0(origin);
                                        } else return run.recur1(origin);
                                } else {
                                        if (depth == 0) {
                                                return run.recur2(origin);
                                        } else return run.recur3(origin);
                                }
                        }

                        struct dfsRALL {
                        public:
                                UDON* self;
                                optional<VERTEX_ID_SET*> knockout;
                                HashMap3* discovered;
                                bool reflect;
                                optional<const Vertex*> target;
                                int* depth;

                                dfsRALL (
                                        UDON* self,
                                        optional<VERTEX_ID_SET*> knockout,
                                        HashMap3* discovered,
                                        bool reflect,
                                        optional<const Vertex*> target,
                                        int* depth) : self(self), knockout(knockout), discovered(discovered),
                                        reflect(reflect), target(target), depth(depth) {}

                                // Target, no depth
                                void recur0(const Vertex* origin) {
                                        EDGES edges = self->incidentEdges(origin, outgoing).value();
                                        for (auto edge: edges) {
                                                const Vertex* vtx = edge->opposite(origin);
                                                if (knockout.has_value() 
                                                        and !reflect
                                                        and (knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}
                                                if (knockout.has_value() 
                                                        and reflect
                                                        and !(knockout.value()->find(vtx->id)
                                                                != knockout.value()->end()))
                                                        {continue;}

                                                if (discovered->find(vtx) == discovered->end()) {
                                                        WAE wae(0, *edge);
                                                        discovered->insert({vtx, wae});
                                                        if (target.value() == vtx) {return;}

                                                        recur0(vtx);
                                                } 
                                        }
                                };

                                // Target, depth
                                void recur1(const Vertex* origin) {
                                        if (*depth <= 0) return;
                                        EDGES edges = self->incidentEdges(origin, outgoing).value();
                                        for (auto edge: edges) {
                                                const Vertex* vtx = edge->opposite(origin);
                                                if (knockout.has_value() 
                                                        and !reflect
                                                        and (knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}
                                                if (knockout.has_value() 
                                                        and reflect
                                                        and !(knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}

                                                if (discovered->find(vtx) == discovered->end()) {
                                                        WAE wae(0, *edge);
                                                        discovered->insert({vtx, wae});
                                                        if (target.value() == vtx) {return;}

                                                        depth -= 1;
                                                        recur1(vtx);
                                                } 
                                        }
                                };

                                // No target, no depth
                                void recur2(const Vertex* origin) {
                                        EDGES edges = self->incidentEdges(origin, outgoing).value();
                                        for (auto edge: edges) {
                                                const Vertex* vtx = edge->opposite(origin);
                                                if (knockout.has_value() 
                                                        and !reflect
                                                        and (knockout.value()->find(vtx->id)
                                                                != knockout.value()->end()))
                                                        {continue;}
                                                if (knockout.has_value() 
                                                        and reflect
                                                        and !(knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}

                                                if (discovered->find(vtx) == discovered->end()) {
                                                        WAE wae(0, *edge);
                                                        discovered->insert({vtx, wae});
                                                        recur2(vtx);
                                                } 
                                        }
                                };
                                // No target, depth
                                void recur3(const Vertex* origin) {
                                        if (*depth <= 0) return;
                                        EDGES edges = self->incidentEdges(origin, outgoing).value();
                                        for (auto edge: edges) {
                                                const Vertex* vtx = edge->opposite(origin);
                                                if (knockout.has_value()
                                                        and !reflect    
                                                        and (knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}
                                                if (knockout.has_value() 
                                                        and reflect
                                                        and !(knockout.value()->find(vtx->id) 
                                                                != knockout.value()->end()))
                                                        {continue;}

                                                if (discovered->find(vtx) == discovered->end()) {
                                                        WAE wae(0, *edge);
                                                        discovered->insert({vtx, wae});

                                                        depth -= 1;
                                                        recur3(vtx);
                                                } 
                                        }
                                };
                };

                // AUX types and containers

                using Degree = std::unordered_map<const Vertex*, u64, std::hash<const Vertex*>>;
                struct WAE2 { decltype(INF) weight; optional<Edge> edge; 
                        bool operator<(const WAE2& other) const {
                                return weight < other.weight;
                        }
                        bool operator>(const WAE2& other) const {
                                return weight > other.weight;
                        }
                };
                using minQ2 = std::priority_queue<WAE2, std::vector<WAE2>, std::greater<WAE2>>;
                using HashMap5 = std::unordered_map<const Vertex*, WAE2, std::hash<const Vertex*>>;

        public:
                struct Topo {
                        public:
                        // `topo` variable is a vector container. 
                        // You can query it by an index or iterate over.
                        Path topo;
                        // Tells if the graph is acyclic.
                        bool acyclic = false;
                };

                // Since this graph is an implementation of an adjacency map, the topological sort is
                // an arbitrary sort, reporting vertices as they sit in the map, of course 
                // preserving the topological constraints that vertices that depend on other vertices
                // cannot be reported before those on which they depend.
                Topo topologicalSort() {
                        codeLine = __LINE__ + 1;
                        if (vertexCount() == 0)
                                {cerr << "GRAPH:ERROR -> empty graph: code line: " << codeLine << endl;}
 
                        static_assert(graphType == directed,
                                      "This implementation of Topological Sort works only on directed graphs!");

                        // possible run time protection
                        codeLine = __LINE__ + 1; 
                        if (graphType != directed)
                                {cerr << "GRAPH:ERROR -> This implementation of Topological Sort works only on directed graphs!" 
                                        << codeLine << endl;
                                exit(EXIT_FAILURE);}

                        // Vector containing the sorting result
                        Path topo;
                        // Linked list, collector for vertices with no incoming edges.
                        Tq freeVtx;
                        // Map, collector of degree numbers for each vertex.
                        Degree degrees;

                        VertexIterator vertices_ = vertices();
                        while (auto vtx = vertices_.next()) {
                                optional<u64> degree_ = degree(vtx.value(), incoming); 
                                degrees.insert({vtx.value(), degree_.value()});

                                // Take vertices with no incoming edges
                                if (degree_ == 0) {
                                        freeVtx.push_back(vtx.value());
                                }
                        }

                        while (freeVtx.size() > 0) {
                                // Release the vertex and add it to the result
                                const Vertex* vtx = freeVtx.back(); freeVtx.pop_back();
                                topo.push_back(vtx);

                                // Take all outgoing neighbors of just outgoing vertex
                                EDGES edges = incidentEdges(vtx, outgoing).value();
                                for (auto edge : edges) {
                                        const Vertex* vtx_ = edge->opposite(vtx);
                                        u64* degree_ = &degrees.find(vtx_)->second;
                                        *degree_ -= 1; // reduce their degree by one

                                        // If it happens to be a vertex with zero degree, put it into queue
                                        if (*degree_ == 0) {
                                                freeVtx.push_back(vtx_);
                                        }
                                }
                        }
                        bool acyclic = (topo.size() != vertexCount()) ? false : true;
                        return { topo, acyclic };
                }
                
                // Minimum Spaning tree struct. See info for each field.
                struct MST {
                        // `.cost` -> cost of the tree
                        decltype(INF) cost;

                        // `.tree` -> a set container where:
                        // ... `key: Endpoints { vertex_origin*, vertex_destination* },`
                        ENDPOINTS_SET tree;

                        // `.len` -> size of the tree
                        u64 len;

                        public:
                        // Query if MST contains the given Pair of vertices by vertex pointers.
                        bool contains(const Vertex* a, const Vertex* b) {
                                Endpoints pair = { a, b };
                                auto query = tree.find(pair);
                                return (query != tree.end()) ? true : false;
                        }
                };
                // Calculate the minimum spanning tree of the graph, using the 
                // Prim-Jarnik Algorithm.
                // The more edges each vertex has, the slover this implementation of the algorithm 
                // runs compared to the Kruskal algorithm.
                MST primJarnikMST() {
                        codeLine = __LINE__ + 1;
                        if (vertexCount() == 0)
                                {cerr << "GRAPH:ERROR -> empty graph: code line: " << codeLine << endl;
                                EXIT_FAILURE;}
                        
                        static_assert(graphType == undirected && graphMode == weighted,
                                      "Prim-Jarnik algorithm works only on " 
                                        "undirected weighted graphs!");

                        codeLine = __LINE__ + 1;
                        if (graphType != undirected && graphMode != weighted)
                                {cerr << "GRAPH:ERROR -> Prim-Jarnik MST algorithm works only on "
                                        "undirected weighted graphs: code line: " << codeLine << endl;
                                EXIT_FAILURE;}

                        // Tree and cost 
                        ENDPOINTS_SET tree;
                        decltype(INF) cost = 0;

                        // aux data
                        HashMap5 discovered;
                        minQ minHeap;

                        // RUN
                        VertexIterator vertices_ = vertices();

                        const Vertex* first = vertices_.next().value();
                        WAE2 wae = { 0, None };
                        discovered.insert( {first, wae} );
                
                        while (auto vtx = vertices_.next()) {
                                discovered.insert({ vtx.value(), { INF, None}});
                        }
                        
                        Element element = { 0, first };
                        minHeap.push(element); 
                        
                        while (!minHeap.empty()) {
                                Element fringe = minHeap.top(); minHeap.pop();

                                if (auto Item = discovered.find(fringe.vertex); 
                                        Item != discovered.end()) {

                                        if (auto edge = Item->second.edge; edge != None) {
                                                Endpoints pair = edge->endpoints();
                                                tree.insert({ pair });
                                                cost += edge->weight;
                                        }
                                        discovered.erase(fringe.vertex);
                                }

                                HashMap2 adjacent = outGoing->find(*fringe.vertex)->second;
                                auto iter = adjacent.begin();                        
                                
                                for (auto iter = adjacent.begin(); 
                                        iter != adjacent.end(); ++iter) {
                                        // adjacent vertex      
                                        const Vertex* vtx = iter->first;
                                        // corresponding edge
                                        Edge* edge = &iter->second;

                                        auto entry = discovered.find(vtx);

                                        if (entry != discovered.end()) {
                                                if (edge->weight < entry->second.weight) {
                                                        WAE2 wae_ = { edge->weight, *edge };
                                                        entry->second = wae_;

                                                        Element element_ { edge->weight, vtx };
                                                        minHeap.push(element_);
                                                }
                                        }
                                }
                        }
                        return  { cost, tree, tree.size() };
                }
                // Calculate the minimum spanning tree of the graph, using the
                // Kruskal Algorithm.          
                MST kruskalMST () {
                        codeLine = __LINE__ + 1;
                        if (vertexCount() == 0)
                                {cerr << "GRAPH:ERROR -> empty graph: code line: " << codeLine << endl;
                                EXIT_FAILURE;}

                        static_assert(graphType == undirected && graphMode == weighted,
                                      "Kruskal algorithm works only on " 
                                        "undirected weighted graphs!");

                        codeLine = __LINE__ + 1;
                        if (graphType != undirected && graphMode != weighted)
                                {cerr << "GRAPH:ERROR -> Kruskal MST algorithm works only on " 
                                        "undirected weighted graphs: code line: " << codeLine << endl;
                                EXIT_FAILURE;}

                        // Type that manages clusters of vertices and finds a 
                        // min-edge in between them.
                        struct Partition {
                                struct Position {
                                        u64 size;
                                        Position* parent;
                                        Position () {}
                                };
                                Position* find(Position* p) {
                                        if (p->parent != p) {
                                                p->parent = find(p->parent);
                                        }
                                        return p->parent;
                                }
                                void union_(Position* p, Position* q) {
                                        if (p != q) {
                                                if (p->size < q->size) {
                                                        q->parent = p;
                                                        p->size += q->size;
                                                } else {
                                                        p->parent = q;
                                                        q->size += p->size;
                                                }
                                        }
                                }
                        };

                        ENDPOINTS_SET tree;
                        minQ2 pq;
                        Partition forest{};
                        std::unordered_map<
                                const Vertex*, typename Partition::Position*, 
                                std::hash<const Vertex*>> position;
                        
                        decltype(INF) cost = 0;

                        // the algorithm
                        VertexIterator vertices_ = vertices();
                        while (auto vtx = vertices_.next()) {
                                auto* pos = new typename Partition::Position;
                                pos->size = 1;
                                pos->parent = pos;

                                position.insert({vtx.value(), pos});

                                auto i_edges = incidentEdges(vtx.value(), outgoing).value();

                                for (auto edge: i_edges) {
                                        WAE2 wae = {edge->weight, *edge};
                                        pq.push(wae);
                                }
                        }

                        auto size = vertexCount();
                        while ((tree.size() != size - 1) and !pq.empty()) {
                                WAE2 item = pq.top(); pq.pop();
                                auto weight = item.weight;
                                Edge edge = item.edge.value();

                                const Vertex* o = edge.origin;
                                const Vertex* d = edge.destination;

                                auto p_ = position.find(o)->second;
                                auto q_ = position.find(d)->second;
                                auto p = forest.find(p_);
                                auto q = forest.find(q_);

                                if (p != q) {
                                        Endpoints pair = edge.endpoints();
                                        tree.insert({pair});
                                        cost += weight;
                                        forest.union_(p, q);
                                }
                        }
                        // we delete allocated Partition::Position
                        for (auto item: position) delete item.second;
                        return { cost, tree, tree.size() };
                };
};
