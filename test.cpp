#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/properties.hpp>

typedef boost::adjacency_list<
    boost::vecS, boost::vecS, boost::undirectedS,
    boost::property<boost::vertex_index_t, int,
        boost::property<boost::vertex_color_t, boost::default_color_type> >,
    boost::property<boost::edge_index_t, int>
> Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor   Edge;
typedef boost::graph_traits<Graph>::edge_iterator     EdgeIter;
typedef std::vector<std::vector<Edge> >               PlanarEmbedding;
typedef boost::iterator_property_map<
    PlanarEmbedding::iterator,
    boost::property_map<Graph, boost::vertex_index_t>::type
> EmbeddingMap;

struct EdgeData {
    int id, source, target;
    double cost, reverse_cost;
};

struct OrderingRow {
    int seq, node_id, prefix_size;
    std::vector<int> neighbors;
};

struct PlanarOrderingResult {
    bool is_planar;
    std::vector<OrderingRow> rows;
    PlanarOrderingResult() : is_planar(false) {}
};

void reindex_edges(Graph& g) {
    int idx = 0; EdgeIter ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei, ++idx)
        put(boost::edge_index, g, *ei, idx);
}

Graph build_graph(const std::vector<EdgeData>& edges,
                  std::map<int, Vertex>& vmap,
                  std::map<Vertex, int>& id_map) {
    Graph g;
    std::set<int> nodes;
    for (size_t i = 0; i < edges.size(); ++i) {
        nodes.insert(edges[i].source);
        nodes.insert(edges[i].target);
    }
    for (std::set<int>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
        Vertex v = add_vertex(g);
        vmap[*it] = v; id_map[v] = *it;
    }
    // deduplicate undirected edges
    std::set<std::pair<int,int> > added;
    for (size_t i = 0; i < edges.size(); ++i) {
        int u = edges[i].source, v = edges[i].target;
        if (u > v) std::swap(u, v);
        if (edges[i].cost > 0 && !added.count(std::make_pair(u,v))) {
            add_edge(vmap[u], vmap[v], g);
            added.insert(std::make_pair(u,v));
        }
    }
    reindex_edges(g);
    return g;
}

bool compute_embedding(Graph& g, PlanarEmbedding& embedding) {
    embedding.assign(num_vertices(g), std::vector<Edge>());
    EmbeddingMap emb(embedding.begin(), get(boost::vertex_index, g));
    return boost::boyer_myrvold_planarity_test(
        boost::boyer_myrvold_params::graph = g,
        boost::boyer_myrvold_params::embedding = emb);
}

std::map<int, std::set<int> > build_adjacency(const std::vector<EdgeData>& edges) {
    std::map<int, std::set<int> > adj;
    for (size_t i = 0; i < edges.size(); ++i)
        if (edges[i].cost > 0) {
            adj[edges[i].source].insert(edges[i].target);
            adj[edges[i].target].insert(edges[i].source);
        }
    return adj;
}

PlanarOrderingResult apply_planar_canonical_ordering(const std::vector<EdgeData>& edges) {
    std::map<int, Vertex> vmap;
    std::map<Vertex, int> id_map;
    Graph g = build_graph(edges, vmap, id_map);

    PlanarOrderingResult result;
    PlanarEmbedding embedding;
    if (!compute_embedding(g, embedding)) return result;
    result.is_planar = true;

    { EmbeddingMap emb(embedding.begin(), get(boost::vertex_index, g));
      boost::make_biconnected_planar(g, emb); }
    reindex_edges(g);
    if (!compute_embedding(g, embedding)) return result;

    { EmbeddingMap emb(embedding.begin(), get(boost::vertex_index, g));
      boost::make_maximal_planar(g, emb); }
    reindex_edges(g);
    if (!compute_embedding(g, embedding)) return result;

    std::vector<Vertex> order;
    { EmbeddingMap emb(embedding.begin(), get(boost::vertex_index, g));
      boost::planar_canonical_ordering(g, emb, std::back_inserter(order)); }

    std::map<int, std::set<int> > orig_adj = build_adjacency(edges);
    std::set<int> placed;
    int prefix_size = 0;

    for (size_t i = 0; i < order.size(); ++i) {
        int nid = id_map[order[i]];
        std::vector<int> nbrs;
        const std::set<int>& adj = orig_adj[nid];
        for (std::set<int>::iterator it = adj.begin(); it != adj.end(); ++it)
            if (placed.count(*it)) nbrs.push_back(*it);
        std::sort(nbrs.begin(), nbrs.end());
        if (nbrs.empty() && !placed.empty()) prefix_size = 0;
        prefix_size++;
        OrderingRow row;
        row.seq = (int)result.rows.size() + 1;
        row.node_id = nid; row.neighbors = nbrs; row.prefix_size = prefix_size;
        result.rows.push_back(row);
        placed.insert(nid);
    }
    return result;
}

int main() {
    std::vector<EdgeData> edges;
    edges.push_back((EdgeData){1,  1,  3,  6,  6});
    edges.push_back((EdgeData){2,  3,  7,  3,  3});
    edges.push_back((EdgeData){3,  2,  4, 17, 17});
    edges.push_back((EdgeData){4,  4,  9, 14, 14});
    edges.push_back((EdgeData){5,  8,  4, 14, 14});
    edges.push_back((EdgeData){6,  7,  8, 10, 10});
    edges.push_back((EdgeData){7,  6,  7,  4,  4});
    edges.push_back((EdgeData){8,  5,  6,  1,  1});
    edges.push_back((EdgeData){9,  7, 11,  8,  8});
    edges.push_back((EdgeData){10,11, 12, 11, 11});
    edges.push_back((EdgeData){11, 8, 12, 12, 12});
    edges.push_back((EdgeData){12,13, 14, 18, 18});
    edges.push_back((EdgeData){13,12, 17, 13, 13});
    edges.push_back((EdgeData){14,11, 16,  9,  9});
    edges.push_back((EdgeData){15,16, 17, 15, 15});
    edges.push_back((EdgeData){16,10, 11,  5,  5});
    edges.push_back((EdgeData){17,10, 15,  3,  3});
    edges.push_back((EdgeData){18,15, 16, 16, 16});
    edges.push_back((EdgeData){19, 6, 10,  2,  2});

    PlanarOrderingResult res = apply_planar_canonical_ordering(edges);

    std::cout << "Planar Canonical Ordering" << std::endl;
    std::cout << "=========================" << std::endl;
    std::cout << "Planar : " << (res.is_planar ? "YES" : "NO") << std::endl;
    if (!res.is_planar) { std::cout << "Not planar — ordering undefined." << std::endl; return 0; }

    // group into components
    std::vector<std::vector<const OrderingRow*> > comps;
    for (size_t i = 0; i < res.rows.size(); ++i) {
        if (res.rows[i].prefix_size == 1)
            comps.push_back(std::vector<const OrderingRow*>());
        comps.back().push_back(&res.rows[i]);
    }

    std::cout << "Total nodes : " << res.rows.size() << std::endl;
    std::cout << "Components  : " << comps.size() << std::endl;

    for (size_t c = 0; c < comps.size(); ++c) {
        const std::vector<const OrderingRow*>& comp = comps[c];
        std::cout << "\nComponent " << c+1 << " (" << comp.size() << " nodes)"
                  << "  anchors: " << comp.front()->node_id
                  << " (v1) -> "   << comp.back()->node_id << " (vn)" << std::endl;
        std::cout << "  sequence : ";
        for (size_t i = 0; i < comp.size(); ++i) { if (i) std::cout << " -> "; std::cout << comp[i]->node_id; }
        std::cout << std::endl;
        std::cout << "  prior neighbors: ";
        for (size_t i = 0; i < comp.size(); ++i) { if (i) std::cout << "  "; std::cout << comp[i]->node_id << ":" << comp[i]->neighbors.size(); }
        std::cout << std::endl;
    }
    return 0;
}