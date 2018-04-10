#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
// #include <chrono>

using namespace std;

// External vector of size 2 - for forward and backward search.
// Internal 2-dimensional vector is vector of adjacency lists for each node.
typedef vector<vector<vector<int>>> Adj;

// Distances can grow out of int type
typedef long long Len;

// Vector of two priority queues - for forward and backward searches.
// Each priority queue stores the closest unprocessed node in its head.
typedef vector<priority_queue<pair<Len, int>,vector<pair<Len,int>>,greater<pair<Len,int>>>> Queue;

const Len INF = numeric_limits<Len>::max() / 4;

class Bidijkstra {
    // Number of nodes
    int n_;
    // Graph adj_[0] and cost_[0] correspond to the initial graph,
    // adj_[1] and cost_[1] correspond to the reversed graph.
    // Graphs are stored as vectors of adjacency lists corresponding
    // to nodes.
    // Adjacency list itself is stored in adj_, and the corresponding
    // edge costs are stored in cost_.
    Adj adj_;
    Adj cost_;
    // distance_[0] stores distances for the forward search,
    // and distance_[1] stores distances for the backward search.
    vector<vector<Len>> distance_;
    // Stores all the nodes visited either by forward or backward search.
    vector<int> processed_;
    // Stores a flag for each node which is True iff the node was visited
    // either by forward or backward search.
    vector<vector<bool>> visited_;

public:
    Bidijkstra(int n, Adj adj, Adj cost)
        : n_(n), adj_(adj), cost_(cost), distance_(2, vector<Len>(n, INF)), visited_(2, vector<bool>(n, false))
    { processed_.reserve(n); }

    // Initialize the data structures before new query,
    // clear the changes made by the previous query.
    void clear() {
        for (int i = 0; i < processed_.size(); ++i) {
            int v = processed_[i];
            distance_[0][v] = distance_[1][v] = INF;
            visited_[0][v] = visited_[1][v] = false;
        }
        processed_.clear();
    }

    // Processes visit of either forward or backward search 
    // (determined by value of side), to node v.
    void visit(Queue& q, int side, int v) {
        for (int i = 0; i < adj_[side][v].size(); ++i) {
            int w = adj_[side][v][i];
            if (!visited_[side][w]) {
                processed_.push_back(w);
                Len new_distance = distance_[side][v]+cost_[side][v][i];
                if (distance_[side][w] > new_distance) {
                    distance_[side][w] = new_distance;
                    q[side].push(pair<Len,int>(new_distance, w));
                }
            }
        }
        visited_[side][v] = true;
    }


    // Returns the shortest path when some node
    // was visited both by forward and backward search
    Len shortest_path() {
        Len distance = INF;
        for (auto u : processed_) {
            Len new_distance = distance_[0][u]+distance_[1][u];
            if (new_distance < distance) { 
                distance = new_distance;
            }
        }
        return distance;
    }


    // Returns the distance from s to t in the graph.
    Len query(int s, int t) {
        clear();
        if (s == t) {
            return 0;
        }
        Queue q(2);
        distance_[0][s] = 0;
        distance_[1][t] = 0;
        visit(q, 0, s);
        visit(q, 1, t);
        processed_.push_back(s);
        processed_.push_back(t);
        while (true) {
            for (int side=0; side <= 1; side++) { 
                int v;
                while (true) {
                    if (q[side].empty()) {
                        return -1;
                    }
                    v = q[side].top().second;
                    q[side].pop();
                    if (!visited_[side][v]) {
                        break;
                    }
                }
                visit(q, side, v);
                if (visited_[int(!bool(side))][v]) {
                    return shortest_path();
                }
            }
        }
        return -1;
    }
};

int main() {

    // auto start = std::chrono::high_resolution_clock::now();
    int n, m;
    scanf("%d%d", &n, &m);
    Adj adj(2, vector<vector<int>>(n));
    Adj cost(2, vector<vector<int>>(n));
    for (int i=0; i<m; ++i) {
        int u, v, c;
        scanf("%d%d%d", &u, &v, &c);
        adj[0][u-1].push_back(v-1);
        cost[0][u-1].push_back(c);
        adj[1][v-1].push_back(u-1);
        cost[1][v-1].push_back(c);
    }

    Bidijkstra bd(n, adj, cost);

    int t;
    scanf("%d", &t);
    for (int i=0; i<t; ++i) {
        int u, v;
        scanf("%d%d", &u, &v);
        printf("%lld\n", bd.query(u-1, v-1));
    }
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> elapsed = end-start;
    // printf("%f ms\n", elapsed.count());
}
