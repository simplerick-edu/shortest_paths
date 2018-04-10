#include <cstdio>
#include <cassert>
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <cmath>
// #include <chrono>


using namespace std;

// External vector of size 2 - for forward and backward search.
// Internal 2-dimensional vector is vector of adjacency lists for each node.
typedef vector<vector<vector<int>>> Adj;

// Distances can grow out of int type
typedef long long Len;

// Vector of two priority queues - for forward and backward searches.
// Each priority queue stores the closest unprocessed node in its head.
typedef vector<priority_queue<pair<double, int>,vector<pair<double,int>>,greater<pair<double,int>>>> Queue;

const Len INF = numeric_limits<Len>::max()/4;

const double pi = 3.14159265359;

class AStar {
    // See the descriptions of these fields in the starter for friend_suggestion
    int n_;
    Adj adj_;
    Adj cost_;
    vector<vector<Len>> distance_;
    vector<int> workset_;
    vector<vector<bool>> visited_;
    // Coordinates of the nodes
    std::vector<std::pair<Len,Len>> xy_;
    vector<int> Landmarks;
    vector<vector<double>> dist;
    // vector<vector<double>> pf; 


public:
    AStar(int n, Adj adj, Adj cost, std::vector<std::pair<Len,Len>> xy)
        : n_(n), adj_(adj), cost_(cost), distance_(2, vector<Len>(n, INF)), visited_(2, vector<bool>(n, false)), xy_(xy), Landmarks(16), dist(16,vector<double>(n))//,  pf(2, vector<double>(n,-1))
    { workset_.reserve(n); }

    //Precomputing

    // Choose landmarks
    void choose_landmarks() {
        Len x_c;     //mass center
        Len y_c;
        for (int i = 0; i < n_; i++) {
            x_c += xy_[i].first;
            y_c += xy_[i].second;
        }
        x_c = x_c/n_;
        y_c = y_c/n_;

        vector<double> tg;
        for (int k =0; k<Landmarks.size()/4; k++) {
            tg.push_back(tan(k*pi/8));
        }
        vector<double> prev_dist (Landmarks.size(), 0);
        for (int i = 0; i < n_; i++) {
            double x = xy_[i].first;
            double y = xy_[i].second;

            int quart = 0; 

            if(x>0) {
                if(y>0) quart = 0;
                else quart = 3;
            } else {
                if(y>0) quart = 1;
                else quart = 2;
            }

            for(int k = tg.size()-1; k >= 0; k--) {
                if (abs(y/x)>tg[k]) {
                    double new_dist =  sqrt( pow((xy_[i].first-x_c),2) + pow((xy_[i].second-y_c),2) );
                    if ( new_dist > prev_dist[quart*4+k] ) {
                        Landmarks[quart*4+k] = i;
                        prev_dist[quart*4+k] = new_dist;
                        break;
                    }

                }
            }
        }

    }

    //Compute distances to landmarks
    void precompute_dist() {
        for (int i = 0; i < n_; i++) {
            for (int j = 0; j< Landmarks.size(); j++) {
                dist[j][i] = sqrt( pow((xy_[i].first-xy_[Landmarks[j]].first),2) + pow((xy_[i].second-xy_[Landmarks[j]].second),2) );
            }
        }
    }



    //Compute potential function for some path
    int compute_landmark(int s, int t) {
        Len p = 0;
        int A = 0;
        for (int j = 0; j< Landmarks.size(); j++)  {
            if (abs(dist[j][s]-dist[j][t]) > p) {
                A = j;
                p = abs(dist[j][s]-dist[j][t]);
            }
        }
        return A;
    }

    // void compute_pf(int A, int v, int t, int s) {
    //     if ( (dist[A][s]-dist[A][t])> 0) {
    //         pf[0][v] = dist[A][v]-dist[A][t];
    //         pf[1][v] = dist[A][s]-dist[A][v];
    //     }
    //     else {
    //         pf[0][v] = dist[A][t]-dist[A][v];
    //         pf[1][v] = dist[A][v]-dist[A][s];
    //     }
    // }

    double pf(int A, int side, int v, int t, int s) {
        if ( (dist[A][s]-dist[A][t])> 0) {
            return (side == 0) ? dist[A][v]-dist[A][t] : dist[A][s]-dist[A][v];
        }
        else {
            return (side == 0) ? dist[A][t]-dist[A][v] : dist[A][v]-dist[A][s];
        }
    }


    // Initialize the data structures before new query,
    // clear the changes made by the previous query.
    void clear() {
        for (int i = 0; i < workset_.size(); ++i) {
            int v = workset_[i];
            distance_[0][v] = distance_[1][v] = INF;
            visited_[0][v] = visited_[1][v] = false;
            // pf[0][v] = pf[1][v] = 0;
        }
        workset_.clear();
    }

    // Processes visit of either forward or backward search 
    // (determined by value of side), to node v.
    void visit(Queue& q, int side, int v, int A,int t, int s) {
        for (int i = 0; i < adj_[side][v].size(); ++i) {
            int w = adj_[side][v][i];
            if (!visited_[side][w]) {
                workset_.push_back(w);
                Len new_distance = distance_[side][v]+cost_[side][v][i];
                if (distance_[side][w] > new_distance) {
                    distance_[side][w] = new_distance;
                    // if (pf[side][w] == -1) compute_pf(A, w,t,s);
                    q[side].push(pair<double,int>(new_distance+pf(A,side,w,t,s), w));
                }
            }
        }
        visited_[side][v] = true;
    }

    // Returns the shortest path when some node
    // was visited both by forward and backward search
    Len shortest_path() {
        Len distance = INF;
        for (auto u : workset_) {
            Len new_distance = distance_[0][u]+distance_[1][u];
            if (new_distance < distance) { 
                distance = new_distance;
            }
        }
        return distance;
    }

    // Returns the distance from s to t in the graph
    Len query(int s, int t) {
        clear();
        if (s == t) {
            return 0;
        }
        int A = compute_landmark(s,t);

        Queue q(2);
        distance_[0][s] = 0;
        distance_[1][t] = 0;
        // compute_pf(A,s,t,s);
        // compute_pf(A,t,t,s);
        visit(q, 0, s, A,t,s);
        visit(q, 1, t, A,t,s);
        workset_.push_back(s);
        workset_.push_back(t);
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
                visit(q, side, v, A,t,s);
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
    std::vector<std::pair<Len,Len>> xy(n);
    for (int i=0;i<n;++i){
        int a, b;
        scanf("%d%d", &a, &b);
        xy[i] = make_pair(a,b);
    }
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

    AStar astar(n, adj, cost, xy);
    astar.choose_landmarks();
    astar.precompute_dist();

    int t;
    scanf("%d", &t);
    for (int i=0; i<t; ++i) {
        int u, v;
        scanf("%d%d", &u, &v);
        printf("%lld\n", astar.query(u-1, v-1));
        // astar.query(u-1, v-1);
        // if ((i+1)%(t/10) == 0) printf("%d %%\n", 100*(i+1)/t);
    }

    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> elapsed = end-start;
    // printf("%f ms\n", elapsed.count());
}
