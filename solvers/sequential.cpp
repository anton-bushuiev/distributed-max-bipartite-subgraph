#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <fstream>

using namespace std;


struct Graph {
    unsigned n_nodes;
    unsigned n_edges;
    vector<vector<int>> edge_list;
    vector<int> edge_weights;

    // Returns coloring with 2 colors if graph is bipartite else null
    // Note: assumes that graph is connected
    bool is_bipartite(vector<int> & colors) const {
        // Construct neighbors list
        vector<vector<int>> neighbors(n_nodes, vector<int>());
        int u, v;
        for (int i = 0; i < n_edges; i++) {
            u = edge_list[i][0];
            v = edge_list[i][1];
            neighbors[u].push_back(v);
            neighbors[v].push_back(u);
        }

        // Prepare to BFS
        vector<bool> visited(n_nodes, false);
        queue<int> q;
        q.push(0);
        vector<int> t(n_nodes, -1);
        t[0] = 0;

        // Run BFS
        int n;
        while (!q.empty()) {
            n = q.front();
            q.pop();
            visited[n] = true;
            for (int neighbor: neighbors[n]) {
                if (visited[neighbor] == false) {
                    q.push(neighbor);
                    t[neighbor] = t[n] + 1;
                }
            }
        }

        // Color layers with 2 colors
        for (int i = 0; i < n_nodes; i++)
            colors[i] = (t[i] % 2 == 0 ? 0 : 1);

        // Check coloring correctness
        for (int i = 0; i < n_nodes; i++)
            for (int neighbor: neighbors[i])
                if (colors[i] == colors[neighbor])
                    return false;
        return true;
    }

    bool is_subgraph_connencted(const vector<bool> & edge_mask) const {
        // Construct node set and neighbors list of subgraph
        unordered_set<int> nodes_sub;
        vector<vector<int>> neighbors_sub(n_nodes, vector<int>());
        int u, v;
        for (int i = 0; i < n_edges; i++)
            if (edge_mask[i] == true){
                u = edge_list[i][0];
                v = edge_list[i][1];
                neighbors_sub[u].push_back(v);
                neighbors_sub[v].push_back(u);
                nodes_sub.insert(u);
                nodes_sub.insert(v);
            }

        if (nodes_sub.empty())
            return true;

        // Prepare to BFS
        vector<bool> visited(n_nodes, false);
        queue<int> q;
        q.push(*begin(nodes_sub));

        // Run BFS
        int n;
        while (!q.empty()) {
            n = q.front();
            q.pop();
            visited[n] = true;
            for (int neighbor: neighbors_sub[n])
                if (visited[neighbor] == false)
                    q.push(neighbor);
        }

        // Check if all nodes were visited
        for (int n: nodes_sub)
            if (!visited[n])
                return false;
        return true;
    }
};


// Note: for sorting only
struct WeightedEdge {
    vector<int> edge;
    int weight;
};

bool operator> (const WeightedEdge & a, const WeightedEdge & b) {
    return (a.weight > b.weight);
}


void bbdfs(Graph & g, vector<vector<bool>> & edge_masks_max,
           vector<vector<int>> & colors_max, int & score_max, int & n_calls,
           int edge, vector<bool> edge_mask, vector<int> colors,
           vector<vector<int>> neighbors)
{
    n_calls += 1;

    // Calculate solution score
    int score = 0;
    for (int i = 0; i < g.n_edges; i++) {
        if (edge_mask[i] == true) {
            score += g.edge_weights[i];
        }
    }

    // ANALYZE TERMINAL STATE
    if (edge >= g.n_edges) {

        // Check connectivity
        if (!g.is_subgraph_connencted(edge_mask))
            return;

        // Update optimal solution
        if (score == score_max) {
            edge_masks_max.push_back(edge_mask);
            colors_max.push_back(colors);
        }
        else if (score > score_max) {
            score_max = score;
            edge_masks_max = {edge_mask};
            colors_max = {colors};
        }
        return;
    }

    // PRUNE
    // Calculate upper bound on child solutions
    int score_future = 0;
    for (int i = edge; i < g.n_edges; i++) {
            score_future += g.edge_weights[i];
    }

    // Prune if better child solution is impossible
    if (score + score_future < score_max)
        return;

    // BRANCH
    vector<int> colors_prev = colors;
    vector<vector<int>> neighbors_prev = neighbors;

    // "Add edge" case
    edge_mask[edge] = true;
    neighbors[g.edge_list[edge][1]].push_back(g.edge_list[edge][0]);
    neighbors[g.edge_list[edge][0]].push_back(g.edge_list[edge][1]);

    // Try to color (u, v) and (v, u) edge directions with (0, 1) colors
    int u, v;
    for (int d : {0, 1}) {
        u = g.edge_list[edge][d];
        v = g.edge_list[edge][d == 0 ? 1 : 0];

        // Not to color if already colored with opposite
        if (colors_prev[u] == 1 ||
            colors_prev[v] == 0)
            continue;

        // Color
        colors[u] = 0;
        colors[v] = 1;

        // Check coloring consistency
        bool is_consistent = true;
        for (int n: {u, v}) {
            for (int neighbor: neighbors[n]) {
                if (colors[n] == colors[neighbor]) {
                    is_consistent = false;
                    break;
                }
            }
        }
        if (is_consistent == true) {
            bbdfs(g, edge_masks_max, colors_max, score_max, n_calls,
                  edge + 1, edge_mask, colors, neighbors);
        }
    }

    // "Ignore edge" case
    edge_mask[edge] = false;
    bbdfs(g, edge_masks_max, colors_max, score_max, n_calls,
          edge + 1, edge_mask, colors_prev, neighbors_prev);
}


int main(int argc, char* argv[]) {

    // Read command line arguments
    // argument 1 - filename
    if (argc < 2) {
        cout << "No input filename specified." << endl;
        return 1;
    }
    ifstream infile(argv[1]);

    // argument 2 - flag for edge sorting by weights
    bool sort_edges = false;
    if (argc >= 3) {
        if ((string)argv[2] == "no-sort")
            sort_edges = false;
        else if ((string)argv[2] == "sort")
            sort_edges = true;
    }

    // Parse input graph
    Graph g;
    infile >> g.n_nodes;

    int a;
    for (int i = 0; i < g.n_nodes; i++) {
        for (int j = 0; j < g.n_nodes; j++) {
            infile >> a;
            if (a != 0 && j < i) {
                g.edge_list.push_back({i, j});
                g.edge_weights.push_back(a);
            }
        }
    }
    g.n_edges = (unsigned) g.edge_list.size();

    // Sort edges by weights
    if (sort_edges) {
        // convert edges and weights to joint structure representation
        vector<WeightedEdge> weighted_edges;
        for (int i = 0; i < g.n_edges; i++)
            weighted_edges.push_back({g.edge_list[i], g.edge_weights[i]});
        // sort
        sort(weighted_edges.begin(), weighted_edges.end(), greater<WeightedEdge>());
        // convert back edges and weights to separate lists
        for (int i = 0; i < g.n_edges; i++) {
            g.edge_list[i] = weighted_edges[i].edge;
            g.edge_weights[i] = weighted_edges[i].weight;
        }
    }

    // Init max solution
    vector<vector<bool>> edge_masks_max;
    vector<vector<int>> colors_max;
    int score_max = -1;
    int n_calls = 0;

    // Check if bipartite
    vector<int> colors(g.n_nodes, -1);
    if (g.is_bipartite(colors)) {
        // If bipartite construct solution
        edge_masks_max.push_back(vector<bool>(g.n_edges, true));
        colors_max.push_back(colors);
        score_max = 0;
        for (int i = 0; i < g.n_edges; i++)
            score_max += g.edge_weights[i];
    }
    else {
        // If not prepare to BB-DFS
        // init empty solution
        vector<bool> edge_mask(g.n_edges, false);
        // init colors
        vector<int> colors(g.n_nodes, -1);
        colors[0] = 0;
        // init empty neighbor lists
        vector<vector<int>> neighbors(g.n_nodes, vector<int>());
        // Run BB-DFS
        bbdfs(g, edge_masks_max, colors_max, score_max, n_calls,
              0, edge_mask, colors, neighbors);
    }

    // Print results
    cout << "Weight: " << score_max << endl;
    cout << "Number of solutions: " << edge_masks_max.size() << endl;
    cout << "Number of calls: " << n_calls << endl;
    cout << "Solutions: " << endl;
    bool is_first = true;
    for (int i = 0; i < edge_masks_max.size(); i++) {
        cout << i + 1 << ")" << endl;
        cout << "U={";
        for (int j = 0 ; j < g.n_nodes; j++)
            if (colors_max[i][j] == 0) {
                cout << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        cout << "}" << endl;
        is_first = true;
        cout << "W={";
        for (int j = 0 ; j < g.n_nodes; j++)
            if (colors_max[i][j] == 1) {
                cout << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        cout << "}" << endl;
        is_first = true;
        cout << "E={";
        for (int j = 0; j < g.n_edges; j++)
            if (edge_masks_max[i][j] == true) {
                cout << (is_first ? "{" : ", {") << g.edge_list[j][0] << ", "
                     << g.edge_list[j][1] << "}";
                 if (is_first)
                     is_first = false;
            }
        cout << "}" << endl;
    }

    return 0;
}
