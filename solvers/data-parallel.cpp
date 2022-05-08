#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <fstream>
#include <omp.h>
#include <cstdlib>

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


// Input graph
Graph G;

// Max solution
vector<vector<bool>> EDGE_MASKS_MAX;
vector<vector<int>> COLORS_MAX;
int SCORE_MAX = -1;

// Number of recursive calls
int N_CALLS = 0;

// Threshold of number of edges left unprocessed to start a task
int TASK_THRESHOLD = 10;


void bbdfs(int edge, vector<bool> edge_mask, vector<int> colors,
           vector<vector<int>> neighbors)
{
    N_CALLS += 1;

    // Calculate solution score
    int score = 0;
    for (int i = 0; i < G.n_edges; i++) {
        if (edge_mask[i] == true) {
            score += G.edge_weights[i];
        }
    }

    // ANALYZE TERMINAL STATE
    if (edge >= G.n_edges) {

        // Check connectivity
        if (!G.is_subgraph_connencted(edge_mask))
            return;

        // Update optimal solution
        if (score >= SCORE_MAX) {
            #pragma omp critical
            {
                if (score == SCORE_MAX) {
                    EDGE_MASKS_MAX.push_back(edge_mask);
                    COLORS_MAX.push_back(colors);
                }
                else if (score > SCORE_MAX) {
                    SCORE_MAX = score;
                    EDGE_MASKS_MAX = {edge_mask};
                    COLORS_MAX = {colors};
                }
            }
        }
        return;
    }

    // PRUNE
    // Calculate upper bound on child solutions
    int score_future = 0;
    for (int i = edge; i < G.n_edges; i++) {
            score_future += G.edge_weights[i];
    }

    // Prune if better child solution is impossible
    if (score + score_future < SCORE_MAX)
        return;

    // BRANCH
    vector<int> colors_prev = colors;
    vector<vector<int>> neighbors_prev = neighbors;

    // "Add edge" case
    edge_mask[edge] = true;
    neighbors[G.edge_list[edge][1]].push_back(G.edge_list[edge][0]);
    neighbors[G.edge_list[edge][0]].push_back(G.edge_list[edge][1]);

    // Try to color (u, v) and (v, u) edge directions with (0, 1) colors
    int u, v;
    for (int d : {0, 1}) {
        u = G.edge_list[edge][d];
        v = G.edge_list[edge][d == 0 ? 1 : 0];

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
            bbdfs(edge + 1, edge_mask, colors, neighbors);
        }
    }

    // "Ignore edge" case
    edge_mask[edge] = false;
    bbdfs(edge + 1, edge_mask, colors_prev, neighbors_prev);
}


struct State {
    int edge;
    vector<bool> edge_mask;
    vector<int> colors;
    vector<vector<int>> neighbors;
};

ostream & operator<< (ostream & os, const State & state) {
    os << "edge: " << state.edge << endl;
    os << "edge_mask:";
    for (int i = 0; i < state.edge_mask.size(); i++)
        os << state.edge_mask[i] << " ";
    os << endl;
    os << "colors:";
    for (int i = 0; i < state.colors.size(); i++)
        os << state.colors[i] << " ";
    os << endl;
    os << "neighbors:";
    for (int i = 0; i < state.neighbors.size(); i++) {
        for (int j = 0; j < state.neighbors[i].size(); j++)
            os << state.neighbors[i][j] << " ";
        os << endl;
    }
    return os;
}


vector<State> bfs(State init_state, int max_n_states)
{
    // Push first state
    queue<State> q;
    q.push(init_state);

    // Generate states with BFS
    State state;
    int edge;
    vector<bool> edge_mask;
    vector<int> colors;
    vector<vector<int>> neighbors;
    while (q.size() <= max_n_states) {
        // Pop
        state = q.front();
        q.pop();
        edge = state.edge;
        edge_mask = state.edge_mask;
        colors = state.colors;
        neighbors = state.neighbors;

        // Check if final
        if (edge >= G.n_edges) {
            continue;
        }

        // Calculate solution score
        int score = 0;
        for (int i = 0; i < G.n_edges; i++) {
            if (edge_mask[i] == true) {
                score += G.edge_weights[i];
            }
        }

        // PRUNE
        // Calculate upper bound on child solutions
        int score_future = 0;
        for (int i = edge; i < G.n_edges; i++) {
                score_future += G.edge_weights[i];
        }

        // Prune if better child solution is impossible
        if (score + score_future < SCORE_MAX)
            continue;

        // BRANCH
        vector<int> colors_prev = colors;
        vector<vector<int>> neighbors_prev = neighbors;

        // "Add edge" case
        edge_mask[edge] = true;
        neighbors[G.edge_list[edge][1]].push_back(G.edge_list[edge][0]);
        neighbors[G.edge_list[edge][0]].push_back(G.edge_list[edge][1]);

        // Try to color (u, v) and (v, u) edge directions with (0, 1) colors
        int u, v;
        for (int d : {0, 1}) {
            u = G.edge_list[edge][d];
            v = G.edge_list[edge][d == 0 ? 1 : 0];

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
                q.push({edge + 1, edge_mask, colors, neighbors});
            }
        }

        // "Ignore edge" case
        edge_mask[edge] = false;
        q.push({edge + 1, edge_mask, colors_prev, neighbors_prev});
    }

    vector<State> states;
    while (!q.empty())
    {
        states.push_back(q.front());
        q.pop();
    }
    return states;
}


// Note: for sorting only
struct WeightedEdge {
    vector<int> edge;
    int weight;
};


bool operator> (const WeightedEdge & a, const WeightedEdge & b) {
    return (a.weight > b.weight);
}


int main(int argc, char* argv[]) {
    int n_threads = omp_get_max_threads();

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
    infile >> G.n_nodes;
    G.edge_list = vector<vector<int>>();
    G.edge_weights = vector<int>();

    int a;
    for (int i = 0; i < G.n_nodes; i++) {
        for (int j = 0; j < G.n_nodes; j++) {
            infile >> a;
            if (a != 0 && j < i) {
                G.edge_list.push_back({i, j});
                G.edge_weights.push_back(a);
            }
        }
    }
    G.n_edges = (unsigned) G.edge_list.size();

    // Sort edges by weights
    if (sort_edges) {
        // convert edges and weights to joint structure representation
        vector<WeightedEdge> weighted_edges;
        for (int i = 0; i < G.n_edges; i++)
            weighted_edges.push_back({G.edge_list[i], G.edge_weights[i]});
        // sort
        sort(weighted_edges.begin(), weighted_edges.end(), greater<WeightedEdge>());
        // convert back edges and weights to separate lists
        for (int i = 0; i < G.n_edges; i++) {
            G.edge_list[i] = weighted_edges[i].edge;
            G.edge_weights[i] = weighted_edges[i].weight;
        }
    }

    // Check if bipartite
    vector<int> colors(G.n_nodes, -1);
    if (G.is_bipartite(colors)) {
        // If bipartite construct solution
        EDGE_MASKS_MAX.push_back(vector<bool>(G.n_edges, true));
        COLORS_MAX.push_back(colors);
        SCORE_MAX = 0;
        for (int i = 0; i < G.n_edges; i++)
            SCORE_MAX += G.edge_weights[i];
    }
    else {
        // If not, prepare to BB-DFS
        // init empty solution
        int edge = 0;
        vector<bool> edge_mask(G.n_edges, false);
        // init colors
        vector<int> colors(G.n_nodes, -1);
        colors[0] = 0;
        // init empty neighbor lists
        vector<vector<int>> neighbors(G.n_nodes, vector<int>());
        State init_state{edge, edge_mask, colors, neighbors};

        // Generate states with BFS
        int max_n_states = n_threads * 30 * 10;
        vector<State> states = bfs(init_state, max_n_states);

        // Print number of threads and cores
        #pragma omp parallel
        {
            #pragma omp single
            {
                cout << "Running " << omp_get_num_threads() << " threads on "
                     << omp_get_num_procs() << " cores." << endl;
            }
        }

        // Run BB-DFS
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < states.size(); i++)
        {
            bbdfs(states[i].edge, states[i].edge_mask,
                  states[i].colors, states[i].neighbors);
        }
    }

    // Print results
    cout << "Result: " << endl;
    cout << "Maximum weight: " << SCORE_MAX << endl;
    cout << "Number of solutions: " << EDGE_MASKS_MAX.size() << endl;
    cout << "Number of recursive calls: " << N_CALLS << endl;
    cout << "Solutions: " << endl;
    bool is_first = true;
    for (int i = 0; i < EDGE_MASKS_MAX.size(); i++) {
        cout << i + 1 << ")" << endl;
        cout << "U={";
        for (int j = 0 ; j < G.n_nodes; j++)
            if (COLORS_MAX[i][j] == 0) {
                cout << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        cout << "}" << endl;
        is_first = true;
        cout << "W={";
        for (int j = 0 ; j < G.n_nodes; j++)
            if (COLORS_MAX[i][j] == 1) {
                cout << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        cout << "}" << endl;
        is_first = true;
        cout << "E={";
        for (int j = 0; j < G.n_edges; j++)
            if (EDGE_MASKS_MAX[i][j] == true) {
                cout << (is_first ? "{" : ", {") << G.edge_list[j][0] << ", "
                     << G.edge_list[j][1] << "}";
                 if (is_first)
                     is_first = false;
            }
        cout << "}" << endl;
    }

    return 0;
}
