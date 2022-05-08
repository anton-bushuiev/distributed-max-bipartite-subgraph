#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <fstream>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>

using namespace std;


// OpenMPI constants
const int TAG_WORK = 1;
const int TAG_DONE = 2;
const int TAG_TERMINATE = 3;


// Struct for input graph
struct Graph
{
    int n_nodes;
    int n_edges;
    vector<vector<int>> edge_list;
    vector<int> edge_weights;

    bool is_bipartite(vector<int> & colors) const
    {
        // Returns coloring with 2 colors if graph is bipartite else null
        // Note: assumes that graph is connected

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

    bool is_subgraph_connencted(const vector<bool> & edge_mask) const
    {
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

    vector<int> to_msg() const
    {
        vector<int> msg;

        // n_nodes, n_edges
        msg.push_back(n_nodes);
        msg.push_back(n_edges);

        // edge_list
        for (int i = 0; i < edge_list.size(); i++) {
            msg.insert(msg.end(), edge_list[i].begin(), edge_list[i].end());
        }

        // edge_weights
        msg.insert(msg.end(), edge_weights.begin(), edge_weights.end());

        return msg;
    }

    static Graph from_msg(vector<int> msg)
    {
        Graph graph;

        // n_nodes, n_edges
        graph.n_nodes = msg[0];
        graph.n_edges = msg[1];

        // edge_list
        vector<int>::const_iterator beg;
        vector<int>::const_iterator end = msg.begin() + 2;
        graph.edge_list = vector<vector<int>>(graph.n_edges, vector<int>());
        for (int i = 0; i < graph.edge_list.size(); i++) {
            beg = end;
            end = beg + 2;
            graph.edge_list[i] = vector<int>(beg, end);
        }

        // edge_weights
        beg = end;
        end = beg + graph.n_edges;
        graph.edge_weights = vector<int>(beg, end);

        return graph;
    }
};

Graph G;


// Struct for maximum solution
struct Solution
{
    int score = -1;
    vector<vector<bool>> edge_masks;
    vector<vector<int>> colors;

    bool update(const int score_other, const vector<bool> & edge_mask_other,
                const vector<int> & colors_other)
    {
        if (score_other == score) {
            edge_masks.push_back(edge_mask_other);
            colors.push_back(colors_other);
            return true;
        }
        else if (score_other > score) {
            score = score_other;
            edge_masks = {edge_mask_other};
            colors = {colors_other};
            return true;
        }
        return false;
    }

    bool update(const Solution other)
    {
        if (other.score == score) {
            edge_masks.insert(edge_masks.end(),
                              other.edge_masks.begin(), other.edge_masks.end());
            colors.insert(colors.end(),
                          other.colors.begin(), other.colors.end());
            return true;
        }
        else if (other.score > score) {
            score = other.score;
            edge_masks = other.edge_masks;
            colors = other.colors;
            return true;
        }
        return false;
    }

    vector<int> to_msg() const
    {
        vector<int> msg;

        // score
        msg.push_back(score);

        // number of solutions
        msg.push_back(edge_masks.size());

        // edge_masks
        vector<int> edge_mask_int;
        for (const vector<bool> & edge_mask: edge_masks) {
            edge_mask_int = vector<int>(edge_mask.begin(), edge_mask.end());
            msg.insert(msg.end(), edge_mask_int.begin(), edge_mask_int.end());
        }

        // colors
        for (const vector<int> & colors_one: colors) {
            msg.insert(msg.end(), colors_one.begin(), colors_one.end());
        }

        return msg;
    }

    static Solution from_msg(const vector<int> & msg)
    {
        Solution solution;

        // score
        solution.score = msg[0];

        // number of solutions
        int n_solutions = msg[1];

        // edge_masks
        vector<int>::const_iterator beg;
        vector<int>::const_iterator end = msg.begin() + 2;
        for (int i = 0; i < n_solutions; i++) {
            beg = end;
            end = beg + G.n_edges;
            solution.edge_masks.push_back(vector<bool>(beg, end));
        }

        // colors
        for (int i = 0; i < n_solutions; i++) {
            beg = end;
            end = beg + G.n_nodes;
            solution.colors.push_back(vector<int>(beg, end));
        }

        return solution;
    }
};

ostream & operator<< (ostream & os, const Solution & solution)
{
    os << "Maximum weight: " << solution.score << endl;
    os << "Number of solutions: " << solution.edge_masks.size() << endl;
    os << "Solutions: " << endl;
    bool is_first;
    for (int i = 0; i < solution.edge_masks.size(); i++) {
        is_first = true;
        os << i + 1 << ")" << endl;
        os << "U={";
        for (int j = 0 ; j < solution.colors[i].size(); j++)
            if (solution.colors[i][j] == 0) {
                os << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        os << "}" << endl;
        is_first = true;
        os << "W={";
        for (int j = 0 ; j < solution.colors[i].size(); j++)
            if (solution.colors[i][j] == 1) {
                os << (is_first ? "" : ", ") << j;
                if (is_first)
                    is_first = false;
            }
        os << "}" << endl;
        is_first = true;
        os << "E={";
        for (int j = 0; j < solution.edge_masks[i].size(); j++)
            if (solution.edge_masks[i][j] == true) {
                os << (is_first ? "{" : ", {") << G.edge_list[j][0] << ", "
                     << G.edge_list[j][1] << "}";
                 if (is_first)
                     is_first = false;
            }
        os << "}";
        os << endl;
    }

    return os;
}

Solution SOLUTION;


// Struct for state in solution space
struct State
{
    int edge;
    vector<bool> edge_mask;
    vector<int> colors;
    vector<vector<int>> neighbors;

    vector<int> to_msg() const {
        vector<int> msg;

        // edge
        msg.push_back(edge);

        // edge_mask
        vector<int> edge_mask_int(edge_mask.begin(), edge_mask.end());
        msg.insert(msg.end(), edge_mask_int.begin(), edge_mask_int.end());

        // colors
        msg.insert(msg.end(), colors.begin(), colors.end());

        // neighbors
        for (int i = 0; i < neighbors.size(); i++) {
            msg.push_back(neighbors[i].size());
            msg.insert(msg.end(), neighbors[i].begin(), neighbors[i].end());
        }

        return msg;
    }

    static State from_msg(vector<int> msg)
    {
        State state;

        // edge
        state.edge = msg[0];

        // edge_mask
        vector<int>::const_iterator beg = msg.begin() + 1;
        vector<int>::const_iterator end = beg + G.n_edges;
        state.edge_mask = vector<bool>(beg, end);

        // colors
        beg = end;
        end = beg + G.n_nodes;
        state.colors = vector<int>(beg, end);

        // neighbors
        int n_neighbors;
        state.neighbors = vector<vector<int>>(G.n_nodes, vector<int>());
        for (int i = 0; i < state.neighbors.size(); i++) {
            beg = end;
            n_neighbors = *(beg++);
            end = beg + n_neighbors;
            state.neighbors[i] = vector<int>(beg, end);
        }

        return state;
    }
};

ostream & operator<< (ostream & os, const State & state)
{
    os << "edge: " << state.edge << endl;
    os << "edge_mask:" << endl;
    for (int i = 0; i < state.edge_mask.size(); i++)
        os << state.edge_mask[i] << " ";
    os << endl;
    os << "colors:" << endl;
    for (int i = 0; i < state.colors.size(); i++)
        os << state.colors[i] << " ";
    os << endl;
    os << "neighbors:" << endl;
    for (int i = 0; i < state.neighbors.size(); i++) {
        for (int j = 0; j < state.neighbors[i].size(); j++)
            os << state.neighbors[i][j] << " ";
        os << endl;
    }
    return os;
}


// Struct for sorting edges by weights
struct WeightedEdge
{
    vector<int> edge;
    int weight;
};

bool operator> (const WeightedEdge & a, const WeightedEdge & b)
{
    return (a.weight > b.weight);
}


// BB-DFS for searching through state-space
void bbdfs(int edge, vector<bool> edge_mask, vector<int> colors,
           vector<vector<int>> neighbors)
{
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
        if (score >= SOLUTION.score) {
            #pragma omp critical
            SOLUTION.update(score, edge_mask, colors);
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
    if (score + score_future < SOLUTION.score)
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


// BFS to construct initial states for BB-DFS
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
    // note: q.size() can increase from 0 to 2 per iteration
    while (q.size() <= max_n_states - 1) {
        // Pop
        if (q.size() == 0)
            return {};
        state = q.front();
        q.pop();
        edge = state.edge;
        edge_mask = state.edge_mask;
        colors = state.colors;
        neighbors = state.neighbors;

        // Stop if state is final (it makes no sense to continue if
        // BFS has alreaday reached the last layer)
        if (edge >= G.n_edges) {
            q.push(state);
            break;
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
        if (score + score_future < SOLUTION.score)
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

    // Convert queue to vector
    vector<State> states;
    while (!q.empty())
    {
        states.push_back(q.front());
        q.pop();
    }
    return states;
}


// TODO Remove
ostream & operator<<(ostream & os, const vector<int> & v)
{
    for (int i = 0; i < v.size(); i++)
        os << v[i] << " ";
    return os;
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank, n_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    int n_threads = omp_get_max_threads();

    // Print number of processes and threads
    if (rank == 0) {
	cout << "Number of processes: " << n_procs << endl;
    }
    cout << "Rank: " << rank << ",number of threads: " << n_threads << endl;

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

    // If bipartite, construct solution immediately
    if (rank == 0) {
        vector<int> colors(G.n_nodes, -1);
        if (G.is_bipartite(colors)) {
            // Construct solution
            SOLUTION.edge_masks.push_back(vector<bool>(G.n_edges, true));
            SOLUTION.colors.push_back(colors);
            SOLUTION.score = 0;
            for (int i = 0; i < G.n_edges; i++)
                SOLUTION.score += G.edge_weights[i];

            // Send terminate signal to all other processes
            for (int dest = 1; dest < n_procs; dest++) {
                MPI_Send(nullptr, 0, MPI_INT,
                         dest, TAG_TERMINATE, MPI_COMM_WORLD);
            }

            // Print solution and exit
            cout << SOLUTION;
            MPI_Finalize();
            return 0;
        }
    }

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


    // Run distributed computation
    int buffer_size = 10000;
    vector<int> msg(buffer_size);
    vector<int> msg_recv(buffer_size);

    int n_jobs_per_process = 10;
    int n_jobs_per_thread = 10;

    int size_received;
    MPI_Status status;
    State state;
    int max_n_states;
    vector<State> states;
    Solution received_solution;

    if (rank == 0) {
        // Construct jobs (states) for processes with BFS
        // init root state
        int edge = 0;
        vector<bool> edge_mask(G.n_edges, false);
        vector<int> colors(G.n_nodes, -1);
        colors[0] = 0;
        vector<vector<int>> neighbors(G.n_nodes, vector<int>());
        State init_state{edge, edge_mask, colors, neighbors};
        // run BFS
        max_n_states = (n_procs - 1) * n_jobs_per_process;
        states = bfs(init_state, max_n_states);
        int i_state = 0;

        // Send jobs to all other processes
        for (int dest = 1; dest < n_procs; dest++) {
            msg = states[i_state++].to_msg();
            MPI_Send(msg.data(), msg.size(), MPI_INT,
                     dest, TAG_WORK, MPI_COMM_WORLD);
        }

        int n_active_workers = n_procs - 1;
        while (n_active_workers > 0) {
            // Recieve done message from other process
            msg_recv.resize(buffer_size);
            MPI_Recv(msg_recv.data(), buffer_size, MPI_INT,
                     MPI_ANY_SOURCE, TAG_DONE, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &size_received);
            msg_recv.resize(size_received);

            // If more jobs left, send one
            if (i_state != states.size()) {
                msg = states[i_state++].to_msg();
                MPI_Send(msg.data(), msg.size(), MPI_INT,
                         status.MPI_SOURCE, TAG_WORK, MPI_COMM_WORLD);
            }
            // If no more jobs left, send terminate signal
            else {
                MPI_Send(nullptr, 0, MPI_INT,
                         status.MPI_SOURCE, TAG_TERMINATE, MPI_COMM_WORLD);
                n_active_workers--;
            }

            // Update solution
            received_solution = Solution::from_msg(msg_recv);
            SOLUTION.update(received_solution);
        }

        // Print result
        cout << SOLUTION;
    }
    else {
        while (true) {
            // Recieve message from master process
            msg.resize(buffer_size);
            MPI_Recv(msg.data(), buffer_size, MPI_INT,
                     0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            // If terminate signal received, stop working
            if (status.MPI_TAG == TAG_TERMINATE) {
                break;
            }
            // If work signal received, keep working
            else if (status.MPI_TAG == TAG_WORK) {
                MPI_Get_count(&status, MPI_INT, &size_received);
                msg.resize(size_received);
                state = State::from_msg(msg);

                // Construct jobs (states) for threads with BFS
                max_n_states = n_threads * n_jobs_per_thread;
                states = bfs(state, max_n_states);

                // Run parallelized BB-DFS
                #pragma omp parallel for schedule(dynamic)
                for (int i = 0; i < states.size(); i++)
                {
                    bbdfs(states[i].edge, states[i].edge_mask,
                          states[i].colors, states[i].neighbors);
                }

                // Send solution
                msg = SOLUTION.to_msg();
                MPI_Send(msg.data(), msg.size(), MPI_INT,
                         0, TAG_DONE, MPI_COMM_WORLD);
                SOLUTION.edge_masks = {};
                SOLUTION.colors = {};
            }
        }
    }

    MPI_Finalize();
    return 0;
}
