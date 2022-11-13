//Group #17
//Dinic Algorithm
//Members : 
//Rohit Kinha - 2020CSB1118
//Arya Phalke - 2020CSB1107
//Abhay Shukralia - 2020CSB1061
//Aakash - 2020CSB1060
//Sushil - 2020CSB1132

#include <iostream>
#include <vector>
#include <queue>
#include <limits.h>
using namespace std;
// Edge class for describing a edge in capacitated graph used for dinic algorithm
class Edge
{
public:
    int to_vertex;    // the vertex where the edge ends
    int reverse_edge; // to store the index of reverse edge
    int flow;         // to get the current flow in the edge
    int Capacity;     // overall capacity of that edge

    Edge(int to, int r, int f, int c)
    {
        this->to_vertex = to;
        this->reverse_edge = r;
        this->flow = f;
        this->Capacity = c;
    }
};
// for creating a graph of n nodes for dinic
class Graph
{
public:
    int nodes;
    vector<vector<Edge>> adj;
    Graph(int n)
    {
        this->nodes = n;
        adj = vector<vector<Edge>>(n);
    }
};
// For the functions to find the maximum flow using dinic's algorithm
class dinic
{
public:
    // addding the normal and backedge to the graph
    void add_edges(int from_vertex, int to_vertex, int Capacity, Graph &g)
    {
        //Edge a is a forward edge storing the vertex to which it is pointing, flow, capacity of that edge and the index of reverse edge.
        //Edge b is back edge pointing to the from vertex, with 0 flow and 0 capacity and index of the reverse edge. 
        Edge a = Edge(to_vertex, (int)g.adj[to_vertex].size(), 0, Capacity);
        Edge b = Edge(from_vertex, (int)g.adj[from_vertex].size(), 0, 0);
        g.adj[from_vertex].push_back(a);
        g.adj[to_vertex].push_back(b);
    }
    //Function to convert a bipartite matching into network flow with a source and a sink all edges with unit capacity
    void construct_network(int m, int n, int source, int sink, Graph &g)
    {
        //Adding edges of left_set to source
        for (int i = 1; i <= m; i++)
        {
            add_edges(source, i, 1, g);
        }
        //Adding edges of right_set to sink
        for (int j = m + 1; j < sink; j++)
        {

            add_edges(j, sink, 1, g);
        }
    }
    // BFS function to get the level graph(if possible to send more flow) from source to sink vertex and labelling levels to all the vertex
    int BFS(int source, int sink, Graph &g, vector<int> &level)
    {
        vector<int> visited(g.nodes, 0);//to keep track of the visited vertex
        for (int i = 0; i < g.nodes; i++)
        {
            level[i] = -1; // level array to store the level of each visited vertex
        }
        queue<int> q; // for BFS traversal
        level[source] = 0; // level of source vertex is zero
        visited[source] = 1;
        q.push(source);
        while (!q.empty())
        {
            int from_vertex = q.front();
            q.pop();
            for (auto i = g.adj[from_vertex].begin(); i != g.adj[from_vertex].end(); i++)
            {
                if (visited[(*i).to_vertex] == 0 && (*i).flow < (*i).Capacity )
                {
                    level[(*i).to_vertex] = level[from_vertex] + 1;
                    q.push((*i).to_vertex);
                    visited[(*i).to_vertex] = 1;
                }
            }
        }
        // if we can't reach the sink vertex that means no more flow can be reached
        if (visited[sink] == 0)
        {
            return 0;
        }
        return 1;
    }
    // DFS function to get the bottleneck flow for a given augumenting path of level graph
    int DFS(int from_vertex, int flow, int sink, vector<int> &next, Graph &g, vector<int> &level)
    {
        if (from_vertex == sink)
        {
            return flow;
        }
        int n = g.adj[from_vertex].size();
        for (int &i = next[from_vertex]; i < n; i++)
        {
            Edge &e = g.adj[from_vertex][i];//checking if this edge can be part of the flow network or not, i.e ensuring it doesn't lead to a dead end
            if (e.flow < e.Capacity && level[e.to_vertex] == level[from_vertex] + 1)
            {
                int bottleneck_flow = DFS(e.to_vertex, min(flow, e.Capacity - e.flow), sink, next, g, level);
                //Augmenting the edges
                if (bottleneck_flow > 0)
                {
                    //adding the bottleneck_flow to forward edge  
                    e.flow += bottleneck_flow;
                    //Subracting the flow from reverse edge
                    g.adj[e.to_vertex][e.reverse_edge].flow -= bottleneck_flow;
                    return bottleneck_flow;
                }
            }
        }
        return 0;
    }
    int dinic_solver(int source, int sink, Graph &g)
    {
        if (source == sink)
            return -1;
        int total_flow = 0;
        vector<int> level(g.nodes);
        // next[i] indicates the next edge index to take in the adjacency list for node i. This is
        // part of the Shimon Even and Alon Itai optimization of pruning deads ends as part of the DFS phase.
        vector<int> next(g.nodes);
        while (BFS(source, sink, g, level) == 1)
        {
            for (int i = 0; i < g.nodes; i++)
            {
                next[i] = 0;
            }
            while (int flow = DFS(source, INT_MAX, sink, next, g, level))
            {
                total_flow += flow;
                for (int i = 0; i < g.nodes;i++)
                {
                    cout << next[i] << " ";
                }
                cout << endl;
            }
        }
        return total_flow;
    }
};
int main()
{
    dinic sol;
    Graph g(6);
    sol.add_edges(0, 1, 5, g);
    sol.add_edges(0, 2, 6, g);
    sol.add_edges(2, 1, 7, g);
    sol.add_edges(1, 3, 6, g);
    sol.add_edges(3, 2, 2, g);
    sol.add_edges(2, 4, 3, g);
    sol.add_edges(3, 4, 4, g);
    sol.add_edges(3, 5, 6, g);
    sol.add_edges(4, 5, 4, g);

    cout << "Maximum flow " << sol.dinic_solver(0, 5, g);

    // More Examples
    // Graph g(6);
    // sol.add_edges(0, 1, 16,g);
    // sol.add_edges(0, 2, 13,g);
    // sol.add_edges(1, 2, 10,g);
    // sol.add_edges(1, 3, 12,g);
    // sol.add_edges(2, 1, 4,g);
    // sol.add_edges(2, 4, 14,g);
    // sol.add_edges(3, 2, 9,g);
    // sol.add_edges(3, 5, 20,g);
    // sol.add_edges(4, 3, 7,g);
    // sol.add_edges(4, 5, 4,g);

    // cout << "Maximum flow " << sol.dinic_solver(0, 5, g);
    // int m, n;
    // m = n = 4;
    // int source = 0;
    // int sink = m + n + 1;
    // Graph g(10);
    // sol.add_edges(1, 6, 1, g);
    // sol.add_edges(1, 7, 1, g);
    // sol.add_edges(2, 5, 1, g);
    // sol.add_edges(3, 6, 1, g);
    // sol.add_edges(4, 8, 1, g);
    // sol.add_edges(4, 6, 1, g);
    // sol.construct_network(m, n, source, sink, g);
    // cout << "maximum answer of the bipartite graph is: " << sol.dinic_solver(0, 9, g);
    return 0;
}