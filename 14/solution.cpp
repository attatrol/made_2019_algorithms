/*
 * 14_2. «MST».
 *
 * Дан неориентированный связный граф.
 * Требуется найти вес минимального остовного дерева в этом графе.
 * С помощью алгоритма Крускала.
 */

#define DEBUG

#include <assert.h>
#ifdef DEBUG
#include <fstream>
#endif
#include <iostream>
#include <queue>

/* Disjoint set of vertices */
struct DisjointSet
{
    //std::size_t count_;
    std::size_t* parent_; // value parent
    std::size_t* rank_;   // value rank

    DisjointSet(std::size_t count):
        /*count_(count), */parent_(new std::size_t[count]), rank_(new std::size_t[count])
    {
        for (std::size_t i = 0; i < count; ++i)
        {
            parent_[i] = i;
            rank_[i] = 1;
        }
    }
    DisjointSet(const DisjointSet&) = delete;
    DisjointSet(DisjointSet&&) = delete;
    DisjointSet& operator=(const DisjointSet&) = delete;
    DisjointSet& operator=(DisjointSet&&) = delete;
    ~DisjointSet()
    {
        delete[] parent_;
        delete[] rank_;
    }
    /* Find characteristic member of a set by its member */
    std::size_t find(std::size_t v)
    {
        return parent_[v] == v ? v : parent_[v] = find(parent_[v]); // recursion depth < 5
    }
    /* Merge 2 sets by its members, return false if they are already belong to one set */
    bool merge(std::size_t v1, std::size_t v2)
    {
        std::size_t parent1 = find(v1);
        std::size_t parent2 = find(v2);

        if (parent1 == parent2)
        {
            return false;
        }

        std::size_t rank1 = rank_[parent1];
        std::size_t rank2 = rank_[parent2];

        if (rank1 < rank2)
        {
            parent_[parent1] = parent2;
        }
        else if (rank2 < rank1)
        {
            parent_[parent2] = parent1;
        }
        else
        {
            parent_[parent1] = parent2;
            rank_[parent2]++;
        }
        return true;
    }
};
/* Weighted edge of an unordered graph */
struct Edge
{
    std::size_t vertex1; // a vertex
    std::size_t vertex2; // a vertex
    std::size_t weight;  // weight

    /* Comparator for the MST queue, lowest weight on top */
    bool operator<(const Edge& other) const &
    {
        return weight > other.weight;
    }
};

/* Calculate total MST edge weight by Kruskal algorithm */
std::size_t findMST(std::priority_queue<Edge>& edges, std::size_t vertexCount)
{
    DisjointSet ds(vertexCount);

    std::size_t result = 0;
    while (!edges.empty())
    {
        Edge topEdge = edges.top();
        if (ds.merge(topEdge.vertex1, topEdge.vertex2))
        {
            result += topEdge.weight;
        }
        edges.pop();
    }

    return result;
}

int main(void) {
#ifdef DEBUG
    std::ifstream cin("input.txt");
    using std::cout;
#else
    using std::cin;
    using std::cout;
#endif

    unsigned vertexCount, edgeCount;
    cin >> vertexCount >> edgeCount;

    std::priority_queue<Edge> queue;

    for (unsigned i = 0; i < edgeCount; ++i)
    {
        std::size_t v1, v2, weight;
        cin >> v1 >> v2 >> weight;
        --v1;
        --v2;
        assert(v1 < vertexCount && v2 < vertexCount);
        queue.push({v1, v2, weight});
    }

    std::size_t result = findMST(queue, vertexCount);
    cout << result;

    return 0;
}
