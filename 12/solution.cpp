/*
 * 12. Мосты.
 * Ребро неориентированного графа называется мостом, если удаление этого ребра из графа увеличивает число компонент связности.
 * Дан неориентированный граф, требуется найти в нем все мосты.
 */

#define DEBUG

#include <algorithm>
#include <assert.h>
#include <fstream>
#ifdef DEBUG
#include <iostream>
#endif
#include <limits.h>
#include <unordered_map>
#include <stack>
#include <vector>

/**
 * Edge
 */
struct Edge
{
    unsigned edgeIndex_; // edge index
    bool repeating_; // edge is repeating

    Edge(): edgeIndex_(std::numeric_limits<unsigned>::max()), repeating_(false)
    {
    }
    Edge(unsigned edgeIndex): edgeIndex_(edgeIndex), repeating_(false)
    {
    }
    Edge(unsigned edgeIndex, bool repeating): edgeIndex_(edgeIndex), repeating_(repeating)
    {
    }
};

/**
 * Vertex state
 */
struct VertexState {
    unsigned parent_;  // dfs parent
    unsigned entry_;   // order of entry
    unsigned low_;     // order of entry propagated from already visited vertici
    std::unordered_map<unsigned, Edge>::const_iterator children_; // children iterator
    bool hasChildren_; // flag of children ownship

    VertexState(unsigned emptyIndex, std::unordered_map<unsigned, Edge>::const_iterator children): parent_(emptyIndex), children_(children), hasChildren_(true)
    {
    }
    VertexState(unsigned emptyIndex): parent_(emptyIndex), hasChildren_(false)
    {
    }
};

/**
 * Find all bridges in a fully connected graph with DFS
 */
void findBridges(const std::unordered_map<unsigned, Edge>** adjastencyList, const unsigned vertexCount,
                 std::vector<VertexState>& state, unsigned startVertex, std::vector<unsigned>& result)
{
    state[startVertex].parent_ = vertexCount; // visited
    assert(state[startVertex].hasChildren_ == (adjastencyList[startVertex] != nullptr));
    if (!state[startVertex].hasChildren_)
    {
        return;
    }
    state[startVertex].entry_ = 0; // is root
    state[startVertex].low_ = 0;

    unsigned order = 1;
    std::stack<unsigned> stack; // dfs stack
    stack.push(startVertex);

    while (!stack.empty())
    {
        auto& vertexState = state[stack.top()];
        assert(vertexState.hasChildren_);
        if (vertexState.children_ == adjastencyList[stack.top()]->cend()) // every child visited
        {
            unsigned vertex = stack.top();
            stack.pop();
            if (!stack.empty())
            {
                auto& parentState = state[stack.top()];
                parentState.low_ = std::min(parentState.low_, vertexState.low_);
                if (vertexState.low_ > parentState.entry_) // bridge condition
                {
                    auto edge = adjastencyList[stack.top()]->find(vertex);
                    if (!(edge->second.repeating_)) // edge is not a repeating one
                    {
                        result.push_back(edge->second.edgeIndex_);
                    }
                }
            }
        }
        else // visit next child
        {
            unsigned nextVertex = vertexState.children_->first;
            auto& nextVertexState = state[nextVertex];
            if (nextVertexState.parent_ == stack.top()) // ignore a loop
            {
                ++vertexState.children_;
                continue;
            }
            if (nextVertexState.parent_ == vertexCount + 1) // an unvisited child
            {
                nextVertexState.parent_ = stack.top();
                nextVertexState.entry_ = nextVertexState.low_ = ++order;
                stack.push(nextVertex);
            }
            else if (vertexState.parent_ != nextVertex) // already visited, but not a direct parent
            {
                vertexState.low_ = std::min(nextVertexState.entry_, vertexState.low_);
            }
            ++vertexState.children_;
        }
    }
}

/**
 * Find all bridges in the graph
 */
std::vector<unsigned> findBridges(const std::unordered_map<unsigned, Edge>** adjastencyList, const unsigned vertexCount)
{
    std::vector<unsigned> result;
    std::vector<VertexState> state; // vertexCount + 1 == unvisited

    for (unsigned i = 0; i < vertexCount; ++i)
    {
        if (adjastencyList[i])
        {
            state.push_back(VertexState(vertexCount + 1, adjastencyList[i]->cbegin()));
        }
        else
        {
            state.push_back(VertexState(vertexCount + 1));
        }
    }
    for (unsigned i = 0; i < vertexCount; ++i)
    {
        if (state[i].parent_ == vertexCount + 1)
        {
            findBridges(adjastencyList, vertexCount, state, i, result);
        }
    }
    std::sort(result.begin(), result.end());
    return result;
}

int main(void) {
#ifdef DEBUG
    std::ifstream cin("input_451.txt");
    using std::cout;
#else
    std::ifstream cin("bridges.in");
    std::ofstream cout("bridges.out");
#endif

    unsigned vertexCount, edgeCount;
    cin >> vertexCount >> edgeCount;

    std::unordered_map<unsigned, Edge>** adjastencyList = new std::unordered_map<unsigned, Edge>*[vertexCount];

    for (unsigned i = 0; i < vertexCount; ++i)
    {
        adjastencyList[i] = nullptr;
    }

    for (unsigned i = 0; i < edgeCount; ++i)
    {
        unsigned p0, p1;
        cin >> p0 >> p1;
        --p0;
        --p1;
        assert(p0 < vertexCount && p1 < vertexCount);

        if (!adjastencyList[p0])
        {
            adjastencyList[p0] = new std::unordered_map<unsigned, Edge>();
            adjastencyList[p0]->insert(std::make_pair(p1, Edge(i)));
        }
        else
        {
            auto edge = adjastencyList[p0]->find(p1);
            if (edge == adjastencyList[p0]->end())
            {
                adjastencyList[p0]->insert(std::make_pair(p1, Edge(i)));
            }
            else
            {
                edge->second = Edge(i, true);
            }
        }

        if (!adjastencyList[p1])
        {
            adjastencyList[p1] = new std::unordered_map<unsigned, Edge>();
            adjastencyList[p1]->insert(std::make_pair(p0, Edge(i)));
        }
        else
        {
            auto edge = adjastencyList[p1]->find(p0);
            if (edge == adjastencyList[p1]->end())
            {
                adjastencyList[p1]->insert(std::make_pair(p0, Edge(i)));
            }
            else
            {
                edge->second = Edge(i, true);
            }
        }
    }

    auto result = findBridges(const_cast<const std::unordered_map<unsigned, Edge>**>(adjastencyList), vertexCount);

    for (unsigned i = 0; i < vertexCount; ++i)
    {
        delete adjastencyList[i];
    }
    delete[] adjastencyList;

    cout << result.size() << '\n';
    for (auto edge: result)
    {
        cout << (edge + 1) << ' ';
    }

    return 0;
}
