/*
 * 11_1. Цикл минимальной длины
 * Дан невзвешенный неориентированный граф. Найдите цикл минимальной длины.
 * Ввод: v: кол-во вершин (макс. 50000), n: кол-во ребер (макс. 200000), n пар реберных вершин.
 * Вывод: одно целое число равное длине минимального цикла. Если цикла нет, то вывести -1.
 */

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <queue>
#include <vector>

/**
 * Find the smallest cycle from startVertex to startVertex by running Dijkstra algorithm
 */
std::size_t findMinCycleLength(std::vector<std::size_t>** adjastencyList, const std::size_t vertexCount, std::size_t startVertex)
{
    // <parent, distance>
    std::vector<std::pair<std::size_t, std::size_t>> distances(vertexCount, std::make_pair(vertexCount, vertexCount + 1));
    std::queue<std::size_t> vertices;

    distances[startVertex].first = startVertex;
    distances[startVertex].second = 0;
    vertices.push(startVertex);

    std::size_t result = vertexCount + 1;
    if (!adjastencyList[startVertex])
    {
        return result;
    }
    while (!vertices.empty())
    {
        std::size_t v = vertices.front();
        const auto& vDistance = distances[v];
        assert(vDistance.second <= vertexCount);
        vertices.pop();
        const std::vector<std::size_t>* adjastent = adjastencyList[v];
        assert(adjastent != nullptr);

        for (auto it = adjastent->begin(); it < adjastent->end(); ++it)
        {
            if (distances[*it].second > vertexCount)
            {
                distances[*it].first = v;
                distances[*it].second = vDistance.second + 1;
                vertices.push(*it);
            }
            else if (distances[*it].first != v && vDistance.first != *it && distances[*it].second + vDistance.second + 1 < result)
            {
                result = distances[*it].second + vDistance.second + 1;
            }
        }
    }
    return result;
}

/**
 * Find the smallest cycle in a graph
 */
std::size_t findMinCycleLength(std::vector<std::size_t>** adjastencyList, const std::size_t vertexCount)
{
    std::size_t result = vertexCount + 1;
    for (std::size_t i = 0; i < vertexCount; ++i)
    {
        std::size_t vertexResult = findMinCycleLength(adjastencyList, vertexCount, i);
        if (vertexResult < result)
        {
            result = vertexResult;
        }
    }
    return result;
}

int main(void) {
    using std::cin;

    std::size_t vertexCount, edgeCount;
    cin >> vertexCount >> edgeCount;

    std::vector<std::size_t>** adjastencyList = new std::vector<std::size_t>*[vertexCount];
    std::fill(adjastencyList, adjastencyList + vertexCount, nullptr);

    for (std::size_t i = 0; i < edgeCount; i++)
    {
        std::size_t p0, p1;
        cin >> p0 >> p1;
        if (!adjastencyList[p0])
        {
            adjastencyList[p0] = new std::vector<std::size_t>();
        }
        adjastencyList[p0]->push_back(p1);
        if (!adjastencyList[p1])
        {
            adjastencyList[p1] = new std::vector<std::size_t>();
        }
        adjastencyList[p1]->push_back(p0);
    }

    std::size_t minCycleLength = findMinCycleLength(adjastencyList, vertexCount);

    for (std::size_t i = 0; i < vertexCount; i++)
    {
        delete adjastencyList[i];
    }
    delete[] adjastencyList;

    if (minCycleLength <= vertexCount)
    {
        std::cout << minCycleLength << std::endl;
    }
    else
    {
        std::cout << "-1" << std::endl;
    }

    return 0;
}
