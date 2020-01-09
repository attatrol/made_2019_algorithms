/*
 * 13. «Пятнашки»
 * Написать алгоритм для решения игры в “пятнашки”. Решением задачи является приведение к виду:
 * [ 1  2  3  4 ]
 * [ 5  6  7  8 ]
 * [ 9  10 11 12]
 * [ 13 14 15 0 ]
 * где 0 задает пустую ячейку.
 * Достаточно найти хотя бы какое-то решение. Число перемещений костяшек не обязано быть минимальным.
 */

//#define NDEBUG // disable asserts

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdint.h>
#include <unordered_map>
#include <queue>
#include <vector>

using packed_t = unsigned int; // type for storing really small integers

/* Move direction */
enum class Move: packed_t
{
    UP = 0,
    RIGHT = 1,
    DOWN = 2,
    LEFT = 3
};

/* Calculate target position need to reach */
constexpr uint64_t targetPosition(uint64_t posValue)
{
    return posValue == 1 ? 1 : (posValue << 4 * (posValue - 1)) | targetPosition(posValue - 1);
}

static const uint64_t POS_MASK = 15; // 4 bits for every puzzle token
static constexpr uint64_t TARGET_POSITION = targetPosition(15); // target position

/** Vertex state, additional info for some puzzle position */
class VertexState
{
private:
    union {
        uint64_t integer; // raw internal state
        struct {
            packed_t zeroPos:4;     // index of the zero in the position
            packed_t movedUp:1;     // move in UP direction happened or forbidden
            packed_t movedLeft:1;   // move in LEFT direction happened or forbidden
            packed_t movedDown:1;   // move in DOWN direction happened or forbidden
            packed_t movedRight:1;  // move in RIGHT direction happened or forbidden
            packed_t moveCounter:3; // total count of moves happened
            packed_t parent:2;      // parent direction
            packed_t heuristics:16; // heuristics mode
            packed_t generation:16; // generation, distance from the starting position
        } bits;
    } data_;

    /* Reset move bitfields, assuming data_.integer was zeroed */
    void initMoves()
    {
        switch (data_.bits.zeroPos)
        {
        case 0:
           data_.bits.movedUp = 1;
           data_.bits.movedLeft = 1;
           data_.bits.moveCounter = 2;
        break;
        case 1:
        case 2:
            data_.bits.movedUp = 1;
            data_.bits.moveCounter = 1;
        break;
        case 3:
            data_.bits.movedUp = 1;
            data_.bits.movedRight = 1;
            data_.bits.moveCounter = 2;
        break;
        case 4:
        case 8:
            data_.bits.movedLeft = 1;
            data_.bits.moveCounter = 1;
        break;
        case 7:
        case 11:
            data_.bits.movedRight = 1;
            data_.bits.moveCounter = 1;
        break;
        case 5:
        case 6:
        case 9:
        case 10:
        break;
        case 12:
            data_.bits.movedDown = 1;
            data_.bits.movedLeft = 1;
            data_.bits.moveCounter = 2;
        break;
        case 13:
        case 14:
            data_.bits.movedDown = 1;
            data_.bits.moveCounter = 1;
        break;
        case 15:
            data_.bits.movedDown = 1;
            data_.bits.movedRight = 1;
            data_.bits.moveCounter = 2;
        break;
        default:
        assert(false);

        }
    }
    /* Reset move bitfields */
    void resetMoves()
    {
        switch (data_.bits.zeroPos)
        {
        case 0:
           data_.bits.movedUp     = 1;
           data_.bits.movedRight  = 0;
           data_.bits.movedDown   = 0;
           data_.bits.movedLeft   = 1;
           data_.bits.moveCounter = 2;
        break;
        case 1:
        case 2:
            data_.bits.movedUp     = 1;
            data_.bits.movedRight  = 0;
            data_.bits.movedDown   = 0;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 1;
        break;
        case 3:
            data_.bits.movedUp     = 1;
            data_.bits.movedRight  = 1;
            data_.bits.movedDown   = 0;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 2;
        break;
        case 4:
        case 8:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 0;
            data_.bits.movedDown   = 0;
            data_.bits.movedLeft   = 1;
            data_.bits.moveCounter = 1;
        break;
        case 7:
        case 11:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 1;
            data_.bits.movedDown   = 0;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 1;
        break;
        case 5:
        case 6:
        case 9:
        case 10:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 0;
            data_.bits.movedDown   = 0;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 0;
        break;
        case 12:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 0;
            data_.bits.movedDown   = 1;
            data_.bits.movedLeft   = 1;
            data_.bits.moveCounter = 2;
        break;
        case 13:
        case 14:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 0;
            data_.bits.movedDown   = 1;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 1;
        break;
        case 15:
            data_.bits.movedUp     = 0;
            data_.bits.movedRight  = 1;
            data_.bits.movedDown   = 1;
            data_.bits.movedLeft   = 0;
            data_.bits.moveCounter = 2;
        break;
        default:
        assert(false);

        }
    }
    /* Generate new position where zero is moved right for TShift positions */
    template <packed_t TShift>
    uint64_t swapZeroRight(const uint64_t position) const
    {
        assert(data_.bits.zeroPos - TShift < data_.bits.zeroPos);
        uint64_t mask = POS_MASK << (data_.bits.zeroPos - TShift) * 4;
        uint64_t swap = (position & mask) << TShift * 4;
        return (position & ~mask) | swap;
    }
    /* Generate new position where zero is moved left for TShift positions */
    template <packed_t TShift>
    uint64_t swapZeroLeft(const uint64_t position) const
    {
        assert(data_.bits.zeroPos + TShift < 16);
        uint64_t mask = POS_MASK << (data_.bits.zeroPos + TShift) * 4;
        uint64_t swap = (position & mask) >> TShift * 4;
        return (position & ~mask) | swap;
    }
public:
    /* Starting position ctor */
    explicit VertexState(packed_t zeroPos, packed_t heuristics): data_ {0}
    {
        assert(zeroPos < 16);
        data_.bits.zeroPos = zeroPos;
        initMoves();
        data_.bits.heuristics = heuristics;
    }
    /* Parent moved ctor, parent is mutated */
    VertexState(VertexState& parent, Move move, packed_t heuristics): data_{0}
    {
        switch (move)
        {
            case Move::UP:
                assert(parent.canMoveUp());
                parent.data_.bits.movedUp = 1;
                data_.bits.zeroPos = parent.data_.bits.zeroPos - 4;
            break;
            case Move::RIGHT:
                assert(parent.canMoveRight());
                parent.data_.bits.movedRight = 1;
                data_.bits.zeroPos = parent.data_.bits.zeroPos + 1;
            break;
            case Move::DOWN:
                assert(parent.canMoveDown());
                parent.data_.bits.movedDown = 1;
                data_.bits.zeroPos = parent.data_.bits.zeroPos + 4;
            break;
            case Move::LEFT:
                assert(parent.canMoveLeft());
                parent.data_.bits.movedLeft = 1;
                data_.bits.zeroPos = parent.data_.bits.zeroPos - 1;
            break;
        }
        initMoves();
        assert(parent.data_.bits.moveCounter < 4);
        ++parent.data_.bits.moveCounter;
        data_.bits.parent = static_cast<packed_t>(move);
        data_.bits.heuristics = heuristics;
        assert(parent.data_.bits.generation < (1 << 22) - 2);
        data_.bits.generation = parent.data_.bits.generation + 1;
    }
    /* Verbose ctor */
    VertexState(packed_t parentZeroPos, packed_t generation, Move move, packed_t heuristics): data_{0}
    {
        switch (move)
        {
            case Move::UP:
                data_.bits.zeroPos = parentZeroPos - 4;
            break;
            case Move::RIGHT:
                data_.bits.zeroPos = parentZeroPos + 1;
            break;
            case Move::DOWN:
                data_.bits.zeroPos = parentZeroPos + 4;
            break;
            case Move::LEFT:
                data_.bits.zeroPos = parentZeroPos - 1;
            break;
        }
        initMoves();
        data_.bits.parent = static_cast<packed_t>(move);
        data_.bits.heuristics = heuristics;
        assert(generation< (1 << 22) - 1);
        data_.bits.generation = generation;
    }
    /* Get child position from current position and move direction */
    uint64_t getMovePosition(const uint64_t position, Move move) const
    {
        switch (move)
        {
            case Move::UP:
                assert(data_.bits.zeroPos >= 4);
                return swapZeroRight<4>(position);
            case Move::RIGHT:
                assert((data_.bits.zeroPos & 3) != 3);
                return swapZeroLeft<1>(position);
            case Move::DOWN:
                assert(data_.bits.zeroPos <= 11);
                return swapZeroLeft<4>(position);
            case Move::LEFT:
                assert(data_.bits.zeroPos & 3);
                return swapZeroRight<1>(position);
        }
        assert(false);
    }
    /* Get parent position from current position */
    uint64_t getParentPosition(const uint64_t position) const
    {
        switch(static_cast<Move>(data_.bits.parent))
        {
            case Move::UP: return getMovePosition(position, Move::DOWN);
            case Move::RIGHT: return getMovePosition(position, Move::LEFT);
            case Move::DOWN: return getMovePosition(position, Move::UP);
            case Move::LEFT: return getMovePosition(position, Move::RIGHT);
        }
        assert(false);
    }
    /* Get parent move direction in output format */
    char getParentMoveChar() const
    {
        // note that they are inverted:
        switch(static_cast<Move>(data_.bits.parent))
        {
            case Move::UP: return 'D';
            case Move::RIGHT: return 'L';
            case Move::DOWN: return 'U';
            case Move::LEFT: return 'R';
        }
        assert(false);
    }
    /* Set new parent, returns if the new candidate parent has lower generation */
    bool setParent(const VertexState& parent, Move move)
    {
        assert(parent.data_.bits.generation < (1 << 22) - 2);
        if (parent.data_.bits.generation + 1 >= data_.bits.generation)
        {
            return false;
        }
        data_.bits.parent = static_cast<packed_t>(move);
        data_.bits.generation = parent.data_.bits.generation + 1;
        data_.bits.movedUp = data_.bits.movedRight = data_.bits.movedDown = data_.bits.movedLeft = 0;
        return true;
    }
    /* Check if a child in UP direction is possible */
    bool canMoveUp() const
    {
        return !data_.bits.movedUp && data_.bits.zeroPos >= 4;
    }
    /* Check if a child in RIGHT direction is possible */
    bool canMoveRight() const
    {
        return !data_.bits.movedRight && (data_.bits.zeroPos & 3) != 3;
    }
    /* Check if a child in DOWN direction is possible */
    bool canMoveDown() const
    {
        return !data_.bits.movedDown && data_.bits.zeroPos <= 11;
    }
    /* Check if a child in LEFT direction is possible */
    bool canMoveLeft() const
    {
        return !data_.bits.movedLeft && (data_.bits.zeroPos & 3);
    }
    /* Check if any child is possible */
    bool canMove() const
    {
        return data_.bits.moveCounter < 4;
    }
    /* Set that all children already exist */
    void setMoved()
    {
        data_.bits.moveCounter = 4;
    }
    /* Get generation */
    packed_t getGeneration() const
    {
        return data_.bits.generation;
    }
    /* Get heuristics */
    packed_t getHeuristics() const
    {
        return data_.bits.heuristics;
    }
    /* Get wiight */
    packed_t getWeight() const
    {
        return data_.bits.generation + data_.bits.heuristics;
    }
    /* Get zero position index */
    packed_t getZeroPos() const
    {
        return data_.bits.zeroPos;
    }
    /* Compare 2 states */
    bool operator<(const VertexState& other) const &
    {
        return getWeight() < other.getWeight();
    }
};

/*
 * Calculate position heuristics.
 * It is a doubled sum of all manhattan distances between all tokens and their target positions, excluding zero token.
 */
packed_t calcHeuristics(uint64_t position)
{
    static const std::pair<packed_t, packed_t> defaultCoords[16] = {
        { 3, 3 }, { 0, 0 }, { 0, 1 }, { 0, 2 }, { 0, 3 }, { 1, 0 }, { 1, 1 }, { 1, 2 }, { 1, 3 },
                  { 2, 0 }, { 2, 1 }, { 2, 2 }, { 2, 3 }, { 3, 0 }, { 3, 1 }, { 3, 2 }
    };

    packed_t result = 0;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            if (!(position & POS_MASK))
            {
                position >>= 4;
                continue;
            }
            const std::pair<int, int>& coords = defaultCoords[position & POS_MASK];
            result += abs(coords.first - i) + abs(coords.second - j);
            position >>= 4;
        }
    }
    return 2 * result;
}

/* Comparator for A* queue */
struct VertexComparator
{
    bool operator() (const std::pair<uint64_t, VertexState>& lhs, const std::pair<uint64_t, VertexState>& rhs)
    {
        return rhs.second < lhs.second;
    }
};

/* Solve puzzle 15 with A* algorithm */
std::vector<char> solve15AStar(const std::pair<uint64_t, packed_t> initialPosition)
{
    std::unordered_map<uint64_t, VertexState> visited;
    std::priority_queue<std::pair<uint64_t, VertexState>, std::vector<std::pair<uint64_t, VertexState>>, VertexComparator> queue;

    auto initVertex = std::make_pair(initialPosition.first, VertexState(initialPosition.second, calcHeuristics(initialPosition.first)));
    queue.emplace(initVertex);
    visited.insert(initVertex);

    while (queue.top().first != TARGET_POSITION)
    {
        auto vertex = queue.top();
        assert(vertex.second.getHeuristics() == calcHeuristics(vertex.first));
        const packed_t generation = vertex.second.getGeneration() + 1;
        const packed_t parentZeroPos = vertex.second.getZeroPos();
        queue.pop();
        auto processChild = [&visited, &queue, &vertex, &generation, &parentZeroPos](uint64_t childPosition, Move move) {
            packed_t heuristic = calcHeuristics(childPosition);
            auto found = visited.find(childPosition);
            if (found == visited.end())
            {
                // VertexState(packed_t parentZeroPos, packed_t generation, Move move, packed_t heuristics)
                auto newVertex = std::make_pair(childPosition, VertexState(parentZeroPos, generation, move, heuristic));
                visited.insert(newVertex);
                queue.push(newVertex);
            }
            else if (found->second.setParent(vertex.second, move))
            {
                queue.push(*found);
            }
        };

        if (vertex.second.canMoveUp())
        {
            processChild(vertex.second.getMovePosition(vertex.first, Move::UP), Move::UP);
        }
        if (vertex.second.canMoveRight())
        {
            processChild(vertex.second.getMovePosition(vertex.first, Move::RIGHT), Move::RIGHT);
        }
        if (vertex.second.canMoveDown())
        {
            processChild(vertex.second.getMovePosition(vertex.first, Move::DOWN), Move::DOWN);
        }
        if (vertex.second.canMoveLeft())
        {
            processChild(vertex.second.getMovePosition(vertex.first, Move::LEFT), Move::LEFT);
        }
    }

    auto head = queue.top();
    while(!queue.empty())
    {
        visited.insert(queue.top());
        queue.pop();
    }

    std::vector<char> path;

    while (head.second.getGeneration())
    {
        path.push_back(head.second.getParentMoveChar());
        head = *visited.find(head.second.getParentPosition(head.first));
    }

    return path;

}

/* Count all inversions among tokens */
std::size_t countInversions(packed_t* position)
{
    std::size_t result = 0;
    for (std::size_t i = 0; i < 15; ++i)
    {
        for (std::size_t j = i + 1; j < 16; ++j)
        {
            if (position[j] && position[i] && position[i] > position[j])
                ++result;
        }
    }
    return result;
}

/* Check if the puzzle starting position is solvable */
bool checkSolvability(packed_t* position)
{
    std::size_t rowIndex = 0;
    for (std::size_t i = 3; i < 4; --i)
    {
        for (std::size_t j = 3; j < 4; --j)
        {
            if (position[4 * i + j] == 0)
            {
                rowIndex = i;
                goto cycleOut;
            }
        }
    }
    cycleOut:

    bool inversionCountParity = (countInversions(position) & 1) == 0;

    return (rowIndex & 1) == inversionCountParity;
}

std::pair<uint64_t, packed_t> calcPosition(const packed_t* rawPosition)
{
    uint64_t position = 0;
    packed_t zeroIndex;
    for (std::size_t i = 0; i < 16; ++i)
    {
        assert(rawPosition[i] < 16);
        position |= static_cast<uint64_t>(rawPosition[i]) << 4 * i;
        if (!rawPosition[i])
        {
            zeroIndex = static_cast<packed_t>(i);
        }
    }
    return std::make_pair(position, zeroIndex);
}

void solve(std::istream& in)
{
    packed_t initialPositionRaw[16];
    for (std::size_t i = 0; i < 16; ++i)
    {
        in >> initialPositionRaw[i];
    }

    if (!checkSolvability(initialPositionRaw))
    {
        std::cout << "-1\n";
        return;
    }

    auto path = solve15AStar(calcPosition(initialPositionRaw));

    std::cout << path.size() << '\n';
    for (auto it = path.crbegin(); it < path.crend(); ++it)
    {
        std::cout << *it;
    }
    std::cout << '\n';
}


int main(void)
{
    assert(sizeof(VertexState) == 8);

    solve(std::cin);

    return 0;
}
