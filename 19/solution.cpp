/* Задача 19. Поиск точек в прямоугольникe
 *
 * На вход подаются точки и прямоугольники.
 * Точка задается двумя координатами (x, y).
 * Прямоугольник - четверкой чисел [left, bottom, right, top].
 * Точка (left, bottom) - принадлежит прямоугольнику, точка (right, top) - нет. (left < right, bottom < top)
 * Для каждого прямоугольника нужно вывести, сколько добавленных точек он содержит.
 *
 * X  в диапазоне [-180, 180)
 * Y -  [-90, 90)
 * Количество точек <= 100000, Количество прямоугольников <= 1000
 * Для решения задачи необходимо реализовать алгоритм “Geohashing”.
 */

#include <assert.h>
#include <iostream>
#include <stack>
#include <vector>

using coord_t = double;

struct Point
{
    coord_t lat_;
    coord_t lon_;
};

constexpr std::size_t calcBoolean(std::size_t searchDepth)
{
    return searchDepth < 1 ? 1 : 2 * calcBoolean(searchDepth - 1);
}
constexpr std::size_t calcScaleSize(std::size_t searchDepth)
{
    return calcBoolean(searchDepth) - 1;
}
constexpr std::size_t calcBinCount(std::size_t searchDepth)
{
    return calcBoolean(searchDepth - 1);
}
//template <std::size_t TSearchDepth>
//std::size_t[] generateInOrderTreeTraverseIndex
template <std::size_t TSearchDepth>
struct TreeInOrderIndex {

    std::size_t index[calcScaleSize(TSearchDepth)];

    constexpr TreeInOrderIndex();
};

template <>
constexpr TreeInOrderIndex<1>::TreeInOrderIndex(): index {}
{
}

template <std::size_t TSearchDepth>
constexpr TreeInOrderIndex<TSearchDepth>::TreeInOrderIndex() : index {}
{
    static_assert(TSearchDepth > 1);
    constexpr TreeInOrderIndex<TSearchDepth - 1> previous;
    const std::size_t expandNodeMask = 1 << (TSearchDepth - 1);

    for (std::size_t i = 0, j = 0; i < calcScaleSize(TSearchDepth - 1); ++i)
    {
        if ((previous.index[i] + 1) && expandNodeMask)
        {
            index[j++] = (previous.index[i] + 1) * 2 - 1;
            index[j++] = previous.index[i];
            index[j++] = (previous.index[i] + 1) * 2;
        }
        else
        {
            index[j++] = previous.index[i];
        }
    }
}

template <std::size_t TSearchDepth, std::size_t TCoupledSearchDepth>
struct GeoHashScale
{
    static constexpr std::size_t TREE_INORDER_INDEX[TSearchDepth] = TreeInOrderIndex<TSearchDepth>();

    coord_t min_;
    coord_t max_;
    coord_t scale_[calcScaleSize(TSearchDepth)]; // this is not a linear array but a BST in form of array
    GeoHashScale<TCoupledSearchDepth, TSearchDepth>* innerScales_[calcBinCount(TSearchDepth)];

    void calcTips(coord_t start, coord_t end, std::size_t treeIndex)
    {
        const coord_t middle = (end - start) / 2;
        scale_[treeIndex] = middle;
        if (2 * treeIndex >= calcScaleSize(TSearchDepth))
        {
            return;
        }
        // recursion depth = TSearchDepth
        calcTips(start, middle, 2 * treeIndex);
        calcTips(middle, end, 2 * treeIndex + 1);
    }

    short searchBinIndex(coord_t coord) const
    {
        short result = 0;
        while (result <= calcScaleSize(TSearchDepth - 1))
        {
            if (scale_[result] < coord)
            {
                result *= 2;
            }
            else
            {
                result = 2 * result + 1;
            }
        }
        result -= calcScaleSize(TSearchDepth - 1);
        assert(result >= 0 && result < calcBinCount(TSearchDepth));
        return result;
    }

    void buildInnerScale(std::size_t index)
    {
        assert(index < calcBinCount(TSearchDepth));
        std::size_t treeNodeIndex = index + calcScaleSize(TSearchDepth - 1);
        assert(index < calcScaleSize(TSearchDepth));
        std::size_t min, max;
        if (!index)
        {
            min = min_;
            max = scale_[treeNodeIndex];
        }
        else if (index == calcBinCount(TSearchDepth) - 1)
        {
            min = scale_[treeNodeIndex];
            max = max_;
        }
        else
        {
            min = scale_[treeNodeIndex];
            max = (treeNodeIndex + 1) / 2;
            while (max & 1)
            {
                max /= 2;
            }
            assert(max);
            assert(max - 1 < calcScaleSize(TSearchDepth));
            max = scale_[max - 1];
        }
        assert(!innerScales_[index]);
        innerScales_[index] = new GeoHashScale<TCoupledSearchDepth, TSearchDepth>(min, max);
    }
    void tryBuildInnerScale(std::size_t index)
    {
        if (!innerScales_[index])
        {
            buildInnerScale(index);
        }
    }
    GeoHashScale(coord_t min, coord_t max) :
        min_(min), max_(max)
    {
        std::fill(innerScales_, innerScales_ + calcBinCount(TSearchDepth), nullptr);
        calcTips(min, max, 0);
    }
    GeoHashScale(const GeoHashScale&) = delete;
    GeoHashScale& operator=(const GeoHashScale&) = delete;
    GeoHashScale(GeoHashScale&&) = delete;
    GeoHashScale& operator=(GeoHashScale&&) = delete;
    ~GeoHashScale()
    {
        delete[] innerScales_;
    }
};

template <std::size_t TSearchDepth, std::size_t TComplementSearchDepth>
struct GeoHashBin
{
    short latIdx_;
    short lonIdx_;
    short depth_;
    std::size_t count_;

    GeoHashScale<TSearchDepth, TComplementSearchDepth>* scale_;
    GeoHashBin<TComplementSearchDepth, TSearchDepth>* bins_;
    std::vector<std::size_t>* points_;

    GeoHashBin(short latIdx, short lonIdx, short depth) :
        latIdx_(latIdx), lonIdx_(lonIdx), depth_(depth), count_(0), bins_(nullptr), points_(new std::vector<std::size_t>)
    {
    }
    GeoHashBin(const GeoHashBin&) = delete;
    GeoHashBin& operator=(const GeoHashBin&) = delete;
    GeoHashBin(GeoHashBin&&) = delete;
    GeoHashBin& operator=(GeoHashBin&&) = delete;
    ~GeoHashBin()
    {
        delete[] bins_;
        delete points_;
    }
};

//class GeoHash
//{
//    GeoHashBin
//};

constexpr TreeInOrderIndex<2> x;

int main(void)
{
    std::cout << x.index[0];

    return 0;
}
