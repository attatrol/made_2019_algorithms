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
//using geo_index_t = short;

struct Point
{
    coord_t x;
    coord_t y;
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
class GeoHashScale
{
public:
    using GeoHashChildScale = GeoHashScale<TCoupledSearchDepth, TSearchDepth>;
    static const std::size_t SEARCH_TERMINATED = 0;
private:
    static constexpr TreeInOrderIndex TREE_INORDER_INDEX = TreeInOrderIndex<TSearchDepth>();

    coord_t min_;
    coord_t max_;
    coord_t scale_[calcScaleSize(TSearchDepth)]; // this is not a linear array but a BST in form of array
    GeoHashChildScale* innerScales_[calcBinCount(TSearchDepth)];

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
        innerScales_[index] = new GeoHashChildScale(min, max);
    }
public:
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
        for (std::size_t i = 0; i < calcBinCount(TSearchDepth); ++i)
        {
            delete innerScales_[i];
        }
    }
    std::size_t searchBinIndex(coord_t coord) const
    {
        assert(coord >= min_ && coord < max_);
        std::size_t result = 0;
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
    GeoHashChildScale* operator[](std::size_t index)
    {
        if (!innerScales_[index])
        {
            buildInnerScale(index);
        }
        return innerScales_[index];
    }
};

template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
struct GeoHashBin
{
    using GeoHashScaleX = GeoHashScale<TXSearchDepth, TYSearchDepth>;
    using GeoHashScaleY = GeoHashScale<TYSearchDepth, TXSearchDepth>;
    using GeoHashChildBin = GeoHashBin<TYSearchDepth, TXSearchDepth>;

    const std::size_t SPLIT_BIN_COUNT = 32;

private:
    const std::size_t xIdx_;
    const std::size_t yIdx_;
    const std::size_t depth_;
    std::size_t count_;

    GeoHashScaleX* scaleX_;
    GeoHashScaleY* scaleY_;
    GeoHashChildBin* bins_;
    std::vector<Point>* points_;

    void addPointNested(Point&& point)
    {
        assert(bins_);
        bins_[scaleX_->searchBinIndex(point.x) * calcBinCount(TXSearchDepth) + scaleY_->searchBinIndex(point.y)]->addPoint(point);
    }
    void splitHashBin()
    {
        assert(points_ && !bins_);
        bins_ = new GeoHashBin<TYSearchDepth, TXSearchDepth>[calcBinCount(TXSearchDepth) * calcBinCount(TYSearchDepth)];
        for (std::size_t i = 0; i < calcBinCount(TXSearchDepth); ++i)
        {
            for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
            {
                bins_[i * calcBinCount(TXSearchDepth) + j] = new GeoHashChildBin(i, j, depth_ + 1, scaleX_[i], scaleY_[j]);
            }
        }
        for (const Point& point : *points_)
        {
            addPointNested(point);
        }
        delete points_;
        points_ = nullptr;
    }
public:
    GeoHashBin(std::size_t xIdx, std::size_t yIdx, std::size_t depth, GeoHashScale<TXSearchDepth, TYSearchDepth>* scaleX, GeoHashScale<TYSearchDepth, TXSearchDepth>* scaleY) :
        xIdx_(xIdx), yIdx_(yIdx), depth_(depth), count_(0), scaleX_(scaleX), scaleY_(scaleY), bins_(nullptr), points_(new std::vector<Point>)
    {
        assert(scaleX && scaleY);
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
    void addPoint(Point&& point)
    {
        ++count_;
        if (points_)
        {
            points_->push_back(point); // forward is unnecessary beacuse of small size
            if (count_ >= SPLIT_BIN_COUNT)
            {
                splitHashBin();
            }
        }
        else
        {
            addPointNested(point);
        }
    }
    std::size_t getCount() const
    {
        return count_;
    }
    //
    //    5_____________4___________3
    //    |                         |
    //    |                         |
    //    |                         |
    //   6|                         | 2
    //    |                         |
    //    |                         |
    //    |                         |
    //    7-------------------------1
    //               0
    //

    /* Get count of all points that are lesser than X coord (4) */
    std::size_t getCount_ltX(coord_t coord) const
    {
        std::size_t result = 0;
        if (points_)
        {
            for (const Point& point : points_)
            {
                if (point.x < coord)
                {
                    ++result;
                }
            }
        }
        else
        {
            std::size_t index = scaleX_->searchBinIndex(coord);
            for (std::size_t i = 0; i < index; ++i)
            {
                for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
                {
                    result += bins_[i * calcBinCount(TXSearchDepth) + j]->getCount();
                }
            }
            for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
            {
                result += bins_[index * calcBinCount(TXSearchDepth) + j]->getCount_ltX(coord);
            }
        }
        return result;
    }
    /* Get count of all points that are greater or equal than X coord (0) */
    std::size_t getCount_gteX(coord_t coord) const
    {
        std::size_t result = 0;
        if (points_)
        {
            for (const Point& point : points_)
            {
                if (point.x >= coord)
                {
                    ++result;
                }
            }
        }
        else
        {
            std::size_t index = scaleX_->searchBinIndex(coord);
            for (std::size_t i = index + 1; i < calcBinCount(TYSearchDepth); ++i)
            {
                for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
                {
                    result += bins_[i * calcBinCount(TXSearchDepth) + j]->getCount();
                }
            }
            for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
            {
                result += bins_[index * calcBinCount(TXSearchDepth) + j]->getCount_gteX(coord);
            }
        }
        return result;
    }
    /* Get count of all points that are lesser than Y coord (2) */
    std::size_t getCount_ltY(coord_t coord) const
    {

    }
    /* Get count of all points that are greater or equal than Y coord (6) */
    std::size_t getCount_gteY(coord_t coord) const
    {

    }
    /* Get count of all points that are lesser than X coord and lesser than Y coord (3) */
    std::size_t getCount_ltX_ltY(const Point& point) const
    {

    }
    /* Get count of all points that are greater or equal than X coord and greater or equal than Y coord (7) */
    std::size_t getCount_gteX_gteY(const Point& point) const
    {

    }
    /* Get count of all points that are lesser than X coord and greater or equal than Y coord (2) */
    std::size_t getCount_ltX_gteY(const Point& point) const
    {

    }
    /* Get count of all points that are greater or equal than X coord and lesser than Y coord (2) */
    std::size_t getCount_gteX_ltY(const Point& point) const
    {

    }
};

class GeoHash
{
private:
    const coord_t minX_;
    const coord_t maxX_;
    const coord_t minY_;
    const coord_t maxY_;
public:
    GeoHash(coord_t minX, coord_t maxX, coord_t minY, coord_t maxY) :
        minX_(minX), maxX_(maxX), minY_(minY), maxY_(maxY)
    {

    }

};

constexpr TreeInOrderIndex<2> x;

int main(void)
{
    std::cout << x.index[0];

    return 0;
}
