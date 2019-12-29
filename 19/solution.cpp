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

//#define DEBUG

#include <algorithm>
#include <assert.h>
#include <iostream>
#include <stack>
#include <vector>

#ifdef DEBUG
#include <fstream>
#include <random>
#endif

using coord_t = double;
//using geo_index_t = short;

/* 2D point */
struct Point
{
    coord_t x;
    coord_t y;

    Point(coord_t _x, coord_t _y) :
        x(_x), y(_y)
    {
    }
};

/* Calculate pow(2, searchDepth) at compile time*/
constexpr std::size_t calcBoolean(std::size_t searchDepth)
{
    return searchDepth < 1 ? 1 : 2 * calcBoolean(searchDepth - 1);
}
/* Calculate number of inner dividing points that split a range into pow(2, searchDepth) distinct ranges at compile time */
constexpr std::size_t calcScaleSize(std::size_t searchDepth)
{
    return calcBoolean(searchDepth) - 1;
}
/* Calculate number of distict results of a search tree of a defined depth at compile time */
constexpr std::size_t calcBinCount(std::size_t searchDepth)
{
    return calcBoolean(searchDepth);
}
/* In-order traversal index for a complete binary search tree */
template <std::size_t TSearchDepth>
struct TreeInOrderIndex {

    std::size_t index[calcScaleSize(TSearchDepth)]; // the index

    constexpr TreeInOrderIndex();
};
template <>
constexpr TreeInOrderIndex<1>::TreeInOrderIndex(): index { 0 }
{
}
template <std::size_t TSearchDepth>
constexpr TreeInOrderIndex<TSearchDepth>::TreeInOrderIndex() : index {}
{
    static_assert(TSearchDepth > 1);
    constexpr TreeInOrderIndex<TSearchDepth - 1> previous;
    const std::size_t expandNodeMask = 1 << (TSearchDepth - 2);

    for (std::size_t i = 0, j = 0; i < calcScaleSize(TSearchDepth - 1); ++i)
    {
        if ((previous.index[i] + 1) & expandNodeMask)
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
/* A scale for a dimension of */
template <std::size_t TSearchDepth, std::size_t TCoupledSearchDepth>
class GeoHashScale
{
public:
    using GeoHashChildScale = GeoHashScale<TCoupledSearchDepth, TSearchDepth>;
//    static const std::size_t SEARCH_TERMINATED = 0;
private:
    static constexpr TreeInOrderIndex<TSearchDepth> TREE_INORDER_INDEX = TreeInOrderIndex<TSearchDepth>();

    coord_t min_;
    coord_t max_;
    coord_t scale_[calcScaleSize(TSearchDepth)]; // this is not a linear array but a BST in form of array
    GeoHashChildScale* innerScales_[calcBinCount(TSearchDepth)];

    void calcTips(coord_t start, coord_t end, std::size_t treeIndex)
    {
        const coord_t middle = start + (end - start) / 2;
        scale_[treeIndex] = middle;
        if (2 * treeIndex + 1 >= calcScaleSize(TSearchDepth))
        {
            return;
        }
        // recursion depth = TSearchDepth
        calcTips(start, middle, 2 * treeIndex + 1);
        calcTips(middle, end, 2 * treeIndex + 2);
    }

    void buildInnerScale(std::size_t index)
    {
        assert(index < calcBinCount(TSearchDepth));
//        std::size_t treeNodeIndex = index + calcScaleSize(TSearchDepth - 1);
//        assert(treeNodeIndex < calcScaleSize(TSearchDepth));
        coord_t min, max;
        if (!index)
        {
            min = min_;
            max = scale_[TREE_INORDER_INDEX.index[0]];
        }
        else if (index == calcBinCount(TSearchDepth) - 1)
        {
            min = scale_[TREE_INORDER_INDEX.index[calcBinCount(TSearchDepth) - 1]];
            max = max_;
        }
        else
        {
            min = scale_[TREE_INORDER_INDEX.index[index - 1]];
            max = scale_[TREE_INORDER_INDEX.index[index]];
//            min = scale_[treeNodeIndex];
//            std::size_t nextTreeNodeIndex = (treeNodeIndex + 1) / 2;
//            while (nextTreeNodeIndex & 1)
//            {
//                nextTreeNodeIndex /= 2;
//            }
//            //assert(nextTreeNodeIndex);
//            assert(nextTreeNodeIndex < calcScaleSize(TSearchDepth));
//            max = scale_[nextTreeNodeIndex];
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
        std::size_t index = 0, result = 0;
        while (index < calcScaleSize(TSearchDepth))
        {
            result <<= 1;
            if (scale_[index] <= coord)
            {
                index = 2 * index + 2;
                ++result;
            }
            else
            {
                index = 2 * index + 1;
            }
        }
        assert(result < calcBinCount(TSearchDepth));
        assert(!innerScales_[result] || innerScales_[result]->getMin() <= coord && innerScales_[result]->getMax() > coord);
        return result;
    }
    GeoHashChildScale* getChildScale(std::size_t index)
    {
        if (!innerScales_[index])
        {
            buildInnerScale(index);
        }
        return innerScales_[index];
    }
    coord_t getMin() const
    {
        return min_;
    }
    coord_t getMax() const
    {
        return max_;
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
    std::size_t xIdx_;
    std::size_t yIdx_;
    std::size_t depth_;
    std::size_t count_;

    GeoHashScaleX* scaleX_;
    GeoHashScaleY* scaleY_;
    GeoHashChildBin* bins_;
    std::vector<Point>* points_;

    void addPointNested(const Point& point)
    {
        assert(bins_);
        std::size_t i = scaleX_->searchBinIndex(point.x);
        std::size_t j = scaleY_->searchBinIndex(point.y);
        bins_[i * calcBinCount(TYSearchDepth) + j].addPoint(point);
    }
    void splitHashBin()
    {
        assert(points_ && !bins_);
        bins_ = new GeoHashChildBin[calcBinCount(TXSearchDepth) * calcBinCount(TYSearchDepth)];
        for (std::size_t i = 0; i < calcBinCount(TXSearchDepth); ++i)
        {
            for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
            {
                bins_[i * calcBinCount(TYSearchDepth) + j] = GeoHashChildBin(i, j, depth_ + 1, scaleX_->getChildScale(i), scaleY_->getChildScale(j));
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
    GeoHashBin():
        bins_(nullptr), points_(nullptr)
    {
    }
    GeoHashBin(std::size_t xIdx, std::size_t yIdx, std::size_t depth, GeoHashScaleX* scaleX, GeoHashScaleY* scaleY) :
        xIdx_(xIdx), yIdx_(yIdx), depth_(depth), count_(0), scaleX_(scaleX), scaleY_(scaleY), bins_(nullptr), points_(new std::vector<Point>)
    {
        assert(scaleX && scaleY);
    }
    GeoHashBin(const GeoHashBin&) = delete;
    GeoHashBin& operator=(const GeoHashBin&) {

    }
    GeoHashBin(GeoHashBin&&) = delete;
    GeoHashBin& operator=(GeoHashBin&& other)
    {
        if (&other != this)
        {
            xIdx_ = other.xIdx_;
            yIdx_ = other.yIdx_;
            depth_ = other.depth_;
            count_ = other.count_;

            scaleX_ = other.scaleX_;
            scaleY_ = other.scaleY_;
            std::swap(bins_, other.bins_);
            std::swap(points_, other.points_);
        }

        return *this;
    }
    ~GeoHashBin()
    {
        delete[] bins_;
        delete points_;

    }
    void addPoint(const Point& point)
    {
        ++count_;
        if (points_)
        {
            points_->push_back(point);
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

    std::size_t getCountInRect(const Point& minPoint, const Point& maxPoint, const bool minXUnderflow, const bool maxXOverflow, const bool minYUnderflow, const bool maxYOverflow) const
    {
        assert (minPoint.x <= maxPoint.x && minPoint.y <= maxPoint.y);
        //assert(!(minXUnderflow && maxXOverflow && minYUnderflow && maxYOverflow));
        if (minXUnderflow && maxXOverflow && minYUnderflow && maxYOverflow)
        {
            return count_;
        }

        std::size_t result = 0;

        if (points_)
        {
            for (const Point& point : *points_)
            {
                if (point.x >= minPoint.x && point.y >= minPoint.y && point.x < maxPoint.x && point.y < maxPoint.y )
                {
                    ++result;
                }
            }
            return result;
        }

        std::size_t minXIndex = minXUnderflow ? 0 : scaleX_->searchBinIndex(minPoint.x);
        std::size_t maxXIndex = maxXOverflow ? calcBinCount(TXSearchDepth) - 1 : scaleX_->searchBinIndex(maxPoint.x);
        std::size_t minYIndex = minYUnderflow ? 0 : scaleY_->searchBinIndex(minPoint.y);
        std::size_t maxYIndex = maxYOverflow ? calcBinCount(TYSearchDepth) - 1: scaleY_->searchBinIndex(maxPoint.y);

        //
        //    5_____________4___________3
        //    |                         |
        //    |                         |
        //    |                         |
        //   6|       inner region      | 2
        //    |                         |
        //    |                         |
        //    |                         |
        //    7-------------0-----------1
        //

        // 0. inner region
        for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i)
        {
            for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j)
            {
                result += bins_[i * calcBinCount(TYSearchDepth) + j].getCount();
            }
        }

        // 1. edges
        if (minXIndex != maxXIndex)
        {
            for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 6
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, minXUnderflow, true, true, true);
            }
            for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 2
            {
                result += bins_[maxXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, true, maxXOverflow, true, true);
            }
        }
        else
        {
            for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 2, 6
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, true, true);
            }
        }

        if (minYIndex != maxYIndex)
        {
            for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 0
            {
                result += bins_[i * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true, true, minYUnderflow, true);
            }
            for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 4
            {
                result += bins_[i * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, true, true, true, maxYOverflow);
            }
        }
        else
        {
            for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 0, 4
            {
                result += bins_[i * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true, true, minYUnderflow, maxYOverflow);
            }
        }

        // 2. vertices
        if (minXIndex != maxXIndex)
        {
            if (minYIndex != maxYIndex)
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         minYUnderflow, true        ); // 7
                result += bins_[maxXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, minYUnderflow, true        ); // 1
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         true,          maxYOverflow); // 5
                result += bins_[maxXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, true,          maxYOverflow); // 3
            }
            else
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         minYUnderflow, maxYOverflow); // 5, 7
                result += bins_[maxXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, minYUnderflow, maxYOverflow); // 1, 3
            }
        }
        else
        {
            if (minYIndex != maxYIndex)
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, minYUnderflow, true        ); // 1, 7
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, true,          maxYOverflow); // 3, 5
            }
            else
            {
                result += bins_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, minYUnderflow, maxYOverflow); // 1, 3, 5, 7
            }
        }
        return result;
    }
};

class GeoHash
{
private:
    static const std::size_t X_SEARCH_DEPTH = 3;
    static const std::size_t Y_SEARCH_DEPTH = 2;
    const coord_t minX_;
    const coord_t maxX_;
    const coord_t minY_;
    const coord_t maxY_;

    GeoHashScale<X_SEARCH_DEPTH, Y_SEARCH_DEPTH> xScale_;
    GeoHashScale<Y_SEARCH_DEPTH, X_SEARCH_DEPTH> yScale_;
    GeoHashBin<X_SEARCH_DEPTH, Y_SEARCH_DEPTH> root_;

public:
    GeoHash(coord_t minX, coord_t maxX, coord_t minY, coord_t maxY) :
        minX_(minX), maxX_(maxX), minY_(minY), maxY_(maxY), xScale_(minX, maxX), yScale_(minY, maxY), root_(0, 0, 0, &xScale_, &yScale_)
    {
    }

    void addPoint(const Point& point)
    {
        assert(point.x >= minX_ && point.x < maxX_ && point.y >= minY_ && point.y < maxY_);
        root_.addPoint(point);
    }

    std::size_t getCountInRect(const Point& minPoint, const Point& maxPoint)
    {
         assert (minPoint.x <= maxPoint.x && minPoint.y <= maxPoint.y);
         if (minPoint.x >= maxX_ || minPoint.y >= maxY_ || maxPoint.x < minX_ || maxPoint.y < minY_)
         {
             return 0;
         }
         bool minXUnderflow = minPoint.x < minX_;
         bool maxXOverflow  = maxPoint.x >= maxX_;
         bool minYUnderflow = minPoint.y < minY_;
         bool maxYOverflow  = maxPoint.y >= maxY_;
         if (minXUnderflow && maxXOverflow && minYUnderflow && maxYOverflow)
         {
             return root_.getCount();
         }
         else
         {
            return root_.getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, minYUnderflow, maxYOverflow);
         }
    }


};

#ifdef DEBUG
coord_t fRand(coord_t fMin, coord_t fMax)
{
    coord_t f = (coord_t)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
Point generatePoint()
{
    return Point(fRand(-180., 180.), fRand(-9, 90));
}
std::pair<Point, Point> generateRect()
{
    coord_t x0 = fRand(-200., 200.);
    coord_t x1 = fRand(-200., 200.);
    coord_t y0 = fRand(-100., 100.);
    coord_t y1 = fRand(-100., 100.);
    return std::make_pair(Point(std::min(x0, x1), std::min(y0, y1)), Point(std::max(x0, x1), std::max(y0, y1)));
}
void generateTestData(std::size_t pointCount, std::size_t rectCount)
{
    std::ofstream out("input_test.txt");
    out << pointCount << '\n';
    for (std::size_t i = 0; i < pointCount; ++i)
    {
        Point point = generatePoint();
        out << point.x << " " << point.y << '\n';
    }
    out << rectCount << '\n';
    for (std::size_t i = 0; i < rectCount; ++i)
    {
        auto rect = generateRect();
        out << rect.first.x << ' ' << rect.first.y << ' ' << rect.second.x << ' ' << rect.second.y << '\n';
    }
}
#endif

int main(void)
{
    GeoHash geoHash(-180, 180, -90, 90);

#ifdef DEBUG
    std::ifstream cin("input_test.txt");
#else
    using std::cin;
#endif
    std::size_t pointCount;
    cin >> pointCount;
#ifdef DEBUG
    std::cout << "calcBinCount(2) " << calcBinCount(2) << '\n';
    std::cout << "calcBinCount(3) " << calcBinCount(3) << '\n';
    std::vector<Point> testPoints;
    testPoints.reserve(pointCount);
#endif
    for (std::size_t i = 0; i < pointCount; ++i)
    {
        coord_t x, y;
        cin >> x >> y;
        geoHash.addPoint(Point(x, y));
#ifdef DEBUG
        testPoints.emplace_back(x, y);
#endif
    }
    std::size_t rectCount;
    cin >> rectCount;
    for (std::size_t i = 0; i < rectCount; ++i)
    {
        coord_t xMin, xMax, yMin, yMax;
        cin >> xMin >> yMin >> xMax >> yMax;
        std::size_t pointCount = geoHash.getCountInRect(Point(xMin, yMin), Point(xMax, yMax));
#ifdef DEBUG
        std::size_t testCount = 0;
        for (const Point& point : testPoints)
        {
            if (point.x >= xMin && point.x < xMax && point.y >= yMin && point.y < yMax)
            {
                ++testCount;
            }
        }
        assert(pointCount == testCount);
#endif
        std::cout << pointCount << '\n';
    }

    //generateTestData(10000, 10000);

//    TreeInOrderIndex<1> t1;
//    TreeInOrderIndex<2> t2;
//    TreeInOrderIndex<3> t3;
//    TreeInOrderIndex<4> t4;

    return 0;
}
