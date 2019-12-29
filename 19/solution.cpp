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

#define DEBUG

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
/*
 * A scale for a linear dimension of geohash.
 * Determines index of geohash cell from sa coordinate.
 * \param TSearchDepth determines number of comparisons needed to determine the index.
 */
template <std::size_t TSearchDepth, std::size_t TCoupledSearchDepth>
class GeoHashScale
{
public:
    using GeoHashChildScale = GeoHashScale<TCoupledSearchDepth, TSearchDepth>;
private:
    /* In-order traversal index fo BST of depth TSearchDepth */
    static constexpr TreeInOrderIndex<TSearchDepth> TREE_INORDER_INDEX = TreeInOrderIndex<TSearchDepth>();

    coord_t min_; // inclusive minimal border of the range
    coord_t max_; // exclusive maximum value of the range
    coord_t tips_[calcScaleSize(TSearchDepth)]; // values that split the range in equal subranges
                                                 // this is not a linear array but a BST in form of array
    GeoHashChildScale* innerScales_[calcBinCount(TSearchDepth)]; // scales of the subranges
    /* Generate tips_ */
    void calcTips(coord_t start, coord_t end, std::size_t treeIndex);
    /* Build a scale for the subrange */
    void buildInnerScale(std::size_t index);
public:
    /* Default ctor */
    GeoHashScale(coord_t min, coord_t max);
    GeoHashScale(const GeoHashScale&) = delete;
    GeoHashScale& operator=(const GeoHashScale&) = delete;
    GeoHashScale(GeoHashScale&&) = delete;
    GeoHashScale& operator=(GeoHashScale&&) = delete;
    ~GeoHashScale();
    /* Finds the subrange index for a coord */
    std::size_t findBinIndex(coord_t coord) const;
    /* Get a scale for the subrange */
    GeoHashChildScale* getChildScale(std::size_t index);
    coord_t getMin() const
    {
        return min_;
    }
    coord_t getMax() const
    {
        return max_;
    }
};
/* Geohash cell, rectangular area, set of poins with the same geohash prefix */
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
class GeoHashCell
{
public:
    using GeoHashScaleX = GeoHashScale<TXSearchDepth, TYSearchDepth>;
    using GeoHashScaleY = GeoHashScale<TYSearchDepth, TXSearchDepth>;
    using GeoHashChildCell = GeoHashCell<TYSearchDepth, TXSearchDepth>;
private:
    /* Determines the number of points that causes division of the cell */
    const std::size_t SPLIT_BIN_COUNT = 32;

    std::size_t xIdx_;  // X coordinate in the parent, part of geocode value
    std::size_t yIdx_;  // Y coordinate in the parent, part of geocode value
    std::size_t depth_; // length of the path to the root cell
    std::size_t count_; // number of values in the cell

    GeoHashScaleX* scaleX_;      // link to the scale for children X coordinate indexing
    GeoHashScaleY* scaleY_;      // link to the scale for children Y coordinate indexing
    GeoHashChildCell* cells_;    // child cells
    std::vector<Point>* points_; // current points, either this or cells_ are present at one time
    /* Add point to the child cells */
    void addPointNested(const Point& point);
    /* Create child cells, pass all points into them */
    void splitHashBin();
public:
    GeoHashCell();
    GeoHashCell(std::size_t xIdx, std::size_t yIdx, std::size_t depth, GeoHashScaleX* scaleX, GeoHashScaleY* scaleY);
    GeoHashCell(const GeoHashCell&) = delete;
    GeoHashCell& operator=(const GeoHashCell&) = delete;
    GeoHashCell(GeoHashCell&&) = delete;
    GeoHashCell& operator=(GeoHashCell&& other);
    ~GeoHashCell();
    /* Add a point into the cell */
    void addPoint(const Point& point);
    /* Calculate points in the rectangle area */
    std::size_t getCountInRect(const Point& minPoint, const Point& maxPoint, const bool minXUnderflow, const bool maxXOverflow, const bool minYUnderflow, const bool maxYOverflow) const;
    std::size_t getCount() const
    {
        return count_;
    }
};

class GeoHash
{
private:
    static const std::size_t X_SEARCH_DEPTH = 3; // root search depth for the X coordinate
    static const std::size_t Y_SEARCH_DEPTH = 2; // root search depth for the Y coordinate
    const coord_t minX_; // root min X
    const coord_t maxX_; // root max X
    const coord_t minY_; // root min Y
    const coord_t maxY_; // root max Y

    GeoHashScale<X_SEARCH_DEPTH, Y_SEARCH_DEPTH> xtips_; // root X coordinate scale
    GeoHashScale<Y_SEARCH_DEPTH, X_SEARCH_DEPTH> ytips_; // root Y coordinate scale
    GeoHashCell<X_SEARCH_DEPTH, Y_SEARCH_DEPTH> root_;   // the root cell, contains the whole addressed area

public:
    GeoHash(coord_t minX, coord_t maxX, coord_t minY, coord_t maxY);
    /* Add a point into the geohash */
    void addPoint(const Point& point);
    /* Calculate points in the rectangle area */
    std::size_t getCountInRect(const Point& minPoint, const Point& maxPoint);
};

//////////////////////////////////////////////////////////////////////////////////////////////////

template<std::size_t TSearchDepth, std::size_t TComplementDepth>
void GeoHashScale<TSearchDepth, TComplementDepth>::calcTips(coord_t start, coord_t end, std::size_t treeIndex)
{
    const coord_t middle = start + (end - start) / 2;
    tips_[treeIndex] = middle;
    if (2 * treeIndex + 1 >= calcScaleSize(TSearchDepth))
    {
        return;
    }
    // recursion depth = TSearchDepth
    calcTips(start, middle, 2 * treeIndex + 1);
    calcTips(middle, end, 2 * treeIndex + 2);
}

template<std::size_t TSearchDepth, std::size_t TComplementDepth>
void GeoHashScale<TSearchDepth, TComplementDepth>::buildInnerScale(std::size_t index)
{
    assert(index < calcBinCount(TSearchDepth));
    coord_t min, max;
    if (!index)
    {
        min = min_;
        max = tips_[TREE_INORDER_INDEX.index[0]];
    }
    else if (index == calcBinCount(TSearchDepth) - 1)
    {
        min = tips_[TREE_INORDER_INDEX.index[calcBinCount(TSearchDepth) - 1]];
        max = max_;
    }
    else
    {
        min = tips_[TREE_INORDER_INDEX.index[index - 1]];
        max = tips_[TREE_INORDER_INDEX.index[index]];
    }
    assert(!innerScales_[index]);
    innerScales_[index] = new GeoHashChildScale(min, max);
}
template<std::size_t TSearchDepth, std::size_t TComplementDepth>
GeoHashScale<TSearchDepth, TComplementDepth>::GeoHashScale(coord_t min, coord_t max) :
    min_(min), max_(max)
{
    std::fill(innerScales_, innerScales_ + calcBinCount(TSearchDepth), nullptr);
    calcTips(min, max, 0);
}
template<std::size_t TSearchDepth, std::size_t TComplementDepth>
GeoHashScale<TSearchDepth, TComplementDepth>::~GeoHashScale()
{
    for (std::size_t i = 0; i < calcBinCount(TSearchDepth); ++i)
    {
        delete innerScales_[i];
    }
}
template<std::size_t TSearchDepth, std::size_t TComplementDepth>
std::size_t GeoHashScale<TSearchDepth, TComplementDepth>::findBinIndex(coord_t coord) const
{
    assert(coord >= min_ && coord < max_);
    std::size_t index = 0, result = 0;
    while (index < calcScaleSize(TSearchDepth))
    {
        result <<= 1;
        if (tips_[index] <= coord)
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
template<std::size_t TSearchDepth, std::size_t TComplementDepth>
typename GeoHashScale<TSearchDepth, TComplementDepth>::GeoHashChildScale* GeoHashScale<TSearchDepth, TComplementDepth>::getChildScale(std::size_t index)
{
    if (!innerScales_[index])
    {
        buildInnerScale(index);
    }
    return innerScales_[index];
}

template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
void GeoHashCell<TXSearchDepth, TYSearchDepth>::addPointNested(const Point& point)
{
    assert(cells_);
    std::size_t i = scaleX_->findBinIndex(point.x);
    std::size_t j = scaleY_->findBinIndex(point.y);
    cells_[i * calcBinCount(TYSearchDepth) + j].addPoint(point);
}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
void GeoHashCell<TXSearchDepth, TYSearchDepth>::splitHashBin()
{
    assert(points_ && !cells_);
    cells_ = new GeoHashChildCell[calcBinCount(TXSearchDepth) * calcBinCount(TYSearchDepth)];
    for (std::size_t i = 0; i < calcBinCount(TXSearchDepth); ++i)
    {
        for (std::size_t j = 0; j < calcBinCount(TYSearchDepth); ++j)
        {
            cells_[i * calcBinCount(TYSearchDepth) + j] = GeoHashChildCell(i, j, depth_ + 1, scaleX_->getChildScale(i), scaleY_->getChildScale(j));
        }
    }
    for (const Point& point : *points_)
    {
        addPointNested(point);
    }
    delete points_;
    points_ = nullptr;
}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
GeoHashCell<TXSearchDepth, TYSearchDepth>::GeoHashCell():
    cells_(nullptr), points_(nullptr)
{
}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
GeoHashCell<TXSearchDepth, TYSearchDepth>::GeoHashCell(std::size_t xIdx, std::size_t yIdx, std::size_t depth, GeoHashScaleX* scaleX, GeoHashScaleY* scaleY) :
    xIdx_(xIdx), yIdx_(yIdx), depth_(depth), count_(0), scaleX_(scaleX), scaleY_(scaleY), cells_(nullptr), points_(new std::vector<Point>)
{
    assert(scaleX && scaleY);
}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
GeoHashCell<TXSearchDepth, TYSearchDepth>& GeoHashCell<TXSearchDepth, TYSearchDepth>::operator=(GeoHashCell&& other)
{
    if (&other != this)
    {
        xIdx_ = other.xIdx_;
        yIdx_ = other.yIdx_;
        depth_ = other.depth_;
        count_ = other.count_;

        scaleX_ = other.scaleX_;
        scaleY_ = other.scaleY_;
        std::swap(cells_, other.cells_);
        std::swap(points_, other.points_);
    }

    return *this;
}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
GeoHashCell<TXSearchDepth, TYSearchDepth>::~GeoHashCell()
{
    delete[] cells_;
    delete points_;

}
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
void GeoHashCell<TXSearchDepth, TYSearchDepth>::addPoint(const Point& point)
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
template <std::size_t TXSearchDepth, std::size_t TYSearchDepth>
std::size_t GeoHashCell<TXSearchDepth, TYSearchDepth>::getCountInRect(const Point& minPoint, const Point& maxPoint,
    const bool minXUnderflow, const bool maxXOverflow, const bool minYUnderflow, const bool maxYOverflow) const
{
    assert (minPoint.x <= maxPoint.x && minPoint.y <= maxPoint.y);
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

    std::size_t minXIndex = minXUnderflow ? 0 : scaleX_->findBinIndex(minPoint.x);
    std::size_t maxXIndex = maxXOverflow ? calcBinCount(TXSearchDepth) - 1 : scaleX_->findBinIndex(maxPoint.x);
    std::size_t minYIndex = minYUnderflow ? 0 : scaleY_->findBinIndex(minPoint.y);
    std::size_t maxYIndex = maxYOverflow ? calcBinCount(TYSearchDepth) - 1: scaleY_->findBinIndex(maxPoint.y);

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
            result += cells_[i * calcBinCount(TYSearchDepth) + j].getCount();
        }
    }

    // 1. edges
    if (minXIndex != maxXIndex)
    {
        for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 6
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, minXUnderflow, true, true, true);
        }
        for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 2
        {
            result += cells_[maxXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, true, maxXOverflow, true, true);
        }
    }
    else
    {
        for (std::size_t j = minYIndex + 1; j < maxYIndex; ++j) // 2, 6
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + j].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, true, true);
        }
    }

    if (minYIndex != maxYIndex)
    {
        for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 0
        {
            result += cells_[i * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true, true, minYUnderflow, true);
        }
        for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 4
        {
            result += cells_[i * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, true, true, true, maxYOverflow);
        }
    }
    else
    {
        for (std::size_t i = minXIndex + 1; i < maxXIndex; ++i) // 0, 4
        {
            result += cells_[i * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true, true, minYUnderflow, maxYOverflow);
        }
    }

    // 2. vertices
    if (minXIndex != maxXIndex)
    {
        if (minYIndex != maxYIndex)
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         minYUnderflow, true        ); // 7
            result += cells_[maxXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, minYUnderflow, true        ); // 1
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         true,          maxYOverflow); // 5
            result += cells_[maxXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, true,          maxYOverflow); // 3
        }
        else
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, true,         minYUnderflow, maxYOverflow); // 5, 7
            result += cells_[maxXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, true,          maxXOverflow, minYUnderflow, maxYOverflow); // 1, 3
        }
    }
    else
    {
        if (minYIndex != maxYIndex)
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, minYUnderflow, true        ); // 1, 7
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + maxYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, true,          maxYOverflow); // 3, 5
        }
        else
        {
            result += cells_[minXIndex * calcBinCount(TYSearchDepth) + minYIndex].getCountInRect(minPoint, maxPoint, minXUnderflow, maxXOverflow, minYUnderflow, maxYOverflow); // 1, 3, 5, 7
        }
    }
    return result;
}

GeoHash::GeoHash(coord_t minX, coord_t maxX, coord_t minY, coord_t maxY) :
    minX_(minX), maxX_(maxX), minY_(minY), maxY_(maxY), xtips_(minX, maxX), ytips_(minY, maxY), root_(0, 0, 0, &xtips_, &ytips_)
{
}
void GeoHash::addPoint(const Point& point)
{
    assert(point.x >= minX_ && point.x < maxX_ && point.y >= minY_ && point.y < maxY_);
    root_.addPoint(point);
}
std::size_t GeoHash::getCountInRect(const Point& minPoint, const Point& maxPoint)
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
