/*
 * Задача 18. Построение выпуклой оболочки
 * Дано множество точек на плоскости (x, y).
 * Постройте выпуклую оболочку этого множества и вычислите ее периметр.
 * 1. С помощью алгоритма Грэхема.
 */

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stack>
#include <vector>

using coord_t = double;
using point_t = std::pair<coord_t, coord_t>;

/* Calculate (pole; p0) x (pole; p1) */
inline double calcCrossProduct(const point_t& pole, const point_t& p0, const point_t& p1)
{
    return (p0.first - pole.first) * (p1.second - pole.second) - (p1.first - pole.first) * (p0.second - pole.second);
}

/* Calculate sqr|(p0, p1)| */
inline double distanceSq(const point_t& p0, const point_t& p1)
{
    const coord_t dX = p0.first - p1.first;
    const coord_t dY = p0.second - p1.second;
    return dX * dX + dY * dY;
}

/*
 * Find the point with min y, place it to the vector top position
 * then sort all other points by angle value (Ox, points[0], points[i]), where Ox is the catresian horizontal axis.
 */
void polarSort(std::vector<point_t>& points)
{
    std::size_t minIndex = 0;
    for (std::size_t i = 1; i < points.size(); ++i)
    {
        if ((points[i].second < points[minIndex].second)
                || (points[i].second == points[minIndex].second && points[i].first < points[minIndex].first))
        {
            minIndex = i;
        }
    }
    std::swap(points[0], points[minIndex]);

    point_t pole = points[0];
    auto sorter = [&pole](const point_t& lhs, const point_t& rhs) -> bool
    {
        const double crossProduct = calcCrossProduct(pole, rhs, lhs);

        bool result;
        if (fabs(crossProduct) == 0.)
        {
            result = distanceSq(pole, lhs) <= distanceSq(pole, rhs);
        }
        else
        {
            result = crossProduct < 0;
        }
        return result;
    };

    std::sort(points.begin() + 1, points.end(), sorter);
}

/* Find convex hull with Graham algorithm */
std::vector<point_t> findConvexHull(const std::vector<point_t>& points)
{
    assert(points.size() > 1);

    std::vector<point_t> hull;
    if (points.size() == 2)
    {
        return hull;
    }
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        const point_t& nextPoint = points[i];

        while (hull.size() > 1 && calcCrossProduct(hull[hull.size() - 2], hull.back(), nextPoint) <= 0)
        {
            hull.pop_back();
        }
        hull.push_back(nextPoint);
    }

    return hull;
}

/* Calculate closed polyline perimiter */
long double calcPerimiter(const std::vector<point_t>& points)
{
    long double result = 0;
    if (points.empty())
    {
        return result;
    }
    for (std::size_t i = 1; i < points.size(); ++i)
    {
        const coord_t dX = points[i].first - points[i - 1].first;
        const coord_t dY = points[i].second - points[i - 1].second;
        result += sqrt(dX * dX + dY * dY);
    }
    const coord_t dX = points.back().first - points.front().first;
    const coord_t dY = points.back().second - points.front().second;
    return result + sqrt(dX * dX + dY * dY);
}


int main(void)
{
    std::size_t count;
    std::cin >> count;

    std::vector<point_t> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        coord_t x, y;
        std::cin >> x >> y;
        points.emplace_back(x, y);
    }

    if (count < 2)
    {
        std::cout << 0;
        return 0;
    }

    polarSort(points);
    auto hull = findConvexHull(points);
    std::cout << std::setprecision(9) << calcPerimiter(hull) << '\n';

    return 0;
}
