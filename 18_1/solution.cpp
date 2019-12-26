/*
 * Задача 18. Построение выпуклой оболочки
 * Дано множество точек на плоскости (x, y).
 * Постройте выпуклую оболочку этого множества и вычислите ее периметр.
 * 1. С помощью алгоритма Грэхема.
 */

#define DEBUG

#include <algorithm>
#include <assert.h>
#include <cmath>
#ifdef DEBUG
#include <fstream>
#endif
#include <iomanip>
#include <iostream>
#include <limits>
#include <stack>
#include <vector>

using coord_t = long double;
using point_t = std::pair<coord_t, coord_t>;
static const double EPS = 0.;

inline double calcCrossProduct(const point_t& pole, const point_t& v0, const point_t& v1)
{
    return (v0.first - pole.first) * (v1.second - pole.second) - (v1.first - pole.first) * (v0.second - pole.second);
}

inline double distanceSq(const point_t& p0, const point_t& p1)
{
    const coord_t dX = p0.first - p1.first;
    const coord_t dY = p0.second - p1.second;
    return dX * dX + dY * dY;
}

template <typename T>
int sign(const T val)
{
    return val < 0 ? -1 : (val > 0 ? 1 : 0);
    //return (T(0) < val) - (val < T(0));
}

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
        if (fabs(crossProduct) <= EPS)
        {
            result = distanceSq(pole, lhs) <= distanceSq(pole, rhs);
        }
        else
        {
            result = crossProduct < 0;
        }
#ifdef DEBUG
        std::cout << '(' << lhs.first << "; " << lhs.second << (result ? ") < " : ") >= (") << rhs.first << "; " << rhs.second <<")\n";
#endif
        return result;
    };
    std::sort(points.begin() + 1, points.end(), sorter);

#ifdef DEBUG
    std::cout << "Sorted\n";
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        std::cout << i << " (" << points[i].first << "; " << points[i].second << ")\n";
    }
    std::cout << std::endl;
#endif
}

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

#ifdef DEBUG
    std::cout << "Hull\n";
    for (std::size_t i = 0; i < hull.size(); ++i)
    {
        std::cout << i << " = (" << hull[i].first << "; " << hull[i].second << ")\n";
    }
    std::cout << std::endl;
#endif

    return hull;
}


int main(void)
{
#ifdef DEBUG
    std::fstream cin("input.txt");
#else
    using std::cin;
#endif

    std::size_t count;
    cin >> count;

    std::vector<point_t> points;
    points.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
    {
        coord_t x, y;
        cin >> x >> y;
        points.emplace_back(x, y);
    }

#ifdef DEBUG
    std::cout << "Input\n";
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        std::cout << i << " (" << points[i].first << "; " << points[i].second << ")\n";
    }
    std::cout << std::endl;
#endif

    if (count < 2)
    {
        std::cout << 0;
        return 0;
    }

    polarSort(points);

    auto hull = findConvexHull(points);


    long double perimiter = 0;
    for (std::size_t i = 1; i < hull.size(); ++i)
    {
        const coord_t dX = hull[i].first - hull[i - 1].first;
        const coord_t dY = hull[i].second - hull[i - 1].second;
        perimiter += sqrt(dX * dX + dY * dY);
    }
    const coord_t dX = hull.back().first - hull.front().first;
    const coord_t dY = hull.back().second - hull.front().second;
    perimiter += sqrt(dX * dX + dY * dY);

    std::cout << std::setprecision(9) << perimiter << '\n';

    return 0;
}
