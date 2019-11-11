/*
 * 4_4. Поиск k-той статистики.
 * Даны неотрицательные целые числа n,k и массив целых чисел из [0..10^9] размера n.
 * Требуется найти k-ю порядковую статистику. т.е. напечатать число, которое бы стояло на позиции с индексом k (0..n-1) в отсортированном массиве.
 * Напишите нерекурсивный алгоритм.
 * Требования к дополнительной памяти: O(n).
 * Требуемое среднее время работы: O(n).
 * Функцию Partition следует реализовывать методом прохода двумя итераторами в одном направлении.
 * Описание для случая прохода от начала массива к концу:
 *
 * - Выбирается опорный элемент. Опорный элемент меняется с последним элементом массива.
 * - Во время работы Partition в начале массива содержатся элементы, не бОльшие опорного.
 *   Затем располагаются элементы, строго бОльшие опорного.
 *   В конце массива лежат нерассмотренные элементы.
 *   Последним элементом лежит опорный.
 * - Итератор (индекс) i указывает на начало группы элементов, строго бОльших опорного.
 * - Итератор j больше i, итератор j указывает на первый нерассмотренный элемент.
 * - Шаг алгоритма. Рассматривается элемент, на который указывает j.
 *   Если он больше опорного, то сдвигаем j.
 *   Если он не больше опорного, то меняем a[i] и a[j] местами, сдвигаем i и сдвигаем j.
 * - В конце работы алгоритма меняем опорный и элемент, на который указывает итератор i.
 *
 * 4. Реализуйте стратегию выбора опорного элемента “случайный элемент”.
 *    Функцию Partition реализуйте методом прохода двумя итераторами от конца массива к началу.
 */

#include <assert.h>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <fstream>

/*
 * Выполнить частичную сортировку относительно опорного элемента.
 * \param array массив значений
 * \param beginIdx начальный индекс сортируемого отрезка
 * \param endIdx конечный индекс сортируемого отрезка
 * \param pivotIdx индекс опорного элемента
 * \param k индекс искомой статистики
 */
template <typename T> bool partition(T* array, std::size_t& beginIdx, std::size_t& endIdx, const std::size_t& pivotIdx, const std::size_t& k)
{
    T pivot = array[pivotIdx];
    std::swap(array[pivotIdx], array[endIdx]);
    std::size_t i = endIdx - 1,
                j = endIdx - 1;
    for (; j >= beginIdx && j <= endIdx; --j)
    {
        if (array[j] > pivot)
        {
            std::swap(array[j], array[i]);
            --i;
        }
    }
    ++i;
    if (k == i)
    {
        return true;
    }
    std::swap(array[i], array[endIdx]);
    if (k < i)
    {
        endIdx = i - 1;
    }
    else
    {
        beginIdx = i + 1;
    }
    return false;
}


/*
 * Вычислить k-тую статистику.
 * \param array массив значений
 * \param size размер массива
 * \param k индекс искомой статистики
 * \param findPivotIdx функция получения индекса опорного элемента
 */
template <typename T> T findKStatistic(T* array, std::size_t size, std::size_t k, std::size_t (*const findPivotIdx)(const T*, std::size_t, std::size_t))
{
    T* temp = new T[size];
    std::copy(array, array + size, temp);
    std::size_t beginIdx = 0;
    std::size_t endIdx = size - 1;
    std::size_t pivotIdx = findPivotIdx(temp, beginIdx, endIdx);


    while (!partition(temp, beginIdx, endIdx, pivotIdx, k))
    {
        pivotIdx = findPivotIdx(temp, beginIdx, endIdx);
    }

    T result = temp[endIdx];
    delete[] temp;
    return result;
}

/*
 * Получить случайный индекс массива на отрезке.
 * \param array массив
 * \param beginIdx начальная точка
 * \param endIdx конечная точка
 */
template <typename T> std::size_t getRandomIdx(const T* array, std::size_t begin, std::size_t end)
{
    return static_cast<std::size_t>(std::rand()) % (end - begin + 1) + begin;
}

int main(void) {
    using std::cin;

    std::size_t k, count;
    cin >> count;
    cin >> k;
    unsigned int* array = new unsigned int[count];
    for (std::size_t i = 0; i < count; i++) {
        cin >> array[i];
    }

    std::srand(std::time(nullptr));
    unsigned int kStat = findKStatistic(array, count, k, getRandomIdx);
    std::cout << kStat<< std::endl;

    delete[] array;

    return 0;
}
