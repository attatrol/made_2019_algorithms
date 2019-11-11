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
 * Вычислить k-тую статистику.
 * \param array массив значений
 * \param size размер массива
 * \param k индекс искомой статистики
 * \param findPivotIdx функция получения индекса опорного элемента
 */
template <typename T> T findKStatistic(T* array, std::size_t size, std::size_t k, std::size_t (*const findPivotIdx)(const T*, std::size_t, std::size_t))
{
    T* temp = new T[size];

    std::size_t pivotIdx = findPivotIdx(array, 0, size - 1);
    T pivot = array[pivotIdx];
    std::size_t i;
    std::size_t j;
    // первый цикл алгоритма совмещен с заполнением вспомогательного массива:
    for (i = size - 1, j = size - 1; j < size; --j)
    {
        if (array[j] >= pivot)
        {
            temp[j] = temp[i];
            temp[i] = array[j];
            --i;
        }
        else
        {
            temp[j] = array[j];
        }
    }
    ++i;
    if (i == k)
    {
        delete[] temp;
        return pivot;
    }
    std::size_t beginIdx;
    std::size_t endIdx;
    if (k < i)
    {
        beginIdx = 0;
        endIdx = i - 1;
    }
    else
    {
        beginIdx = i;
        endIdx = size - 1;
    }

    // Основной рабочий цикл
    while (k != i)
    {
        pivotIdx = findPivotIdx(array, beginIdx, endIdx);
        pivot = temp[pivotIdx];
        std::swap(temp[pivotIdx], temp[endIdx]);
        for (i = endIdx - 1, j = endIdx - 1; j >= beginIdx && j <= endIdx; --j)
        {
            if (temp[j] > pivot)
            {
                std::swap(temp[j], temp[i]);
                --i;
            }
        }
        ++i;
        std::swap(temp[i], temp[endIdx]);
        if (k < i)
        {
            endIdx = i - 1;
        }
        else
        {
            beginIdx = i + 1;
        }
    }

    delete[] temp;

    return pivot;
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
