/*
 * 3_2. Сортировка почти упорядоченной последовательности.
 * Дана последовательность целых чисел a1...an и натуральное число k,
 * такое что для любых i, j: если j >= i + k, то a[i] <= a[j].
 * Требуется отсортировать последовательность.
 * Последовательность может быть очень длинной.
 * Время работы O(n * log(k)). Память O(k).
 * Использовать слияние.
 */

#include <assert.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>

/* Вычислить ближайшую большую чем значение степень двойки. */
std::size_t calcNextPow2(std::size_t value) {
    if (!value) {
        return 1;
    }
    std::size_t result  = value - 1;
    // за основу взято http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
    // поскольку размер size_t не известен, то пришлось сделать цикл,
    // но есть надежда, что компилятор развернет цикл в оригинальную последовательность.
    for (std::size_t shift = 1; shift < sizeof(size_t); shift *= 2) {
        result |= result >> shift;
    }
    return ++result;
}

/*
 * Выполнить слияние двух рядом лежащих массивов 
 * \param array два массива, лежащие рядом: [0; leftSize - 1], [leftSize; rightSize - 1]
 * \param tempArray вспомогательный буфер
 * \param leftSize правая граница левого массива
 * \param rightSize правая граница правого массива
 */
template <typename T> void mergeSortIteration(T* array, T* tempArray, std::size_t leftSize, std::size_t rightSize)
{
    std::size_t leftIdx = 0;
    std::size_t rightIdx = 0;
    std::size_t idx = 0;
    T* rightArray = array + leftSize;
    while (leftIdx < leftSize && rightIdx < rightSize)
    {
        if (array[leftIdx] < rightArray[rightIdx])
        {
            tempArray[idx++] = std::move(array[leftIdx++]);
        }
        else
        {
            tempArray[idx++] = std::move(rightArray[rightIdx++]);
        }
    }
    while (leftIdx < leftSize)
    {
        tempArray[idx++] = std::move(array[leftIdx++]);
    }
//    while (rightIdx < rightSize)
//    {
//        tempArray[idx++] = std::move(rightArray[rightIdx++]);
//    }
    for (std::size_t i = 0; i < leftSize + rightIdx; i++)
    {
        array[i] = std::move(tempArray[i]);
    }
}

/*
 * Произвести сортировку частично упорядоченного массива
 * \param array сортируемый массив
 * \param size размер массива
 * \param k коэффициент в правиле частичного порядка j >= i + k, то array[i] <= array[j].
 */
template <typename T> void mergeSortPartialOrder(T* array, std::size_t size, std::size_t k)
{
    // Для простоты совмещения границ удобно перейти от k к ближайшей степени двойки при вычислении границ интервалов
    // Пространственная сложность остается O(k)
    std::size_t kPow2 = calcNextPow2(k); 
    std::size_t fullSortLimit = kPow2 > size ? size : kPow2;
    std::size_t limit = 4 * kPow2 > size ? size : 4 * kPow2;
    T* temp = new T[limit];

    // Сначала - независимый mergeSort для каждого из подмассивов длины fullSortLimit ~ k
    std::size_t frameSize = 1;
    for (; frameSize <= fullSortLimit; frameSize *= 2)
    {
        std::size_t leftOffset = 0;
        for (; size >= 2 * frameSize + leftOffset; leftOffset += 2 * frameSize)
        {
            mergeSortIteration(array + leftOffset, temp, frameSize, frameSize);
        }
        if (size > frameSize + leftOffset)
        {
            mergeSortIteration(array + leftOffset, temp, frameSize, size - frameSize - leftOffset);
        }
    }

    // затем - последовательный mergeSort двух соседних подмассивов длины fullSortLimit
    // так, что левый подмассив уже участвовал в сортировке ранее
    std::size_t leftOffset = 0;
    for (; size >= 2 * frameSize + leftOffset; leftOffset += frameSize)
    {
        mergeSortIteration(array + leftOffset, temp, frameSize, frameSize);
    }
    if (size > leftOffset + frameSize)
    {
        mergeSortIteration(array + leftOffset, temp, frameSize, size - frameSize - leftOffset);
    }

    delete[] temp;
}

int main(void) {
    using std::cin;

    std::size_t k, count;
    cin >> count;
    cin >> k;
    int* array = new int[count];
    for (std::size_t i = 0; i < count; i++) {
        cin >> array[i];
    }
    if (count && k)
    {
        mergeSortPartialOrder(array, count, k);
    }
    for (std::size_t i = 0; i < count; i++) {
        std::cout << array[i] << ' ';
    }

    delete[] array;

    return 0;
}
