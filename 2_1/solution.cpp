/*
 * 2_1. Жадина.
 *
 * Вовочка ест фрукты из бабушкиной корзины.
 * В корзине лежат фрукты разной массы.
 * Вовочка может поднять не более K грамм.
 * Каждый фрукт весит не более K грамм.
 * За раз он выбирает несколько самых тяжелых фруктов, которые может поднять одновременно,
 * откусывает от каждого половину и кладет огрызки обратно в корзину.
 * Если фрукт весит нечетное число грамм, он откусывает большую половину.
 * Фрукт массы 1гр он съедает полностью.
 *
 * Определить за сколько подходов Вовочка съест все фрукты в корзине.
 *
 * Формат входных данных.
 * Вначале вводится n - количество фруктов и n строк с массами фруктов.
 * n ≤ 50000.
 * Затем K - "грузоподъемность". K ≤ 1000.
*/

#include <assert.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>


/* Начальная емкость кучи */
constexpr unsigned int HEAP_DEFAULT_CAPACITY = 16;

/* Бинарная куча, отсортированная по убыванию */
template <typename T> class MaxBinaryHeap {
    private:
        T* m_data;              /* Бинарное дерево кучи в виде массива*/
        std::size_t m_capacity; /* Емкость дерева, размер m_data */
        std::size_t m_count;    /* Количество элементов в куче */
        /**
         * Увеличить емкость до ближайшей большей степени двойки.
         * Так гарантируется амортизированное время добавления О(1).
         */
        void grow();
        /** Восстановить порядок узла, проталкивать его в направлении вершины дерева */
        void siftUp(std::size_t idx);
        /** Восстановить порядок узла, проталкивать его в направлении листьев дерева */
        void siftDown(std::size_t idx);
        /**
         * Восстановить порядок кучи на m_data.
         * Текущий порядок на m_data не определен.
         */
        void rebuildHeap();
    protected:
        /**
         * Получить значение, возвращаемое из пустой кучи.
         * Можно выбросить здесь исключение.
         */
        T getEmptyValue();
        /**
         * Получить ссылку на значение, возвращаемое из пустой кучи.
         * Можно выбросить здесь исключение.
         */
        const T& peekEmptyValue();
        /**
         * Признак заполненного дерева.
         * Дальнейшее добавление элемента приведет к увеличению размера массива данных.
         */
        bool isFull() { return m_count == m_capacity; };
    public:
        /** Ктор по умолчанию */
        MaxBinaryHeap(std::size_t capacity = HEAP_DEFAULT_CAPACITY);
        /** Ктор на основе внешних данных. */
        MaxBinaryHeap(const T* data, std::size_t capacity, std::size_t count);
        MaxBinaryHeap(const T* data, std::size_t capacity) : MaxBinaryHeap(data, capacity, capacity) {};
        MaxBinaryHeap(const MaxBinaryHeap&) = delete;
        MaxBinaryHeap(MaxBinaryHeap&&) = delete;
        MaxBinaryHeap&operator=(const MaxBinaryHeap&) = delete;
        MaxBinaryHeap&operator=(MaxBinaryHeap&&) = delete;
        /** Деструктор */
        ~MaxBinaryHeap() { delete[] m_data; };
        /** Признак пустой кучи */
        bool isEmpty() { return m_count == 0; };
        /** Количество элементов в куче */
        std::size_t getCount() { return m_count; };
        /** Добавить элемент в кучу */
        void add(T&& value);
        /** Извлечь элемент из кучи */
        T extractMax();
        /** Посмотреть элемент на вершине кучи */
        const T& peekMax() { return isEmpty() ? peekEmptyValue() : m_data[0]; };
};

/**
 * Вычислить ближайшую большую чем значение степень двойки. */
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
    auto res = ++result == value ? result << 1 : result;
    return res;
}

template <typename T> MaxBinaryHeap<T>::MaxBinaryHeap(std::size_t capacity): m_data(new T[capacity]), m_capacity(capacity), m_count(0) {
}

template <typename T> MaxBinaryHeap<T>::MaxBinaryHeap(const T* data, std::size_t capacity, std::size_t count): m_data(new T[capacity]), m_capacity(capacity), m_count(count) {
    for (std::size_t i = 0; i < m_count; i++) {
        m_data[i] = data[i];
    }
    rebuildHeap();
    //test();
}

template <typename T> void MaxBinaryHeap<T>::grow() {
    m_capacity = calcNextPow2(m_capacity);
    T* newData = new T[m_capacity];
    for (std::size_t i = 0; i < m_count; i++) {
        newData[i] = std::move(m_data[i]);
    }
    delete[] m_data;
    m_data = newData;
}

template <typename T> void MaxBinaryHeap<T>::siftUp(std::size_t idx) {
    std::size_t parentIdx = (idx - 1) / 2;
    while (idx > parentIdx && m_data[idx] > m_data[parentIdx]) {
        std::swap(m_data[idx], m_data[parentIdx]);
        idx = parentIdx;
        parentIdx = (parentIdx - 1) / 2;
    }
}

template <typename T> void MaxBinaryHeap<T>::siftDown(std::size_t idx) {
    std::size_t minChildIdx = 2 * idx + 2;
    while (minChildIdx < m_count) {
        if (m_data[minChildIdx] <= m_data[minChildIdx - 1]) {
            minChildIdx--;
        }
        if (m_data[minChildIdx] > m_data[idx]) {
            std::swap(m_data[idx], m_data[minChildIdx]);
            idx = minChildIdx;
        }
        minChildIdx = 2 * minChildIdx + 2;
    }
    if (minChildIdx == m_count && m_data[--minChildIdx] > m_data[idx]) {
        std::swap(m_data[idx], m_data[minChildIdx]);
    }
}

template <typename T> void MaxBinaryHeap<T>::rebuildHeap() {
    for (std::size_t i = m_count / 2; i > 0; --i) {
        siftDown(i - 1);
    }
}

template <typename T> void MaxBinaryHeap<T>::add(T&& value) {
    if (isFull()) {
        grow();
    }
    m_data[m_count] = value;
    siftUp(m_count);
    m_count++;
}

template <typename T> T MaxBinaryHeap<T>::extractMax() {
    if (isEmpty()) {
        return getEmptyValue();
    }
    T result = m_data[0];

    m_count--;
    if (m_count) {
        m_data[0] = std::move(m_data[m_count]);
        siftDown(0);
    }
    return result;
}

template <> int MaxBinaryHeap<int>::getEmptyValue() {
    throw std::runtime_error("Failed to read value from the empty heap.");
}

template <> const int& MaxBinaryHeap<int>::peekEmptyValue() {
    throw std::runtime_error("Failed to read value from the empty heap.");
}

int countFruitEatingSessions(MaxBinaryHeap<int>& heap, const int maxSum) {
    int count = 0;
    int* usedValues = new int[maxSum];
    while (!heap.isEmpty()) {
        int currentSum = 0;
        int usedValuesCount = 0;
        while (!heap.isEmpty() && maxSum - currentSum >= heap.peekMax()) {
            int value = heap.extractMax();
            currentSum += value;
            if (value > 1) {
                usedValues[usedValuesCount++] = value / 2;
            }
        }
        while(usedValuesCount) {
            heap.add(std::move(usedValues[--usedValuesCount]));
        }
        count++;
    }
    delete[] usedValues;
    return count;
}

int main(void) {
    using std::cin;

    int fruitCount;
    cin >> fruitCount;
    int* fruits = new int[fruitCount];
    for (int i = 0; i < fruitCount; i++) {
        cin >> fruits[i];
    }
    int maxWeight;
    cin >> maxWeight;

    MaxBinaryHeap<int>heap(fruits, fruitCount);

    std::cout << countFruitEatingSessions(heap, maxWeight);

    delete[] fruits;

    return 0;
}
