/*
* Во всех задачах из следующего списка следует написать структуру данных, обрабатывающую команды push* и pop*.
* Формат входных данных.
* В первой строке количество команд n. n ≤ 1000000.
* Каждая команда задаётся как 2 целых числа: a b.
* a = 1 - push front
* a = 2 - pop front
* a = 3 - push back
* a = 4 - pop back
* Для очереди используются команды 2 и 3. Для дека используются все четыре команды.
* Если дана команда pop*, то число b - ожидаемое значение.Если команда pop вызвана для пустой структуры данных, то ожидается “-1”.
* Формат выходных данных.
* Требуется напечатать YES - если все ожидаемые значения совпали. Иначе, если хотя бы одно ожидание не оправдалось, то напечатать NO.
*
* 1_2. Реализовать дек с динамическим зацикленным буфером.
*/

#include <assert.h>
#include <limits>
#include <iostream>
//#include <fstream>

/* Log2 от начальной емкости дека - 16 */
constexpr unsigned int DEQUE_DEFAULT_CAPACITY = 27;

/* Двухсторонний список */
template <typename T> class Deque {
private:
    T* m_data;              /* Динамический массив с данными */
    std::size_t m_capacity; /* Емкость дека, размер m_data */
    std::size_t m_idxMask;  /* Маска, используется для усечения вычисляемых индексов массива:
                            * Чтобы не использовать дорогое взятие остатка, введено требование -
                            * m_capacity должно всегда быть степенью двойки.
                            * Тогда взятие остатка превращается в побитовое И с маской.
                            */
    std::size_t m_frontIdx; /* Индекс на первый элемент дека */
    std::size_t m_count;    /* Количество элементов в деке */
    /*
         * Увеличить емкость в 2 раза.
         * Так гарантируется амортизированное время добавления О(1)
         */
    void grow();
protected:
    /*
    * Получить значение, возвращаемое из пустого стека.
    * Можно выбросить здесь исключение.
    */
    T&& getEmptyValue();
    /* Признак заполненного стека */
    bool isFull() { return m_count == m_capacity; };
public:
    /* Ктор по умолчанию */
    Deque(std::size_t capacity = DEQUE_DEFAULT_CAPACITY);
    /* Деструктор */
    ~Deque() { delete[] m_data; };
    /* Признак пустого дека */
    bool isEmpty() { return m_count == 0; };
    /* Количество элементов в деке */
    std::size_t getCount() { return m_count; };
    /* Положить в начало дека */
    void pushFront(T&& value);
    /* Положить в конец дека */
    void pushBack(T&& value);
    /* Вытолкнуть из начала дека */
    T&& popFront();
    /* Вытолкнуть из концв дека */
    T&& popBack();
    /* Значение в начале дека */
    const T& front() { return isEmpty() ? getEmptyValue() : m_data[m_frontIdx]; };
    /* Значение в конце дека */
    const T& back() { return isEmpty() ? getEmptyValue() : m_data[(m_frontIdx + m_count) & m_idxMask]; };
};

std::size_t calcNextPow2(std::size_t value) {
    --value;
    for (std::size_t shift = 1; shift < sizeof(size_t); shift *= 2) {
        value |= value >> shift;
    }
    return ++value;
}

template <typename T> Deque<T>::Deque(std::size_t capacity): m_data(nullptr), m_frontIdx(0), m_count(0) {
    assert(capacity > 0);
    m_capacity = calcNextPow2(capacity);
    m_idxMask = m_capacity - 1;
    m_data = new T[m_capacity];
}

template <typename T> void Deque<T>::grow() {
    m_capacity <<= 1;
    m_idxMask = m_capacity - 1;
    T* newData = new T[m_capacity];
    std::size_t destIdx = 0;
    for (std::size_t sourceIdx = m_frontIdx; sourceIdx < m_count; sourceIdx++) {
        newData[destIdx++] = std::move(m_data[sourceIdx]);
    }
    for (std::size_t sourceIdx = 0; sourceIdx < m_frontIdx; sourceIdx++) {
        newData[destIdx++] = std::move(m_data[sourceIdx]);
    }
    assert(destIdx == m_capacity / 2);
    delete[] m_data;
    m_data = newData;
    m_frontIdx = 0;
}

template <typename T> void Deque<T>::pushFront(T&& value) {
    if (isFull()) {
        grow();
    }
    --m_frontIdx;
    m_frontIdx &= m_idxMask;
    m_count++;
    m_data[m_frontIdx] = value;
}

template <typename T> void Deque<T>::pushBack(T&& value) {
    if (isFull()) {
        grow();
    }
    m_data[(m_frontIdx + m_count++) & m_idxMask] = value;
}

template <typename T> T&& Deque<T>::popFront() {
    if (isEmpty()) {
        return getEmptyValue();
    }
    std::size_t oldFrontIdx = m_frontIdx;
    ++m_frontIdx;
    m_frontIdx &= m_idxMask;
    m_count--;
    return std::move(m_data[oldFrontIdx]);
}

template <typename T> T&& Deque<T>::popBack() {
    if (isEmpty()) {
        return getEmptyValue();
    }
    return std::move(m_data[(m_frontIdx + --m_count) & m_idxMask]);
}

template <> int&& Deque<int>::getEmptyValue() {
    return std::move(-1);
}

bool test(Deque<int>& deque, int& opcode, int& value) {
    switch (opcode) {
    case 1:
        deque.pushFront(std::move(value));
        return true;
    case 2: {
        int popped = deque.popFront();
        return popped == value;
    }
    case 3:
        deque.pushBack(std::move(value));
        return true;
    case 4: {
        int popped = deque.popBack();
        return popped == value;
    }
    }
    assert(false);
}

int main(void) {
    using std::cin;
    //std::ifstream cin("input.txt");
    std::size_t lineCount;
    cin >> lineCount;
    Deque<int> deque;
    for (std::size_t i = 0; i < lineCount; i++) {
        int opcode, value;
        cin >> opcode >> value;
        if (!test(deque, opcode, value)) {
            std::cout<<"NO";
            return 0;
        }
    }
    std::cout<<"YES";
    return 0;
}
