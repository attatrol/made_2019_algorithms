/*
 * 8_1. Хеш-таблица
 * Реализуйте структуру данных типа “множество строк” на основе динамической хеш-таблицы с открытой адресацией.
 * Хранимые строки непустые и состоят из строчных латинских букв.
 * Хеш-функция строки должна быть реализована с помощью вычисления значения многочлена методом Горнера.
 * Начальный размер таблицы должен быть равным 8-ми. Перехеширование выполняйте при добавлении элементов в случае, когда коэффициент заполнения таблицы достигает 3/4.
 * Структура данных должна поддерживать операции добавления строки в множество, удаления строки из множества и проверки принадлежности данной строки множеству.
 *
 * При добавлении элемента в множество НЕ ГАРАНТИРУЕТСЯ, что он отсутствует в этом множестве.
 * При удалении элемента из множества НЕ ГАРАНТИРУЕТСЯ, что он присутствует в этом множестве.
 *
 * 1. Для разрешения коллизий используйте квадратичное пробирование.
 * i-ая проба
 * g(k, i)=g(k, i-1) + i (mod m). m - степень двойки.
 */

#include <assert.h>
#include <iostream>
#include <limits>
#include <stack>
#include <vector>

#include <fstream>

/* value -> 2^ceil(log2(value)) */
std::size_t ceilPow2(std::size_t value)
{
    std::size_t result = 1;
    while(result - value > value)
    {
        result <<= 1;
    }
    return result;
}

/* Hash set based on hash table */
template<typename T> class HashSet
{
    using hashFn_t = std::size_t (* const)(const std::string& value, const std::size_t mask); /* Hash function type */
private:
    static const std::size_t NON_EXISTING_INDEX = std::numeric_limits<std::size_t>::max(); /* Result of a failed index search. Greater than any power of 2. */
    /* Describes state of the hash table cell */
    enum state_t {
        virgin,   /* empty and never used cell */
        empty,    /* empty used cell */
        occupied, /* cell with value */
    };

    const hashFn_t m_hashFn;       /* hash function */
    const float m_rehashRatio;     /* determines max share of occupied cells before rehashing */
    state_t* m_state;              /* array of the table cell states */
    T** m_table;                   /* the hash table */
    std::size_t m_size;            /* the hash table size, always a power of 2 */
    std::size_t m_count;           /* number of occupied cells */
    std::size_t m_sizeMask;        /* m_size - 1, used to trunkate indices */
    std::size_t m_resizeThreshold; /* threshold used to determine when to rehash */

    /* Expand size and rearrange elements to the new cells */
    void rehash();
    /* Find current position of the value in m_table */
    std::size_t findIndex(const T& value) const;
    /* Find position where the value may be inserted */
    std::size_t findFreeIndex(const T& value) const;
    /* Unsafe private ctor */
    enum unchecked_t { unchecked };
    HashSet(hashFn_t hashFn, std::size_t initialSize, float rehashRatio, unchecked_t);
public:
    /* Ctor */
    HashSet(hashFn_t hashFn, std::size_t initialSize = 8, float rehashRatio = 0.75f);
    /* Dtor */
    ~HashSet();
    HashSet(const HashSet&) = delete;
    HashSet(HashSet&&) = delete;
    HashSet& operator=(const HashSet&) = delete;
    HashSet& operator=(HashSet&&) = delete;
    /* Load factor */
    float loadFactor() const
    {
        return static_cast<float>(m_count) / m_size;
    }
    /* Hash set size */
    std::size_t size() const
    {
        return m_size;
    }
    /* Hash set number of stored elements */
    std::size_t count() const
    {
        return m_count;
    }
    /* Test if the value is stored within the set */
    bool exists(const T& value) const
    {
        return NON_EXISTING_INDEX != findIndex(value);
    }
    /*
     * Move element to the set.
     * \return if the element was successfully added.
     */
    bool add(T&& value);
    /*
     * Copy and add element to the set.
     * \return if the element was successfully added.
     */
    bool add(const T& value);
    /*
     * Delete element from the set.
     * \return if the element was found.
     */
    bool remove(const T& value);
};

template<typename T> HashSet<T>::HashSet(hashFn_t hashFn, std::size_t initialSize, float rehashRatio, unchecked_t):
    m_hashFn(hashFn), m_rehashRatio(rehashRatio), m_state(new state_t[initialSize]), m_table(new T*[initialSize]), m_size(initialSize),
    m_count(0), m_sizeMask(initialSize - 1), m_resizeThreshold(initialSize * m_rehashRatio)
{
    std::fill(m_table, m_table + initialSize, nullptr);
    std::fill(m_state, m_state + initialSize, state_t::virgin);
}

template<typename T> HashSet<T>::HashSet(hashFn_t hashFn, std::size_t initialSize, float rehashRatio):
    HashSet(hashFn, initialSize < 8 ? 8 : ceilPow2(initialSize), rehashRatio > 0.75f ? 0.75f : (rehashRatio < 0.5f ? 0.5f : rehashRatio), unchecked)
{
}

template<typename T> HashSet<T>::~HashSet()
{
    for (std::size_t i = 0; i < m_size; ++i)
    {
        if (m_state[i] == state_t::occupied)
        {
            delete m_table[i];
        }
    }
    delete[] m_table;
    delete[] m_state;
}

template<typename T> void HashSet<T>::rehash()
{
    const std::size_t size = m_size;
    m_size *= 2;
    m_sizeMask = m_size - 1;
    m_resizeThreshold = m_rehashRatio * m_size;

    T** table = m_table;
    m_table = new T*[m_size];
    std::fill(m_table, m_table + m_size, nullptr);
    state_t* state = m_state;
    m_state = new state_t[m_size];
    std::fill(m_state, m_state + m_size, state_t::virgin);

    for (std::size_t i = 0; i < size; ++i)
    {
        if (state[i] == state_t::occupied)
        {
            std::size_t index = findFreeIndex(*table[i]);
            assert(NON_EXISTING_INDEX != index);
            m_table[index] = table[i];
            m_state[index] = state_t::occupied;
        }
    }

    delete[] table;
    delete[] state;
}

template<typename T> std::size_t HashSet<T>::findIndex(const T& value) const
{
    std::size_t index = m_hashFn(value, m_sizeMask);
    for (std::size_t i = 0; i < m_size; index = (index + ++i) & m_sizeMask)
    {
        if (m_state[index] == state_t::occupied)
        {
            if (value == *m_table[index])
            {
                return index;
            }
        }
        else if (m_state[index] == state_t::virgin)
        {
            break;
        }
    }
    return NON_EXISTING_INDEX;
}

template<typename T> std::size_t HashSet<T>::findFreeIndex(const T& value) const
{
    std::size_t index = m_hashFn(value, m_sizeMask);
    std::size_t result = NON_EXISTING_INDEX;
    for (std::size_t i = 0; i < m_size; index = (index + ++i) & m_sizeMask)
    {
        if (m_state[index] == state_t::empty)
        {
            if (result == NON_EXISTING_INDEX)
            {
                result = index;
            }
        }
        else if (m_state[index] == state_t::virgin)
        {
            if (result == NON_EXISTING_INDEX)
            {
                result = index;
            }
            break;
        }
        else if (value == *m_table[index])
        {
            return NON_EXISTING_INDEX;
        }
    }
    return result;
}

template<typename T> bool HashSet<T>::add(const T& value)
{
    if (m_count == m_resizeThreshold)
    {
        rehash();
    }
    std::size_t index = findFreeIndex(value);
    if (NON_EXISTING_INDEX == index)
    {
        return false;
    }
    m_table[index] = new T(value);
    m_state[index] = state_t::occupied;
    ++m_count;
    return true;
}

template<typename T> bool HashSet<T>::add(T&& value)
{
    if (m_count == m_resizeThreshold)
    {
        rehash();
    }
    std::size_t index = findFreeIndex(value);
    if (NON_EXISTING_INDEX == index)
    {
        return false;
    }
    m_table[index] = new T(value);
    m_state[index] = state_t::occupied;
    ++m_count;
    return true;
}

template<typename T> bool HashSet<T>::remove(const T& value)
{
    std::size_t index = findIndex(value);
    if (NON_EXISTING_INDEX == index)
    {
        return false;
    }
    assert(m_state[index] == state_t::occupied);
    m_state[index] = state_t::empty;
    delete m_table[index];
    m_table[index] = nullptr;
    --m_count;
    return true;
}

std::size_t hashString(const std::string& value, const std::size_t mask)
{
    static const std::size_t FACTOR = 31415;
    std::size_t result = 0;
    for(const char& c : value) {
        result = (result * FACTOR + static_cast<std::size_t>(c)) & mask;
    }
    return result;
}

int main(void) {
  HashSet<std::string>set(hashString);
  char command = ' ';
  std::string value;
  while (std::cin >> command >> value) {
    switch (command) {
      case '?':
        std::cout << (set.exists(value) ? "OK" : "FAIL") << std::endl;
        break;
      case '+':
        std::cout << (set.add(value) ? "OK" : "FAIL") << std::endl;
        break;
      case '-':
        std::cout << (set.remove(value) ? "OK" : "FAIL") << std::endl;
        break;
    }
  }
  return 0;
}
