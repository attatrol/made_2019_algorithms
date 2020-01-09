/*
 * 6_4. Порядок обхода.
 *
 * Дано число N < 10^6 и последовательность целых чисел из [-2^31..2^31] длиной N.
 * Требуется построить бинарное дерево, заданное наивным порядком вставки.
 * Т.е., при добавлении очередного числа K в дерево с корнем root, если root→Key ≤ K, то узел K добавляется в правое поддерево root; иначе в левое поддерево root.
 * Рекурсия запрещена.
 *
 * 4. Выведите элементы в порядке level-order (по слоям, “в ширину”).
 */


#include <iostream>
#include <stack>
#include <queue>

/*
 * Узел дерева
 */
template <typename T>
struct Node
{
    Node* right; /* Правое поддерево */
    Node* left;  /* Левое поддерево */
    T value;     /* Значение */
    Node(T&& value): right(nullptr), left(nullptr), value(value)
    {
    }
};

/*
 * Бинарное дерево поиска с наивной вставкой
 */
template <typename T>
class NaiveTree {
private:
    Node<T>* m_head;
public:
    NaiveTree(): m_head(nullptr)
    {
    }
    ~NaiveTree();
    /* Вставить элемент */
    void insert(T&& value);
    /*
     * Обойти дерево послойно
     * \param visitor колбэк, вызываемый для значения узла
     */
    void traverseLevelOrder(void (*const visitor)(const T& value));
};

template <typename T>
NaiveTree<T>::~NaiveTree()
{
    if (!m_head)
    {
        return;
    }
    std::stack<Node<T>*> stack;
    stack.push(m_head);
    // post-order
    while (!stack.empty())
    {
        Node<T>* top = stack.top();
        Node<T>* leftChild = top->left;
        if (leftChild)
        {
            top->left = nullptr;
            stack.push(leftChild);
        }
        else
        {
            Node<T>* rightChild = stack.top()->right;
            if (rightChild)
            {
                top->right = nullptr;
                stack.push(rightChild);
            }
            else
            {
                delete top;
                stack.pop();
            }
        }
    }
}
/*
 * Вставить элемент
 */
template <typename T>
void NaiveTree<T>::insert(T&& value)
{
    if (!m_head)
    {
        m_head = new Node<T>(std::forward<T>(value));
        return;
    }
    Node<T>* currentNode = m_head;
    while (true)
    {
        if (value <= currentNode->value)
        {
            if (currentNode->left == nullptr)
            {
                currentNode->left = new Node<T>(std::forward<T>(value));
                return;
            }
            else
            {
                currentNode = currentNode->left;
            }
        }
        else
        {
            if (currentNode->right == nullptr)
            {
                currentNode->right = new Node<T>(std::forward<T>(value));
                return;
            }
            else
            {
                currentNode = currentNode->right;
            }
        }
    }
}

template <typename T>
void NaiveTree<T>::traverseLevelOrder(void (*const visitor)(const T& value))
{
    if (m_head == nullptr)
    {
        return;
    }
    std::queue<Node<T>*> queue;
    queue.push(m_head);
    do
    {
        Node<T>* current = queue.front();
        visitor(current->value);
        if (current->left)
        {
            queue.push(current->left);
        }
        if (current->right)
        {
            queue.push(current->right);
        }
        queue.pop();
    }
    while(!queue.empty());
}

/*
 * Вывести целое число в стандартный поток вывода
 */
void printOutValue(const int& value)
{
    std::cout << value << ' ';
}

int main(void) {
    using std::cin;


    std::size_t count;
    cin >> count;

    NaiveTree<int> tree;
    for (std::size_t i = 0; i < count; i++) {
        int value;
        cin >> value;
        tree.insert(std::move(value));
    }

    tree.traverseLevelOrder(printOutValue);

    return 0;
}
