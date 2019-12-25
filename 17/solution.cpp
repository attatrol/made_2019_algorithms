/* Задача 17. Шаблон с ?
 *
 * Шаблон поиска задан строкой длины m, в которой кроме обычных символов могут встречаться символы “?”.
 * Найти позиции всех вхождений шаблона в тексте длины n.
 * Каждое вхождение шаблона предполагает, что все обычные символы совпадают с соответствующими из текста,
 * а вместо символа “?” в тексте встречается произвольный символ.
 *
 * Время работы - O(n + m + Z), где Z - общее -число вхождений подстрок шаблона “между вопросиками” в исходном тексте. (Ахо-Корасик)
 * m ≤ 5000, n ≤ 2000000. Время: 10с, память 32Mb.
 */

#define NDEBUG

#include <algorithm>
#include <assert.h>
#ifdef DEBUG
#include <fstream>
#endif
#include <iostream>
#include <stack>
#include <unordered_map>
#include <queue>
#include <vector>

using index_t = int;

/* Trie that supports Aho-Corasick pattern matching search */
class Trie
{
private:
    /* Trie node */
    struct TrieNode
    {
        static const index_t NOT_CALCULATED = -1;

        char value_;                                  // value
        index_t parent_;                              // index of the parent node
        std::unordered_map<char, index_t>* children_; // map value -> child index
        index_t suffixLink_;                          // link to the node which path has the largest common suffix with this node path
        index_t entrySuffixLink_;                     // link in the sequence of suffix links that leads to the first node that has en entry
        index_t entry_;                               // entry, index of a pattern that ends with this node

        TrieNode():
            children_(nullptr), suffixLink_(NOT_CALCULATED), entrySuffixLink_(NOT_CALCULATED), entry_(NOT_CALCULATED)
        {
        }
        TrieNode(char value, index_t parent):
            value_(value), parent_(parent), children_(nullptr), suffixLink_(NOT_CALCULATED), entrySuffixLink_(NOT_CALCULATED), entry_(NOT_CALCULATED)
        {
        }
        TrieNode(const TrieNode&) = delete;
        TrieNode& operator=(const TrieNode&) = delete;
        TrieNode(TrieNode&& other) noexcept :
            value_(other.value_), parent_(other.parent_), children_(nullptr), suffixLink_(other.suffixLink_), entrySuffixLink_(other.entrySuffixLink_), entry_(other.entry_)
        {
            std::swap(other.children_, children_);
        }
        TrieNode& operator=(TrieNode&& other) noexcept
        {
            if (&other != this)
            {
                value_ = other.value_;
                parent_ = other.parent_;
                std::swap(other.children_, children_);
                suffixLink_ = other.suffixLink_;
                entrySuffixLink_ = other.entrySuffixLink_;
                entry_ = other.entry_;
            }
            return *this;
        }
        ~TrieNode()
        {
            delete children_;
        }
    };

    std::vector<TrieNode> nodes_;                  // tree nodes
    index_t entryCount_;
    index_t currentNode_;

    /* Get index of a child node by its value and the parent node index */
    index_t getChildIndex(index_t nodeIdx, char childValue)
    {
        assert(nodeIdx >= 0 && nodeIdx < nodes_.size());
        TrieNode& node = nodes_[static_cast<std::size_t>(nodeIdx)];
        if (!node.children_)
        {
            return -1l;
        }
        auto found = node.children_->find(childValue);
        return found == node.children_->end() ? -1l : found->second;
    }
    /* Add child node */
    index_t addChild(index_t parentNodeIdx, char childValue)
    {
        assert(parentNodeIdx >= 0 && parentNodeIdx < nodes_.size());
        assert(getChildIndex(parentNodeIdx, childValue) < 0);
        TrieNode& parentNode = nodes_[static_cast<std::size_t>(parentNodeIdx)];
        if (!parentNode.children_)
        {
            parentNode.children_ = new std::unordered_map<char, index_t>();
        }
        index_t result = static_cast<index_t>(nodes_.size());
        parentNode.children_->insert({childValue, nodes_.size()});
        nodes_.emplace_back(childValue, parentNodeIdx);
        return result;
    }
    /* Get the next trie node that supports new incoming value */
    index_t getNext(index_t nodeIdx, char value)
    {
        index_t childIdx = getChildIndex(nodeIdx, value);
        if (childIdx != TrieNode::NOT_CALCULATED)
        {
            return childIdx;
        }
        if (!nodeIdx)
        {
            return 0;
        }
        return getNext(getSuffixLink(nodeIdx), value);
    }
    /* Get node suffix link */
    index_t getSuffixLink(index_t nodeIdx)
    {
        assert(nodeIdx >= 0 && nodeIdx < nodes_.size());
        TrieNode& node = nodes_[static_cast<std::size_t>(nodeIdx)];
        if (node.suffixLink_ != TrieNode::NOT_CALCULATED)
        {
            return node.suffixLink_;
        }
        else
        {
            return node.suffixLink_ = node.parent_ ? getNext(getSuffixLink(node.parent_), node.value_) : 0;
        }
    }
    /* Get node entry suffix link */
    index_t getEntrySuffixLink(index_t nodeIdx)
    {
        assert(nodeIdx >= 0 && nodeIdx < nodes_.size());
        TrieNode& node = nodes_[static_cast<std::size_t>(nodeIdx)];
        if (node.entrySuffixLink_ != TrieNode::NOT_CALCULATED)
        {
            return node.entrySuffixLink_;
        }
        else
        {
            index_t suffixLink = getSuffixLink(nodeIdx);
            return node.entrySuffixLink_ = nodes_[static_cast<std::size_t>(suffixLink)].entry_ >= 0 ? suffixLink : getEntrySuffixLink(suffixLink);
        }
    }
public:
    explicit Trie() :
        entryCount_(0), currentNode_(0)
    {
        // init root node
        nodes_.emplace_back();
        nodes_[0].suffixLink_ = 0;
        nodes_[0].entrySuffixLink_ = 0;
    }
    /* Add new pattern to the trie, pattern is a substring */
    index_t addPattern(const std::string& pattern, std::size_t beginIdx, std::size_t endIdx)
    {
        index_t currentNode = 0;
        for (std::size_t i = beginIdx; i < endIdx; ++i)
        {
            char ch = pattern[i];
            index_t child = getChildIndex(currentNode, ch);
            currentNode = child >= 0 ? child : addChild(currentNode, ch);
        }
        if (nodes_[static_cast<std::size_t>(currentNode)].entry_ == TrieNode::NOT_CALCULATED)
        {
            nodes_[static_cast<std::size_t>(currentNode)].entry_ = entryCount_++;
        }
#ifdef DEBUG
        std::cout << "Added trie pattern " << nodes_[static_cast<std::size_t>(currentNode)].entry_ << '\n';
        for (std::size_t i = beginIdx; i < endIdx; ++i)
        {
            std::cout << pattern[i];
        }
        std::cout << std::endl;
#endif
        return nodes_[static_cast<std::size_t>(currentNode)].entry_;
    }
    /* Add new pattern to the trie, pattern is a whole string */
    index_t addPattern(const std::string& pattern)
    {
        return addPattern(pattern, 0, pattern.size());
    }
    /* Reset the trie before new pattern matching */
    void init()
    {
        currentNode_ = 0;
    }
    /**
      * Push next value to the heap
      * \param value next char from the text
      * \param entries bool array where encountered entries are registered
      * \param hasEntries flag, is true if any entry was encountered
      */
    void next(char value, bool* entries, bool& hasEntries)
    {
        hasEntries = false;
        index_t childIdx;
        assert(currentNode_ >= 0 && currentNode_ < nodes_.size());
#ifdef DEBUG
            std::cout << "previous node " << currentNode_ << " (" << nodes_[static_cast<std::size_t>(currentNode_)].value_ << ")\n";
#endif
        while (-1 == (childIdx = getChildIndex(currentNode_, value)) && currentNode_)
        {
            currentNode_ = getSuffixLink(currentNode_);
#ifdef DEBUG
            std::cout << "suffix link move to " << currentNode_ << " (" << nodes_[static_cast<std::size_t>(currentNode_)].value_ << ')';
#endif
            if (nodes_[static_cast<std::size_t>(currentNode_)].entry_ != TrieNode::NOT_CALCULATED)
            {
#ifdef DEBUG
                std::cout << " entry found " << nodes_[static_cast<std::size_t>(currentNode_)].entry_;
#endif
                //assert(!entries[nodes_[static_cast<std::size_t>(currentNode_)].entry_]);
                entries[nodes_[static_cast<std::size_t>(currentNode_)].entry_] = hasEntries = true;
            }
#ifdef DEBUG
            std::cout << std::endl;
#endif
        }
        // child was found
        assert(childIdx >= 0 || !currentNode_);
        if (childIdx >= 0)
        {
            currentNode_ = childIdx;
        }
#ifdef DEBUG
            std::cout << "child move " << currentNode_ << " (" << nodes_[static_cast<std::size_t>(currentNode_)].value_ << ')';
#endif
        if (currentNode_)
        {
            if (nodes_[static_cast<std::size_t>(currentNode_)].entry_ != TrieNode::NOT_CALCULATED)
            {
#ifdef DEBUG
                std::cout << " entry found " << nodes_[static_cast<std::size_t>(currentNode_)].entry_;
#endif
                //assert(!entries[nodes_[static_cast<std::size_t>(currentNode_)].entry_]);
                entries[nodes_[static_cast<std::size_t>(currentNode_)].entry_] = hasEntries = true;
            }
            while ((childIdx = getEntrySuffixLink(childIdx)))
            {
                if (nodes_[static_cast<std::size_t>(childIdx)].entry_ != TrieNode::NOT_CALCULATED)
                {
    #ifdef DEBUG
                    std::cout << " entry found " << nodes_[static_cast<std::size_t>(childIdx)].entry_;
    #endif
                    //assert(!entries[nodes_[static_cast<std::size_t>(childIdx)].entry_]);
                    entries[nodes_[static_cast<std::size_t>(childIdx)].entry_] = hasEntries = true;
                }
            }
        }
#ifdef DEBUG
        else
        {
            std::cout << "root encountered";
        }
        std::cout << std::endl;
#endif

    }
    const TrieNode* getRoot() const
    {
        return &nodes_[0];
    }
};

/* Token that describes pattern registered in the trie as a part of the wildcard pattern */
struct SearcherEntry
{
    std::size_t offset_; // number of characters from end of the previous pattern to the last character of this pattern
    std::size_t length_; // total number of characters in this pattern
    index_t entry_;      // index of trie entry

    SearcherEntry()
    {
    }
    SearcherEntry(std::size_t offset, std::size_t length, index_t entry) :
        offset_(offset), length_(length), entry_(entry)
    {
    }
};

/* Token used by priority queue of the pending patterns */
struct TimerToken
{
    std::size_t stringIndex_;   // the text index where the pattern expected to be
    index_t entry_;             // pattern index
    std::size_t sequenceIndex_; // wildcard sequence token index

    TimerToken()
    {
    }
    TimerToken(std::size_t stringIndex, index_t entry, std::size_t sequenceIndex) :
        stringIndex_(stringIndex), entry_(entry), sequenceIndex_(sequenceIndex)
    {
    }

    friend bool operator<(const TimerToken& lhs, const TimerToken& rhs)
    {
        return lhs.stringIndex_ > rhs.stringIndex_;
    }
};

/* Searches a text for the pattern with wildcards */
class WildcardSearcher
{
public:
    static const char WILDCARD_SYMBOL = '?'; // wildcard character
private:
    std::size_t size_;                    // whole pattern size
    Trie trie_;                           // trie that stores partial patterns
    std::size_t wildcardPrefix_;          // number of wildcard characters before any definite characters
    std::size_t wildcardPostfix_;         // number of wildcard characters after any definite characters
    std::vector<SearcherEntry> sequence_; // sequence of tokens to trie patterns
public:
    /* Default ctor, splits wildcard pattern into definite character patters, whish are stored in the trie */
    explicit WildcardSearcher(const std::string& pattern):
        size_(pattern.size())
    {
        std::size_t i, offsetStart, patternStart;
        bool readingPattern = pattern[0] != WILDCARD_SYMBOL;
        offsetStart = 0;
        if (readingPattern)
        {
            patternStart = 0;
        }
        for (i = 0; i < size_; ++i)
        {
            bool isPatternChar = pattern[i] != WILDCARD_SYMBOL;
            if (isPatternChar && !readingPattern)
            {
                patternStart = i;
                readingPattern = true;
            }
            else if (!isPatternChar && readingPattern)
            {
                long trieEntryIndex = trie_.addPattern(pattern, patternStart, i);
                sequence_.emplace_back(i - offsetStart, i - patternStart, trieEntryIndex);
                offsetStart = i;
                readingPattern = false;
            }
        }
        if (readingPattern)
        {
            long trieEntryIndex = trie_.addPattern(pattern, patternStart, size_);
            sequence_.emplace_back(size_ - offsetStart, size_ - patternStart, trieEntryIndex);
            wildcardPrefix_ = sequence_[0].offset_ - sequence_[0].length_;
            wildcardPostfix_ = 0;
        }
        else if (sequence_.size())
        {
            wildcardPrefix_ =  sequence_[0].offset_ - sequence_[0].length_;
            wildcardPostfix_ = size_ - offsetStart;
        }
        else
        {
            wildcardPrefix_ =  size_;
            wildcardPostfix_ = 0;
        }
#ifdef DEBUG
        std::cout << "Sequence:\n";
        std::cout << "wildcardPrefix = " << wildcardPrefix_ << std::endl;
        for (const auto& token : sequence_)
        {
            std::cout << "length = " << token.length_ << " offset = " << token.offset_ << " entry = " << token.entry_ << '\n';
        }
        std::cout << "wildcardPostfix = " << wildcardPostfix_ << '\n' << std::endl;
#endif
    }

    /* Search a text for the wildcard pattern */
    std::vector<std::size_t>search(const std::string& text)
    {
        std::vector<std::size_t> result;
        if (text.size() < size_)
        {
            return result; // text is smaller than pattern
        }
        else if (!sequence_.size())
        {
            for (std::size_t i = 0; i < text.size() - wildcardPrefix_ + 1; ++i)
            {
                result.push_back(i);
            }
            return result; // pattern contains only wildcards
        }

        trie_.init();
        auto initialEntry = sequence_[0];
        bool hasEntries = false;
        bool* entries = new bool[sequence_.size()];
        std::fill(entries, entries + sequence_.size(), false);

        if (sequence_.size() == 1)
        {
            for (std::size_t i = wildcardPrefix_; i < text.size() - wildcardPostfix_; ++i)
            {
                trie_.next(text[i], entries, hasEntries);
                if (hasEntries)
                {
                    result.push_back(i + 1 - initialEntry.length_ - wildcardPrefix_);
                    entries[0] = false;
                }
            }

            delete[] entries;

            return result; // there is a single pattern in the trie
        }

        std::priority_queue<TimerToken> timers;
        for (std::size_t i = wildcardPrefix_; i < text.size() - wildcardPostfix_; ++i)
        {
#ifdef DEBUG
            std::cout << "////// " << i << " //////\n";
#endif
            std::stack<std::size_t> newTokens;
            while (!timers.empty() && timers.top().stringIndex_ == i)
            {
                if (entries[timers.top().entry_])
                {
                    std::size_t sequenceIndex = timers.top().sequenceIndex_;
                    if (sequenceIndex == sequence_.size() - 1)
                    {
                        assert(i >= size_ - wildcardPostfix_);
                        result.push_back(i + wildcardPostfix_ - size_);
#ifdef DEBUG
                        std::cout << "full sequence" << std::endl;
#endif
                    }
                    else
                    {
                        ++sequenceIndex;
                        newTokens.emplace(sequenceIndex);
#ifdef DEBUG
                        std::cout << "finished entry = " << timers.top().entry_ << " sequenceIndex = " << (sequenceIndex - 1) << std::endl;
#endif
                    }
                }
#ifdef DEBUG
                else
                {
                    std::cout << "rejected entry " << timers.top().entry_ << " sequenceIndex = " << timers.top().sequenceIndex_ << std::endl;
                }
#endif
                timers.pop();
            }
            while(!newTokens.empty())
            {
               timers.emplace(i + sequence_[newTokens.top()].offset_, sequence_[newTokens.top()].entry_, newTokens.top());
#ifdef DEBUG
               std::cout << "sequence registered index = " << i + sequence_[newTokens.top()].offset_ <<
                            " entry = " << sequence_[newTokens.top()].entry_ << " sequenceIndex = " << newTokens.top() << std::endl;
#endif
               newTokens.pop();

            }
            if (hasEntries)
            {
                std::fill(entries, entries + sequence_.size(), false);
            }
            trie_.next(text[i], entries, hasEntries);
#ifdef DEBUG
            if (hasEntries)
            {
                for (std::size_t i = 0; i < sequence_.size(); ++i)
                {
                    if (entries[i])
                    {
                        std::cout << "got entry " << i << '\n';
                    }
                }
            }
#endif
            if (entries[0])
            {
                timers.emplace(i + sequence_[1].offset_ + 1, sequence_[1].entry_, 1);
#ifdef DEBUG
                std::cout << "Start found, sequence registered index = " << i + sequence_[1].offset_ + 1 <<
                             " entry = " << sequence_[1].entry_ << " sequenceIndex = " << 1 << std::endl;
#endif
            }
        }

        if (entries[timers.top().entry_])
        {
            while (!timers.empty() && timers.top().stringIndex_ == text.size() - wildcardPostfix_)
            {
                if (timers.top().sequenceIndex_ == sequence_.size() - 1)
                {
                    result.push_back(text.size() - wildcardPostfix_ + wildcardPostfix_ - size_);
                    break;
                }
                timers.pop();
            }
        }

        delete[] entries;

        return result;
    }
};

int main(void)
{
    std::string pattern, text;
#ifdef DEBUG
    std::ifstream in("input.txt");
    std::getline(in, pattern);
    std::getline(in, text);
#else
    std::iostream::sync_with_stdio(false);
    std::cin >> pattern;
    std::cin >> text;
#endif

    WildcardSearcher searcher(pattern);
    auto result = searcher.search(text);

    for (auto index : result)
    {
        std::cout << index << ' ';
    }

    return 0;
}
