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

        char value_;                                 // value
        index_t parent_;                             // index of the parent node
        std::unordered_map<char, index_t> children_; // map value -> child index
        index_t suffixLink_;                         // link to the node which path has the largest common suffix with this node path
        index_t entrySuffixLink_;                    // link in the sequence of suffix links that leads to the first node that has en entry
        index_t entry_;                              // entry, index of a pattern that ends with this node

        TrieNode():
            suffixLink_(NOT_CALCULATED), entrySuffixLink_(NOT_CALCULATED), entry_(NOT_CALCULATED)
        {
        }
        TrieNode(char value, index_t parent):
            value_(value), parent_(parent), suffixLink_(NOT_CALCULATED), entrySuffixLink_(NOT_CALCULATED), entry_(NOT_CALCULATED)
        {
        }
    };

    std::vector<TrieNode> nodes_; // tree nodes
    index_t entryCount_;          // number of registered patterns
    index_t currentNode_;         // index of the current node

    /* Get index of a child node by its value and the parent node index */
    index_t getChildIndex(index_t nodeIdx, char childValue) const;
    /* Add child node */
    index_t addChild(index_t parentNodeIdx, char childValue);
    /* Get the next trie node that supports new incoming value */
    index_t getNext(index_t nodeIdx, char value);
    /* Get node suffix link */
    index_t getSuffixLink(index_t nodeIdx);
    /* Get node entry suffix link */
    index_t getEntrySuffixLink(index_t nodeIdx);
public:
    explicit Trie();
    /* Add new pattern to the trie, pattern is a substring */
    index_t addPattern(const std::string& pattern, std::size_t beginIdx, std::size_t endIdx);
    /* Add new pattern to the trie, pattern is a whole string */
    index_t addPattern(const std::string& pattern);
    /* Reset the trie before new pattern matching */
    void init()
    {
        currentNode_ = 0;
    }
    /**
      * Push next value to the heap
      * \param value next char from the text
      * \param matches bool array where encountered matches are registered
      * \param hasMatches flag, is true if any entry was encountered
      */
    void next(char value, bool* matches, bool& hasMatches);
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
    explicit WildcardSearcher(const std::string& pattern);
    /* Search a text for the wildcard pattern */
    std::vector<std::size_t> search(const std::string& text);

};

///////////////////////////////////////////////////////////////////////

index_t Trie::getChildIndex(index_t nodeIdx, char childValue) const
{
    assert(nodeIdx >= 0 && nodeIdx < nodes_.size());
    const TrieNode& node = nodes_[static_cast<std::size_t>(nodeIdx)];
    auto found = node.children_.find(childValue);
    return found == node.children_.end() ? -1l : found->second;
}
index_t Trie::addChild(index_t parentNodeIdx, char childValue)
{
    assert(parentNodeIdx >= 0 && parentNodeIdx < nodes_.size());
    assert(getChildIndex(parentNodeIdx, childValue) < 0);
    TrieNode& parentNode = nodes_[static_cast<std::size_t>(parentNodeIdx)];
    index_t result = static_cast<index_t>(nodes_.size());
    parentNode.children_.insert({childValue, nodes_.size()});
    nodes_.emplace_back(childValue, parentNodeIdx);
    return result;
}
index_t Trie::getNext(index_t nodeIdx, char value)
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
    index_t result = getNext(getSuffixLink(nodeIdx), value);
    addChild(result, value);
    return result;
}
index_t Trie::getSuffixLink(index_t nodeIdx)
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
index_t Trie::getEntrySuffixLink(index_t nodeIdx)
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
Trie::Trie() :
    entryCount_(0), currentNode_(0)
{
    // init root node
    nodes_.emplace_back();
    nodes_[0].suffixLink_ = 0;
    nodes_[0].entrySuffixLink_ = 0;
}
index_t Trie::addPattern(const std::string& pattern, std::size_t beginIdx, std::size_t endIdx)
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
    return nodes_[static_cast<std::size_t>(currentNode)].entry_;
}
index_t Trie::addPattern(const std::string& pattern)
{
    return addPattern(pattern, 0, pattern.size());
}
void Trie::next(char value, bool* matches, bool& hasMatches)
{
    index_t childIdx;
    assert(currentNode_ >= 0 && currentNode_ < nodes_.size());
    while (-1 == (childIdx = getChildIndex(currentNode_, value)) && currentNode_)
    {
        currentNode_ = getSuffixLink(currentNode_);
        if (nodes_[static_cast<std::size_t>(currentNode_)].entry_ != TrieNode::NOT_CALCULATED)
        {
            matches[nodes_[static_cast<std::size_t>(currentNode_)].entry_] = hasMatches = true;
        }
    }
    // child was found
    assert(childIdx >= 0 || !currentNode_);
    if (childIdx >= 0)
    {
        currentNode_ = childIdx;
    }
    if (currentNode_)
    {
        if (nodes_[static_cast<std::size_t>(currentNode_)].entry_ != TrieNode::NOT_CALCULATED)
        {
            matches[nodes_[static_cast<std::size_t>(currentNode_)].entry_] = hasMatches = true;
        }
        while ((childIdx = getEntrySuffixLink(childIdx)))
        {
            if (nodes_[static_cast<std::size_t>(childIdx)].entry_ != TrieNode::NOT_CALCULATED)
            {
                matches[nodes_[static_cast<std::size_t>(childIdx)].entry_] = hasMatches = true;
            }
        }
    }
}

WildcardSearcher::WildcardSearcher(const std::string& pattern):
    size_(pattern.size())
{
    std::size_t offsetStart = 0;  // end index of the previous pattern
    std::size_t patternStart;     // start index of current pattern
    bool readingPattern = pattern[0] != WILDCARD_SYMBOL; // true, if the last symbol is not a wildcard
    if (readingPattern)
    {
        patternStart = 0;
    }
    std::size_t i;
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
    // all characters were processed, there are 3 possible outcomes
    if (readingPattern) // in the end is a pattern, put it into the trie
    {
        long trieEntryIndex = trie_.addPattern(pattern, patternStart, size_);
        sequence_.emplace_back(size_ - offsetStart, size_ - patternStart, trieEntryIndex);
        wildcardPrefix_ = sequence_[0].offset_ - sequence_[0].length_;
        wildcardPostfix_ = 0;
    }
    else if (sequence_.size()) // in the end are some wildcards
    {
        wildcardPrefix_ =  sequence_[0].offset_ - sequence_[0].length_;
        wildcardPostfix_ = size_ - offsetStart;
    }
    else // all pattern consists of wildcards
    {
        wildcardPrefix_ =  size_;
        wildcardPostfix_ = 0;
    }
}

/* Search a text for the wildcard pattern */
std::vector<std::size_t> WildcardSearcher::search(const std::string& text)
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
    bool hasMatches = false;
    bool* matches = new bool[sequence_.size()];
    std::fill(matches, matches + sequence_.size(), false);

    if (sequence_.size() == 1)
    {
        for (std::size_t i = wildcardPrefix_; i < text.size() - wildcardPostfix_; ++i)
        {
            hasMatches = false;
            trie_.next(text[i], matches, hasMatches);
            if (hasMatches)
            {
                result.push_back(i + 1 - initialEntry.length_ - wildcardPrefix_);
                matches[0] = false;
            }
        }

        delete[] matches;

        return result; // there is a single pattern in the trie
    }

    std::priority_queue<TimerToken> timers; // priority queue for patterns, expected to be finished at certain position in future
    for (std::size_t i = wildcardPrefix_; i < text.size() - wildcardPostfix_; ++i)
    {
        // 0. check if there are any expected patterns at symbol (i - 1)
        std::stack<std::size_t> newTimers; // new items for the timer queue
        while (!timers.empty() && timers.top().stringIndex_ == i)
        {
            if (matches[timers.top().entry_])
            {
                std::size_t sequenceIndex = timers.top().sequenceIndex_;
                if (sequenceIndex == sequence_.size() - 1)
                {
                    assert(i >= size_ - wildcardPostfix_);
                    result.push_back(i + wildcardPostfix_ - size_);
                }
                else
                {
                    ++sequenceIndex;
                    newTimers.emplace(sequenceIndex);
                }
            }
            timers.pop();
        }
        // 1. put next expected patterns in the queue
        while(!newTimers.empty())
        {
           timers.emplace(i + sequence_[newTimers.top()].offset_, sequence_[newTimers.top()].entry_, newTimers.top());
           newTimers.pop();

        }
        // 2. check the trie if there are any matches
        if (hasMatches)
        {
            std::fill(matches, matches + sequence_.size(), false);
            hasMatches = false;
        }
        trie_.next(text[i], matches, hasMatches);
        if (matches[0])
        {
            timers.emplace(i + sequence_[1].offset_ + 1, sequence_[1].entry_, 1);
        }
    }

    // process the terminal matches
    if (matches[timers.top().entry_])
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

    delete[] matches;

    return result;
}

int main(void)
{
    std::string pattern, text;
    std::iostream::sync_with_stdio(false);

    std::cin >> pattern >> text;

    WildcardSearcher searcher(pattern);
    auto result = searcher.search(text);
    for (auto index : result)
    {
        std::cout << index << ' ';
    }

    return 0;
}
