/**
 * Задача 9. Алгоритм сжатия данных Хаффмана
 *
 * Реализуйте алгоритм построения оптимального префиксного кода Хаффмана.
 * При помощи алгоритма реализуйте две функции для создания архива из одного файла и извлечения файла из архива.
 */


#include <algorithm>
#include <assert.h>
#include <limits>
#include <stack>
#include <queue>
#include <vector>

#include <iomanip>
#include <iostream>

#include "Huffman.h"

////////////////////////////////////////////////////////////////////////////////
/// Bit managing routines
////////////////////////////////////////////////////////////////////////////////

/* Increment the unsigned number sored in array */
void incCode(std::size_t* code, std::size_t codeLength)
{
    for (std::size_t i = 0; i < codeLength; ++i)
    {
        if (++code[i])
        {
            return;
        }
    }
    assert(false);
}
/* Shift left the unsigned number stored in array */
void leftShiftCode(std::size_t* code, std::size_t codeLength, std::size_t shift)
{
    const std::size_t interElementShift = shift / (sizeof(std::size_t) * 8);
    const std::size_t intraElementShift = shift % (sizeof(std::size_t) * 8);
    const std::size_t complementaryShift = sizeof(std::size_t) * 8 - intraElementShift;
    assert(!(code[codeLength - 1 - interElementShift] >> complementaryShift));
    std::size_t carry = 0;
    for (std::size_t i = 0; i < codeLength; ++i)
    {
        std::size_t carry1 = code[i] >> complementaryShift;
        code[i] = (code[i] << intraElementShift) + carry;
        carry = carry1;
    }
    if (interElementShift)
    {
        for (std::size_t i = codeLength - 1; i - interElementShift < codeLength; --i)
        {
            code[i] = code[i - interElementShift];
        }
        for (std::size_t i = 0; i < interElementShift; ++i)
        {
            code[i] = 0;
        }
    }
}
/* Add the summand to the number stored in array */
void addCode(std::size_t* code, std::size_t codeLength, std::size_t summand)
{
    bool carry = code[0] + summand < code[0];
    code[0] += summand;
    if (carry)
    {
        for (std::size_t i = 1; i < codeLength; ++i)
        {
            if (++code[i])
            {
                return;
            }
        }
        assert(false);
    }
}
/*
 * Write the number stored in array to the output stream.
 * The number is denormalized what means its most significat digit aligned with the left border of the array.
 */
void writeDenormalizedCode(const std::size_t* code, std::size_t codeLength, const std::size_t size, IOutputStream& output, byte& lastByte, std::size_t& lastUnusedBitIndex)
{
    assert(size < codeLength * sizeof(std::size_t) * 8);
    assert(lastUnusedBitIndex <= sizeof(byte) * 8);
    assert(sizeof(std::size_t) > sizeof(byte));

    std::size_t bitIdx = sizeof(std::size_t) * 8;
    std::size_t elementIdx = codeLength - 1;
    std::size_t idx = 0;

    const std::size_t complement = sizeof(byte) * 8 - lastUnusedBitIndex;

    if (complement)
    {
        bitIdx -= lastUnusedBitIndex;
        byte b = static_cast<byte>(code[elementIdx] >> (sizeof(std::size_t) * 8 - lastUnusedBitIndex));
        lastByte += b;
        if (size < lastUnusedBitIndex)
        {
            lastUnusedBitIndex -= size;
            return;
        }
        else
        {
            output.Write(lastByte);
            idx += lastUnusedBitIndex;
        }
    }

    while (size - idx >= sizeof(byte) * 8)
    {
        if (bitIdx > sizeof(byte) * 8)
        {
            bitIdx -= sizeof(byte) * 8;
            byte b = static_cast<byte>(code[elementIdx] >> bitIdx);
            output.Write(b);
        }
        else if (bitIdx == sizeof(byte) * 8)
        {
            bitIdx = sizeof(std::size_t) * 8;
            byte b = static_cast<byte>(code[elementIdx]);
            output.Write(b);
        }
        else
        {
            byte b = static_cast<byte>(code[elementIdx] << (sizeof(byte) * 8 - bitIdx));
            bitIdx = sizeof(std::size_t) * 8 - bitIdx;
            --elementIdx;
            assert(elementIdx < codeLength);
            b += code[elementIdx] >> bitIdx;
            output.Write(b);
        }
        idx += sizeof(byte) * 8;
    }

    std::size_t residue = size - idx;
    if (!residue)
    {
        lastUnusedBitIndex = sizeof(byte) * 8;
    }
    else
    {
        lastUnusedBitIndex = sizeof(byte) * 8 - residue;
        if (bitIdx > sizeof(byte) * 8)
        {
            bitIdx -= sizeof(byte) * 8;
            lastByte = static_cast<byte>(code[elementIdx] >> bitIdx);
        }
        else if (bitIdx == sizeof(byte) * 8)
        {
            bitIdx = sizeof(std::size_t) * 8;
            lastByte = static_cast<byte>(code[elementIdx]);
        }
        else
        {
            lastByte = static_cast<byte>(code[elementIdx] << (sizeof(byte) * 8 - bitIdx));
            bitIdx = sizeof(std::size_t) * 8 - bitIdx;
            --elementIdx;
            assert(elementIdx < codeLength);
            lastByte += code[elementIdx] >> bitIdx;
        }
    }
}
/* Write to the stream number that repsedsents file size, position in file, etc. (must be in range of the unsigned int) */
void encodeSize(IOutputStream& output, std::size_t size)
{
    assert(size < std::numeric_limits<unsigned int>::max());
    for (std::size_t  i = 0; i < sizeof(unsigned int) / sizeof(byte); ++i)
    {
        byte b = static_cast<byte>(size);
        output.Write(b);
        size >>= sizeof(byte) * 8;
    }
    assert(!size);
}

////////////////////////////////////////////////////////////////////////////////
/// Entities
////////////////////////////////////////////////////////////////////////////////

/**
 * Basic class for Haffman coder and decoder.
 * \param TFeatureSize defines size of features that can be encoded in bytes.
 */
template<std::size_t TFeatureSize>
class Coder
{
public:
    /* Cardinality of set of all features of size TFeatureSize */
    const static std::size_t FEATURE_SPACE_SIZE = 1 << sizeof(byte) * 8 * TFeatureSize;
protected:
    std::size_t m_codeSize; /* How many elements of array is needed to store any code for a feature */
protected:
    /* Check if 2 features are equal */
    bool compareFeatures(const byte* feature1, const byte* feature2)
    {
        for (std::size_t i = 0; i < TFeatureSize; ++i)
        {
            if (feature1[i] != feature2[i])
            {
                return false;
            }
        }
        return true;
    }
    /* 
     * Generate a set of the first codes for each level of the tree
     * based on the numbers of leaves on each depth level of the tree.
     * New codes should be generated by incrementing these.
     */
    std::size_t* generateInitialCode(std::size_t* levelCounts, std::size_t maxLevel)
    {
        m_codeSize = maxLevel / (sizeof(std::size_t) * 8) + 1;
        assert(maxLevel < sizeof(std::size_t) * 8 * m_codeSize);
        std::size_t* levelCodes = new std::size_t[m_codeSize * (maxLevel + 1)];
        std::fill(levelCodes, levelCodes + m_codeSize * (maxLevel + 1), 0);
        for (std::size_t i = 1; i <= maxLevel; ++i)
        {
            std::size_t* previousCode = levelCodes + (i - 1) * m_codeSize;
            std::size_t* currentCode = levelCodes + i * m_codeSize;
            for (std::size_t i = 0; i < m_codeSize; ++i)
            {
                currentCode[i] = previousCode[i];
            }
            addCode(currentCode, m_codeSize, levelCounts[i - 1]);
            leftShiftCode(currentCode, m_codeSize, 1);
        }
        return levelCodes;
    }
#if defined(VERBOSE_CODER) || defined(STATS_CODER)
    /** Debug output of the feature value */
    void printFeature(const byte* feature) const
    {
        for (std::size_t i = 0; i < TFeatureSize - 1; ++i)
        {
            std::cout << feature[i] << ":";
        }
        std::cout << feature[TFeatureSize - 1];
    }
#endif
};

/** Description of some feature of the input data*/
template<std::size_t TFeatureSize>
struct CodeBookEntry
{
    unsigned int m_size;          /** Occurence count */
    unsigned int m_level;         /** Depth level in the Huffman tree */
    std::size_t* m_code;          /** Denormalied code prepared to be writthen out */
    byte m_feature[TFeatureSize]; /** Sequence of bytes that are the feature */
};

/** Comparator of features for the priority queue */
template<std::size_t TFeatureSize, std::size_t TCodeSize>
bool compareCodebookEntries(const CodeBookEntry<TFeatureSize>* lhs, const CodeBookEntry<TFeatureSize>* rhs)
{
  return lhs->m_size > rhs->m_size;
}

/** Comparator of features for the priority queue */
struct FeatureTreeNode
{
    std::size_t m_count;      /* Occurence count of features in the tree */
    FeatureTreeNode* m_left;  /* Left subtree */
    FeatureTreeNode* m_right; /* Right subtree */

    FeatureTreeNode (std::size_t count, FeatureTreeNode* left, FeatureTreeNode* right): m_count(count), m_left(left), m_right(right)
    {
    }
    FeatureTreeNode (FeatureTreeNode* left, FeatureTreeNode* right): m_count(left->m_count + right->m_count), m_left(left), m_right(right)
    {
    }
};

/** Leaf node of the Haffman tree that is linked to the the feature */
template<std::size_t TFeatureSize>
struct FeatureTreeLeafNode: FeatureTreeNode
{
    CodeBookEntry<TFeatureSize>* m_entry;

    FeatureTreeLeafNode(CodeBookEntry<TFeatureSize>* entry): FeatureTreeNode(entry->m_size, nullptr, nullptr), m_entry(entry)
    {
    }
};

/** Translates stream of bytes into Haffman code */
template<std::size_t TFeatureSize>
class Encoder: public Coder<TFeatureSize>
{
    using SizedCodeBookEntry = CodeBookEntry<TFeatureSize>;
    using SizedFeatureTreeLeafNode = FeatureTreeLeafNode<TFeatureSize>;
protected:
    const std::vector<byte>& m_data;    /** Readonly storage of incoming data */
    SizedCodeBookEntry* m_featureSpace; /** Ordered array of all possible features of size TFeatureSize */
    std::size_t m_nonZeroFeatureCount;  /** Number of features that were encountered at least once */
    byte m_inputBuffer[TFeatureSize];   /** Bufferer where byte by byte the next feature from the stream is formed */
    std::size_t m_inputBufferIndex;     /** Index of the used element in m_inputBuffer */
    std::size_t m_totalFeatureCount;    /** Total number of feature encounters, size of the data in features */
    std::size_t* m_levelCounts;         /** Number of leaves on each depth level of the Haffman tree */
    std::size_t m_minNonZeroLevel;      /** Index of the first nonzero m_levelCounts */
    std::size_t m_maxLevel;             /** Maximal index in m_levelCounts */

    /** Set initial values for all features */
    void createFeatureSpace()
    {
        m_featureSpace = new SizedCodeBookEntry[Coder<TFeatureSize>::FEATURE_SPACE_SIZE];
        for (std::size_t i = 0; i < Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            m_featureSpace[i].m_size = 0;
            m_featureSpace[i].m_level = 0;
            m_featureSpace[i].m_code = nullptr;
            for (std::size_t j = 0; j < TFeatureSize; ++j)
            {
                m_featureSpace[i].m_feature[j] = static_cast<byte>(i >> sizeof(byte) * 8 * j);
            }
        }
    }
    /** Set initial values for all features */
    std::size_t getFeatureIndex(byte* feature)
    {
        std::size_t index = 0;
        for (std::size_t i = TFeatureSize - 1; i < TFeatureSize; --i)
        {
            index <<= sizeof(byte) * 8;
            index += feature[i];
        }
        assert(this->compareFeatures(feature, m_featureSpace[index].m_feature));
        return index;
    }
    /** Register feature occurence */
    virtual SizedCodeBookEntry* addFeature()
    {
        auto result = &m_featureSpace[getFeatureIndex(m_inputBuffer)];
        ++result->m_size;
        ++m_totalFeatureCount;
        return result;
    }
    /** Build Haffman tree */
    FeatureTreeNode* buildFeatureTree()
    {
        auto comparator = [](FeatureTreeNode* x, FeatureTreeNode* y) { return x->m_count > y->m_count; };
        std::priority_queue<FeatureTreeNode*, std::vector<FeatureTreeNode*>, decltype(comparator)> queue(comparator);
        for (std::size_t i = 0; i < Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            if (m_featureSpace[i].m_size)
            {
                queue.push(new SizedFeatureTreeLeafNode(&m_featureSpace[i]));
            }
        }
        m_nonZeroFeatureCount = queue.size();

        while (queue.size() > 1)
        {
            FeatureTreeNode* first = queue.top();
            queue.pop();
            FeatureTreeNode* second = queue.top();
            queue.pop();
            FeatureTreeNode* merge = new FeatureTreeNode(first, second);
            queue.push(merge);
        }
        return queue.top();
    }
    /** Calculate depth level for features in Haffman tree */
    void calculateLevels()
    {
        // 0. build tree
        FeatureTreeNode *node = this->buildFeatureTree();

        // 1. traverse tree, set m_featureSpace[i]->m_space from encounter counts to size of a code
        m_levelCounts = new std::size_t[this->m_nonZeroFeatureCount + 1];
        std::fill(m_levelCounts, m_levelCounts + this->m_nonZeroFeatureCount + 1, 0);
        std::stack<std::pair<FeatureTreeNode*, std::size_t>> stack;
        unsigned int depth = 0;
        m_maxLevel = 0;

        while (node || !stack.empty())
        {
            while (node)
            {
                stack.push(std::make_pair(node, depth++));
                node = node->m_left;
            }
            auto current = stack.top();
            depth = current.second;
            #ifdef  VERBOSE_CODER
            std::cout << "size: " << current.first->m_count << " depth: " << depth << std::endl;
            #endif
            stack.pop();

            if (!current.first->m_right && !current.first->m_left)
            {
                SizedFeatureTreeLeafNode* leafNode = static_cast<SizedFeatureTreeLeafNode*>(current.first);
                SizedCodeBookEntry* entry = leafNode->m_entry;
                ++m_levelCounts[entry->m_level = depth];
                if (depth > m_maxLevel)
                {
                    m_maxLevel = depth;
                }
                #ifdef VERBOSE_CODER
                std::cout << "leaf ";
                this->printFeature(leafNode->m_entry->m_feature);
                std::cout << std::endl;
                #endif
            }

            node = current.first->m_right;
            ++depth;
            delete current.first;
        }

        // 2. calculate minimal nonzero level
        for (std::size_t i = 0; i <= m_maxLevel; ++i)
        {
            if (m_levelCounts[i])
            {
                m_minNonZeroLevel = i;
                break;
            }
        }
    }
    /** Generate codes for each feature, then store features as linearized Haffman tree in the output */
    void encodeCodebook(IOutputStream& output)
    {
        // 0. generate initial codes for each level
        std::size_t* levelCodes = this->generateInitialCode(m_levelCounts, m_maxLevel);

        // 1. build codes
        for (std::size_t i = 0; i < Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            auto& entry = this->m_featureSpace[i];
            const std::size_t level = entry.m_level;
            if (level)
            {
                #ifdef VERBOSE_CODER
                std::cout << "building code for ";
                this->printFeature(entry.m_feature);
                std::cout << std::endl;
                #endif
                std::size_t* levelCode = levelCodes + level * this->m_codeSize;
                entry.m_code = new std::size_t[this->m_codeSize];
                for (std::size_t j = 0; j < this->m_codeSize; ++j)
                {
                    entry.m_code[j] = levelCode[j];
                }
                leftShiftCode(entry.m_code, this->m_codeSize, (sizeof(std::size_t) * 8 * this->m_codeSize - level));
                #ifdef VERBOSE_CODER
                for (std::size_t j = this->m_codeSize - 1; j < this->m_codeSize; --j)
                {
                    std::cout << std::hex /*<< std::setw(sizeof(std::size_t))*/ << std::setfill('0') << entry.m_code[j] << ":";
                }
                std::cout <<std::dec << std::endl;
                #endif
                incCode(levelCode, this->m_codeSize);
            }
        }
        delete[] levelCodes;

        // 2. Encode level sizes
        assert(m_minNonZeroLevel > 0 && m_minNonZeroLevel <= std::numeric_limits<byte>::max());
        output.Write(static_cast<byte>(m_minNonZeroLevel));
        assert(m_maxLevel <= std::numeric_limits<byte>::max());
        output.Write(static_cast<byte>(m_maxLevel));

        #ifdef STATS_CODER
        std::cout << "encode m_minNonZeroLevel: " << m_minNonZeroLevel << std::endl;
        std::cout << "encode m_maxLevel: " << m_maxLevel << std::endl;
        for (std::size_t i = 0; i <= m_maxLevel; ++i)
        {
            std::cout << "encode level count i = " << i << " : " << m_levelCounts[i] << std::endl;
        }
        #endif

        for (std::size_t i = m_minNonZeroLevel; i <= m_maxLevel; ++i)
        {
            assert(m_levelCounts[i] <= 1 << sizeof(byte) * 8 * TFeatureSize);
            std::size_t level = m_levelCounts[i];
            for (std::size_t j = 0; j < TFeatureSize; ++j)
            {
                output.Write(static_cast<byte>(level));
                level >>= sizeof(byte) * 8;
            }
            m_levelCounts[i] += m_levelCounts[i - 1];
        }

        // 3. Calc feature codes
        byte* features = new byte[this->m_nonZeroFeatureCount * TFeatureSize];
        for (std::size_t i = 0; i < this->Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            if (this->m_featureSpace[i].m_size)
            {
                #ifdef VERBOSE_CODER
                std::cout << "prepare feature ";
                this->printFeature(this->m_featureSpace[i].m_feature);
                std::cout << " by index " << m_levelCounts[this->m_featureSpace[i].m_level - 1] << std::endl;
                #endif
                std::size_t offset = m_levelCounts[this->m_featureSpace[i].m_level - 1]++;
                byte* feature = this->m_featureSpace[i].m_feature;
                for (std::size_t j = 0; j < TFeatureSize; ++j)
                {
                    features[offset * TFeatureSize + j] = feature[j];
                }
            }
        }

        // 4. Encode features
        for (std::size_t i = 0; i < this->m_nonZeroFeatureCount * TFeatureSize; ++i)
        {
            #ifdef VERBOSE_CODER
            if (i % TFeatureSize == 0)
            {
                std::cout << "encode feature ";
                this->printFeature(features + i);
                std::cout << std::endl;
            }
            #endif
            output.Write(features[i]);
        }

        delete[] features;
    }
    /** Encode number of features in the data */
    void encodeDataSize(IOutputStream& output)
    {
        #ifdef STATS_CODER
        std::cout << "encode size: " << m_totalFeatureCount << std::endl;
        #endif
        encodeSize(output, m_totalFeatureCount);

    }
    /** Translate data into Haffman codes and export it */
    virtual void encodeData(IOutputStream& output)
    {
        assert(m_totalFeatureCount * TFeatureSize + m_inputBufferIndex == m_data.size());
        std::size_t lastUnusedBitIndex = sizeof(byte) * 8;
        byte feature[TFeatureSize];
        byte buffer;
        for (std::size_t i = 0; i < m_totalFeatureCount * TFeatureSize; i += TFeatureSize)
        {
            for (std::size_t j = 0; j < TFeatureSize; ++j)
            {
                feature[j] = this->m_data[i + j];
            }
            const SizedCodeBookEntry& entry = m_featureSpace[this->getFeatureIndex(feature)];
            #ifdef VERBOSE_CODER
            std::cout << "writing ";
            this->printFeature(feature);
            std::cout << std::endl;
            #endif
            writeDenormalizedCode(entry.m_code, this->m_codeSize, entry.m_level, output, buffer, lastUnusedBitIndex);
        }
        if (lastUnusedBitIndex != sizeof(byte) * 8)
        {
            output.Write(buffer);
        }
    }
    /** Export incomplete feature as is without encoding */
    void encodeResidue(IOutputStream& output)
    {
        if (TFeatureSize == 1)
        {
            return;
        }
        assert(m_inputBufferIndex < TFeatureSize);
        output.Write(static_cast<byte>(m_inputBufferIndex));
        for (std::size_t i = 0; i < m_inputBufferIndex; ++i)
        {
            output.Write(m_inputBuffer[i]);
        }
        #ifdef VERBOSE_CODER
        std::cout << "encode residue size: " << m_inputBufferIndex << std::endl << "encode residue: ";
        for (std::size_t i = 0; i < m_inputBufferIndex; ++i)
        {
            std::cout << m_inputBuffer[i] << ":";
        }
        std::cout << std::endl;
        #endif
    }
public:
    Encoder(const std::vector<byte>& data):
        Coder<TFeatureSize>(), m_data(data), m_featureSpace(nullptr), m_inputBufferIndex(0), m_totalFeatureCount(0), m_levelCounts(nullptr)
    {
        createFeatureSpace();
    }
    virtual ~Encoder()
    {
        for (std::size_t i = 0; i < Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            delete[] m_featureSpace[i].m_code;
        }
        delete[] m_featureSpace;
        delete[] m_levelCounts;
    }
    /** Register a byte from the input data */
    void addByte(byte value)
    {
        m_inputBuffer[m_inputBufferIndex++] = value;
        if (m_inputBufferIndex == TFeatureSize)
        {
            m_inputBufferIndex = 0;
            addFeature();
        }
    }
    /** Claculate precise number of bits in the compressed sequence before starting its creation */
    virtual std::size_t calculateCompressedSize()
    {
        calculateLevels();

        std::size_t size = 2 * 8 + // size bytes
                  (1 + m_inputBufferIndex) * 8 + // residue
                4 * 8 + // lever array sizes
                (m_maxLevel - m_minNonZeroLevel) * 8 * TFeatureSize; // level array
        //  aggregate encoded data size
        for (std::size_t i = 0; i < Coder<TFeatureSize>::FEATURE_SPACE_SIZE; ++i)
        {
            if (m_featureSpace[i].m_size)
            {
                size += m_featureSpace[i].m_size * m_featureSpace[i].m_level + 8 * TFeatureSize;
            }
        }
        return size;
    }
    /** Stream data as an encoded sequence */
    void execute(IOutputStream& output)
    {
        encodeResidue(output);
        encodeCodebook(output);
        encodeDataSize(output);
        encodeData(output);
    }
};

/** Leaf node of the Haffman tree used in decoding */
template<std::size_t TFeatureSize>
struct DecodeFeatureTreeLeafNode: FeatureTreeNode
{
    byte m_feature[TFeatureSize]; /* The feature value */

    DecodeFeatureTreeLeafNode(byte* feature): FeatureTreeNode(0, nullptr, nullptr)
    {
        std::copy(feature, feature + TFeatureSize, m_feature);
    }
};

/** Transforms encoded sequence back into the original */
template<std::size_t TFeatureSize>
class Decoder: public Coder<TFeatureSize>
{
    using SizedCodeBookEntry = CodeBookEntry<TFeatureSize>;
    using SizedDecodeFeatureTreeLeafNode = DecodeFeatureTreeLeafNode<TFeatureSize>;
protected:
    IInputStream& m_input;        /* input stream */
    IOutputStream& m_output;      /* output stream */
    byte m_residue[TFeatureSize]; /* residual ending bytes */
    byte m_residueSize;           /* number of residual bytes */
    /** Read unencoded residual bytes */
    void decodeResidue()
    {
        if (TFeatureSize == 1)
        {
            m_residueSize = 0;
            return;
        }
        this->m_input.Read(m_residueSize);
        for (std::size_t i = 0; i < m_residueSize; ++i)
        {
            this->m_input.Read(m_residue[i]);
        }
        #ifdef VERBOSE_CODER
        std::cout << "decode residue size: " << static_cast<std::size_t>(m_residueSize) << std::endl << "decode residue: ";
        for (std::size_t i = 0; i < m_residueSize; ++i)
        {
            std::cout << m_residue[i] << ":";
        }
        std::cout << std::endl;
        #endif
    }
    /** 
     * Insert leaf to the Haffman tree.
     * It uses code as path in the tree.
     */
    void addNode(FeatureTreeNode* root, byte* feature, const std::size_t* code, std::size_t codeLength)
    {
        std::size_t elementIndex = codeLength / (sizeof(std::size_t) * 8);
        std::size_t bitIndex = codeLength % (sizeof(std::size_t) * 8);
        std::size_t mask = 1 << (bitIndex - 1);
        while (codeLength > 1)
        {
            if (code[elementIndex] & mask)
            {
                if (!root->m_left)
                {
                    root->m_left = new FeatureTreeNode(0, nullptr, nullptr);
                }
                root = root->m_left;
            }
            else
            {
                if (!root->m_right)
                {
                    root->m_right = new FeatureTreeNode(0, nullptr, nullptr);
                }
                root = root->m_right;
            }
            --bitIndex;
            --codeLength;
            if (bitIndex < sizeof(std::size_t) * 8)
            {
                mask >>= 1;
            }
            else
            {
                bitIndex = sizeof(std::size_t) * 8 - 1;
                mask = 1 << bitIndex;
            }
        }
        if (code[elementIndex] & mask)
        {
            assert(!root->m_left);
            root->m_left = new SizedDecodeFeatureTreeLeafNode(feature);
        }
        else
        {
            assert(!root->m_right);
            root->m_right = new SizedDecodeFeatureTreeLeafNode(feature);
        }
    }
    /** Decode features and their codes, build the tree */
    FeatureTreeNode* decodeCodebook()
    {
        byte minNonZeroLevel, maxLevel;
        this->m_input.Read(minNonZeroLevel);
        #ifdef STATS_CODER
        std::cout << "decode minNonZeroLevel= " << static_cast<int>(minNonZeroLevel) << std::endl;
        #endif
        this->m_input.Read(maxLevel);
        #ifdef STATS_CODER
        std::cout << "decode maxLevel= " << static_cast<int>(maxLevel) <<std::endl;
        #endif

        std::size_t* levelCounts = new std::size_t[maxLevel + 1];
        for (std::size_t i = 0; i <= maxLevel; ++i)
        {
            levelCounts[i] = 0;
        }
        std::size_t featureCount = 0;
        for (std::size_t i = minNonZeroLevel; i <= maxLevel; ++i)
        {
            for (std::size_t j = 0; j < TFeatureSize; ++j)
            {
                byte b;
                this->m_input.Read(b);
                std::size_t buffer = b;
                buffer <<= sizeof(byte) * 8 * j;
                levelCounts[i] += buffer;
            }
            featureCount += levelCounts[i];
        }
        #ifdef STATS_CODER
        for (std::size_t i = 0; i <= maxLevel; ++i)
        {
            std::cout << "decode level counts i = " << i <<" : " << levelCounts[i] <<std::endl;
        }
        #endif
        assert(featureCount <= 1 << sizeof(byte) * 8 * TFeatureSize);

        std::size_t* codes = this->generateInitialCode(levelCounts, maxLevel);
        FeatureTreeNode* root = new FeatureTreeNode(0, nullptr, nullptr);

        for (std::size_t i = minNonZeroLevel; i <= maxLevel; ++i)
        {
            for (std::size_t j = 0; j < levelCounts[i]; ++j)
            {
                byte feature[TFeatureSize];
                for (std::size_t i = 0; i < TFeatureSize; ++i)
                {
                    this->m_input.Read(feature[i]);
                }
                std::size_t* code = codes + i * this->m_codeSize;
                #ifdef VERBOSE_CODER
                std::cout << "building code for ";
                this->printFeature(feature);
                std::cout << std::endl;
                for (std::size_t k = this->m_codeSize - 1; k < this->m_codeSize; --k)
                {
                    std::cout << std::hex << std::setfill('0') << code[k] << ":";
                }
                std::cout << std::endl;
                #endif
                addNode(root, feature, code, i);
                incCode(code, this->m_codeSize);
            }
        }

        delete[] codes;
        delete[] levelCounts;

        return root;
    }
    /** Read data size */
    unsigned int decodeSize()
    {
        unsigned int size = 0;
        for (std::size_t i = 0; i < sizeof(unsigned int) / sizeof(byte); ++i)
        {
            byte byte;
            this->m_input.Read(byte);
            unsigned int c = byte;
            c <<= sizeof(byte) * 8 * i;
            size += c;
        }
        return size;
    }
    /** Decode data and stream it */
    virtual void decodeData(const FeatureTreeNode* tree, unsigned int dataSize)
    {
        byte b;
        const FeatureTreeNode* node = tree;
        while(this->m_input.Read(b))
        {
            std::size_t index = sizeof(byte) * 8 - 1;
            byte mask = static_cast<byte>(1 << index);
            while (index < sizeof(byte) * 8)
            {
                if (mask & b)
                {
                    node = node->m_left;
                }
                else
                {
                    node = node->m_right;
                }
                assert(node);
                if (!node->m_left && !node->m_right)
                {
                    const SizedDecodeFeatureTreeLeafNode* leafNode = static_cast<const SizedDecodeFeatureTreeLeafNode*>(node);
                    #ifdef VERBOSE_CODER
                    std::cout << "Reading: ";
                    this->printFeature(leafNode->m_feature);
                    std::cout << std::endl;
                    #endif
                    for (std::size_t i = 0; i < TFeatureSize; ++i)
                    {
                        m_output.Write(leafNode->m_feature[i]);
                    }
                    if (!--dataSize)
                    {
                        for (std::size_t i = 0; i < m_residueSize; ++i)
                        {
                           m_output.Write(m_residue[i]);
                        }
                        return;
                    }
                    node = tree;
                }
                --index;
                mask >>= 1;
            }
        }
        assert(node == tree);
    }
    /** Delete the tree without recursion */
    void deleteTree(FeatureTreeNode* root)
    {
        std::stack<FeatureTreeNode*> stack;
        while (root || !stack.empty())
        {
            while (root)
            {
                stack.push(root);
                root = root->m_left;
            }
            root = stack.top();
            stack.pop();

            FeatureTreeNode* previousNode = root;
            root = root->m_right;
            delete previousNode;
        }
    }
public:
    Decoder(IInputStream& input, IOutputStream& output): Coder<TFeatureSize>(), m_input(input), m_output(output)
    {
    }
     /** Translate compressed sequence into the */
    virtual void execute()
    {
        decodeResidue();
        FeatureTreeNode* tree = decodeCodebook();
        unsigned int dataSize = decodeSize();
        #ifdef STATS_CODER
        std::cout << "Read size = " << std::dec << dataSize << std::endl;
        #endif
        decodeData(tree, dataSize);
        deleteTree(tree);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////
/// Encoding and decoding routines, select the best encoder possible
////////////////////////////////////////////////////////////////////////////////////////////

void Encode(IInputStream& original, IOutputStream& compressed)
{
    std::vector<byte> data;
    Encoder<1> encoder(data);

    byte b;
    while (original.Read(b))
    {
        data.push_back(b);
        encoder.addByte(b);
    }

    encoder.calculateCompressedSize();
    encoder.execute(compressed);

}

void Decode(IInputStream& compressed, IOutputStream& original)
{
    Decoder<1> decoder(compressed, original);
    decoder.execute();
}