/*Задача 16. Поиск подстроки
 * Найдите все вхождения шаблона в строку.
 * Длина шаблона – p, длина строки ­– n. Время O(n + p), доп. память – O(p).
 * Вариант 2. С помощью z-функции.
 */


#include <algorithm>
#include <iostream>
#include <vector>

/**
 * Find indices of substring occurencies in the input string with the Z function (Knuth-Morris-Pratt algorithm).
 * Algorithm slightly modified:
 * instead of concatenation of strings it is split into 2 distinct loops.
 */
std::vector<std::size_t> findSubstrings(const std::string& substring, const std::string& input)
{
    const std::size_t substringLength = substring.length();
    const std::size_t inputLength  = input.length();
    std::vector<std::size_t> result;
    if (substringLength > inputLength)
    {
        return result;
    }

    std::vector<std::size_t> zFn(substringLength + inputLength);
    zFn[0] = 0;
    std::size_t left = 0, right = 0;

    // 1. classic z function calculation for the substring
    for (std::size_t i = 1; i < substringLength; ++i)
    {
        std::size_t zFnI;
        if (right <= i)
        {
            zFnI = 0;
        }
        else
        {
            zFnI = std::min(right - i, zFn[i - left]);
        }
        while (i + zFnI < substringLength && substring[zFnI] == substring[i + zFnI])
        {
            ++zFnI;
        }
        if (i + zFnI > right)
        {
            left = i;
            right = i + zFnI;
        }
        zFn[i] = zFnI;
    }

    // 2. z function calculation for the input string:
    // keep in mind that the substring precedes the input thus indexes are modified:
    for (std::size_t i = 0; i < inputLength; ++i)
    {
        std::size_t zFnI;
        if (right <= i + substringLength)
        {
            zFnI = 0;
        }
        else
        {
            zFnI = std::min(right - i - substringLength, zFn[i + substringLength - left]);
        }
        while (zFnI < substringLength && i + zFnI < inputLength && substring[zFnI] == input[i + zFnI])
        {
            ++zFnI;
        }
        if (i + substringLength + zFnI > right)
        {
            left = i + substringLength;
            right = i + substringLength + zFnI;
        }
        zFn[i + substringLength] = zFnI;
        //std::cout << i << " : " << zFnI << '\n';
        if (zFnI == substringLength)
        {
            result.push_back(i);
        }
    }
    return result;
}


int main(void)
{
    std::string substring, input;
    std::cin >> substring >> input;

    auto substringIndices = findSubstrings(substring, input);

    for (auto index: substringIndices)
    {
        std::cout << index << ' ';
    }

    return 0;
}
