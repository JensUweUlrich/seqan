// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author:  Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_

// --------------------------------------------------------------------------
// Class
// --------------------------------------------------------------------------
namespace seqan{

template<typename TValue>
struct BDHash<TValue, Normal>
{
public:
    Shape<TValue, SimpleShape> kmerShape;
    uint16_t kmerSize{0};

    inline void resize(TKmerSize newKmerSize)
    {
        kmerSize = newKmerSize;
        seqan::resize(kmerShape, kmerSize);
    }

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline uint32_t get_threshold(uint32_t t, uint16_t e)
    {
        uint32_t threshold = (t - kmerSize * (1+e) + 1);
        return threshold > 0 ? threshold : 0;
    }

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            uint32_t possible = seqan::length(text) - kmerSize + 1;

            std::vector<uint64_t> kmerHashes(possible, 0);

            hashInit(begin(text));
            auto it = begin(text);

            for (uint32_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = hashNext(it);
                ++it;
            }
            // std::sort(std::begin(kmerHashes), std::end(kmerHashes));
            kmerHashes.erase(std::unique(std::begin(kmerHashes), std::end(kmerHashes)), std::end(kmerHashes));
            return kmerHashes;
        }
    }
};

template<typename TValue, uint16_t o>
struct BDHash<TValue, Offset<o>>
{
public:
    Shape<TValue, SimpleShape> kmerShape;
    uint16_t offset = o;
    uint16_t kmerSize{0};

    inline void resize(TKmerSize newKmerSize)
    {
        kmerSize = newKmerSize;
        seqan::resize(kmerShape, kmerSize);
    }

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline uint32_t get_threshold(uint32_t t, uint16_t e)
    {
        uint32_t threshold = std::floor((t - kmerSize * (1 + e) + 1)/o);
        return threshold > 0 ? threshold : 0;
    }

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            // how many test positions are left if we take every offset'th kmer
            uint16_t x = (seqan::length(text) - kmerSize) % offset;
            // how many kmers are there when we take every offset'th kmer
            // possible = something left (1/0) + how many fit in the text
            // if something is left, we add a kmer that covers these positions
            uint32_t possible = bool(x) + (seqan::length(text) - kmerSize + offset - x) / offset;

            std::vector<uint64_t> kmerHashes(possible, 0);

            hashInit(begin(text));
            auto it = begin(text);

            uint32_t positions = seqan::length(text) - kmerSize + 1;

            for (uint32_t i = 0, j = 0; i < positions; ++i)
            {
                uint64_t kmerHash = hashNext(it);
                if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
                {
                    kmerHashes[j] = kmerHash;
                    break;
                }
                if (i - j * offset == 0) // we found the j'th kmer with offset
                {
                    kmerHashes[j] = kmerHash;
                    ++j;
                }
                ++it;
            }
            // std::sort(std::begin(kmerHashes), std::end(kmerHashes));
            kmerHashes.erase(std::unique(std::begin(kmerHashes), std::end(kmerHashes)), std::end(kmerHashes));
            return kmerHashes;
        }
    }
};

template<typename TValue, uint16_t k, uint32_t w>
struct BDHash<TValue, Minimizer<k, w>>
{
public:
    uint64_t static constexpr seed{0x8F3F73B5CF1C9ADE};
    uint16_t kmerSize{k};
    uint32_t windowSize{w};

    // All positions are inclusive, i.e. kmer starts at b and ends in e => [b,e]
    std::vector<uint64_t> minBegin;
    std::vector<uint64_t> minEnd;
    std::vector<uint32_t> coverage;
    std::vector<uint64_t> coverageBegin;
    std::vector<uint64_t> coverageEnd;

    Shape<TValue, SimpleShape> kmerShape;
    Shape<TValue, SimpleShape> revCompShape;

    inline void resize(uint16_t newKmerSize, uint32_t newWindowSize = w)
    {
        kmerSize = newKmerSize;
        windowSize = newKmerSize > newWindowSize ? newKmerSize : newWindowSize;
        seqan::resize(kmerShape, kmerSize);
        seqan::resize(revCompShape, kmerSize);
    }

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    template<typename TIt>
    inline void revHashInit(TIt it)
    {
        seqan::hashInit(revCompShape, it);
    }

    template<typename TIt>
    inline auto revHashNext(TIt it)
    {
        return seqan::hashNext(revCompShape, it);
    }

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString & text) // TODO cannot be const for ModifiedString
    {
        if (kmerSize > seqan::length(text))
            return std::vector<uint64_t> {};
        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) > windowSize ? seqan::length(text) - windowSize + 1 : 1;
        uint32_t windowKmers = windowSize - kmerSize + 1;

        std::vector<uint64_t> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);

        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            uint64_t distance = std::distance(begin(text), it);
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(std::get<0>(*min));
        minBegin.push_back(std::get<1>(*min));
        minEnd.push_back(std::get<2>(*min));

        for (uint32_t i = 1; i < possible; ++i)
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            uint64_t distance = std::distance(begin(text), it);
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;

            if (std::get<0>(windowValues.back()) < std::get<0>(*min))
                min = std::end(windowValues) - 1;

            kmerHashes.push_back(std::get<0>(*min));
            minBegin.push_back(std::get<1>(*min));
            minEnd.push_back(std::get<2>(*min));
        }
        // std::sort(std::begin(kmerHashes), std::end(kmerHashes));
        kmerHashes.erase(std::unique(std::begin(kmerHashes), std::end(kmerHashes)), std::end(kmerHashes));
        return kmerHashes;
    }

    inline uint32_t get_threshold(uint32_t t, uint16_t e)
    {
        (void) t;
        uint32_t destroyed{0};
        minBegin.erase(std::unique(std::begin(minBegin), std::end(minBegin)), std::end(minBegin));
        uint32_t available{static_cast<uint32_t>(minBegin.size())};
        // destroyed += e;

        for (uint16_t i = 0; i < e; ++i)
        {
            get_coverage();
            auto max = std::max_element(coverage.begin(), coverage.end());
            if (i == e - 1)
                destroyed += *max;
            else
            {
                destroyed += *max;
                std::vector<uint64_t> newBegin;
                std::vector<uint64_t> newEnd;
                newBegin.reserve(available - destroyed);
                newEnd.reserve(available - destroyed);

                auto idx = std::distance(coverage.begin(), max);
                auto cb = coverageBegin[idx];
                auto ce =  coverageEnd[idx];
                for (uint64_t i = 0; i < minBegin.size(); ++i)
                {
                    auto mb = minBegin[i];
                    auto me = minEnd[i];
                    if ((mb >= cb && mb <= ce) || (me >= cb && me <= ce))
                        continue;
                    newBegin.push_back(mb);
                    newEnd.push_back(me);
                }
                minBegin = std::move(newBegin);
                minEnd = std::move(newEnd);
                coverageBegin.clear();
                coverageEnd.clear();
                coverage.clear();
            }
        }
        return destroyed > available ? 0 : available - destroyed;
    }

    inline void get_coverage()
    {
        uint64_t bIndex{1};
        uint64_t eIndex{0};
        coverageBegin.push_back(minBegin[0]);
        coverage.push_back(1);

        minBegin.erase(std::unique(std::begin(minBegin), std::end(minBegin)), std::end(minBegin));
        minEnd.erase(std::unique(std::begin(minEnd), std::end(minEnd)), std::end(minEnd));

        while ((bIndex < minBegin.size() ) || (eIndex < minEnd.size()))
        {
            uint64_t begin = bIndex < minBegin.size() ? minBegin[bIndex] : 0xFFFFFFFFFFFFFFFFULL;
            uint64_t end   = minEnd[eIndex];
            // Overlap
            if (begin < end)
            {
                coverageEnd.push_back(begin-1);
                coverageBegin.push_back(begin);
                coverage.push_back(coverage.back()+1);
                ++bIndex;
            }
            // Flatten consecutive positions, where one kmer ends and other one starts
            if (begin == end)
            {
                coverageEnd.push_back(begin-1);
                coverageBegin.push_back(begin);
                coverage.push_back(coverage.back()+1);
                while (minBegin[bIndex] == minEnd[eIndex])
                {
                    ++bIndex;
                    ++eIndex;
                }
                --eIndex;
            }
            // Kmer ends
            if (end < begin)
            {
                coverageEnd.push_back(end);
                coverageBegin.push_back(end+1);
                coverage.push_back(coverage.back()-1);
                ++eIndex;
            }
        }
        coverageBegin.pop_back();
        coverage.pop_back();
    }
};


}   // namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_
