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

template<typename TValue, typename TChunks>
struct BDHashBase
{
public:
    uint8_t chunks{TChunks::VALUE};
    std::vector<uint8_t> chunkMap{0};
    Shape<TValue, SimpleShape> kmerShape;
    uint16_t kmerSize{0};
    uint16_t effectiveChunks{1};
    uint8_t significantBits{0};
    uint8_t significantPositions{0};
    uint64_t chunkOffset{0};

    inline void setChunkOffset(uint64_t chunkOffset_)
    {
        chunkOffset = chunkOffset_;
    }

    inline void setEffective(uint16_t effectiveChunks_)
    {
        effectiveChunks = effectiveChunks_;
    }

    inline void setPos(uint8_t significantPositions_)
    {
        significantPositions = significantPositions_;
    }

    inline void setBits(uint8_t significantBits_)
    {
        significantBits = significantBits_;
    }

    inline void setMap(std::vector<uint8_t> chunkMap_)
    {
        chunkMap = chunkMap_;
    }

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
};

template<uint16_t k, typename TChunks>
struct BDHash<Dna, Normal<k>, TChunks> : BDHashBase<Dna, TChunks>
{
public:
    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            uint32_t possible = seqan::length(text) - this->kmerSize + 1;

            std::vector<uint64_t> kmerHashes(possible);

            this->hashInit(begin(text));
            auto it = begin(text);

            for (uint32_t i = 0; i < possible; ++i)
            {
                uint64_t kmerHash = this->hashNext(it);

                kmerHash = (kmerHash >> this->significantBits) + this->chunkMap[(kmerHash & (this->effectiveChunks - 1))] * this->chunkOffset;
                kmerHashes[i] = kmerHash;

                ++it;
            }

            return kmerHashes;
        }
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        }
        else
        {
            uint32_t possible = seqan::length(text) - this->kmerSize + 1;

            std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes(possible);

            this->hashInit(begin(text));
            auto it = begin(text);

            for (uint32_t i = 0; i < possible; ++i)
            {
                uint64_t kmerHash = this->hashNext(it);

                kmerHashes[i] = std::make_tuple(kmerHash >> this->significantBits, this->chunkMap[(kmerHash & (this->effectiveChunks - 1))]);

                ++it;
            }

            return kmerHashes;
        }
    }
};

template<typename TValue, uint16_t k, typename TChunks>
struct BDHash<TValue, Normal<k>, TChunks> : BDHashBase<TValue, TChunks>
{
    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            uint32_t possible = seqan::length(text) - this->kmerSize + 1;

            std::vector<uint64_t> kmerHashes(possible);

            Shape<TValue, SimpleShape> chunkShape;
            seqan::resize(chunkShape, this->significantPositions);
            uint16_t cacheKmerSize = this->kmerSize;
            this->resize(this->kmerSize - this->significantPositions);
            auto it = begin(text);
            this->hashInit(it);
            auto itChunk = begin(text) + this->kmerSize;
            if (this->significantPositions > 1)
                seqan::hashInit(chunkShape, itChunk);

            for (uint32_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = this->hashNext(it) + this->chunkMap[seqan::hashNext(chunkShape, itChunk)] * this->chunkOffset;
                ++it;
                ++itChunk;
            }
            this->resize(cacheKmerSize);
            return kmerHashes;
        }
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        }
        else
        {
            uint32_t possible = seqan::length(text) - this->kmerSize + 1;

            std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes(possible);

            Shape<TValue, SimpleShape> chunkShape;
            seqan::resize(chunkShape, this->significantPositions);
            uint16_t cacheKmerSize = this->kmerSize;
            this->resize(this->kmerSize - this->significantPositions);
            auto it = begin(text);
            this->hashInit(it);
            auto itChunk = begin(text) + this->kmerSize;
            if (this->significantPositions > 1)
                seqan::hashInit(chunkShape, itChunk);

            for (uint32_t i = 0; i < possible; ++i)
            {
                kmerHashes[i] = std::make_tuple(this->hashNext(it), this->chunkMap[seqan::hashNext(chunkShape, itChunk)]);
                ++it;
                ++itChunk;
            }
            this->resize(cacheKmerSize);
            return kmerHashes;
        }
    }
};

template<uint16_t k, uint16_t o, typename TChunks>
struct BDHash<Dna, Offset<k, o>, TChunks> : BDHashBase<Dna, TChunks>
{
public:
    uint16_t offset{o};

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            // how many test positions are left if we take every offset'th kmer
            uint16_t x = (seqan::length(text) - this->kmerSize) % offset;
            // how many kmers are there when we take every offset'th kmer
            // possible = something left (1/0) + how many fit in the text
            // if something is left, we add a kmer that covers these positions
            uint32_t possible = bool(x) + (seqan::length(text) - this->kmerSize + offset - x) / offset;

            std::vector<uint64_t> kmerHashes(possible);

            this->hashInit(begin(text));
            auto it = begin(text);

            uint32_t positions = seqan::length(text) - this->kmerSize + 1;

            for (uint32_t i = 0, j = 0; i < positions; ++i)
            {
                // std::cerr << "OFFSET\n";
                uint64_t kmerHash = this->hashNext(it);
                if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
                {
                    kmerHash = (kmerHash >> this->significantBits) + this->chunkMap[(kmerHash & (this->effectiveChunks - 1))] * this->chunkOffset;
                    kmerHashes[j] = kmerHash;
                    break;
                }
                if (i - j * offset == 0) // we found the j'th kmer with offset
                {
                    kmerHash = (kmerHash >> this->significantBits) + this->chunkMap[(kmerHash & (this->effectiveChunks - 1))] * this->chunkOffset;
                    kmerHashes[j] = kmerHash;
                    ++j;
                }
                ++it;
            }

            return kmerHashes;
        }
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        }
        else
        {
            // how many test positions are left if we take every offset'th kmer
            uint16_t x = (seqan::length(text) - this->kmerSize) % offset;
            // how many kmers are there when we take every offset'th kmer
            // possible = something left (1/0) + how many fit in the text
            // if something is left, we add a kmer that covers these positions
            uint32_t possible = bool(x) + (seqan::length(text) - this->kmerSize + offset - x) / offset;

            std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes(possible);

            this->hashInit(begin(text));
            auto it = begin(text);

            uint32_t positions = seqan::length(text) - this->kmerSize + 1;

            for (uint32_t i = 0, j = 0; i < positions; ++i)
            {
                // std::cerr << "OFFSET\n";
                uint64_t kmerHash = this->hashNext(it);
                if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
                {
                    kmerHashes[j] = std::make_tuple(kmerHash >> this->significantBits, this->chunkMap[(kmerHash & (this->effectiveChunks - 1))]);
                    break;
                }
                if (i - j * offset == 0) // we found the j'th kmer with offset
                {
                    kmerHashes[j] = std::make_tuple(kmerHash >> this->significantBits, this->chunkMap[(kmerHash & (this->effectiveChunks - 1))]);
                    ++j;
                }
                ++it;
            }

            return kmerHashes;
        }
    }
};

template<typename TValue, uint16_t k, uint16_t o, typename TChunks>
struct BDHash<TValue, Offset<k, o>, TChunks> : BDHashBase<TValue, TChunks>
{
    uint16_t offset{o};

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<uint64_t> {};
        }
        else
        {
            // how many test positions are left if we take every offset'th kmer
            uint16_t x = (seqan::length(text) - this->kmerSize) % offset;
            // how many kmers are there when we take every offset'th kmer
            // possible = something left (1/0) + how many fit in the text
            // if something is left, we add a kmer that covers these positions
            uint32_t possible = bool(x) + (seqan::length(text) - this->kmerSize + offset - x) / offset;
            uint32_t positions = seqan::length(text) - this->kmerSize + 1;

            std::vector<uint64_t> kmerHashes(possible);

            Shape<TValue, SimpleShape> chunkShape;
            seqan::resize(chunkShape, this->significantPositions);
            uint16_t cacheKmerSize = this->kmerSize;
            this->resize(this->kmerSize - this->significantPositions);
            auto it = begin(text);
            this->hashInit(it);
            auto itChunk = begin(text) + this->kmerSize;
            if (this->significantPositions > 1)
                seqan::hashInit(chunkShape, itChunk);

            for (uint32_t i = 0, j = 0; i < positions; ++i)
            {
                // std::cerr << "OFFSET\n";
                uint64_t kmerHash = this->hashNext(it);
                uint16_t  chunkIdentifier = seqan::hashNext(chunkShape, itChunk);
                if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
                {
                    kmerHashes[j] = kmerHash + this->chunkMap[chunkIdentifier] * this->chunkOffset;
                    break;
                }
                if (i - j * offset == 0) // we found the j'th kmer with offset
                {
                    kmerHashes[j] = kmerHash + this->chunkMap[chunkIdentifier] * this->chunkOffset;
                    ++j;
                }
                ++it;
                ++itChunk;
            }

            return kmerHashes;
        }
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
    {
        if (this->kmerSize > seqan::length(text))
        {
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        }
        else
        {
            // how many test positions are left if we take every offset'th kmer
            uint16_t x = (seqan::length(text) - this->kmerSize) % offset;
            // how many kmers are there when we take every offset'th kmer
            // possible = something left (1/0) + how many fit in the text
            // if something is left, we add a kmer that covers these positions
            uint32_t possible = bool(x) + (seqan::length(text) - this->kmerSize + offset - x) / offset;
            uint32_t positions = seqan::length(text) - this->kmerSize + 1;

            std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes(possible);

            Shape<TValue, SimpleShape> chunkShape;
            seqan::resize(chunkShape, this->significantPositions);
            uint16_t cacheKmerSize = this->kmerSize;
            this->resize(this->kmerSize - this->significantPositions);
            auto it = begin(text);
            this->hashInit(it);
            auto itChunk = begin(text) + this->kmerSize;
            if (this->significantPositions > 1)
                seqan::hashInit(chunkShape, itChunk);

            for (uint32_t i = 0, j = 0; i < positions; ++i)
            {
                // std::cerr << "OFFSET\n";
                uint64_t kmerHash = this->hashNext(it);
                uint16_t  chunkIdentifier = seqan::hashNext(chunkShape, itChunk);
                if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
                {
                    kmerHashes[j] = std::make_tuple(kmerHash, this->chunkMap[chunkIdentifier]);
                    break;
                }
                if (i - j * offset == 0) // we found the j'th kmer with offset
                {
                    kmerHashes[j] = std::make_tuple(kmerHash, this->chunkMap[chunkIdentifier]);
                    ++j;
                }
                ++it;
                ++itChunk;
            }

            return kmerHashes;
        }
    }
};

template<uint16_t k, uint32_t w, typename TChunks>
struct BDHash<Dna, Minimizer<k, w>, TChunks> :BDHashBase<Dna, TChunks>
{
public:

    // All positions are inclusive, i.e. kmer starts at b and ends in e => [b,e]
    std::vector<uint64_t> minBegin;
    std::vector<uint64_t> minEnd;
    std::vector<uint32_t> coverage;
    std::vector<uint64_t> coverageBegin;
    std::vector<uint64_t> coverageEnd;
    uint32_t windowSize;
    decltype(BDHash::kmerShape) revCompShape{BDHash::kmerShape};

    inline void resize(TKmerSize newKmerSize = k, uint32_t newWindowSize = w)
    {
        this->kmerSize = newKmerSize;
        windowSize = newWindowSize;
        seqan::resize(this->kmerShape, this->kmerSize);
        seqan::resize(revCompShape, this->kmerSize);
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

    inline uint64_t hash(uint64_t & h)
    {
        return (h >> this->significantBits) + this->chunkMap[(h & (this->effectiveChunks - 1))] * this->chunkOffset;
    }

    inline std::tuple<uint64_t, uint8_t> rawHash(uint64_t & h)
    {
        return std::make_tuple(h >> this->significantBits, this->chunkMap[(h & (this->effectiveChunks - 1))]);
    }

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString & text) // TODO cannot be const for ModifiedString
    {
        if (this->kmerSize > seqan::length(text))
            return std::vector<uint64_t> {};
        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) > windowSize ? seqan::length(text) - windowSize + 1 : 1;
        uint32_t windowKmers = windowSize - this->kmerSize + 1;

        std::vector<uint64_t> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        this->hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(hash(std::get<0>(max)));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        for (uint32_t i = 1; i < possible; ++i)
        {
            windowValues.pop_front();
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
            kmerHashes.push_back(hash(std::get<0>(max)));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }
        return kmerHashes;
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString & text)
    {
        if (this->kmerSize > seqan::length(text))
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) > windowSize ? seqan::length(text) - windowSize + 1 : 1;
        uint32_t windowKmers = windowSize - this->kmerSize + 1;

        std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        this->hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(rawHash(std::get<0>(max)));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        for (uint32_t i = 1; i < possible; ++i)
        {
            windowValues.pop_front();
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
            kmerHashes.push_back(rawHash(std::get<0>(max)));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }
        return kmerHashes;
    }

    inline uint32_t maxCoverage()
    {
        get_coverage();
        return *std::max_element(coverage.begin(), coverage.end());
    }

    inline void get_coverage()
    {
        uint64_t bIndex{1};
        uint64_t eIndex{0};
        auto uMinEnd = minEnd;
        uMinEnd.erase(unique(uMinEnd.begin(), uMinEnd.end()), uMinEnd.end());
        auto uMinBegin = minBegin;
        uMinBegin.erase(unique(uMinBegin.begin(), uMinBegin.end()), uMinBegin.end());
        coverageBegin.push_back(uMinBegin[0]);
        coverage.push_back(1);

        while ((bIndex < uMinBegin.size() ) || (eIndex < uMinEnd.size()))
        {
            uint64_t begin = bIndex < uMinBegin.size() ? uMinBegin[bIndex] : 0xFFFFFFFFFFFFFFFFULL;
            uint64_t end   = uMinEnd[eIndex];
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
                while (uMinBegin[bIndex] == uMinEnd[eIndex])
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

template<typename TValue, uint16_t k, uint32_t w, typename TChunks>
struct BDHash<TValue, Minimizer<k, w>, TChunks> :BDHashBase<TValue, TChunks>
{
public:

    // All positions are inclusive, i.e. kmer starts at b and ends in e => [b,e]
    std::vector<uint64_t> minBegin;
    std::vector<uint64_t> minEnd;
    std::vector<uint32_t> coverage;
    std::vector<uint64_t> coverageBegin;
    std::vector<uint64_t> coverageEnd;
    uint32_t windowSize;
    decltype(BDHash::kmerShape) revCompShape{BDHash::kmerShape};

    inline void resize(TKmerSize newKmerSize = k, uint32_t newWindowSize = w)
    {
        this->kmerSize = newKmerSize;
        windowSize = newWindowSize;
        seqan::resize(this->kmerShape, this->kmerSize);
        seqan::resize(revCompShape, this->kmerSize);
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

    inline uint64_t hash(std::tuple<uint64_t, uint8_t> & h)
    {
        return (std::get<0>(h) + this->chunkMap[std::get<1>(h)] * this->chunkOffset);
    }

    inline std::tuple<uint64_t, uint8_t> rawHash(std::tuple<uint64_t, uint8_t> & h)
    {
        return std::make_tuple(std::get<0>(h), this->chunkMap[std::get<1>(h)]);
    }

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString & text) // TODO cannot be const for ModifiedString
    {
        if (this->kmerSize > seqan::length(text))
            return std::vector<uint64_t> {};
        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) > windowSize ? seqan::length(text) - windowSize + 1 : 1;
        uint32_t windowKmers = windowSize - this->kmerSize + 1;

        decltype(BDHash::kmerShape) chunkShape;
        decltype(BDHash::kmerShape) revChunkShape;
        seqan::resize(chunkShape, this->significantPositions);
        seqan::resize(revChunkShape, this->significantPositions);

        uint16_t cacheKmerSize = this->kmerSize;
        this->resize(this->kmerSize - this->significantPositions, windowSize);

        std::vector<uint64_t> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        this->hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);
        auto itChunk = begin(text) + this->kmerSize;
        auto itRevChunk = begin(revComp) + this->kmerSize;
        if (this->significantPositions > 1)
            seqan::hashInit(chunkShape, itChunk);
        if (this->significantPositions > 1)
            seqan::hashInit(revChunkShape, itRevChunk);
        std::deque<std::tuple<std::tuple<uint64_t, uint8_t>, uint64_t, uint64_t>> windowValues;

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(std::make_tuple(kmerHash, hashNext(itChunk)), distance, distance + this->kmerSize - 1));
                hashNext(itRevChunk);
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(std::make_tuple(revcHash, hashNext(itRevChunk)), distance, distance + this->kmerSize - 1));
                hashNext(itChunk);
            }
            ++it;
            ++rcit;
            ++itChunk;
            ++itRevChunk;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(hash(std::get<0>(max)));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        for (uint32_t i = 1; i < possible; ++i)
        {
            windowValues.pop_front();
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(std::make_tuple(kmerHash, hashNext(itChunk)), distance, distance + this->kmerSize - 1));
                hashNext(itRevChunk);
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(std::make_tuple(revcHash, hashNext(itRevChunk)), distance, distance + this->kmerSize - 1));
                hashNext(itChunk);
            }
            ++it;
            ++rcit;
            ++itChunk;
            ++itRevChunk;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
            kmerHashes.push_back(hash(std::get<0>(max)));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }
        this->resize(cacheKmerSize);
        return kmerHashes;
    }

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString & text)
    {
        if (this->kmerSize > seqan::length(text))
            return std::vector<std::tuple<uint64_t, uint8_t>> {};
        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) > windowSize ? seqan::length(text) - windowSize + 1 : 1;
        uint32_t windowKmers = windowSize - this->kmerSize + 1;

        decltype(BDHash::kmerShape) chunkShape;
        decltype(BDHash::kmerShape) revChunkShape;
        seqan::resize(chunkShape, this->significantPositions);
        seqan::resize(revChunkShape, this->significantPositions);

        uint16_t cacheKmerSize = this->kmerSize;
        this->resize(this->kmerSize - this->significantPositions, windowSize);

        std::vector<std::tuple<uint64_t, uint8_t>> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        this->hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);
        auto itChunk = begin(text) + this->kmerSize;
        auto itRevChunk = begin(revComp) + this->kmerSize;
        if (this->significantPositions > 1)
            seqan::hashInit(chunkShape, itChunk);
        if (this->significantPositions > 1)
            seqan::hashInit(revChunkShape, itRevChunk);
        std::deque<std::tuple<std::tuple<uint64_t, uint8_t>, uint64_t, uint64_t>> windowValues;

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(std::make_tuple(kmerHash, hashNext(chunkShape, itChunk)), distance, distance + this->kmerSize - 1));
                hashNext(revChunkShape, itRevChunk);
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(std::make_tuple(revcHash, hashNext(revChunkShape, itRevChunk)), distance, distance + this->kmerSize - 1));
                hashNext(chunkShape, itChunk);
            }
            ++it;
            ++rcit;
            ++itChunk;
            ++itRevChunk;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(rawHash(std::get<0>(max)));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        for (uint32_t i = 1; i < possible; ++i)
        {
            windowValues.pop_front();
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(std::make_tuple(kmerHash, hashNext(chunkShape, itChunk)), distance, distance + this->kmerSize - 1));
                hashNext(revChunkShape, itRevChunk);
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(std::make_tuple(revcHash, hashNext(revChunkShape, itRevChunk)), distance, distance + this->kmerSize - 1));
                hashNext(chunkShape, itChunk);
            }
            ++it;
            ++rcit;
            ++itChunk;
            ++itRevChunk;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues));
            kmerHashes.push_back(rawHash(std::get<0>(max)));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }
        this->resize(cacheKmerSize);
        return kmerHashes;
    }

    inline uint32_t maxCoverage()
    {
        get_coverage();
        return *std::max_element(coverage.begin(), coverage.end());
    }

    inline void get_coverage()
    {
        uint64_t bIndex{1};
        uint64_t eIndex{0};
        auto uMinEnd = minEnd;
        uMinEnd.erase(unique(uMinEnd.begin(), uMinEnd.end()), uMinEnd.end());
        auto uMinBegin = minBegin;
        uMinBegin.erase(unique(uMinBegin.begin(), uMinBegin.end()), uMinBegin.end());
        coverageBegin.push_back(uMinBegin[0]);
        coverage.push_back(1);

        while ((bIndex < uMinBegin.size() ) || (eIndex < uMinEnd.size()))
        {
            uint64_t begin = bIndex < uMinBegin.size() ? uMinBegin[bIndex] : 0xFFFFFFFFFFFFFFFFULL;
            uint64_t end   = uMinEnd[eIndex];
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
                while (uMinBegin[bIndex] == uMinEnd[eIndex])
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
