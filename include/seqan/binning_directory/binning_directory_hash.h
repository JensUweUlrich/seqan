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

template<typename TChunks>
struct BDHash<Dna, Normal, TChunks> : BDHashBase<Dna, TChunks>
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

template<typename TValue, typename TChunks>
struct BDHash<TValue, Normal, TChunks> : BDHashBase<TValue, TChunks>
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

template<uint16_t o, typename TChunks>
struct BDHash<Dna, Offset<o>, TChunks> : BDHashBase<Dna, TChunks>
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

template<typename TValue, uint16_t o, typename TChunks>
struct BDHash<TValue, Offset<o>, TChunks> : BDHashBase<TValue, TChunks>
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

inline bool windowCompare(std::tuple<uint64_t, uint64_t, uint64_t> const & a,
                          std::tuple<uint64_t, uint64_t, uint64_t> const & b)
{
    return (std::get<0>(a) < std::get<0>(b));
}

template<uint16_t k, uint32_t w, typename TChunks>
struct BDHash<Dna, Minimizer<k, w>, TChunks> :BDHashBase<Dna, TChunks>
{
public:

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

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString & text) // TODO cannot be const for ModifiedString
    {
        if (this->kmerSize > seqan::length(text))
            return std::vector<uint64_t> {};

        typedef ModifiedString<ModifiedString<TString, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);
        uint32_t possible = seqan::length(text) - windowSize + 1;
        uint32_t windowKmers = windowSize - this->kmerSize + 1;

        std::vector<uint64_t> kmerHashes;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        this->hashInit(begin(text));
        revHashInit(begin(revComp));
        auto it = begin(text);
        auto rcit = begin(revComp);
        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;
        windowValues.reserve(windowKmers);

        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(hash(kmerHash), distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(hash(revcHash), distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;
        }

        auto max = *std::min_element(std::begin(windowValues), std::end(windowValues), windowCompare);
        kmerHashes.push_back(std::get<0>(max));
        minBegin.push_back(std::get<1>(max));
        minEnd.push_back(std::get<2>(max));

        for (uint32_t i = 1; i < possible; ++i)
        {
            windowValues.erase(std::begin(windowValues));
            uint64_t kmerHash = this->hashNext(it);
            uint64_t revcHash = revHashNext(rcit);
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(hash(kmerHash), distance, distance + this->kmerSize - 1));
            }
            else
            {
                uint64_t distance = std::distance(begin(revComp), rcit);
                windowValues.push_back(std::make_tuple(hash(revcHash), distance, distance + this->kmerSize - 1));
            }
            ++it;
            ++rcit;

            auto max = *std::min_element(std::begin(windowValues), std::end(windowValues), windowCompare);
            kmerHashes.push_back(std::get<0>(max));
            minBegin.push_back(std::get<1>(max));
            minEnd.push_back(std::get<2>(max));
        }
        // DEBUG
        // std::cerr << "inserted:\n";
        // std::cerr << "Hashes: " << kmerHashes.size() << '\n';
        // std::cerr << "Begin\tEnd\n";
        // for (uint64_t i = 0; i < kmerHashes.size(); ++i)
        // {
        //     std::cerr << minBegin[i] << '\t' << minEnd[i] << '\n';
        // }
        get_coverage();
        std::cerr << "Begin\tEnd\tCoverage\n";
        for (uint64_t i = 0; i < coverageBegin.size(); ++i)
        {
            std::cerr << coverageBegin[i] << '\t' << coverageEnd[i] << '\t' << coverage[i] << '\n';
        }
        return kmerHashes;
    }

    inline uint32_t maxCoverage()
    {
        get_coverage();
        return std::max_element(coverage.begin(), coverage.end());
    }

    inline void get_coverage()
    {
        // Index of minBegin
        uint64_t bIndex{1};
        // Index of minEnd
        uint64_t eIndex{0};
        auto uMinEnd = minEnd;
        uMinEnd.erase( unique( uMinEnd.begin(), uMinEnd.end() ), uMinEnd.end() );
        auto uMinBegin = minBegin;
        uMinBegin.erase( unique( uMinBegin.begin(), uMinBegin.end() ), uMinBegin.end() );
        for (uint64_t i = 0; i < uMinBegin.size(); ++i)
        {
            std::cerr << uMinBegin[i] << '\t' << uMinEnd[i] << '\n';
        }
        coverageBegin.push_back(uMinBegin[0]);
        coverage.push_back(1);
        bool skipped{false};

        while ((bIndex < uMinBegin.size() ) || (eIndex < uMinEnd.size()))
        {
            if (skipped)
            {
                skipped = false;
                // std::cerr << "skip pushes " << uMinEnd[eIndex] << '\n';
                // coverageEnd.push_back(uMinEnd[eIndex]);
                ++eIndex;
                ++coverage.back();
                continue;
            }
            // Overlapping kmers
            if (bIndex < uMinBegin.size() && uMinBegin[bIndex] < uMinEnd[eIndex])
            {
                std::cerr << "Case 1\n";

                std::cerr << "1 pushes " << uMinBegin[bIndex] << '\n';
                coverageEnd.push_back(uMinBegin[bIndex]);
                coverageBegin.push_back(uMinBegin[bIndex]);
                coverage.push_back(coverage.back()+1);
                ++bIndex;
                continue;
            }
            // kmer ends
            if ((eIndex < uMinEnd.size() && uMinBegin[bIndex] > uMinEnd[eIndex]) || (bIndex >= uMinBegin.size()))
            {
                std::cerr << "Case 2\n";
                coverageBegin.push_back(uMinEnd[eIndex]);
                std::cerr << "2 pushes " << uMinEnd[eIndex] << '\n';
                coverageEnd.push_back(uMinEnd[eIndex]);
                coverage.push_back(coverage.back()-1);
                ++eIndex;
                continue;
            }
            // Another kmer on current position
            std::cerr << "Case 3\n";
            std::cerr << "3 pushes " << uMinEnd[eIndex] << '\n';
            coverageEnd.push_back(uMinEnd[eIndex]);
            coverageBegin.push_back(uMinEnd[eIndex]);
            coverage.push_back(coverage.back());
            while (uMinBegin[bIndex] == uMinEnd[eIndex])
            {
                ++eIndex;
                ++bIndex;
            }
            skipped = true;
            continue;
        }
    }
};

}   // namespace seqan







































#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_
