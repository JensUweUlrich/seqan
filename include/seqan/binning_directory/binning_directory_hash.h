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

uint64_t ipow(uint64_t base, uint64_t exp)
{
    uint64_t result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

template<typename TValue, uint8_t chunks>
struct BDHash<TValue, Normal, chunks>
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

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        // std::cerr << "NORMAL\n";
        uint32_t possible = seqan::length(text) - kmerSize + 1;

        std::vector<uint64_t> kmerHashes(possible, 0);

        uint8_t significantPos = std::ceil((double) ValueSize<TValue>::VALUE / chunks);

        Shape<TValue, SimpleShape> chunkShape;
        seqan::resize(chunkShape, significantPos);
        uint16_t cacheKmerSize = kmerSize;
        resize(kmerSize - significantPos);
        auto it = begin(text) + significantPos;
        hashInit(it);
        auto itChunk = begin(text);
        seqan::hashInit(chunkShape, itChunk);

        for (uint32_t i = 0; i < possible; ++i)
        {
            uint8_t chunkIdentifier = seqan::hashNext(chunkShape, itChunk);

            uint8_t chunkId = std::ceil((double) chunkIdentifier / chunks);

            uint8_t chunk = std::max(0, chunkId - 1);

            uint64_t chunkOffset = ipow(ValueSize<TValue>::VALUE, kmerSize - chunk) * chunk;
            uint64_t temp = hashNext(it);
            kmerHashes[i] = chunkOffset + temp; // hashNext(it);
            ++it;
            ++itChunk;
        }
        resize(cacheKmerSize);
        return kmerHashes;
    }
};

template<typename TValue, uint16_t o, uint8_t chunks>
struct BDHash<TValue, Offset<o>, chunks>
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

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
    {
        // how many test positions are left if we take every offset'th kmer
        uint16_t x = (seqan::length(text) - kmerSize) % offset;
        // how many kmers are there when we take every offset'th kmer
        // possible = something left (1/0) + how many fit in the text
        // if something is left, we add a kmer that covers these positions
        uint32_t possible = bool(x) + (seqan::length(text) - kmerSize + offset - x) / offset;
        uint32_t positions = seqan::length(text) - kmerSize + 1;

        std::vector<uint64_t> kmerHashes(possible, 0);

        uint8_t significantPos = std::ceil((double) ValueSize<TValue>::VALUE / chunks);

        Shape<TValue, SimpleShape> chunkShape;
        seqan::resize(chunkShape, significantPos);
        uint16_t cacheKmerSize = kmerSize;
        resize(kmerSize - significantPos);
        auto it = begin(text) + significantPos;
        hashInit(it);
        auto itChunk = begin(text);
        seqan::hashInit(chunkShape, itChunk);

        for (uint32_t i = 0, j = 0; i < positions; ++i)
        {
            // std::cerr << "OFFSET\n";
            uint64_t kmerHash = hashNext(it);
            uint8_t  chunkIdentifier = seqan::hashNext(chunkShape, itChunk);
            if (x && i == positions - 1) // we take the last kmer that covers otherwise uncovered positions
            {
                uint8_t chunkId = std::ceil((double) chunkIdentifier / chunks);

                uint8_t chunk = std::max(0, chunkId - 1);

                uint64_t chunkOffset = ipow(ValueSize<TValue>::VALUE, kmerSize - chunk) * chunk;

                kmerHashes[j] = chunkOffset + kmerHash;
                break;
            }
            if (i - j * offset == 0) // we found the j'th kmer with offset
            {
                uint8_t chunkId = std::ceil((double) chunkIdentifier / chunks);

                uint8_t chunk = std::max(0, chunkId - 1);

                uint64_t chunkOffset = ipow(ValueSize<TValue>::VALUE, kmerSize - chunk) * chunk;

                kmerHashes[j] = chunkOffset + kmerHash;
                ++j;
            }
            ++it;
            ++itChunk;
        }

        return kmerHashes;
    }
};

}   // namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_
