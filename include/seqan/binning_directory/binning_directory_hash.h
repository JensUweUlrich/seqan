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

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
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
};

template<typename TValue, typename TChunks>
struct BDHash<TValue, Normal, TChunks> : BDHashBase<TValue, TChunks>
{
    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
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

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
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
};

template<uint16_t o, typename TChunks>
struct BDHash<Dna, Offset<o>, TChunks> : BDHashBase<Dna, TChunks>
{
public:
    uint16_t offset{o};

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
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

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
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
};

template<typename TValue, uint16_t o, typename TChunks>
struct BDHash<TValue, Offset<o>, TChunks> : BDHashBase<TValue, TChunks>
{
    uint16_t offset{o};

    template<typename TString>
    inline std::vector<uint64_t> getHash(TString const & text)
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

    template<typename TString>
    inline std::vector<std::tuple<uint64_t, uint8_t>> getRawHash(TString const & text)
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
};

}   // namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_HASH_H_
