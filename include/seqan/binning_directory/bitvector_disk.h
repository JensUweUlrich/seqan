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
#include <random>

#if __has_include(<filesystem>)
#include <filesystem>
#else
#include <experimental/filesystem>
#endif

#if __has_include(<filesystem>)
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif

namespace seqan {

template<>
struct Bitvector<CompressedDisk> : BitvectorBase
{
    // static const uint64_t MAX_VEC = 1ULL<<32; //512 MB, 36 -> 8GB, 39 -> 64

    CharString PREFIX{random_string()};
    TNoOfBins noOfBins;
    TNoOfBits noOfBits;
    TBinWidth binWidth;
    TBlockBitSize blockBitSize;
    TNoOfBlocks noOfBlocks;
    uint8_t noOfChunks;
    uint64_t chunkSize;

    std::vector<std::tuple<bool,std::unique_ptr<sdsl::bit_vector>, std::unique_ptr<sdsl::sd_vector<> > > > filterVector;

    double size_in_mega_bytes()
    {
        double size{0};
        for (uint8_t j = 0; j < noOfChunks; j++)
        {
            size += sdsl::size_in_mega_bytes(*std::get<2>(filterVector[j]));
        }
        return size;
    }

    inline void decompress(uint8_t chunk)
    {
        if (std::get<0>(filterVector[chunk]))
        {
            for (uint8_t c = 0; c < noOfChunks; ++c)
            {
                compress(c);
            }
            sdsl::load_from_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
            std::get<0>(filterVector[chunk]) = false;
        }
    }

    inline void compress(uint8_t chunk)
    {
        if (!std::get<0>(filterVector[chunk]))
        {
            std::get<2>(filterVector[chunk]) = std::make_unique<sdsl::sd_vector<> >(*std::get<1>(filterVector[chunk]));
            sdsl::store_to_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
            std::get<1>(filterVector[chunk]) = std::make_unique<sdsl::bit_vector>(0,0);
            std::get<0>(filterVector[chunk]) = true;
        }
    }

    Bitvector() {}

    Bitvector(uint32_t bins, uint64_t bits, uint8_t chunks):
        noOfBins(bins),
        noOfBits(bits),
        noOfChunks(chunks)
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((float)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits / blockBitSize;

        // If we split, we need to split at the end of a block.
        chunkSize = (noOfBlocks / noOfChunks) * blockBitSize;
        for (uint64_t i = 0; i < noOfChunks; ++i)
        {
            if (i == noOfChunks -1)
            {
                filterVector.emplace_back(std::make_tuple(false,
                                                       std::make_unique<sdsl::bit_vector>(chunkSize + FILTER_METADATA_SIZE, 0),
                                                       std::make_unique<sdsl::sd_vector<> >()));
            }
            else
            {
                filterVector.emplace_back(std::make_tuple(false,
                                                       std::make_unique<sdsl::bit_vector>(chunkSize, 0),
                                                       std::make_unique<sdsl::sd_vector<> >()));
            }
            compress(i);
        }
    }

    Bitvector<CompressedDisk> & operator=(Bitvector<CompressedDisk> & other)
    {
        for (const auto& element : other.filterVector)
            filterVector.emplace_back(std::make_tuple(std::get<0>(element), std::make_unique<sdsl::bit_vector>(*std::get<1>(element)), std::make_unique<sdsl::sd_vector<> >(*std::get<2>(element))));

        noOfBins = other.noOfBins;
        noOfBits = other.noOfBits;
        binWidth = other.binWidth;
        blockBitSize = other.blockBitSize;
        noOfBlocks = other.noOfBlocks;
        chunkSize = other.chunkSize;
        noOfChunks = other.noOfChunks;
        for (uint8_t j = 0; j < noOfChunks; j++)
        {
            filesystem::copy_file(toCString(other.PREFIX)+std::to_string(j), toCString(PREFIX)+std::to_string(j));
        }

        return *this;
    }

    Bitvector<CompressedDisk> & operator=(Bitvector<CompressedDisk> && other)
    {
        filterVector = std::move(other.filterVector);
        noOfBins = std::move(other.noOfBins);
        noOfBits = std::move(other.noOfBits);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);
        chunkSize = std::move(other.chunkSize);
        noOfChunks = std::move(other.noOfChunks);
        PREFIX = std::move(other.PREFIX);

        return *this;
    }

    ~Bitvector()
    {
        for (uint8_t chunk = 0; chunk < noOfChunks; ++chunk)
        {
            std::remove((toCString(PREFIX)+std::to_string(chunk)).c_str());
        }
    }

    Bitvector(CharString fileName)
    {
        uint8_t chunk = 0;
        while (true)
        {
            filterVector.emplace_back(
                std::make_tuple(
                    false,
                    std::make_unique<sdsl::bit_vector>(0,0),
                    std::make_unique<sdsl::sd_vector<> >()));
            if (sdsl::load_from_file(*std::get<1>(filterVector[chunk]), toCString(fileName)+std::to_string(chunk)))
            {
                compress(chunk);
                ++chunk;
            }
            else
            {
                break;
            }
        }
        chunkSize = std::get<2>(filterVector[0])->size();
        noOfChunks = chunk;
        noOfBits = chunkSize;
        for (uint8_t c = 1; c < noOfChunks; ++c)
        {
            noOfBits += std::get<2>(filterVector[c])->size();
        }
        noOfBits -= FILTER_METADATA_SIZE;
        noOfBins = get_int(noOfBits, noOfChunks);
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits  / blockBitSize;
    }

    uint64_t get_int(uint64_t idx, /*uint64_t len = 1ULL<<6,*/ uint8_t chunk)
    {
        uint8_t chunkNo = idx / chunkSize;
        if (chunkNo != chunk)
            return 0;
        if (chunk == noOfChunks)
            --chunkNo;
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        compress(chunkNo);
        return std::get<2>(filterVector[chunkNo])->get_int(chunkPos/*, len*/);
    }

    uint64_t get_pos(uint64_t vecIndex, uint8_t chunk)
    {
        uint8_t chunkNo = vecIndex / chunkSize;
        if (chunkNo != chunk)
            return 0;
        if (chunk == noOfChunks)
            --chunkNo;
        uint64_t chunkPos = vecIndex - chunkNo * chunkSize;
        compress(chunk);
        return (*std::get<2>(filterVector[chunkNo]))[chunkPos];
    }

    void set_int(uint64_t idx, uint64_t val, uint8_t chunk)
    {
        uint8_t chunkNo = idx / chunkSize;
        if (chunkNo != chunk)
            return;
        if (chunk == noOfChunks)
            --chunkNo;
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        decompress(chunkNo);
        (*std::get<1>(filterVector[chunkNo])).set_int(chunkPos, val);
        // compress(chunkNo);
    }

    void set_pos(uint64_t idx, uint8_t chunk)
    {
        uint8_t chunkNo = idx / chunkSize;
        if (chunkNo != chunk)
            return;
        if (chunk == noOfChunks)
            --chunkNo;
        decompress(chunkNo);
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        (*std::get<1>(filterVector[chunkNo]))[chunkPos] = true;
    }

    void unset_pos(uint64_t idx, uint8_t chunk)
    {
        uint8_t chunkNo = idx / chunkSize;
        if (chunkNo != chunk)
            return;
        if (chunk == noOfChunks)
            --chunkNo;
        decompress(chunkNo);
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        (*std::get<1>(filterVector[chunkNo]))[chunkPos] = false;
    }

    bool store(CharString fileName)
    {
        bool res = true;
        for (uint8_t chunk = 0; chunk < noOfChunks; ++chunk)
        {
            decompress(chunk);
            res && sdsl::store_to_file(*std::get<1>(filterVector[chunk]), toCString(fileName)+std::to_string(chunk));
            compress(chunk);
        }
        return res;
    }

    void retrieve(CharString fileName)
    {
        *this = Bitvector(fileName);
    }
};

}
