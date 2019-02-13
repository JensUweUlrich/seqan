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

#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSEDDISK_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSEDDISK_H_

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
struct Bitvector<UncompressedDisk> : BitvectorBase
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
    int16_t currentChunk;

    std::vector<std::unique_ptr<sdsl::bit_vector> > filterVector;

    double size_in_mega_bytes()
    {
        double size{0};
        for (uint8_t j = 0; j < noOfChunks; j++)
        {
            decompress(j);
            size += sdsl::size_in_mega_bytes(*filterVector[j]);
        }
        return size;
    }

    inline void decompress(uint8_t chunk)
    {
        if (currentChunk < 0 || currentChunk != static_cast<int16_t>(chunk))
        {
            compress(currentChunk);
            sdsl::load_from_file(*filterVector[chunk], toCString(PREFIX)+std::to_string(chunk));
            currentChunk = chunk;
        }
    }

    inline void compress(uint8_t chunk)
    {
        if (currentChunk == static_cast<int16_t>(chunk))
        {
            sdsl::store_to_file(*filterVector[chunk], toCString(PREFIX)+std::to_string(chunk));
            filterVector[chunk] = std::make_unique<sdsl::bit_vector>(0,0);
            currentChunk = -1;
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
        for (uint8_t i = 0; i < noOfChunks; ++i)
        {
            currentChunk = i;
            if (i == noOfChunks - 1)
            {
                filterVector.emplace_back(std::make_unique<sdsl::bit_vector>(chunkSize + FILTER_METADATA_SIZE, 0));
            }
            else
            {
                filterVector.emplace_back(std::make_unique<sdsl::bit_vector>(chunkSize, 0));
            }
            compress(i);
        }
    }

    Bitvector<UncompressedDisk> & operator=(Bitvector<UncompressedDisk> & other)
    {
        for (const auto& element : other.filterVector)
            filterVector.emplace_back(std::make_unique<sdsl::bit_vector>(*element));

        noOfBins = other.noOfBins;
        noOfBits = other.noOfBits;
        binWidth = other.binWidth;
        blockBitSize = other.blockBitSize;
        noOfBlocks = other.noOfBlocks;
        chunkSize = other.chunkSize;
        noOfChunks = other.noOfChunks;
        currentChunk = other.currentChunk;
        for (uint8_t j = 0; j < noOfChunks; j++)
        {
            filesystem::copy_file(toCString(other.PREFIX)+std::to_string(j), toCString(PREFIX)+std::to_string(j), filesystem::copy_options::overwrite_existing);
        }

        return *this;
    }

    Bitvector<UncompressedDisk> & operator=(Bitvector<UncompressedDisk> && other)
    {
        filterVector = std::move(other.filterVector);
        noOfBins = std::move(other.noOfBins);
        noOfBits = std::move(other.noOfBits);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);
        chunkSize = std::move(other.chunkSize);
        noOfChunks = std::move(other.noOfChunks);
        currentChunk = std::move(other.currentChunk);
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
            auto tmp = sdsl::bit_vector(0,0);
            if (sdsl::load_from_file(tmp, toCString(fileName)+std::string(".chunk_")+std::to_string(chunk)))
            {
                filterVector.emplace_back(std::make_unique<sdsl::bit_vector>(std::move(tmp)));
                currentChunk = chunk;
                chunkSize = filterVector[currentChunk]->size();
                compress(chunk);
                ++chunk;
            }
            else
            {
                currentChunk--;
                break;
            }
        }
        noOfChunks = chunk;
        chunkSize -= FILTER_METADATA_SIZE;
        noOfBits = chunkSize * chunk;
        noOfBins = get_int(noOfBits, noOfChunks - 1);
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits  / blockBitSize;
    }

    uint64_t get_int(uint64_t idx, uint8_t chunk, uint8_t len = 64) //const
    {
        uint8_t chunkNo = idx / (chunkSize + (chunk == noOfChunks - 1 ? FILTER_METADATA_SIZE : 0));
        if (chunkNo != chunk)
            return 0;
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        decompress(chunkNo);
        return filterVector[chunkNo]->get_int(chunkPos, len);
    }

    uint64_t get_pos(uint64_t idx, uint8_t chunk) //const
    {
        uint8_t chunkNo = idx / (chunkSize + (chunk == noOfChunks - 1 ? FILTER_METADATA_SIZE : 0));
        if (chunkNo != chunk)
            return 0;
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        decompress(chunkNo);
        return (*filterVector[chunkNo])[chunkPos];
    }

    void set_int(uint64_t idx, uint64_t val, uint8_t chunk, uint8_t len = 64)
    {
        uint8_t chunkNo = idx / (chunkSize + (chunk == noOfChunks - 1 ? FILTER_METADATA_SIZE : 0));
        if (chunkNo != chunk)
            return;
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        decompress(chunkNo);
        filterVector[chunkNo]->set_int(chunkPos, val, len);
    }

    void set_pos(uint64_t idx, uint8_t chunk)
    {
        uint8_t chunkNo = idx / (chunkSize + (chunk == noOfChunks - 1 ? FILTER_METADATA_SIZE : 0));
        if (chunkNo != chunk)
            return;
        decompress(chunkNo);
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        (*filterVector[chunkNo])[chunkPos] = true;
    }

    void unset_pos(uint64_t idx, uint8_t chunk)
    {
        uint8_t chunkNo = idx / (chunkSize + (chunk == noOfChunks - 1 ? FILTER_METADATA_SIZE : 0));
        if (chunkNo != chunk)
            return;
        decompress(chunkNo);
        uint64_t chunkPos = idx - chunkNo * chunkSize;
        (*filterVector[chunkNo])[chunkPos] = false;
    }

    bool store(CharString fileName)
    {
        bool res = true;
        for (uint8_t chunk = 0; chunk < noOfChunks; ++chunk)
        {
            res && filesystem::copy_file(toCString(PREFIX)+std::to_string(chunk), toCString(fileName)+std::string(".chunk_")+std::to_string(chunk), filesystem::copy_options::overwrite_existing);
        }
        return res;
    }

    void retrieve(CharString fileName)
    {
        *this = Bitvector(fileName);
    }
};

}

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSEDDISK_H_
