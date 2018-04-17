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

using namespace seqan;
struct FilterVector
{
    static const uint64_t FILTER_METADATA_SIZE = 256;
    static const uint64_t BUFFER_SIZE = 10000000;
    static const uint64_t MAX_VEC = 1ULL<<32;
    static const uint64_t INT_SIZE = 0x40;

    std::string random_string()
    {
         std::string str("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");

         std::random_device rd;
         std::mt19937 generator(rd());

         std::shuffle(str.begin(), str.end(), generator);

         return str.substr(0, 32);
    }

    CharString PREFIX{random_string()};

    uint64_t noOfBins;
    uint64_t noOfBits;
    uint64_t binWidth;
    uint64_t blockBitSize;
    uint64_t noOfBlocks;
    uint64_t noOfChunks;
    uint64_t chunkSize;

    std::vector<std::tuple<bool,std::unique_ptr<sdsl::bit_vector>, std::unique_ptr<sdsl::sd_vector<> > > > filterVector;

    double size_in_mega_bytes()
    {
        double size{0};
        for (uint64_t j = 0; j < noOfChunks; j++)
        {
            size += sdsl::size_in_mega_bytes(*std::get<2>(filterVector[j]));
        }
        return size;
    }

    inline void decompress(uint64_t chunk)
    {
        if (std::get<0>(filterVector[chunk]))
        {
            sdsl::load_from_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
            std::get<0>(filterVector[chunk]) = false;
        }
    }

    inline void compress(uint64_t chunk)
    {
        if (!std::get<0>(filterVector[chunk]))
        {
            std::get<2>(filterVector[chunk]) = std::make_unique<sdsl::sd_vector<> >(*std::get<1>(filterVector[chunk]));
            sdsl::store_to_file(*std::get<1>(filterVector[chunk]), toCString(PREFIX)+std::to_string(chunk));
            std::get<1>(filterVector[chunk]) = std::make_unique<sdsl::bit_vector>(0,0);
            std::get<0>(filterVector[chunk]) = true;
        }
    }

    FilterVector() {}

    FilterVector(uint64_t bins, uint64_t bits):
        noOfBins(bins),
        noOfBits(bits)
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((float)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = (noOfBits + FILTER_METADATA_SIZE) / blockBitSize;

        // If we split, we need to split at the end of a block.
        chunkSize = std::min((double)noOfBlocks * blockBitSize, std::ceil((double)MAX_VEC/blockBitSize) * blockBitSize);
        // This is how many chunks we need.
        noOfChunks = std::ceil((double)noOfBits/chunkSize);
        if (chunkSize * noOfChunks < noOfBits + FILTER_METADATA_SIZE)
            ++noOfChunks;
        uint64_t size = noOfBits + FILTER_METADATA_SIZE;
        for (uint64_t i = 0; i < noOfChunks; ++i)
        {
            filterVector.emplace_back(std::make_tuple(false,
                                                   std::make_unique<sdsl::bit_vector>(std::min(chunkSize, size), 0),
                                                   std::make_unique<sdsl::sd_vector<> >()));
            compress(i);
            size -= chunkSize;
        }
    }

    FilterVector & operator=(FilterVector & other)
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

        return *this;
    }

    FilterVector & operator=(FilterVector && other)
    {
        filterVector = std::move(other.filterVector);
        noOfBins = std::move(other.noOfBins);
        noOfBits = std::move(other.noOfBits);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);
        chunkSize = std::move(other.chunkSize);
        noOfChunks = std::move(other.noOfChunks);

        return *this;
    }

    ~FilterVector()
    {
        for (uint64_t chunk = 0; chunk < noOfChunks; ++chunk)
        {
            std::remove((toCString(PREFIX)+std::to_string(chunk)).c_str());
        }
    }

    FilterVector(CharString fileName)
    {
        uint64_t chunk = 0;
        while (true)
        {
            filterVector.emplace_back(
                std::make_tuple(
                    true,
                    std::make_unique<sdsl::bit_vector>(0,0),
                    std::make_unique<sdsl::sd_vector<> >()));
            if (sdsl::load_from_file(*std::get<2>(filterVector[chunk]), toCString(fileName)+std::to_string(chunk)))
            {
                ++chunk;
            }
            else
            {
                break;
            }
        }
        noOfChunks = chunk;
        noOfBits = 0;
        for (uint64_t c = 0; c < noOfChunks; ++c)
        {
            noOfBits += std::get<2>(filterVector[c])->size();
        }
        noOfBits -= FILTER_METADATA_SIZE;
        noOfBins = get_int(noOfBits);
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = (noOfBits + FILTER_METADATA_SIZE) / blockBitSize;

        // We need to split at the end of a block.
        chunkSize = std::min((double)noOfBlocks * blockBitSize, std::ceil((double)MAX_VEC/blockBitSize) * blockBitSize);
    }

    uint64_t get_int(uint64_t idx, uint64_t len = 1ULL<<6)
    {
        uint64_t access = idx;
        uint64_t chunkNo = access / MAX_VEC;
        uint64_t chunkPos = access - chunkNo * MAX_VEC;
        return std::get<2>(filterVector[chunkNo])->get_int(chunkPos, len);
    }

    uint64_t get_pos(uint64_t vecIndex)
    {
        uint64_t access = vecIndex;
        uint64_t chunkNo = access / MAX_VEC;
        uint64_t chunkPos = access - chunkNo * MAX_VEC;
        return (*std::get<2>(filterVector[chunkNo]))[chunkPos];
    }

    void set_int(uint64_t idx, uint64_t val)
    {
        uint64_t access = idx;
        uint64_t chunkNo = access / MAX_VEC;
        uint64_t chunkPos = access - chunkNo * MAX_VEC;
        decompress(chunkNo);
        (*std::get<1>(filterVector[chunkNo])).set_int(chunkPos, val);
        compress(chunkNo);
    }

    void set_pos(uint64_t chunkNo, uint64_t chunkPos)
    {
        (*std::get<1>(filterVector[chunkNo]))[chunkPos] = true;
    }

    void unset_pos(uint64_t chunkNo, uint64_t chunkPos)
    {
        (*std::get<1>(filterVector[chunkNo]))[chunkPos] = false;
    }

    bool store(CharString fileName)
    {
        bool res = true;
        for (uint64_t chunk = 0; chunk < noOfChunks; ++chunk)
        {
            res && sdsl::store_to_file(*std::get<2>(filterVector[chunk]), toCString(fileName)+std::to_string(chunk));
        }
        return res;
    }

    void retrieve(CharString fileName)
    {
        *this = FilterVector(fileName);
    }
};
