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

#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSED_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSED_H_

namespace seqan {

template<>
struct Bitvector<Uncompressed> : BitvectorBase
{
    std::unique_ptr<sdsl::bit_vector> uncompressed_vector;
    TNoOfBins noOfBins;
    TNoOfBits noOfBits;
    TBinWidth binWidth;
    TBlockBitSize blockBitSize;
    TNoOfBlocks noOfBlocks;

    double size_in_mega_bytes()
    {
        return sdsl::size_in_mega_bytes(*uncompressed_vector);
    }

    uint64_t size()
    {
        return uncompressed_vector->size();
    }

    Bitvector() {}

    Bitvector(uint32_t bins, uint64_t bits):
        noOfBins(bins),
        noOfBits(bits)
    {
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((float)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits / blockBitSize;

        uncompressed_vector = std::make_unique<sdsl::bit_vector>(noOfBits+FILTER_METADATA_SIZE,0);
    }

    Bitvector & operator=(Bitvector & other)
    {
        uncompressed_vector = std::make_unique<sdsl::bit_vector>(*other.uncompressed_vector);
        noOfBins = other.noOfBins;
        noOfBits = other.noOfBits;
        binWidth = other.binWidth;
        blockBitSize = other.blockBitSize;
        noOfBlocks = other.noOfBlocks;

        return *this;
    }

    Bitvector & operator=(Bitvector && other)
    {
        uncompressed_vector = std::move(other.uncompressed_vector);
        noOfBins = std::move(other.noOfBins);
        noOfBits = std::move(other.noOfBits);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);

        return *this;
    }

    ~Bitvector() = default;

    Bitvector(CharString fileName)
    {
        uncompressed_vector = std::make_unique<sdsl::bit_vector>(0,0);

        sdsl::load_from_file(*uncompressed_vector, toCString(fileName));

        noOfBits = uncompressed_vector->size();
        noOfBits -= FILTER_METADATA_SIZE;
        noOfBins = get_int(noOfBits);
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits / blockBitSize;
    }

    uint64_t get_int(uint64_t idx, uint64_t len = 1ULL<<6)
    {
        return uncompressed_vector->get_int(idx, len);
    }

    uint64_t get_pos(uint64_t vecIndex)
    {
        return (*uncompressed_vector)[vecIndex];
    }

    void set_int(uint64_t idx, uint64_t val)
    {
        uncompressed_vector->set_int(idx, val);
    }

    inline void set_pos(uint64_t idx)
    {
        (*uncompressed_vector)[idx] = true;
    }

    void unset_pos(uint64_t idx)
    {
        (*uncompressed_vector)[idx] = false;
    }

    bool store(CharString fileName)
    {
        return sdsl::store_to_file(*uncompressed_vector, toCString(fileName));
    }

    void retrieve(CharString fileName)
    {
        *this = Bitvector(fileName);
    }

    auto resize(uint32_t bins)
    {
        CharString file{random_string()};
        store(file);
        TBlockBitSize newBlockBitSize = std::ceil((double)bins / INT_SIZE) * INT_SIZE;
        TBlockBitSize delta = newBlockBitSize - blockBitSize + 1;
        if (delta == 1)
        {
            return std::make_tuple(noOfBits, bins, binWidth, blockBitSize);
        }
        TNoOfBits newNoOfBits = noOfBlocks * newBlockBitSize;
        uncompressed_vector.reset(new sdsl::bit_vector(newNoOfBits+FILTER_METADATA_SIZE,0));
        sdsl::int_vector_buffer<1> buffered_vector(toCString(file));
        TNoOfBits pos{0};
        TNoOfBits posBuff{0};
        for (auto it = buffered_vector.begin(); it != buffered_vector.end() && pos != newNoOfBits; ++it)
        {
            if (*it)
            {
              set_pos(pos);
            }
            if (posBuff == blockBitSize -1)
            {
                posBuff = 0;
                pos += delta;
            }
            else
            {
                ++pos;
                ++posBuff;
            }
        }
        sdsl::remove(toCString(file));
        noOfBits = newNoOfBits;
        noOfBins = bins;
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        blockBitSize = newBlockBitSize;
        return std::make_tuple(noOfBits, noOfBins, binWidth, blockBitSize);
    }
};
}   //  namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSED_H_
