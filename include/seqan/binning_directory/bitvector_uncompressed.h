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

    double size_in_mega_bytes() const
    {
        return sdsl::size_in_mega_bytes(*uncompressed_vector);
    }

    uint64_t size() const
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

    uint64_t get_int(uint64_t idx, uint64_t len = 1ULL<<6) const
    {
        return uncompressed_vector->get_int(idx, len);
    }

    uint64_t get_pos(uint64_t vecIndex) const
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

    bool store(CharString fileName) const
    {
        return sdsl::store_to_file(*uncompressed_vector, toCString(fileName));
    }

    void retrieve(CharString fileName)
    {
        *this = Bitvector(fileName);
    }

    void resize(TNoOfBins bins, TNoOfBits newNoOfBits, TBlockBitSize newBlockBitSize, TBinWidth newBinWidth)
    {
        TNoOfBits idxNew{newNoOfBits}, idxOld{noOfBits};
        TBlockBitSize delta = newBlockBitSize - blockBitSize + 64;

        uncompressed_vector->resize(newNoOfBits + FILTER_METADATA_SIZE);

        // Shift MetaData
        for (TNoOfBits i = idxNew, j = idxOld, c = 0; c < FILTER_METADATA_SIZE / 64; ++c, i += 64, j += 64)
        {
            uint64_t old = get_int(j);
            set_int(j, 0);
            set_int(i, old);
        }
        // Shift Data
        for (TNoOfBits i = idxNew, j = idxOld; j > 0; i -= newBlockBitSize, j -= blockBitSize)
        {
            TNoOfBits stop = i - newBlockBitSize;
            for (TNoOfBits ii = i - delta, jj = j - 64; stop && ii >= stop; ii -= 64, jj -= 64)
            {
                uint64_t old = get_int(jj);
                set_int(jj, 0);
                set_int(ii, old);
            }
        }
        noOfBins = bins;
        binWidth = newBinWidth;
        blockBitSize = newBlockBitSize;
        noOfBits = newNoOfBits;
    }
};
}   //  namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_UNCOMPRESSED_H_
