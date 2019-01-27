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

#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_COMPRESSED_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_COMPRESSED_H_

namespace seqan {

template<>
struct Bitvector<Compressed> : BitvectorBase
{
    CharString PREFIX{random_string()};
    TNoOfBins noOfBins;
    TNoOfBits noOfBits;
    TBinWidth binWidth;
    TBlockBitSize blockBitSize;
    TNoOfBlocks noOfBlocks;

    std::unique_ptr<sdsl::bit_vector> uncompressed_vector;
    std::unique_ptr<sdsl::sd_vector<> > compressed_vector;

    double size_in_mega_bytes() const
    {
        return sdsl::size_in_mega_bytes(*compressed_vector);
    }

    inline void decompress(uint8_t = 0)
    {
        if (compressed)
        {
            sdsl::load_from_file(*uncompressed_vector, toCString(PREFIX));
            compressed = false;
        }
    }

    inline void compress(uint8_t = 0)
    {
        if (!compressed)
        {
            compressed_vector = std::make_unique<sdsl::sd_vector<> >(*uncompressed_vector);
            sdsl::store_to_file(*uncompressed_vector, toCString(PREFIX));
            uncompressed_vector = std::make_unique<sdsl::bit_vector>(0,0);
            compressed = true;
        }
    }

    Bitvector() {}

    Bitvector(uint32_t bins, uint64_t bits, uint8_t = 0):
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
        compressed_vector = std::make_unique<sdsl::sd_vector<> >();
    }

    Bitvector & operator=(Bitvector & other)
    {
        uncompressed_vector = std::make_unique<sdsl::bit_vector>(*other.uncompressed_vector);
        compressed_vector = std::make_unique<sdsl::sd_vector<> >(*other.compressed_vector);
        compressed = other.compressed;
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
        compressed_vector = std::move(other.compressed_vector);
        compressed = std::move(other.compressed);
        noOfBins = std::move(other.noOfBins);
        noOfBits = std::move(other.noOfBits);
        binWidth = std::move(other.binWidth);
        blockBitSize = std::move(other.blockBitSize);
        noOfBlocks = std::move(other.noOfBlocks);

        return *this;
    }

    ~Bitvector()
    {
        std::remove(toCString(PREFIX));
    }

    Bitvector(CharString fileName)
    {
        /*
        uncompressed_vector = std::make_unique<sdsl::bit_vector>(0,0);
        sdsl::load_from_file(*uncompressed_vector, toCString(fileName));

        noOfBits = uncompressed_vector->size();
        */
        compressed_vector = std::make_unique<sdsl::sd_vector<>>();
        sdsl::load_from_file(*compressed_vector, toCString(fileName));

        noOfBits = compressed_vector->size();
        noOfBits -= FILTER_METADATA_SIZE;
        // How many blocks of 64 bit do we need to represent our noOfBins
        binWidth = std::ceil((double)noOfBins / INT_SIZE);
        // How big is then a block (multiple of 64 bit)
        blockBitSize = binWidth * INT_SIZE;
        // How many hash values can we represent
        noOfBlocks = noOfBits / blockBitSize;
        compressed = true;
    }

    inline uint64_t get_int(uint64_t idx, uint8_t = 0, uint8_t len = 64) const
    {
        if (compressed)
            return compressed_vector->get_int(idx, len);
        else
            return uncompressed_vector->get_int(idx, len);
    }

    inline uint64_t get_pos(uint64_t vecIndex, uint8_t = 0) const
    {
        if (compressed)
            return (*compressed_vector)[vecIndex];
        else
            return (*uncompressed_vector)[vecIndex];
    }

    void set_int(uint64_t idx, uint64_t val, uint8_t = 0, uint8_t len = 64)
    {
        decompress();
        uncompressed_vector->set_int(idx, val, len);
    }

    void set_pos(uint64_t idx, uint8_t = 0)
    {
        decompress();
        (*uncompressed_vector)[idx] = true;
    }

    void unset_pos(uint64_t idx, uint8_t = 0)
    {
        decompress();
        (*uncompressed_vector)[idx] = false;
    }

    bool store(CharString fileName)
    {
        compress();
        return sdsl::store_to_file(*compressed_vector, toCString(fileName));
        // bool store_compressed = sdsl::store_to_file(*compressed_vector, toCString(fileName));
        // decompress();
        // return store_compressed && sdsl::store_to_file(*uncompressed_vector, toCString(fileName)+std::string(".decomp"));
    }

    void retrieve(CharString fileName)
    {
        *this = Bitvector(fileName);
    }
};
}   // namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BITVECTOR_COMPRESSED_H_
