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
// Author:      Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#include <algorithm>
#undef NDEBUG
#include <cassert>
#include <stdio.h>
#include <vector>

#include <seqan/seq_io.h>
#include <seqan/kmer.h>

static const uint32_t filterMetadataSize = 256;
static const uint8_t INT_WIDTH = 0x40;

using namespace seqan;

int main()
{
    //SeqAnBloomFilter<> filter (10, 1, 12, 16777216);
    uint64_t threads{3};
    uint64_t noBins{10};
    uint64_t kmerSize{14};
    [[maybe_unused]]
    uint64_t hashFunc{3};
    [[maybe_unused]]
    uint64_t bits{1ULL<<33};


    // ==========================================================================
    // Test constructors
    // ==========================================================================
    std::cout << "Testing ctors" << '\n';

    typedef InterleavedBloomFilter TSpec;
    // typedef DirectAddressing       TSpec;

    // typedef Uncompressed TFilter;
    // typedef CompressedSimple TFilter;
    typedef CompressedArray TFilter;

    // Empty default constructor
    KmerFilter<Dna, TSpec, TFilter> ctor_empty;
    // Default constructor

    KmerFilter<Dna, TSpec, TFilter> ctor_default (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec, TFilter> ctor_default_helper1 (noBins, hashFunc, kmerSize, bits);
    KmerFilter<Dna, TSpec, TFilter> ctor_default_helper2 (noBins, hashFunc, kmerSize, bits);

    // KmerFilter<Dna, TSpec, TFilter> ctor_default (noBins, kmerSize);
    // KmerFilter<Dna, TSpec, TFilter> ctor_default_helper1 (noBins, kmerSize);
    // KmerFilter<Dna, TSpec, TFilter> ctor_default_helper2 (noBins, kmerSize);

    // Copy constructor
    KmerFilter<Dna, TSpec, TFilter> ctor_copy (ctor_default);
    // Copy assignment
    KmerFilter<Dna, TSpec, TFilter> assignment_copy;
    assignment_copy = ctor_default;
    // Move constructor
    KmerFilter<Dna, TSpec, TFilter> ctor_move(std::move(ctor_default_helper1));
    // Move assignment
    KmerFilter<Dna, TSpec, TFilter> assignment_move;
    assignment_move = std::move(ctor_default_helper2);


    // ==========================================================================
    // Test insertKmer()
    // ==========================================================================
    std::cout << "Testing insertKmer" << '\n';

        CharString fasta("../../unit/test.fasta");
        insertKmer(ctor_default, toCString(fasta), 1);
        insertKmer(ctor_default, toCString(fasta), 5);
        insertKmer(ctor_default, toCString(fasta), 8);
        insertKmer(ctor_copy, toCString(fasta), 2);
        insertKmer(assignment_copy, toCString(fasta), 3);
        insertKmer(ctor_move, toCString(fasta), 0);
        insertKmer(assignment_move, toCString(fasta), 9);


    // ==========================================================================
    // Test storing and retrieving
    // ==========================================================================
    std::cout << "Testing store/retrieve" << '\n';

    CharString store1("file");
    store(ctor_default, store1);
    retrieve(ctor_empty, store1);



    // ==========================================================================
    // Test select()
    // ==========================================================================
    std::cout << "Testing select" << '\n';

    std::vector<uint16_t> ctor_default_set;

    std::vector<bool> which = select(ctor_default, DnaString(std::string(kmerSize, 'A')), 1);
    (void) select(ctor_default, DnaString(std::string(kmerSize, 'T')));
    for (uint16_t i = 0; i < which.size(); ++i)
    {
        if (i == 1 || i == 5 || i == 8)
        {
            assert(which[i]);
            ctor_default_set.push_back(i);
        }
        else
            assert(!which[i]);
    }


    // ==========================================================================
    // Test clear()
    // ==========================================================================
    std::cout << "Testing clear" << '\n';

    // Check if any elements are set in the filters.
    bool ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector.get_pos(i))
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == true);

    // Reset the filter vectors.
    clear(ctor_default, ctor_default_set, threads);

    // Check if filter Vectors are empty.
    ctor_default_any = false;
    for (uint64_t i = 0; i < ctor_default.noOfBlocks * ctor_default.noOfBins; ++i)
    {
        if (ctor_default.filterVector.get_pos(i))
        {
            ctor_default_any = true;
            break;
        }
    }
    assert(ctor_default_any == false);

    return 0;
}
#define NDEBUG
