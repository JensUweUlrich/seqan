// ==========================================================================
//                             binning_directory
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
// Author: Your Name <your.email@example.net>
// ==========================================================================

#ifndef TESTS_BINNING_DIRECTORY_TEST_BINNING_DIRECTORY_H_
#define TESTS_BINNING_DIRECTORY_TEST_BINNING_DIRECTORY_H_

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/binning_directory.h>

using namespace seqan;

// // A test for strings.
// SEQAN_DEFINE_TEST(test_binning_directory_strings_example1)
// {
//     using namespace seqan;
//
//     // Define some constant test data for comparison...
//     CharString const STRING1 = "test 1";
//     CharString const STRING2 = "test 2";
//
//     // Append to a string and make equality assertion on the result.
//     CharString myStr = "test ";
//     append(myStr, "1");
//     SEQAN_ASSERT_EQ(STRING1, myStr);
//
//     // Demonstration of other assertions.
//     SEQAN_ASSERT_GT(STRING2, myStr);
//     SEQAN_ASSERT_GEQ(STRING2, myStr);
//     SEQAN_ASSERT_LT(myStr, STRING2);
//     SEQAN_ASSERT_LEQ(STRING2, STRING2);
// }

// A test for strings.
typedef
    TagList<BinningDirectory<Dna,   Shape<Dna, SimpleShape>,    InterleavedBloomFilter,     Uncompressed>,
    TagList<BinningDirectory<Dna,   Shape<Dna, SimpleShape>,    InterleavedBloomFilter,     Compressed>
    > >
    BinningDirectoriesIBF;

typedef
    TagList<BinningDirectory<Dna,   Shape<Dna, SimpleShape>,    DirectAddressing,           Uncompressed>,
    TagList<BinningDirectory<Dna,   Shape<Dna, SimpleShape>,    DirectAddressing,           Compressed>
    > >
    BinningDirectoriesDA;

template <typename TBinning_>
class BinningDirectoryIBFTest : public Test
{
public:
    typedef TBinning_ TBinning;
};

template <typename TBinning_>
class BinningDirectoryDATest : public Test
{
public:
    typedef TBinning_ TBinning;
};

SEQAN_TEST(BinningDirectoryIBFTest, literals)
{
    SEQAN_ASSERT_EQ(0_m,    0ULL);
    SEQAN_ASSERT_EQ(0_g,    0ULL);
    SEQAN_ASSERT_EQ(1_m,    8388608ULL);
    SEQAN_ASSERT_EQ(2_m,    16777216ULL);
    SEQAN_ASSERT_EQ(4_m,    33554432ULL);
    SEQAN_ASSERT_EQ(8_m,    67108864ULL);
    SEQAN_ASSERT_EQ(16_m,   134217728ULL);
    SEQAN_ASSERT_EQ(32_m,   268435456ULL);
    SEQAN_ASSERT_EQ(64_m,   536870912ULL);
    SEQAN_ASSERT_EQ(128_m,  1073741824ULL);
    SEQAN_ASSERT_EQ(256_m,  2147483648ULL);
    SEQAN_ASSERT_EQ(512_m,  4294967296ULL);
    SEQAN_ASSERT_EQ(1_g,    8589934592ULL);
    SEQAN_ASSERT_EQ(2_g,    17179869184ULL);
    SEQAN_ASSERT_EQ(4_g,    34359738368ULL);
    SEQAN_ASSERT_EQ(8_g,    68719476736ULL);
    SEQAN_ASSERT_EQ(16_g,   137438953472ULL);
    SEQAN_ASSERT_EQ(32_g,   274877906944ULL);
    SEQAN_ASSERT_EQ(64_g,   549755813888ULL);
    SEQAN_ASSERT_EQ(128_g,  1099511627776ULL);
    SEQAN_ASSERT_EQ(256_g,  2199023255552ULL);
    SEQAN_ASSERT_EQ(512_g,  4398046511104ULL);
    SEQAN_ASSERT_EQ(1024_g, 8796093022208ULL);
}

SEQAN_TEST(BinningDirectoryDATest, literals)
{
    SEQAN_ASSERT_EQ(0_m,    0ULL);
    SEQAN_ASSERT_EQ(0_g,    0ULL);
    SEQAN_ASSERT_EQ(1_m,    8388608ULL);
    SEQAN_ASSERT_EQ(2_m,    16777216ULL);
    SEQAN_ASSERT_EQ(4_m,    33554432ULL);
    SEQAN_ASSERT_EQ(8_m,    67108864ULL);
    SEQAN_ASSERT_EQ(16_m,   134217728ULL);
    SEQAN_ASSERT_EQ(32_m,   268435456ULL);
    SEQAN_ASSERT_EQ(64_m,   536870912ULL);
    SEQAN_ASSERT_EQ(128_m,  1073741824ULL);
    SEQAN_ASSERT_EQ(256_m,  2147483648ULL);
    SEQAN_ASSERT_EQ(512_m,  4294967296ULL);
    SEQAN_ASSERT_EQ(1_g,    8589934592ULL);
    SEQAN_ASSERT_EQ(2_g,    17179869184ULL);
    SEQAN_ASSERT_EQ(4_g,    34359738368ULL);
    SEQAN_ASSERT_EQ(8_g,    68719476736ULL);
    SEQAN_ASSERT_EQ(16_g,   137438953472ULL);
    SEQAN_ASSERT_EQ(32_g,   274877906944ULL);
    SEQAN_ASSERT_EQ(64_g,   549755813888ULL);
    SEQAN_ASSERT_EQ(128_g,  1099511627776ULL);
    SEQAN_ASSERT_EQ(256_g,  2199023255552ULL);
    SEQAN_ASSERT_EQ(512_g,  4398046511104ULL);
    SEQAN_ASSERT_EQ(1024_g, 8796093022208ULL);
}

SEQAN_TYPED_TEST_CASE(BinningDirectoryIBFTest, BinningDirectoriesIBF);
SEQAN_TYPED_TEST_CASE(BinningDirectoryDATest, BinningDirectoriesDA);

// Empty constructor
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, empty_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd;

    SEQAN_ASSERT_EQ(bd.kmerSize, 0u);
    SEQAN_ASSERT_EQ(bd.noOfBins, 0u);
    SEQAN_ASSERT_EQ(bd.noOfHashFunc, 0u);
    SEQAN_ASSERT_EQ(bd.noOfBits, 0u);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, empty_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd;

    SEQAN_ASSERT_EQ(bd.kmerSize, 0u);
    SEQAN_ASSERT_EQ(bd.noOfBins, 0u);
}

// Default constructor
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, default_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);

    SEQAN_ASSERT_EQ(bd.kmerSize, 12u);
    SEQAN_ASSERT_EQ(bd.noOfBins, 64u);
    SEQAN_ASSERT_EQ(bd.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, default_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);

    SEQAN_ASSERT_EQ(bd.kmerSize, 4u);
    SEQAN_ASSERT_EQ(bd.noOfBins, 64u);
}

// Copy constructor
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, copy_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);
    TBinning bd2(bd);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 12u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
    SEQAN_ASSERT_EQ(bd2.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd2.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, copy_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);
    TBinning bd2(bd);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 4u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
}

// Copy assignment
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, copy_assignment)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);
    TBinning bd2;
    bd2 = bd;

    SEQAN_ASSERT_EQ(bd2.kmerSize, 12u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
    SEQAN_ASSERT_EQ(bd2.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd2.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, copy_assignment)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);
    TBinning bd2;
    bd2 = bd;

    SEQAN_ASSERT_EQ(bd2.kmerSize, 4u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
}

// Move constructor
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, move_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);
    TBinning bd2(std::move(bd));

    SEQAN_ASSERT_EQ(bd2.kmerSize, 12u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
    SEQAN_ASSERT_EQ(bd2.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd2.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, move_constructor)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);
    TBinning bd2(std::move(bd));

    SEQAN_ASSERT_EQ(bd2.kmerSize, 4u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
}

// Move assignment
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, move_assignment)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);
    TBinning bd2;
    bd2 = std::move(bd);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 12u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
    SEQAN_ASSERT_EQ(bd2.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd2.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, move_assignment)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);
    TBinning bd2;
    bd2 = std::move(bd);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 4u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 64u);
}

// insertKmer from text
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, insertKmerText)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 4, 32_m);

    insertKmer(bd, DnaString{""}, 1);
    insertKmer(bd, DnaString{"ACG"}, 63);
    insertKmer(bd, DnaString{"ACGA"}, 0);
    insertKmer(bd, DnaString{"ACGATGCTAGCTAGCTGAC"}, 5);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, insertKmerText)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);

    insertKmer(bd, DnaString{""}, 1);
    insertKmer(bd, DnaString{"A"}, 63);
    insertKmer(bd, DnaString{"ACGA"}, 0);
    insertKmer(bd, DnaString{"ACGATGCTAGCTAGCTGAC"}, 5);
}

// insertKmer from file
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, insertKmerFile)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);

    insertKmer(bd, getAbsolutePath("tests/binning_directory/test.fasta").c_str(), 5);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, insertKmerFile)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);

    insertKmer(bd, getAbsolutePath("tests/binning_directory/test.fasta").c_str(), 5);
}

// count
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, count)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);
    insertKmer(bd, getAbsolutePath("tests/binning_directory/test.fasta").c_str(), 0);

    auto result = count(bd, DnaString{"TAACTTTTTTAT"});

    SEQAN_ASSERT_NEQ(result[0], 0u);

    for (uint16_t i = 1; i < 64; ++i)
    {
        SEQAN_ASSERT_EQ(result[i], 0u);
    }
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, count)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);
    insertKmer(bd, getAbsolutePath("tests/binning_directory/test.fasta").c_str(), 0);

    auto result = count(bd, DnaString{"TAAC"});

    SEQAN_ASSERT_NEQ(result[0], 0u);

    for (uint16_t i = 1; i < 64; ++i)
    {
        SEQAN_ASSERT_EQ(result[i], 0u);
    }
}

// store
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, store)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 3, 12, 32_m);

    store(bd, CharString(SEQAN_TEMP_FILENAME()));
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, store)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(64, 4);

    store(bd, CharString(SEQAN_TEMP_FILENAME()));
}

// retrieve
SEQAN_TYPED_TEST(BinningDirectoryIBFTest, retrieve)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(32, 3, 11, 32_m);
    TBinning bd2(64, 3, 12, 16_m);

    auto tmp = CharString(SEQAN_TEMP_FILENAME());
    store(bd, tmp);
    retrieve(bd2, tmp);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 11u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 32u);
    SEQAN_ASSERT_EQ(bd2.noOfHashFunc, 3u);
    SEQAN_ASSERT_EQ(bd2.noOfBits, 32_m);
}

SEQAN_TYPED_TEST(BinningDirectoryDATest, retrieve)
{
    typedef typename TestFixture::TBinning      TBinning;

    TBinning bd(32, 2);
    TBinning bd2(64, 4);

    auto tmp = CharString(SEQAN_TEMP_FILENAME());
    store(bd, tmp);
    retrieve(bd2, tmp);

    SEQAN_ASSERT_EQ(bd2.kmerSize, 2u);
    SEQAN_ASSERT_EQ(bd2.noOfBins, 32u);
}

#endif  // TESTS_BINNING_DIRECTORY_TEST_BINNING_DIRECTORY_H_
