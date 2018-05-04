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
#include <random>
#include <benchmark/benchmark.h>
#include <seqan/kmer.h>

using namespace seqan;

struct kmers{
    StringSet<String<Dna> > input;

    kmers()
    {
        std::mt19937 RandomNumber;

        reserve(input, 1000000);

        for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
        {
            String<Dna> tmp;
            for (int32_t i = 0; i < 19; ++i)
            {
                appendValue(tmp, Dna(RandomNumber() % ValueSize<Dna>::VALUE));
            }
            appendValue(input, tmp);
        }
    }

    String<Dna> operator[] (const int index)
    {
        return input[index];
    }
};

kmers input;

template <typename TAlphabet, typename TFilter>
static void insertKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf (bins, hash, k, (1ULL<<bits));

    uint64_t i{0};

    for (auto _ : state)
    {
        auto in = input[i % 1000000];
        uint64_t b = i % bins;
        auto start = std::chrono::high_resolution_clock::now();
        insertKmer(ibf, in, b,0);
        auto end   = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        state.SetIterationTime(elapsed_seconds.count());

        ++i;
    }
}

template <typename TAlphabet, typename TFilter>
static void select_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf (bins, hash, k, (1ULL<<bits));

    uint64_t i{0};
    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        insertKmer(ibf, input[seqNo], i%bins,0);
        ++i;
    }

    for (auto _ : state)
    {
        auto in = input[i % 1000000];
        auto start = std::chrono::high_resolution_clock::now();
        select(ibf, in, 1);
        auto end   = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        state.SetIterationTime(elapsed_seconds.count());

        ++i;
    }
}

template <typename TAlphabet, typename TFilter>
static void insertKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter> da (bins, k);

    uint64_t i{0};

    for (auto _ : state)
    {
        auto in = input[i % 1000000];
        uint64_t b = i % bins;
        auto start = std::chrono::high_resolution_clock::now();
        insertKmer(da, in, b, 0);
        auto end   = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        state.SetIterationTime(elapsed_seconds.count());

        ++i;
    }
}

template <typename TAlphabet, typename TFilter>
static void select_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter> da (bins, k);

    uint64_t i{0};
    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        insertKmer(da, input[seqNo], i%bins,0);
        ++i;
    }

    for (auto _ : state)
    {
        auto in = input[i % 1000000];
        auto start = std::chrono::high_resolution_clock::now();
        select(da, in, 1);
        auto end   = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
        state.SetIterationTime(elapsed_seconds.count());

        ++i;
    }
}

static void IBFArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 17; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 35; bits <= 37; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {

                    b->Args({binNo, k, bits, hashNo});
                }
            }
        }
    }
}

static void DAArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 13; k <= 15; ++k)
        {
            b->Args({binNo, k});
        }
    }
}

BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed)->Apply(IBFArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedArray)->Apply(IBFAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedArray)->Apply(DAAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed)->Apply(IBFArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedArray)->Apply(DAWhichArguments)->UseManualTime();

BENCHMARK_MAIN();
