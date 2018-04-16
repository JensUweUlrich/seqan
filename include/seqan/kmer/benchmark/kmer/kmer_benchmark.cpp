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

template <typename TAlphabet>
static void addKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter> ibf (bins, hash, k, (1ULL<<bits)+256);
    std::mt19937 RandomNumber;
    StringSet<String<TAlphabet> > input;
    reserve(input, 1000000);

    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        String<TAlphabet> tmp;
        for (int32_t i = 0; i < k; ++i)
        {
            appendValue(tmp, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
        }
        appendValue(input, tmp);
    }

    uint64_t i{0};
    for (uint64_t chunk = 0; chunk < ibf.filterVector.noOfChunks; ++chunk)
    {
        auto start = std::chrono::high_resolution_clock::now();
        ibf.filterVector.decompress(chunk);
        auto end   = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.counters["decompress"] = elapsed_seconds.count();

        for (auto _ : state)
        {
            addKmer(ibf, input[i % 1000000], i % bins, chunk);
            ++i;
        }

        start = std::chrono::high_resolution_clock::now();
        ibf.filterVector.compress(chunk);
        end   = std::chrono::high_resolution_clock::now();
        elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        state.counters["compress"] = elapsed_seconds.count();
    }

    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet>
static void whichBins_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    auto occ  = state.range(4);
    KmerFilter<TAlphabet, InterleavedBloomFilter> ibf(bins, hash, k, (1ULL<<bits)+256);
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");

    auto vecPos = (1ULL<<bits) - occ;
    while (vecPos > 0)
    {
        uint64_t chunk = vecPos / ibf.filterVector.chunkSize;
        uint64_t chunkPos = vecPos - chunk * ibf.filterVector.chunkSize;
        ibf.filterVector.set_pos(chunk, chunkPos);
        vecPos -= occ;
    }
    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();

    StringSet<String<TAlphabet> > input;
    reserve(input, 1000000);

    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        String<TAlphabet> tmp;
        for (int32_t i = 0; i < k; ++i)
        {
            appendValue(tmp, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
        }
        appendValue(input, tmp);
    }

    uint64_t i{0};

    for (auto _ : state)
    {
        whichBins(ibf, input[i % 1000000], 0);
        ++i;
    }
}
/*
template <typename TAlphabet>
static void addKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing> da (bins, k);
    std::mt19937 RandomNumber;
    StringSet<String<TAlphabet> > input;
    reserve(input, 1000000);

    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        String<TAlphabet> tmp;
        for (int32_t i = 0; i < k; ++i)
        {
            appendValue(tmp, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
        }
        appendValue(input, tmp);
    }

    uint64_t i{0};

    for (auto _ : state)
    {
        addKmer(da, input[i % 1000000], i % bins);
        ++i;
    }

    auto start = std::chrono::high_resolution_clock::now();
    da.filterVector.unload();
    auto end   = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

    state.counters["Size"] = da.filterVector.size_in_mega_bytes();
    state.counters["unload"] = elapsed_seconds.count();
}

template <typename TAlphabet>
static void whichBins_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto occ = state.range(2);
    KmerFilter<TAlphabet, DirectAddressing> da (bins, k);
    std::mt19937 RandomNumber;
    String<TAlphabet> kmer("");

    auto vecPos = da.noOfBits - da.filterMetadataSize - occ;
    while (vecPos > 0)
    {
        da.filterVector.set_pos(vecPos);
        // ibf.filterVector.set_pos(vecPos);
        vecPos -= occ;
    }

    da.filterVector.unload();
    state.counters["Size"] = da.filterVector.size_in_mega_bytes();

    StringSet<String<TAlphabet> > input;
    reserve(input, 1000000);

    for (uint64_t seqNo = 0; seqNo < 1000000; ++seqNo)
    {
        String<TAlphabet> tmp;
        for (int32_t i = 0; i < k; ++i)
        {
            appendValue(tmp, TAlphabet(RandomNumber() % ValueSize<TAlphabet>::VALUE));
        }
        appendValue(input, tmp);
    }

    uint64_t i{0};

    for (auto _ : state)
    {
        whichBins(da, input[i % 1000000], 0);
        ++i;
    }
}
*/
static void IBFAddArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 20; k <= 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 32; bits < 37; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {
                    b->Args({binNo, k, bits, hashNo});
                }
            }
        }
    }
}

static void IBFWhichArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 20; k <= 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 32; bits < 37; ++bits )
            {
                for (int32_t hashNo = 3; hashNo < 4; ++hashNo)
                {
                    for (int32_t occ = 1; occ <= 8192; occ *= 2)
                    {
                        if (occ==4 || (occ >= 16 && occ < 64) || occ==256 || occ==512 || occ==2048 || occ==4096)
                            continue;
                        b->Args({binNo, k, bits, hashNo, occ});
                    }
                }
            }
        }
    }
}
/*
static void DAAddArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 1; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 10; k <= 13; ++k)
        {
            b->Args({binNo, k});
        }
    }
}

static void DAWhichArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 10; k <= 13; ++k)
        {
            for (int32_t occ = 1; occ <= 8192; occ *= 2)
            {
                if (occ==4 || (occ >= 16 && occ < 64) || occ==256 || occ==512 || occ==2048 || occ==4096)
                    continue;
                b->Args({binNo, k, occ});
            }
        }
    }
}
*/
BENCHMARK_TEMPLATE(addKmer_IBF, Dna)->Apply(IBFAddArguments);
// BENCHMARK_TEMPLATE(addKmer_DA, Dna)->Apply(DAAddArguments);
BENCHMARK_TEMPLATE(whichBins_IBF, Dna)->Apply(IBFWhichArguments);
// BENCHMARK_TEMPLATE(whichBins_DA, Dna)->Apply(DAWhichArguments);

BENCHMARK_MAIN();
