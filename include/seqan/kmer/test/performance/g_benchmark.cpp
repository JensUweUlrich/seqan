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

CharString baseDir{"/Users/enricoseiler/dev/eval/64/"};
uint64_t e{2};
uint64_t readNo{1000000};

template <typename TAlphabet, typename TFilter>
static void addKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf (bins, hash, k, (1ULL<<bits)+256);

    for (auto _ : state)
    {
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString{"bins/"});
            if (i < 10)
            {
                append(file, CharString("bin_0"));
            }
            else
            {
                append(file, CharString("bin_"));
            }
            append(file, CharString(std::to_string(i)));
            append(file, CharString(".fasta"));

            auto start = std::chrono::high_resolution_clock::now();
            addFastaFile(ibf, toCString(file), i);
            auto end   = std::chrono::high_resolution_clock::now();
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
            state.SetIterationTime(elapsed_seconds.count());
            std::cerr << "IBF Bin " << i << " done." << '\n';
        }
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));
    append(storage, CharString("_ibf.filter"));
    store(ibf, storage);

    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter>
static void whichBins_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf(bins, hash, k, (1ULL<<bits)+256);

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));
    append(storage, CharString("_ibf.filter"));
    retrieve(ibf, storage);

    uint64_t verifications{0};

    for (auto _ : state)
    {
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString{"reads/"});

            if (i < 10)
            {
                append(file, CharString("bin_0"));
            }
            else
            {
                append(file, CharString("bin_"));
            }
            append(file, CharString(std::to_string(i)));
            append(file, CharString(".fastq"));

            CharString id;
            String<TAlphabet> seq;
            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(file)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(file));
                throw toCString(msg);
            }
            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                auto start = std::chrono::high_resolution_clock::now();
                auto res = whichBins(ibf, seq, 100-k+1 - k*e);
                auto end   = std::chrono::high_resolution_clock::now();
                auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
                state.SetIterationTime(elapsed_seconds.count());

                verifications += count(res.begin(), res.end(), true);
            }
            ++i;
        }
        state.counters["Verifications"] = verifications/readNo;
    }
}

template <typename TAlphabet, typename TFilter>
static void addKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing> da (bins, k);

    for (auto _ : state)
    {
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString{"bins/"});
            if (i < 10)
            {
                append(file, CharString("bin_0"));
            }
            else
            {
                append(file, CharString("bin_"));
            }
            append(file, CharString(std::to_string(i)));
            append(file, CharString(".fasta"));

            auto start = std::chrono::high_resolution_clock::now();
            addFastaFile(da, toCString(file), i);
            auto end   = std::chrono::high_resolution_clock::now();
            auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
            state.SetIterationTime(elapsed_seconds.count());
            // std::cerr << "DA Iteration " << r << " Bin " << i << " done." << '\n';
        }
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_da.filter"));
    store(da, storage);

    state.counters["Size"] = da.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter>
static void whichBins_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing> da (bins, k);

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_da.filter"));
    retrieve(da, storage);

    uint64_t verifications{0};

    for (auto _ : state)
    {
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString{"reads/"});

            if (i < 10)
            {
                append(file, CharString("bin_0"));
            }
            else
            {
                append(file, CharString("bin_"));
            }
            append(file, CharString(std::to_string(i)));
            append(file, CharString(".fastq"));

            CharString id;
            String<TAlphabet> seq;
            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(file)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(file));
                throw toCString(msg);
            }
            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                auto start = std::chrono::high_resolution_clock::now();
                auto res = whichBins(da, seq, 100-k+1 - k*e);
                auto end   = std::chrono::high_resolution_clock::now();
                auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double> >(end - start);
                state.SetIterationTime(elapsed_seconds.count());

                verifications += count(res.begin(), res.end(), true);
            }
            ++i;
        }
        state.counters["Verifications"] = verifications/readNo;
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

BENCHMARK_TEMPLATE(addKmer_IBF, Dna, Uncompressed)->Apply(IBFArguments)->UseManualTime();
BENCHMARK_TEMPLATE(addKmer_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(addKmer_IBF, Dna, CompressedArray)->Apply(IBFAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(addKmer_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
BENCHMARK_TEMPLATE(addKmer_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(addKmer_DA, Dna, CompressedArray)->Apply(DAAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(whichBins_IBF, Dna, Uncompressed)->Apply(IBFArguments)->UseManualTime();
BENCHMARK_TEMPLATE(whichBins_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(whichBins_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
BENCHMARK_TEMPLATE(whichBins_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
BENCHMARK_TEMPLATE(whichBins_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(whichBins_DA, Dna, CompressedArray)->Apply(DAWhichArguments)->UseManualTime();

BENCHMARK_MAIN();
