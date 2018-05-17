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

// CharString baseDir{"/srv/public/enricoseiler/benchmark/"}; // lncrna
CharString baseDir{"/group/ag_abi/seiler/benchmark/"}; // redwood
uint64_t e{2};
uint64_t readNo{1000000};

template <typename TAlphabet, typename TFilter>
static void insertKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf (bins, hash, k, (1ULL<<bits));

    for (auto _ : state)
    {
        insertKmerDir(ibf, toCString(baseDir), 8);
        // double elapsed_seconds{0.0};
        // for(int32_t i = 0; i < bins; ++i)
        // {
        //     CharString file(baseDir);
        //     append(file, CharString(std::to_string(bins)));
        //     append(file, CharString{"/bins/bin_"});
        //     append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
        //     append(file, CharString(".fasta"));
        //
        //     auto start = std::chrono::high_resolution_clock::now();
        //     insertKmer(ibf, toCString(file), i);
        //     auto end   = std::chrono::high_resolution_clock::now();
        //     elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();
        //     std::cerr << "IBF Bin " << i << " done." << '\n';
        // }
        // state.SetIterationTime(elapsed_seconds);
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));
    if constexpr (std::is_same_v<TFilter, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_ibf.filter"));
    store(ibf, storage);

    state.counters["Size"] = ibf.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter>
static void select_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    KmerFilter<TAlphabet, InterleavedBloomFilter, TFilter> ibf(bins, hash, k, (1ULL<<bits));

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));
    if constexpr (std::is_same_v<TFilter, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_ibf.filter"));
    retrieve(ibf, storage);
    ibf.filterVector.compress(0);

    uint64_t verifications{0};
    uint64_t tp{0};

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
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
                auto res = select(ibf, seq, 100-k+1 - k*e);
                auto end   = std::chrono::high_resolution_clock::now();
                elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();
                if (res[i])
                {
                    ++tp;
                }
                verifications += count(res.begin(), res.end(), true);
            }
        }
        state.SetIterationTime(elapsed_seconds);
        state.counters["Verifications"] = static_cast<double>(verifications)/readNo;
        state.counters["Sensitivity"] = static_cast<double>(tp)/verifications;
    }
}

template <typename TAlphabet, typename TFilter>
static void insertKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter> da (bins, k);

    for (auto _ : state)
    {
        insertKmerDir(da, toCString(baseDir), 8);
        // double elapsed_seconds{0.0};
        // for(int32_t i = 0; i < bins; ++i)
        // {
        //     CharString file(baseDir);
        //     append(file, CharString(std::to_string(bins)));
        //     append(file, CharString{"/bins/bin_"});
        //     append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
        //     append(file, CharString(".fasta"));
        //
        //     auto start = std::chrono::high_resolution_clock::now();
        //     insertKmer(da, toCString(file), i);
        //     auto end   = std::chrono::high_resolution_clock::now();
        //     elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();
        //     std::cerr << "DA Bin " << i << " done." << '\n';
        // }
        // state.SetIterationTime(elapsed_seconds);
    }

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    if constexpr (std::is_same_v<TFilter, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_da.filter"));
    store(da, storage);

    state.counters["Size"] = da.filterVector.size_in_mega_bytes();
}

template <typename TAlphabet, typename TFilter>
static void select_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    KmerFilter<TAlphabet, DirectAddressing, TFilter> da (bins, k);

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    if constexpr (std::is_same_v<TFilter, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_da.filter"));
    retrieve(da, storage);
    da.filterVector.compress(0);

    uint64_t verifications{0};
    uint64_t tp{0};

    for (auto _ : state)
    {
        double elapsed_seconds{0.0};
        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
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
                auto res = select(da, seq, 100-k+1 - k*e);
                auto end   = std::chrono::high_resolution_clock::now();
                elapsed_seconds += (std::chrono::duration_cast<std::chrono::duration<double> >(end - start)).count();

                if (res[i])
                {
                    ++tp;
                }
                verifications += count(res.begin(), res.end(), true);
            }
        }
        state.SetIterationTime(elapsed_seconds);
        state.counters["Verifications"] = verifications/readNo;
        state.counters["Sensitivity"] = tp/verifications;
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

[[maybe_unused]]
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

// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedArray)->Apply(IBFAddArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedArray)->Apply(DAAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedArray)->Apply(DAWhichArguments)->UseManualTime();

BENCHMARK_MAIN();
