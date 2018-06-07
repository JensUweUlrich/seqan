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
#include <atomic>

using namespace seqan;

// CharString baseDir{"/srv/public/enricoseiler/benchmark/"}; // lncrna
CharString baseDir{"/group/ag_abi/seiler/small_benchmark/"}; // redwood
uint64_t e{2};

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
    if constexpr (std::is_same_v<TFilter, Uncompressed> || std::is_same_v<TFilter, CompressedSimple>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        throw "Invalid template";
    }
    append(storage, CharString("_ibf.filter"));

    auto fullTime = std::chrono::high_resolution_clock::now();

    double loadingTime{0.0};
    auto start = std::chrono::high_resolution_clock::now();
    retrieve(ibf, storage);
    auto end   = std::chrono::high_resolution_clock::now();
    loadingTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();


    double compressionTime{0.0};
    start = std::chrono::high_resolution_clock::now();
    ibf.filterVector.compress(0); // Loading automatically compresses
    end   = std::chrono::high_resolution_clock::now();
    compressionTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();

    std::atomic<uint64_t> verifications{0};
    std::atomic<uint64_t> tp{0};
    std::atomic<uint64_t> p{0};
    std::atomic<uint64_t> fp{0};
    std::atomic<uint64_t> fn{0};
    std::atomic<uint64_t> readNo{0};

    for (auto _ : state)
    {
        double selectTime{0.0};
        double ioTime{0.0};
        Semaphore thread_limiter(8);
        std::mutex mtx;
        std::mutex mtx2;
        std::vector<std::future<void>> tasks;

        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            tasks.emplace_back(
                std::async(std::launch::async, [&, file, i] {
                    Critical_section _(thread_limiter);
                    CharString id;
                    String<TAlphabet> seq;
                    SeqFileIn seqFileIn;
                    uint64_t c{0};
                    if (!open(seqFileIn, toCString(file)))
                    {
                        CharString msg = "Unable to open contigs file: ";
                        append(msg, CharString(file));
                        throw toCString(msg);
                    }
                    while(!atEnd(seqFileIn))
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        readRecord(id, seq, seqFileIn);
                        auto end   = std::chrono::high_resolution_clock::now();
                        mtx2.lock();
                        ioTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx2.unlock();

                        start = std::chrono::high_resolution_clock::now();
                        auto res = select(ibf, seq, 100-k+1 - k*e);
                        end   = std::chrono::high_resolution_clock::now();
                        ++readNo;
                        mtx.lock();
                        selectTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx.unlock();

                        if (res[i])
                            ++tp;
                        else
                            ++fn;
                        c = count(res.begin(), res.end(), true);
                        verifications += c;
                        if (c > 1)
                        {
                            if (res[i])
                                fp += c - 1;
                            else
                                fp += c;
                        }
                        p += c;
                    }
                })
            );
        }
        for (auto &&task : tasks){
            task.get();
        }

        auto fullTime2   = std::chrono::high_resolution_clock::now();
        state.counters["5_TP"] = tp.load();
        state.counters["6_FN"] = fn.load();
        state.counters["7_FP"] = fp.load();
        state.counters["8_P"] = p.load();
        state.counters["99_readNo"] = readNo.load();
        state.counters["9_verifications"] = verifications.load();
        state.counters["0_Verifications"] = static_cast<double>(verifications.load())/readNo.load();
        state.counters["1_Sensitivity"] = static_cast<double>(tp.load())/readNo.load();
        state.counters["2_Precision"] = static_cast<double>(tp.load())/p.load();
        state.counters["3_FNR"] = static_cast<double>(fn.load())/readNo.load();
        state.counters["4_FDR"] = static_cast<double>(fp.load())/p.load();
        state.counters["compressionTime"] = compressionTime;
        state.counters["loadingTime"] = loadingTime;
        state.counters["ioTime"] = ioTime;
        state.counters["selectTime"] = selectTime;
        state.counters["vectorSize"] = ibf.filterVector.size_in_mega_bytes();
        state.counters["fullTime"] = std::chrono::duration_cast<std::chrono::duration<double> >(fullTime2 - fullTime).count();
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
    if constexpr (std::is_same_v<TFilter, Uncompressed> || std::is_same_v<TFilter, CompressedSimple>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        throw "Invalid template";
    }
    append(storage, CharString("_da.filter"));

    auto fullTime = std::chrono::high_resolution_clock::now();

    double loadingTime{0.0};
    auto start = std::chrono::high_resolution_clock::now();
    retrieve(da, storage);
    auto end   = std::chrono::high_resolution_clock::now();
    loadingTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();


    double compressionTime{0.0};
    start = std::chrono::high_resolution_clock::now();
    da.filterVector.compress(0); // Loading automatically compresses
    end   = std::chrono::high_resolution_clock::now();
    compressionTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();

    std::atomic<uint64_t> verifications{0};
    std::atomic<uint64_t> tp{0};
    std::atomic<uint64_t> p{0};
    std::atomic<uint64_t> fp{0};
    std::atomic<uint64_t> fn{0};
    std::atomic<uint64_t> readNo{0};

    for (auto _ : state)
    {
        double selectTime{0.0};
        double ioTime{0.0};
        Semaphore thread_limiter(8);
        std::mutex mtx;
        std::mutex mtx2;
        std::vector<std::future<void>> tasks;

        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            tasks.emplace_back(
                std::async(std::launch::async, [&, file, i] {
                    Critical_section _(thread_limiter);
                    CharString id;
                    String<TAlphabet> seq;
                    SeqFileIn seqFileIn;
                    uint64_t c{0};
                    if (!open(seqFileIn, toCString(file)))
                    {
                        CharString msg = "Unable to open contigs file: ";
                        append(msg, CharString(file));
                        throw toCString(msg);
                    }
                    while(!atEnd(seqFileIn))
                    {
                        auto start = std::chrono::high_resolution_clock::now();
                        readRecord(id, seq, seqFileIn);
                        auto end   = std::chrono::high_resolution_clock::now();
                        mtx2.lock();
                        ioTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx2.unlock();

                        start = std::chrono::high_resolution_clock::now();
                        auto res = select(da, seq, 100-k+1 - k*e);
                        end   = std::chrono::high_resolution_clock::now();
                        ++readNo;
                        mtx.lock();
                        selectTime += std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();
                        mtx.unlock();

                        if (res[i])
                            ++tp;
                        else
                            ++fn;
                        c = count(res.begin(), res.end(), true);
                        verifications += c;
                        if (c > 1)
                        {
                            if (res[i])
                                fp += c - 1;
                            else
                                fp += c;
                        }
                        p += c;
                    }
                })
            );
        }
        for (auto &&task : tasks){
            task.get();
        }

        auto fullTime2   = std::chrono::high_resolution_clock::now();
        state.counters["5_TP"] = tp.load();
        state.counters["6_FN"] = fn.load();
        state.counters["7_FP"] = fp.load();
        state.counters["8_P"] = p.load();
        state.counters["99_readNo"] = readNo.load();
        state.counters["9_verifications"] = verifications.load();
        state.counters["0_Verifications"] = static_cast<double>(verifications.load())/readNo.load();
        state.counters["1_Sensitivity"] = static_cast<double>(tp.load())/readNo.load();
        state.counters["2_Precision"] = static_cast<double>(tp.load())/p.load();
        state.counters["3_FNR"] = static_cast<double>(fn.load())/readNo.load();
        state.counters["4_FDR"] = static_cast<double>(fp.load())/p.load();
        state.counters["compressionTime"] = compressionTime;
        state.counters["loadingTime"] = loadingTime;
        state.counters["ioTime"] = ioTime;
        state.counters["selectTime"] = selectTime;
        state.counters["vectorSize"] = da.filterVector.size_in_mega_bytes();
        state.counters["fullTime"] = std::chrono::duration_cast<std::chrono::duration<double> >(fullTime2 - fullTime).count();
        break;
    }
}

[[maybe_unused]]
static void IBFArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 14; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 20; bits <= 24; ++bits )
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
        for (int32_t k = 3; k <= 4; ++k)
        {
            if (binNo == 8192 && k == 14)
                continue;
            b->Args({binNo, k});
        }
    }
}

BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Uncompressed)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedSimple)->Apply(IBFArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedArray)->Apply(IBFAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Uncompressed)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedSimple)->Apply(DAArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedArray)->Apply(DAAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Uncompressed)->Apply(IBFArguments);
BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedSimple)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Uncompressed)->Apply(DAArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, CompressedSimple)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedArray)->Apply(DAWhichArguments));

BENCHMARK_MAIN();
