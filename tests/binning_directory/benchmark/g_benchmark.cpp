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
#include <seqan/binning_directory.h>
#include <atomic>

using namespace seqan;

// CharString baseDir{"/srv/public/enricoseiler/benchmark/"}; // lncrna
CharString baseDir{"/srv/hdd/seiler/IBF/"}; // redwood
uint64_t e{3};

template <typename TValue, typename THash, typename TBitvector>
static void insertKmer_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    auto hash = state.range(3);
    BinningDirectory<InterleavedBloomFilter, BDConfig<TValue, THash, TBitvector> > ibf (bins, hash, k, (1ULL<<bits));

    for (auto _ : state)
    {
        insertKmerDir(ibf, toCString(baseDir), 32);
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
    if constexpr (std::is_same_v<TBitvector, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_ibf.filter"));
    store(ibf, storage);

    state.counters["Size"] = size(ibf);
}

template <typename TValue, typename THash, typename TBitvector>
static void select_IBF(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    auto bits = state.range(2);
    // auto hash = state.range(3);
    // BinningDirectory<InterleavedBloomFilter, BDConfig<TValue, THash, TBitvector> > ibf (bins, hash, k, (1ULL<<bits));

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(bits)));
    if constexpr (std::is_same_v<TBitvector, Uncompressed>)
    {
        append(storage, CharString("_Uncompressed"));
    }
    else if constexpr (std::is_same_v<TBitvector, Compressed>)
    {
        append(storage, CharString("_Compressed"));
    }
    else
    {
        throw "Invalid template";
    }
    append(storage, CharString("_ibf.filter"));

    auto fullTime = std::chrono::high_resolution_clock::now();

    double loadingTime{0.0};
    auto start = std::chrono::high_resolution_clock::now();
    BinningDirectory<InterleavedBloomFilter, BDConfig<TValue, THash, TBitvector> > ibf (storage);
    // retrieve(ibf, storage);
    auto end   = std::chrono::high_resolution_clock::now();
    loadingTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();


    double compressionTime{0.0};
    start = std::chrono::high_resolution_clock::now();
    ibf.bitvector.compress(0); // Loading automatically compresses
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
        Semaphore thread_limiter(32);
        std::mutex mtx;
        std::mutex mtx2;
        std::vector<std::future<void>> tasks;

        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file{"/dev/shm/IBF/"};
            // CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            tasks.emplace_back(
                std::async(std::launch::async, [&, file, i] {
                    Critical_section _(thread_limiter);
                    CharString id;
                    String<TValue> seq;
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
                        auto res = select(ibf, seq, e);
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
        state.counters["TP"] = tp.load();
        state.counters["FN"] = fn.load();
        state.counters["FP"] = fp.load();
        state.counters["P"] = p.load();
        state.counters["readNo"] = readNo.load();
        state.counters["verifications"] = verifications.load();
        state.counters["Verifications"] = static_cast<double>(verifications.load())/readNo.load();
        state.counters["Sensitivity"] = static_cast<double>(tp.load())/readNo.load();
        state.counters["Precision"] = static_cast<double>(tp.load())/p.load();
        state.counters["FNR"] = static_cast<double>(fn.load())/readNo.load();
        state.counters["FDR"] = static_cast<double>(fp.load())/p.load();
        state.counters["compressionTime"] = compressionTime;
        state.counters["loadingTime"] = loadingTime;
        state.counters["ioTime"] = ioTime;
        state.counters["selectTime"] = selectTime;
        state.counters["vectorSize"] = size(ibf);
        state.counters["fullTime"] = std::chrono::duration_cast<std::chrono::duration<double> >(fullTime2 - fullTime).count();
    }
}

template <typename TValue, typename THash, typename TBitvector>
static void insertKmer_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    BinningDirectory<DirectAddressing, BDConfig<TValue, THash, TBitvector> > da (bins, k);

    for (auto _ : state)
    {
        insertKmerDir(da, toCString(baseDir), 32);
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
    if constexpr (std::is_same_v<TBitvector, Uncompressed>) {
        append(storage, CharString("_Uncompressed"));
    }
    else
    {
        append(storage, CharString("_Compressed"));
    }
    append(storage, CharString("_da.filter"));
    store(da, storage);

    state.counters["Size"] = size(da);
}

template <typename TValue, typename THash, typename TBitvector>
static void select_DA(benchmark::State& state)
{
    auto bins = state.range(0);
    auto k = state.range(1);
    // BinningDirectory<DirectAddressing, BDConfig<TValue, THash, TBitvector> > da (bins, k);

    CharString storage("");
    append(storage, CharString(std::to_string(bins)));
    append(storage, CharString("_"));
    append(storage, CharString(std::to_string(k)));
    if constexpr (std::is_same_v<TBitvector, Uncompressed>)
    {
        append(storage, CharString("_Uncompressed"));
    }
    else if constexpr (std::is_same_v<TBitvector, Compressed>)
    {
        append(storage, CharString("_Compressed"));
    }
    else
    {
        throw "Invalid template";
    }
    append(storage, CharString("_da.filter"));

    auto fullTime = std::chrono::high_resolution_clock::now();

    double loadingTime{0.0};
    auto start = std::chrono::high_resolution_clock::now();
    // retrieve(da, storage);
    BinningDirectory<DirectAddressing, BDConfig<TValue, THash, TBitvector> > da (storage);
    auto end   = std::chrono::high_resolution_clock::now();
    loadingTime = std::chrono::duration_cast<std::chrono::duration<double> >(end - start).count();


    double compressionTime{0.0};
    start = std::chrono::high_resolution_clock::now();
    da.bitvector.compress(0); // Loading automatically compresses
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
        Semaphore thread_limiter(32);
        std::mutex mtx;
        std::mutex mtx2;
        std::vector<std::future<void>> tasks;

        for(int32_t i = 0; i < bins; ++i)
        {
            CharString file{"/dev/shm/IBF/"};
            // CharString file(baseDir);
            append(file, CharString(std::to_string(bins)));
            append(file, CharString{"/reads/bin_"});
            append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
            append(file, CharString(".fastq"));

            tasks.emplace_back(
                std::async(std::launch::async, [&, file, i] {
                    Critical_section _(thread_limiter);
                    CharString id;
                    String<TValue> seq;
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
                        auto res = select(da, seq, e);
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
        state.counters["TP"] = tp.load();
        state.counters["FN"] = fn.load();
        state.counters["FP"] = fp.load();
        state.counters["P"] = p.load();
        state.counters["readNo"] = readNo.load();
        state.counters["verifications"] = verifications.load();
        state.counters["Verifications"] = static_cast<double>(verifications.load())/readNo.load();
        state.counters["Sensitivity"] = static_cast<double>(tp.load())/readNo.load();
        state.counters["Precision"] = static_cast<double>(tp.load())/p.load();
        state.counters["FNR"] = static_cast<double>(fn.load())/readNo.load();
        state.counters["FDR"] = static_cast<double>(fp.load())/p.load();
        state.counters["compressionTime"] = compressionTime;
        state.counters["loadingTime"] = loadingTime;
        state.counters["ioTime"] = ioTime;
        state.counters["selectTime"] = selectTime;
        state.counters["vectorSize"] = size(da);
        state.counters["fullTime"] = std::chrono::duration_cast<std::chrono::duration<double> >(fullTime2 - fullTime).count();
    }
}

[[maybe_unused]]
static void IBFArguments(benchmark::internal::Benchmark* b)
{
    for (int32_t binNo = 64; binNo <= 8192; binNo *= 2)
    // for (int32_t binNo = 256; binNo <= 256; binNo *= 2)
    {
        if ((binNo > 1 && binNo < 64) || binNo==128 || binNo==512 || binNo==2048 || binNo==4096)
            continue;
        for (int32_t k = 19; k < 20; ++k)
        {
            // 35 = 4GiB, 36 = 8GiB, 37 = 16GiB
            for (int32_t bits = 33; bits <= 37; ++bits )
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
        for (int32_t k = 10; k <= 15; ++k)
        {
            if (binNo == 8192 && k >= 14)
                continue;
            b->Args({binNo, k});
        }
    }
}
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Normal<5>, Uncompressed)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, Normal<5>, Compressed)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(insertKmer_IBF, Dna, CompressedArray)->Apply(IBFAddArguments)->UseManualTime();
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Normal<5>, Uncompressed)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, Normal<5>, Compressed)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(insertKmer_DA, Dna, CompressedArray)->Apply(DAAddArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_IBF, Dna, Normal<5>, Uncompressed)->Apply(IBFArguments);
BENCHMARK_TEMPLATE(select_IBF, Dna, Normal<5>, Compressed)->Apply(IBFArguments);
// BENCHMARK_TEMPLATE(select_IBF, Dna, CompressedArray)->Apply(IBFWhichArguments)->UseManualTime();
BENCHMARK_TEMPLATE(select_DA, Dna, Normal<5>, Uncompressed)->Apply(DAArguments);
BENCHMARK_TEMPLATE(select_DA, Dna, Normal<5>, Compressed)->Apply(DAArguments);
// BENCHMARK_TEMPLATE(select_DA, Dna, CompressedArray)->Apply(DAWhichArguments));

BENCHMARK_MAIN();
