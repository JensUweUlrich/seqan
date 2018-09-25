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
// Author:  Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
//          Enrico Seiler <enrico.seiler@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_BASE_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_BASE_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/seq_io.h>
#include <valarray>
#include <algorithm>
#include <future>
#include <mutex>
#include <type_traits>

namespace seqan {

// ==========================================================================
// Forwards
// ==========================================================================

class Semaphore
{
    std::mutex m;
    std::condition_variable cv;
    int count;

public:
    Semaphore(int n) : count{n} {}
    void notify()
    {
        std::unique_lock<std::mutex> l(m);
        ++count;
        cv.notify_one();
    }
    void wait()
    {
        std::unique_lock<std::mutex> l(m);
        cv.wait(l, [this]{ return count!=0; });
        --count;
    }
};

class Critical_section
{
    Semaphore &s;
public:
    Critical_section(Semaphore &ss) : s{ss} { s.wait(); }
    ~Critical_section() { s.notify(); }
};

template <typename T>
int numDigits(T number)
{
    int digits = 0;
    if (number <= 0) digits = 1; // remove this line if '-' counts as a digit
    while (number) {
        number /= 10;
        digits++;
    }
    return digits;
}

// ==========================================================================
// Tags, Classes, Enums
// ==========================================================================

// --------------------------------------------------------------------------
// Tag k-mer Filter Tags
// --------------------------------------------------------------------------

//!\brief A tag for the IBF.
struct InterleavedBloomFilter_;
typedef Tag<InterleavedBloomFilter_> InterleavedBloomFilter;

//!\brief A tag for direct addressing.
struct DirectAddressing_;
typedef Tag<DirectAddressing_> DirectAddressing;

//!\brief A tag for the uncompressed Bitvector.
struct Uncompressed_;
typedef Tag<Uncompressed_> Uncompressed;

//!\brief A tag for the compressed Bitvector.
struct Compressed_;
typedef Tag<Compressed_> Compressed;

//!\brief A tag for the compressed array Bitvector.
struct CompressedArray_;
typedef Tag<CompressedArray_> CompressedArray;

struct Normal_;
typedef Tag<Normal_> Normal;

template<typename TSpec>
struct Bitvector;

template<typename TValue, typename TSpec, typename TChunks>
struct BDHash;

template<uint16_t>
struct Offset;

template<typename>
struct is_offset : std::false_type {};

template<uint16_t o>
struct is_offset<Offset<o>> : std::true_type {};

// --------------------------------------------------------------------------
// Class BinningDirectory
// --------------------------------------------------------------------------

template<uint8_t chunks>
struct Chunks
{
    static const uint8_t VALUE = chunks;
};

template<typename TValue_ = Dna, typename THash_ = Normal, typename TBitvector_ = Uncompressed, typename TChunks_ = Chunks<4>>
struct BDConfig
{
    typedef TValue_      TValue;
    typedef THash_       THash;
    typedef TBitvector_  TBitvector;
    typedef TChunks_     TChunks;
};

//!\brief The BinningDirectory class.
template<typename TSpec, typename TConfig>
class BinningDirectory;

// ==========================================================================
// Metafunctions
// ==========================================================================

//!\brief Type definition for variables.
template<typename TSpec, typename TConfig>
struct Value<BinningDirectory<TSpec, TConfig> >
{
    typedef uint32_t noOfBins;
    typedef uint16_t kmerSize;
    typedef uint64_t noOfBits;
    typedef uint64_t noOfBlocks;
    typedef uint32_t binWidth;
    typedef uint32_t blockBitSize;
    typedef uint8_t  intSize;
    typedef uint16_t filterMetadataSize;
    typedef uint8_t  noOfHashFunc;
    typedef uint8_t  shiftValue;
    typedef uint64_t preCalcValues;
    typedef uint64_t seedValue;
    typedef uint8_t  TNoOfChunks;
};

typedef uint32_t TNoOfBins;
typedef uint16_t TKmerSize;
typedef uint64_t TNoOfBits;
typedef uint64_t TNoOfBlocks;
typedef uint32_t TBinWidth;
typedef uint32_t TBlockBitSize;
typedef uint8_t  TIntSize;
typedef uint16_t TFilterMetadataSize;
typedef uint8_t  TNoOfHashFunc;
typedef uint8_t  TShiftValue;
typedef uint64_t TPreCalcValues;
typedef uint64_t TSeedValue;
typedef uint8_t  TNoOfChunks;

// --------------------------------------------------------------------------
// Metafunction MetafunctionName
// --------------------------------------------------------------------------

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function functionName()
// --------------------------------------------------------------------------

/*!
 * \brief Adds k-mers from a text to a bin in a given filter.
 * \param me The BinningDirectory instance.
 * \param text The text from which the k-mers are to be added.
 * \param binNo The bin to add the k-mers to.
 */
template<typename TSpec, typename TConfig>
inline void insertKmer(BinningDirectory<TSpec, TConfig> & me, String<typename TConfig::TValue> const & text, TNoOfBins binNo)
{
    if (!me.chunkMapSet)
        configureChunkMap(me);
    typedef typename TConfig::THash THash;
    if(length(text) >= me.kmerSize)
    {
        if (std::is_same<THash, Normal>::value || is_offset<THash>::value)
        {
            me.template insertKmer<Normal>(text, binNo);
        }
        else
        {
            me.template insertKmer<THash>(text, binNo);
        }
    }
}

/*!
 * \brief Sets the vectors for given bins to 0.
 * \param me The BinningDirectory instance.
 * \param bins A vector containing the bin numbers.
 * \param threads The number of threads to use.
 */
template<typename TSpec, typename TConfig, typename TInt>
inline void clear(BinningDirectory<TSpec, TConfig> &  me, std::vector<TNoOfBins> & bins, TInt&& threads)
{
    // me.clear(bins, static_cast<uint64_t>(threads));
    me.clear(bins, threads);
}

/*!
 * \brief Adds all k-mers from a fasta file to a bin of a given BinningDirectory.
 * \param me The BinningDirectory instance.
 * \param fastaFile The fasta file to process.
 * \param binNo The bin to add the k-mers to.
 * \param batch Parallel batch insertion.
 * \param batchChunkNo Current chunk for CompressedArray in batch mode.
 *
 * If the BinningDirectory's bitvector specialisation is CompressedArray, we need to iterate bitvector.noOfChunks times
 * over the input files. In each iteration the kmers are inserted into one chunk.
 */
template<typename TSpec, typename TConfig>
inline void insertKmer(BinningDirectory<TSpec, TConfig> &  me, const char * fastaFile, TNoOfBins binNo)
{
    if (!me.chunkMapSet)
        configureChunkMap(me);
    typedef typename TConfig::TValue TValue;
    CharString id;
    String<TValue> seq;
    SeqFileIn seqFileIn;
    if (!open(seqFileIn, fastaFile))
    {
        CharString msg = "Unable to open contigs file: ";
        append(msg, CharString(fastaFile));
        std::cerr << msg << std::endl;
        throw toCString(msg);
    }
    while(!atEnd(seqFileIn))
    {
        readRecord(id, seq, seqFileIn);
        if(length(seq) < me.kmerSize)
            continue;
        insertKmer(me, seq, binNo);
    }
    close(seqFileIn);
}

inline uint64_t nextPow4(uint64_t x)
{
    auto clz = sdsl::bits::hi(x);
    if (clz & 1)
        return 1 << (clz + 1);
    else
        return 1 << clz;
}

template<typename TSpec, typename TConfig>
inline void configureChunkMap(BinningDirectory<TSpec, TConfig> & me)
{
    me.chunkMapSet = true;
    typedef typename TConfig::TValue TValue;
    TNoOfChunks chunks = TConfig::TChunks::VALUE;
    TKmerSize kmerSize{me.kmerSize};
    TNoOfChunks effectiveChunks = nextPow4(chunks);
    TNoOfChunks significantBits = sdsl::bits::hi(effectiveChunks);
    TNoOfChunks significantPositions = significantBits >> 1;
    const double alpha = 1.1;
    const double threshold = (1.0 / chunks) * alpha;
    std::vector<uint8_t> chunkMap;
    chunkMap.resize(effectiveChunks, 0);
    std::vector<double> frequencyTable;
    frequencyTable.resize(effectiveChunks, 1.0/effectiveChunks);
    double f{0.0};
    TNoOfChunks chunk{0};
    for (uint8_t index = 0; index < effectiveChunks; ++index)
    {
        double freq = frequencyTable[index];
        if (f + freq >= threshold)
        {
            f = freq;
            if (chunk != chunks -1)
                ++chunk;
            chunkMap[index] = chunk;
        }
        else
        {
            f += freq;
            chunkMap[index] = chunk;
        }
    }
    me.chunkMap = chunkMap;
    me.effectiveChunks = effectiveChunks;
    me.significantBits = significantBits;
    me.significantPositions = significantPositions;
}

template<typename TSpec, typename TConfig>
inline void insertKmer(BinningDirectory<TSpec, TConfig> & me, StringSet<String<typename TConfig::TValue>> & text, std::vector<TNoOfBins> & bins)
{
    // std::cerr << "===========INIT===========\n";
    assert(length(text) == bins.size());
    if (!me.chunkMapSet)
    {
        me.chunkMapSet = true;
        typedef typename TConfig::TValue TValue;
        TNoOfChunks chunks = TConfig::TChunks::VALUE;
        TKmerSize kmerSize{me.kmerSize};
        TNoOfChunks effectiveChunks = nextPow4(chunks);
        TNoOfChunks significantBits = sdsl::bits::hi(effectiveChunks);
        TNoOfChunks significantPositions = significantBits >> 1;
        const double alpha = 1.25;
        const double threshold = (1.0 / chunks) * alpha;

        Shape<TValue, SimpleShape> kmerShape;

        std::vector<uint8_t> chunkMap;

        chunkMap.resize(effectiveChunks, 0);
        std::vector<uint64_t> countTable;
        countTable.resize(effectiveChunks, 0);
        Shape<TValue, SimpleShape> freqShape;
        resize(freqShape, significantPositions);
        uint64_t totalCount{0};

        for (auto & t : text)
        {
            if ( length(t) >= kmerSize)
            {
                // String<TValue> reconstructedText{};
                // std::cerr << "Text = " << t << '\n';
                auto it = begin(t) + (kmerSize - significantPositions);
                if (significantPositions > 1)
                    hashInit(freqShape, it);
                while (it != end(t) - significantPositions + 1)
                {
                    countTable[hashNext(freqShape, it)]++;
                    ++totalCount;
                    // appendValue(reconstructedText, *it);
                    ++it;
                }
                // std::cerr << "RecT = " << reconstructedText << '\n';
                // std::cerr << "new text\n";
            }
        }
        // std::cerr << "counting done\n";
        std::vector<double> frequencyTable;
        frequencyTable.resize(effectiveChunks, 0.0);
        std::transform(countTable.begin(), countTable.end(), frequencyTable.begin(),
                       [&](auto & count) { return (double) count / totalCount;});
        // std::cerr << "transform done\n";
        double f{0.0};
        TNoOfChunks chunk{0};
        for (uint8_t index = 0; index < effectiveChunks; ++index)
        {
            double freq = frequencyTable[index];
            if (f + freq >= threshold)
            {
                f = freq;
                if (chunk != chunks -1)
                    ++chunk;
                chunkMap[index] = chunk;
            }
            else
            {
                f += freq;
                chunkMap[index] = chunk;
            }
        }
        // std::cerr << "==================================\nchunks = " << (int)chunks << "\teffectiveChunks = " << (int)effectiveChunks << "\nalpha = " << alpha << "\t threshold = " << threshold << "\nsigPos = " << (int) significantPositions << "\tkmerSize = " << kmerSize << "\nCount Table\n";
        // for (uint8_t index = 0; index < effectiveChunks; ++index)
        //     std::cerr << (int)index << '\t' << countTable[index] << '\n';
        // std::cerr << "Total Count = " << totalCount << "\nfrequencyTable\n";
        // for (uint8_t index = 0; index < effectiveChunks; ++index)
        //     std::cerr << (int)index << '\t' << frequencyTable[index] << '\n';
        // std::cerr << "Resulting mapping\n";
        // for (uint8_t index = 0; index < effectiveChunks; ++index)
        //     std::cerr << (int)index << '\t' << (int)chunkMap[index] << '\n';

        me.chunkMap = chunkMap;
        me.effectiveChunks = effectiveChunks;
        me.significantBits = significantBits;
        me.significantPositions = significantPositions;
    }
    for (uint32_t i = 0; i < length(text); ++i)
        insertKmer(me, text[i], bins[i]);
}

/*!
 * \brief Adds all fasta files from a directory to the respective bins.
 * \param me The BinningDirectory instance.
 * \param baseDir The directory containing the fasta files in a "bins" subdirectory.
 * \param threads Number of threads to use.
 *
 * The fasta files are expected to follow the pattern <baseDir>/bins/bin_xxxx.fasta, where xxxx stands for the bin
 * number. All bin numbers must have the same number of digits as the total number of bins. E.g. for 8192 bins, the
 * bins are expected to be named bin_0000.fasta, bin_0001.fasta, ..., bin_8191.fasta; or for 64 bins: bin_00.fasta,
 * bin_01.fasta, ..., bin_63.fasta.
 * Up to <threads> fasta files are added to the bitvector at the same time.
 */
// template<typename TValue, typename TSpec, typename TBitvector>
// inline void insertKmerDir(BinningDirectory<TValue, TSpec, TBitvector> &  me, const char * baseDir, uint8_t threads)
// {
//     Semaphore thread_limiter(threads);
//     // std::mutex mtx;
//     std::vector<std::future<void>> tasks;
//
//     uint32_t bins = me.noOfBins;
//     for (uint8_t c = 0; c < me.bitvector.noOfChunks; ++c)
//     {
//         me.bitvector.decompress(c);
//         for(uint32_t i = 0; i < bins; ++i)
//         {
//             CharString file(baseDir);
//             append(file, CharString(std::to_string(bins)));
//             append(file, CharString{"/bins/bin_"});
//             append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
//             append(file, CharString(".fasta"));
//             tasks.emplace_back(
//                 std::async(std::launch::async, [=, &thread_limiter, &me] { // &mtx
//                     Critical_section _(thread_limiter);
//                     insertKmer(me, toCString(file), i, true, c);
//                     // mtx.lock();
//                     // std::cerr << "IBF Bin " << i << " done." << '\n';
//                     // mtx.unlock();
//                 })
//             );
//         }
//
//         for (auto &&task : tasks){
//             task.get();
//         }
//         me.bitvector.compress(c);
//     }
// }

/*!
 * \brief Calculates the k-mer counts of a given text.
 * \param me The BinningDirectory instance.
 * \param counts Vector of length binNo to save counts to.
 * \param text A single text to count all contained k-mers for.
 */
template<typename THashCount, typename TSpec, typename TConfig>
inline void count(BinningDirectory<TSpec, TConfig> &  me, std::vector<TNoOfBins> & counts, String<typename TConfig::TValue> const & text)
{
    typedef typename TConfig::THash THash;
    if (std::is_same<THash, Normal>::value || is_offset<THash>::value)
    {
        me.template count<THashCount>(counts, text);
    }
    else
    {
        me.template count<THash>(counts, text);
    }
}

template<typename TSpec, typename TConfig>
inline void count(BinningDirectory<TSpec, TConfig> &  me, std::vector<TNoOfBins> & counts, String<typename TConfig::TValue> const & text)
{
    typedef typename TConfig::THash THash;
    me.template count<THash>(counts, text);
}

/*!
 * \brief Returns the k-mer counts of a given text.
 * \param me The BinningDirectory instance.
 * \param text A single text to count all contained k-mers for.
 * \returns std::vector<uint64_t> of size binNo containing counts.
 */
template<typename THashCount, typename TSpec, typename TConfig>
inline std::vector<TNoOfBins> count(BinningDirectory<TSpec, TConfig> &  me, String<typename TConfig::TValue> const & text)
{
    typedef typename TConfig::THash THash;
    std::vector<TNoOfBins> counts(me.noOfBins, 0);

    if (std::is_same<THash, Normal>::value || is_offset<THash>::value)
    {
        count<THashCount>(me, counts, text);
    }
    else
    {
        count<THash>(me, counts, text);
    }

    return counts;
}

template<typename TSpec, typename TConfig>
inline std::vector<TNoOfBins> count(BinningDirectory<TSpec, TConfig> &  me, String<typename TConfig::TValue> const & text)
{
    std::vector<TNoOfBins> counts(me.noOfBins, 0);
    count(me, counts, text);
    return counts;
}

/*!
 * \brief Returns the number of bins.
 * \param me The BinningDirectory instance.
 * \returns Value<BinningDirectory<TValue, TSpec> >::Type Number of bins.
 */
template<typename TSpec, typename TConfig>
inline TNoOfBins getNumberOfBins(BinningDirectory<TSpec, TConfig> &  me)
{
    return me.noOfBins;
}

/*!
 * \brief Returns the k-mer size.
 * \param me The BinningDirectory instance.
 * \returns Value<BinningDirectory<TValue, TSpec> >::Type k-mer size.
 */
template<typename TSpec, typename TConfig>
inline TKmerSize getKmerSize(BinningDirectory<TSpec, TConfig> &  me)
{
    return me.kmerSize;
}

/*!
 * \brief Reads the metadata.
 * \param me The BinningDirectory instance.
 */
template<typename TSpec, typename TConfig>
inline void getMetadata(BinningDirectory<TSpec, TConfig> &  me)
{
    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------
    me.noOfBits = me.bitvector.noOfBits;

    TNoOfBits metadataStart = me.bitvector.noOfBits;
    me.noOfBins = me.bitvector.get_int(metadataStart);
    me.noOfHashFunc = me.bitvector.get_int(metadataStart+64);
    me.kmerSize = me.bitvector.get_int(metadataStart+128);
}

/*!
 * \brief Writes the metadata.
 * \param me The BinningDirectory instance.
 */
template<typename TSpec, typename TConfig>
inline void setMetadata(BinningDirectory<TSpec, TConfig> &  me)
{
    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------

    TNoOfBits metadataStart = me.noOfBits;

    me.bitvector.set_int(metadataStart, me.noOfBins);
    me.bitvector.set_int(metadataStart + 64, me.noOfHashFunc);
    me.bitvector.set_int(metadataStart+128, me.kmerSize);
}

/*!
 * \brief Returns the of the filter vector in MB.
 * \param me The BinningDirectory instance.
 * \returns double filter vector size in MB.
 */
template<typename TSpec, typename TConfig>
inline double size(BinningDirectory<TSpec, TConfig> &  me)
{
    return me.bitvector.size_in_mega_bytes();
}

/*!
 * \brief Writes the filter vector to a file.
 * \param me The BinningDirectory instance.
 * \param fileName Name of the file to write to.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TSpec, typename TConfig>
inline bool store(BinningDirectory<TSpec, TConfig> &  me, CharString fileName)
{
    setMetadata(me);
    return me.bitvector.store(fileName);
}

/*!
 * \brief Loads the filter vector from a file.
 * \param me The BinningDirectory instance.
 * \param fileName Name of the file to read from.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TSpec, typename TConfig>
inline bool retrieve(BinningDirectory<TSpec, TConfig> &  me, CharString fileName)
{
    me.bitvector.retrieve(fileName);
    getMetadata(me);
    me.init();
    return true;
}

constexpr unsigned long long int operator""_g ( unsigned long long int g )
{
    return g*8*1024*1024*1024;
}

constexpr unsigned long long int operator""_m ( unsigned long long int m )
{
    return m*8*1024*1024;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_BASE_H_
