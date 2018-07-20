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

template<typename TSpec>
struct Bitvector;

// --------------------------------------------------------------------------
// Class BinningDirectory
// --------------------------------------------------------------------------

//!\brief The BinningDirectory class.
template<typename TValue = Dna, typename TShape = Shape<TValue, SimpleShape>, typename TSpec = InterleavedBloomFilter, typename TBitvector = Uncompressed>
class BinningDirectory;

// ==========================================================================
// Metafunctions
// ==========================================================================

//!\brief Type definition for variables.
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
struct Value<BinningDirectory<TValue, TShape, TSpec, TBitvector> >
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
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline void insertKmer(BinningDirectory<TValue, TShape, TSpec, TBitvector> & me, String<TValue> const & text, TNoOfBins binNo)
{
    me.insertKmer(text, binNo);
}

/*!
 * \brief Sets the vectors for given bins to 0.
 * \param me The BinningDirectory instance.
 * \param bins A vector containing the bin numbers.
 * \param threads The number of threads to use.
 */
template<typename TValue, typename TSpec, typename TBitvector, typename TInt>
inline void clear(BinningDirectory<TValue, TSpec, TBitvector> &  me, std::vector<TNoOfBins> & bins, TInt&& threads)
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
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline void insertKmer(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me, const char * fastaFile, TNoOfBins binNo)
{
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
template<typename TValue, typename TShape, typename TSpec,  typename TBitvector>
inline void count(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me, std::vector<TNoOfBins> & counts, String<TValue> const & text)
{
    me.count(counts, text);
}

/*!
 * \brief Returns the k-mer counts of a given text.
 * \param me The BinningDirectory instance.
 * \param text A single text to count all contained k-mers for.
 * \returns std::vector<uint64_t> of size binNo containing counts.
 */
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline std::vector<TNoOfBins> count(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me, String<TValue> const & text)
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
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline TNoOfBins getNumberOfBins(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me)
{
    return me.noOfBins;
}

/*!
 * \brief Returns the k-mer size.
 * \param me The BinningDirectory instance.
 * \returns Value<BinningDirectory<TValue, TSpec> >::Type k-mer size.
 */
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline TKmerSize getKmerSize(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me)
{
    return me.kmerSize;
}

/*!
 * \brief Reads the metadata.
 * \param me The BinningDirectory instance.
 */
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline void getMetadata(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me)
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
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline void setMetadata(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me)
{
    //-------------------------------------------------------------------
    //| kmer_size | n_hash_func | n_bins |              bf              |
    //-------------------------------------------------------------------

    TNoOfBits metadataStart = me.noOfBits;

    // TODO also store TValue (alphabet)
    me.bitvector.set_int(metadataStart, me.noOfBins);
    me.bitvector.set_int(metadataStart + 64, me.noOfHashFunc);
    me.bitvector.set_int(metadataStart+128, me.kmerSize);
}

/*!
 * \brief Returns the of the filter vector in MB.
 * \param me The BinningDirectory instance.
 * \returns double filter vector size in MB.
 */
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline double size(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me)
{
    return me.bitvector.size_in_mega_bytes();
}

/*!
 * \brief Writes the filter vector to a file.
 * \param me The BinningDirectory instance.
 * \param fileName Name of the file to write to.
 * \returns bool Indicates if the operation was successful.
 */
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline bool store(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me, CharString fileName)
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
template<typename TValue, typename TShape, typename TSpec, typename TBitvector>
inline bool retrieve(BinningDirectory<TValue, TShape, TSpec, TBitvector> &  me, CharString fileName)
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
