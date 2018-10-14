#ifndef INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_BASE_H_
#define INCLUDE_SEQAN_BINNING_DIRECTORY_BINNING_DIRECTORY_BASE_H_

#include <sdsl/bit_vectors.hpp>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>
#include <valarray>
#include <algorithm>
#include <future>
#include <mutex>
#include <type_traits>

namespace seqan {

// ====== Binning Directory Implementations ======
template<typename TSpec, typename TConfig>
class BinningDirectory;

struct InterleavedBloomFilter_; // TSpec - Default
typedef Tag<InterleavedBloomFilter_> InterleavedBloomFilter;

struct DirectAddressing_; // TSpec
typedef Tag<DirectAddressing_> DirectAddressing;

// ====== BitVector Implementations ======
struct Uncompressed_; // TSpec - Default
typedef Tag<Uncompressed_> Uncompressed;

struct Compressed_; // TSpec
typedef Tag<Compressed_> Compressed;

struct OnDisk_; // TMemory
typedef Tag<OnDisk_> OnDisk;

struct InMemory_; // TMemory - Default
typedef Tag<InMemory_> InMemory;

template<typename TSpec = Uncompressed, typename TMemory = InMemory>
struct Bitvector;

// ====== Hash Implementations ======
struct Normal_;
typedef Tag<Normal_> Normal;

struct Offset_;
typedef Tag<Offset_> Offset;

struct Minimizer_;
typedef Tag<Minimizer_> Minimizer;

template<typename TValue, typename TSpec>
struct BDHash;

// ====== Config struct ======
template<typename TValue_ = Dna, typename THash_ = Normal, typename TBitvector_ = Bitvector<>>
struct BDConfig
{
    typedef TValue_      TValue;
    typedef THash_       THash;
    typedef TBitvector_  TBitvector;
};

// ====== Metafunctions ======
// template<typename TSpec, typename TConfig>
// struct noOfChunks<BinningDirectory<TSpec, TConfig> >
// {
//     typedef uint8_t  Type;
// };
// template<typename TSpec, typename TConfig>
// struct kmerSize<BinningDirectory<TSpec, TConfig> >
// {
//     typedef uint16_t  Type;
// };
// template<typename TSpec, typename TConfig>
// struct noOfBins<BinningDirectory<TSpec, TConfig> >
// {
//     typedef uint32_t  Type;
// };

// ====== Typedefs ======
typedef uint8_t  TNoOfChunks;
typedef uint16_t TKmerSize;
typedef uint32_t TNoOfBins;
typedef std::vector<std::tuple<uint64_t, TNoOfChunks>> THashVector;

// ====== Functions ======

// TODO How to handle chunks?
// TODO Should chunkMap creation be preprocessing?
// TODO be able to specify different Hashing approach, i.e. Offset and Normal are interchangeable for insert, count, select
// TODO Should k, o, w be metadata?

// Inserts a given text into a bin
void insertKmer(BinningDirectory<TSpec, TConfig> & bd, TString const & text, TNoOfBins bin) {}

// Inserts a given text into a bin, only access one chunk
void insertKmer(BinningDirectory<TSpec, TConfig> & bd, TString const & text, TNoOfBins bin, TNoOfChunks chunk)

// Inserts all kmers from a file into a bin. Will iterate over all chunks.
void insertKmer(BinningDirectory<TSpec, TConfig> & bd, const char * file, TNoOfBins bin) // TODO const char * ?

// Inserts all kmers from StringSet into corresponding bins. Used for chunkMap
void insertKmer(BinningDirectory<TSpec, TConfig> & bd, StringSet<TString> const & text, std::vector<TNoOfBins> & bin)

// Parallel clearing (setting to 0) of specified bins
void clear(BinningDirectory<TSpec, TConfig> & bd, std::vector<TNoOfBins> & bins, uint16_t threads)

// Configures a chunk map
void configureChunkMap(BinningDirectory<TSpec, TConfig> & bd)

// Returns a count vector
std::vector<uint64_t> count(BinningDirectory<TSpec, TConfig> const & bd, TString const & text)

// Returns a count vector for a chunk
std::vector<uint64_t> count(BinningDirectory<TSpec, TConfig> const & bd, TString const & text, TNoOfChunks chunk)

// Counts for given count vector
void count(BinningDirectory<TSpec, TConfig> const & bd, std::vector<uint64_t> & counts, TString const & text)

// Counts for given count vector in a chunk
void count(BinningDirectory<TSpec, TConfig> const & bd, std::vector<uint64_t> & counts, TString const & text, TNoOfChunks chunk)

// Returns a select vector
std::vector<bool> select(BinningDirectory<TSpec, TConfig> const & bd, TString const & text, uint64_t threshold)

// Returns a select vector for a chunk
std::vector<bool> select(BinningDirectory<TSpec, TConfig> const & bd, TString const & text, uint64_t threshold, TNoOfChunks chunk)

// Selects for given select vector
void select(BinningDirectory<TSpec, TConfig> const & bd, std::vector<bool> & selected, TString const & text, uint64_t threshold)

// Selects for given select vector in a chunk
void select(BinningDirectory<TSpec, TConfig> const & bd, std::vector<bool> & selected, TString const & text, uint64_t threshold, TNoOfChunks chunk)

// Set Metadata
void setMetadata(BinningDirectory<TSpec, TConfig> & bd)

// Get Metadata
void getMetadata(BinningDirectory<TSpec, TConfig> & bd)

// Returns size
double size(BinningDirectory<TSpec, TConfig> & bd)

// Stores bitvector and chunkMap
bool store(BinningDirectory<TSpec, TConfig> & bd, const char * file)

// Loads bitvector and chunkMap
bool retrieve(BinningDirectory<TSpec, TConfig> & bd, const char * file) // TODO may actually be called load if we use const char * instead of CharString

// Literal for GiB -> Bit
constexpr unsigned long long int operator""_g(unsigned long long int GiB)
// Literal for MiB -> Bit
constexpr unsigned long long int operator""_m(unsigned long long int MiB)
// Literal for KiB -> Bit
constexpr unsigned long long int operator""_k(unsigned long long int KiB)
