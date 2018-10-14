template<>
class BDHash<Dna, Normal> : BDHashBase<Dna>
{
public:

    BDHash() {}

    BDHash(TKmerSize k):
        this->k(k)
    {
        this->resize(this->k);
    }

    BDHash(TKmerSize k, TNoOfChunks chunks):
        this->k(k),
        this->chunks(chunks)
    {
        this->resize(this->k);
        chunkMap.resize(chunks);
        std::iota(std::begin(chunkMap), std::end(chunkMap), 0);
    }

    BDHash(TKmerSize k, TNoOfChunks chunks, std::vector<TNoOfChunks> chunkMap)
        this->k(k),
        this->chunks(chunks),
        this->chunkMap(chunkMap)
    {
        this->resize(this->k);
    }

    BDHash(BDHash<Dna> & other)  = default;
    BDHash(BDHash<Dna> && other) = default;
    BDHash<Dna> & operator=(BDHash<Dna> & other)  = default;
    BDHash<Dna> & operator=(BDHash<Dna> && other) = default;

    template<typename TString>
    inline uint64_t threshold (TString const & text, uint64_t e)
    {
        return std::max(0ULL, seqan::length(text) - this->k * (e + 1) - 1);
    }

    template<typename TString>
    THashVector getHash(TString const & text)
    {
        if (this->k > seqan::length(text))
            return THashVector {};

        uint64_t possible = seqan::length(text) - this->k - 1;

        THashVector kmerHashes;
        kmerHashes.reserve(possible);

        auto it = begin(text);
        this->hashInit(it);

        for (uint64_t i = 0; i < possible; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);

            kmerHashes.push_back(std::make_tuple(kmerHash >> significantBits,
                                                 this->chunkMap[(kmerHash & (this->effectiveChunks - 1))]));

            ++it;
        }

        return kmerHashes;
    }
}

template<typename TString>
class BDHash<TString, Normal> : BDHashBase<TString>
{
public:

    BDHash() {}

    BDHash(TKmerSize k):
        this->k(k)
    {
        this->resize(this->k);
    }

    BDHash(TKmerSize k, TNoOfChunks chunks):
        this->k(k),
        this->chunks(chunks)
    {
        this->resize(this->k);
        chunkMap.resize(chunks);
        std::iota(std::begin(chunkMap), std::end(chunkMap), 0);
    }

    BDHash(TKmerSize k, TNoOfChunks chunks, std::vector<TNoOfChunks> chunkMap)
        this->k(k),
        this->chunks(chunks),
        this->chunkMap(chunkMap)
    {
        this->resize(this->k);
    }

    BDHash(BDHash<Dna> & other)  = default;
    BDHash(BDHash<Dna> && other) = default;
    BDHash<Dna> & operator=(BDHash<Dna> & other)  = default;
    BDHash<Dna> & operator=(BDHash<Dna> && other) = default;

    template<typename TString>
    inline uint64_t threshold (TString const & text, uint64_t e)
    {
        return std::max(0ULL, seqan::length(text) - this->k * (e + 1) - 1);
    }

    template<typename TString>
    THashVector getHash(TString const & text)
    {
        if (this->k > seqan::length(text))
            return THashVector {};

        uint64_t possible = seqan::length(text) - this->k - 1;

        THashVector kmerHashes;
        kmerHashes.reserve(possible);

        auto it = begin(text);
        this->hashInit(it);

        for (uint64_t i = 0; i < possible; ++i)
        {
            uint64_t kmerHash = this->hashNext(it);

            kmerHashes.push_back(std::make_tuple(kmerHash >> significantBits,
                                                 this->chunkMap[(kmerHash & (this->effectiveChunks - 1))]));

            ++it;
        }

        return kmerHashes;
    }
}
