template<typename TValue>
class BDHashBase
{
public:

    TNoOfChunks chunks{1};
    std::vector<TNoOfChunks> chunkMap{0};
    Shape<TValue, SimpleShape> kmerShape;
    TKmerSize k;

    inline void resize(TKmerSize newK)
    {
        k = newK;
        seqan::resize(kmerShape, k);
    }

    template<typename TIterator>
    inline void hashInit(TIterator it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIterator>
    inline auto hashNext(TIterator it)
    {
        seqan::hashNext(kmerShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }
}
