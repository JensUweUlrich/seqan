#include <seqan/kmer.h>
#include <iterator>
using namespace seqan;

int main()
{
    KmerFilter<Dna, InterleavedBloomFilter, Uncompressed> ibf (64, 3, 17, (1ULL<<35));
    retrieve(ibf, CharString{"/group/ag_abi/seiler/benchmark/filter/256_17_35_Uncompressed_ibf.filter"});
    uint16_t k{17};
    uint16_t e{2};
    uint16_t bins{256};
    CharString baseDir{"/group/ag_abi/seiler/benchmark/"};

    for(int32_t i = 0; i < bins; ++i)
    {
        CharString file(baseDir);
        append(file, CharString(std::to_string(bins)));
        append(file, CharString{"/reads/bin_"});
        append(file, CharString(std::string(numDigits(bins)-numDigits(i), '0') + (std::to_string(i))));
        append(file, CharString(".fastq"));

        CharString id;
        String<Dna> seq;
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
            readRecord(id, seq, seqFileIn);

            auto res = select(ibf, seq, 100-k+1 - k*e);

            c = count(res.begin(), res.end(), true);
            if (c > 1)
            {
                std::cout << "False positive. Read from bin " << i << " matched in bins ";
                auto iter = res.begin();
                while((iter = std::find_if(iter, res.end(), true)) != res.end())
                    std::cout << std::distance(iter, res.begin()) << ' ';
                std::cout << '\n';
            }
        }
    }
}
