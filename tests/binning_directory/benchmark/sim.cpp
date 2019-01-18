#include <seqan/seq_io.h>
#include <random>
using namespace seqan;

CharString baseDir{"/srv/hdd/seiler/IBF/"};
uint8_t maxErrors{3};
uint16_t readLength{100};
uint32_t noOfReads{1UL<<20};
uint8_t noOfHaplotypes{16};
std::random_device rd;
std::mt19937 rng(rd());
std::uniform_real_distribution<double> dist(0, 1);

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

int main()
{
    for(uint16_t noOfBins : {64, 256, 1024, 8192})
    // for(uint16_t noOfBins : {64})
    {
        // Since numerator and denominator are powers of two, each bin should get an equal number of reads
        uint32_t readsPerBin = noOfReads / noOfBins;
        // Since numerator and denominator are powers of two, each haplotype should get an equal number of reads
        uint32_t readsPerHaplotype = readsPerBin / noOfHaplotypes;

        for(int32_t bin = 0; bin < noOfBins; ++bin)
        {
            CharString fileIn(baseDir);
            append(fileIn, CharString(std::to_string(noOfBins)));
            append(fileIn, CharString{"/bins/bin_"});
            append(fileIn, CharString(std::string(numDigits(noOfBins)-numDigits(bin), '0') + (std::to_string(bin))));
            append(fileIn, CharString(".fasta"));

            CharString fileOut(baseDir);
            append(fileOut, CharString(std::to_string(noOfBins)));
            append(fileOut, CharString{"/reads/bin_"});
            append(fileOut, CharString(std::string(numDigits(noOfBins)-numDigits(bin), '0') + (std::to_string(bin))));
            append(fileOut, CharString(".fastq"));

            SeqFileIn seqFileIn;
            if (!open(seqFileIn, toCString(fileIn)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(fileIn));
                std::cerr << msg << '\n';
                throw toCString(msg);
            }

            SeqFileOut seqFileOut;
            if (!open(seqFileOut, toCString(fileOut)))
            {
                CharString msg = "Unable to open contigs file: ";
                append(msg, CharString(fileOut));
                std::cerr << msg << '\n';
                throw toCString(msg);
            }

            CharString id;
            String<Dna> seq;

            while(!atEnd(seqFileIn))
            {
                readRecord(id, seq, seqFileIn);
                uint32_t refLength = length(seq);
                std::uniform_int_distribution<uint32_t> readPos(0, refLength - readLength);
                for(uint32_t r = 0; r < readsPerHaplotype; ++r)
                {
                    uint32_t pos = readPos(rng);
                    DnaString segment = infixWithLength(seq, pos, readLength);
                    std::uniform_int_distribution<uint16_t> errorPos(0, readLength-1);
                    for(uint8_t e = 0; e < maxErrors; ++e)
                    {
                        uint32_t epos = errorPos(rng);
                        Dna currentBase = segment[epos];
                        Dna newBase = currentBase;
                        while (newBase == currentBase)
                            newBase = Dna(static_cast<int>(dist(rng) / 0.25));
                        segment[epos] = newBase;
                        if (dist(rng) < 1.0/static_cast<double>(maxErrors))
                            break;
                    }
                    writeRecord(seqFileOut, id, segment);
                }
            }
        }
    }
}
