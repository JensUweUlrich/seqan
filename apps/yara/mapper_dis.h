// ==========================================================================
//                                 classify-reads
// ==========================================================================
// Copyright (c) 2017-2022, Temesgen H. Dadi, FU Berlin
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
//     * Neither the name of Temesgen H. Dadi or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL TEMESGEN H. DADI OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#ifndef APP_YARA_DIS_MAPPER_H_
#define APP_YARA_DIS_MAPPER_H_

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class DisOptions
// --------------------------------------------------------------------------

struct DisOptions : public Options
{
    CharString          superContigsIndicesFile;
    CharString          superOutputFile;
    uint32_t            NUM_OF_BINS = 5;
//    uint32_t            NUM_OF_BINS = 64;

//    Pair<CharString>    superReadsFile;

//    uint32_t            max_errors = 3;
//    CharString          kmer_index_file;
//    uint32_t error_diff;
//
//    void reload()
//    {
//        error_diff =  K_MER_LENGTH - 1 + (max_errors * K_MER_LENGTH);
//    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function set_current_index_file()
// ----------------------------------------------------------------------------

inline void set_current_index_file(Options & yaraOptions, DisOptions const & options, uint32_t const file_no)
{
    // Get the current file name
    yaraOptions.contigsIndexFile = options.superContigsIndicesFile;
//    if (file_no < 10)
//        append(yaraOptions.contigsIndexFile, "0");
    append(yaraOptions.contigsIndexFile, std::to_string(file_no));
//    append(yaraOptions.contigsIndexFile, "-genomes");
}

template <typename TOptions>
inline void set_current_index_file(TOptions & options, uint32_t const file_no)
{
    // Get the current file name
    options.contigsIndexFile = options.superContigsIndicesFile;
//    if (file_no < 10)
//        append(options.contigsIndexFile, "0");
    append(options.contigsIndexFile, std::to_string(file_no));
//    append(options.contigsIndexFile, "-genomes");
}

inline void set_output_file(Options & yaraOptions, DisOptions const & options, uint32_t const file_no)
{
    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(options.superOutputFile); ++first_dot_pos)
    {
        if(options.superOutputFile[first_dot_pos] == '.')
            break;
    }
    // Get the current file name
    yaraOptions.outputFile = prefix(options.superOutputFile, first_dot_pos);
    append(yaraOptions.outputFile, "_");
    if (file_no < 10)
        append(yaraOptions.outputFile, "0");
    append(yaraOptions.outputFile, std::to_string(file_no));
    append(yaraOptions.outputFile, suffix(options.superOutputFile, first_dot_pos));
}

template <typename TOptions>
inline void set_output_file(TOptions & options, uint32_t const file_no)
{
    uint32_t first_dot_pos = 0;
    for (; first_dot_pos < length(options.superOutputFile); ++first_dot_pos)
    {
        if(options.superOutputFile[first_dot_pos] == '.')
            break;
    }
    // Get the current file name
    options.outputFile = prefix(options.superOutputFile, first_dot_pos);
    append(options.outputFile, "_");
    if (file_no < 10)
        append(options.outputFile, "0");
    append(options.outputFile, std::to_string(file_no));
    append(options.outputFile, suffix(options.superOutputFile, first_dot_pos));
}

// ----------------------------------------------------------------------------
// Function copyMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copyMatches(Mapper<TSpec, TMainConfig> & target, Mapper<TSpec, TConfig> & source, uint32_t const & contigOffset)
{
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatch       TMatch;

    TMatch currentMatch;

    uint32_t matchCount = length(source.matchesByCoord);
    for (uint32_t i = 0; i<matchCount; ++i)
    {
        currentMatch.readId        =source.matchesByCoord[i].readId;
        currentMatch.contigId      =source.matchesByCoord[i].contigId + contigOffset;
        currentMatch.isRev         =source.matchesByCoord[i].isRev;
        currentMatch.contigBegin   =source.matchesByCoord[i].contigBegin;
        currentMatch.contigEnd     =source.matchesByCoord[i].contigEnd;
        currentMatch.errors        =source.matchesByCoord[i].errors;
        appendValue(target.matchesByCoord, currentMatch);
        setSeedErrors(target.ctx, currentMatch.readId, currentMatch.errors);
        setMinErrors(target.ctx, currentMatch.readId, currentMatch.errors);
        setMapped(target.ctx, currentMatch.readId);
        setPaired(target.ctx, currentMatch.readId);
    }
}
// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void mapReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, uint32_t const & contigOffset)
{
    swap(me.reads.seqs, mainMapper.reads.seqs);
    swap(me.reads.names, mainMapper.reads.names);
    me.stats.loadedReads += getReadsCount(me.reads.seqs);
    _mapReadsImpl(me, mainMapper, me.reads.seqs, contigOffset);
    swap(me.reads.seqs, mainMapper.reads.seqs);
    swap(me.reads.names, mainMapper.reads.names);
}

template <typename TSpec, typename TConfig, typename TMainConfig, typename TReadSeqs>
inline void _mapReadsImpl(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, TReadSeqs & readSeqs, uint32_t const & contigOffset)
{
    initReadsContext(me, readSeqs);
    initSeeds(me, readSeqs);

    collectSeeds<0>(me, readSeqs);
    findSeeds<0>(me, 0);
    classifyReads(me);
    collectSeeds<1>(me, readSeqs);
    collectSeeds<2>(me, readSeqs);
    findSeeds<0>(me, 1);
    findSeeds<0>(me, 2);
    rankSeeds(me);
    reserveMatches(me);
    extendHits<0>(me, 0);
    extendHits<0>(me, 1);
    extendHits<0>(me, 2);
    clearSeeds(me);
    clearHits(me);

    initSeeds(me, readSeqs);
    collectSeeds<1>(me, readSeqs);
    findSeeds<1>(me, 1);
    collectSeeds<2>(me, readSeqs);
    findSeeds<1>(me, 2);
    rankSeeds(me);
    // TODO(esiragusa): filter out hits with distance < 1.
    extendHits<1>(me, 1);
    extendHits<1>(me, 2);
    clearSeeds(me);
    clearHits(me);

    if (me.options.sensitivity > LOW)
    {
        initSeeds(me, readSeqs);
        collectSeeds<2>(me, readSeqs);
        findSeeds<2>(me, 2);
        rankSeeds(me);
        // TODO(esiragusa): filter out hits with distance < 2.
        extendHits<2>(me, 2);
        clearHits(me);
        clearSeeds(me);
    }
    aggregateMatches(me, readSeqs);
    rankMatches(me, readSeqs);
    if (me.options.verifyMatches)
        verifyMatches(me);
    alignMatches(me);
//    writeMatches(mainMapper);
    copyMatches(mainMapper, me, contigOffset);

//    clearMatches(me);
//    clearAlignments(me);
}

template <typename TSpec, typename TConfig, typename TMainConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig> & mainMapper, uint32_t const & contigOffset)
{
    configureThreads(me);
    loadContigs(me);
    loadContigsIndex(me);
    // Write on the main output file.
    mapReads(me, mainMapper, contigOffset);
}


// ----------------------------------------------------------------------------
// Function spawnMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TContigsSum,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void spawnMapper(Options const & options,
                 Mapper<TSpec, TMainConfig> & mainMapper,
                 uint32_t const & contigOffset,
                 TThreading const & /*threading*/,
                 TSequencing const & /*sequencing*/,
                 TSeedsDistance const & /*distance*/)
{

    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TConfig;
    Mapper<void, TConfig> mapper(options);
    runMapper(mapper, mainMapper, contigOffset);
}
// ----------------------------------------------------------------------------
// Function configureMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize,
          typename TContigsLen,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     uint32_t const & contigOffset,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsSum <= MaxValue<uint32_t>::VALUE)
    {
        spawnMapper<TContigsSize, TContigsLen, uint32_t>(options, mainMapper, contigOffset, threading, sequencing, distance);
    }
    else
    {
        spawnMapper<TContigsSize, TContigsLen, uint64_t>(options, mainMapper, contigOffset, threading, sequencing, distance);
    }
}

template <typename TContigsSize,
          typename TSpec,
          typename TMainConfig,
          typename TThreading,
          typename TSequencing,
          typename TSeedsDistance>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     uint32_t const & contigOffset,
                     TThreading const & threading,
                     TSequencing const & sequencing,
                     TSeedsDistance const & distance)
{
    if (options.contigsMaxLength <= MaxValue<uint32_t>::VALUE)
    {
        configureMapper<TContigsSize, uint32_t>(options, mainMapper, contigOffset, threading, sequencing, distance);
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        configureMapper<TContigsSize, uint64_t>(options, mainMapper, contigOffset, threading, sequencing, distance);
#else
        throw RuntimeError("Maximum contig length exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}

template <typename TSpec, typename TMainConfig>
void configureMapper(Options const & options,
                     Mapper<TSpec, TMainConfig> & mainMapper,
                     uint32_t const & contigOffset)
{
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading       TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSequencing      TSequencing;
    typedef typename MapperTraits<TSpec, TMainConfig>::TSeedsDistance   TSeedsDistance;

    if (options.contigsSize <= MaxValue<uint8_t>::VALUE)
    {
        configureMapper<uint8_t>(options, mainMapper, contigOffset, TThreading(), TSequencing(), TSeedsDistance());
    }
    else if (options.contigsSize <= MaxValue<uint16_t>::VALUE)
    {
        configureMapper<uint16_t>(options, mainMapper, contigOffset, TThreading(), TSequencing(), TSeedsDistance());
    }
    else
    {
#ifdef YARA_LARGE_CONTIGS
        configureMapper<uint32_t>(options, mainMapper, contigOffset, TThreading(), TSequencing(), TSeedsDistance());
#else
        throw RuntimeError("Maximum number of contigs exceeded. Recompile with -DYARA_LARGE_CONTIGS=ON.");
#endif
    }
}

// ----------------------------------------------------------------------------
// Function runMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TContigs>
inline void runDisMapper(Mapper<TSpec, TConfig> & me, DisOptions & disOptions, TContigs & allContigs, std::vector<uint32_t> const & contigOffsets)
{

    Timer<double> timer;

    start(timer);

    configureThreads(me);

    if (me.options.verbose > 1) printRuler(std::cerr);


    assign(me.contigs.names, allContigs.names);
    assign(me.contigs.seqs, allContigs.seqs);

    // Open output file and write header.
    openOutputFile(me);
    openReads(me);

    // Process reads in blocks.
    // load reads here

    // classify reads here
    // create mappers and run them on subsets
    // Process reads in blocks.
    while (true)
    {
        if (me.options.verbose > 1) printRuler(std::cerr);
        loadReads(me);
        if (empty(me.reads.seqs)) break;
        initReadsContext(me, me.reads.seqs);
        for (uint32_t i=0; i < disOptions.NUM_OF_BINS; ++i)
        {
            Options options = me.options;
            set_current_index_file(options, disOptions, i);
            if (!openContigsLimits(options))
                throw RuntimeError("Error while opening reference file.");
            configureMapper<TSpec, TConfig>(options, me, contigOffsets[i]);
        }

        aggregateMatches(me, me.reads.seqs);
        rankMatches(me, me.reads.seqs);
        if (me.options.verifyMatches)
            verifyMatches(me);
        uint32_t matchCount = length(me.matchesByCoord);
        uint32_t cigar_count = length(me.cigars);
        uint32_t cigarSet_count = length(me.cigarSet);
        alignMatches(me);
        matchCount = length(me.matchesByCoord);
        cigar_count = length(me.cigars);
        cigarSet_count = length(me.cigarSet);
        writeMatches(me);
        clearMatches(me);
        clearAlignments(me);

        clearReads(me);
    }
    closeReads(me);
    closeOutputFile(me);
    stop(timer);
    if (me.options.verbose > 0)
        printStats(me, timer);

 }

// ----------------------------------------------------------------------------
// Function spawnDisMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
typename TThreading, typename TSequencing, typename TSeedsDistance, typename TContigs>
inline void spawnDisMapper(DisOptions & disOptions,
                           TContigs & allContigs,
                           std::vector<uint32_t> const & contigOffsets,
                           TThreading const & /* tag */,
                           TSequencing const & /* tag */,
                           TSeedsDistance const & /* tag */)
{
    disOptions.outputFile = disOptions.superOutputFile;
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TMainConfig;
    Mapper<void, TMainConfig> disMapper(disOptions);

    runDisMapper(disMapper, disOptions, allContigs, contigOffsets);
}

#endif  // #ifndef APP_YARA_MAPPER_H_
