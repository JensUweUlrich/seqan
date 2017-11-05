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

    CharString              superContigsIndicesFile;
    CharString              superOutputFile;
    uint32_t                NUM_OF_BINS = 5;
    std::vector<uint32_t>   contigOffsets;
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
// Function appendStats()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void appendStats(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper)
{
    mainMapper.stats.loadContigs    += childMapper.stats.loadContigs;
    mainMapper.stats.loadReads      += childMapper.stats.loadReads;
    mainMapper.stats.collectSeeds   += childMapper.stats.collectSeeds;
    mainMapper.stats.findSeeds      += childMapper.stats.findSeeds;
    mainMapper.stats.classifyReads  += childMapper.stats.classifyReads;
    mainMapper.stats.rankSeeds      += childMapper.stats.rankSeeds;
    mainMapper.stats.extendHits     += childMapper.stats.extendHits;
    mainMapper.stats.sortMatches    += childMapper.stats.sortMatches;
    mainMapper.stats.compactMatches += childMapper.stats.compactMatches;
    mainMapper.stats.selectPairs    += childMapper.stats.selectPairs;
    mainMapper.stats.verifyMatches  += childMapper.stats.verifyMatches;
    mainMapper.stats.alignMatches   += childMapper.stats.alignMatches;
    mainMapper.stats.writeMatches   += childMapper.stats.writeMatches;
}

// ----------------------------------------------------------------------------
// Function copyMatches()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void copyMatches(Mapper<TSpec, TMainConfig> & mainMapper, Mapper<TSpec, TConfig> & childMapper, uint32_t const & contigOffset)
{
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatch             TMatch;
    typedef typename MapperTraits<TSpec, TMainConfig>::TThreading         TThreading;
    typedef typename MapperTraits<TSpec, TMainConfig>::TMatchesAppender   TMatchesAppender;

    TMatchesAppender appender(mainMapper.matchesByCoord);

    uint32_t matchCount = length(childMapper.matchesByCoord);
    for (uint32_t i = 0; i < matchCount; ++i)
    {
        TMatch currentMatch;
        currentMatch.readId        = childMapper.matchesByCoord[i].readId;
        currentMatch.contigId      = childMapper.matchesByCoord[i].contigId + contigOffset;
        currentMatch.isRev         = childMapper.matchesByCoord[i].isRev;
        currentMatch.contigBegin   = childMapper.matchesByCoord[i].contigBegin;
        currentMatch.contigEnd     = childMapper.matchesByCoord[i].contigEnd;
        currentMatch.errors        = childMapper.matchesByCoord[i].errors;
        appendValue(appender, currentMatch, Generous(), TThreading());

        if(getMinErrors(mainMapper.ctx, currentMatch.readId) > currentMatch.errors)
        {
            setMinErrors(mainMapper.ctx, currentMatch.readId, currentMatch.errors);
        }

        if (!isPaired(mainMapper.ctx, currentMatch.readId) && isPaired(childMapper.ctx, currentMatch.readId))
        {
            setPaired(mainMapper.ctx, currentMatch.readId);
        }
        setMapped(mainMapper.ctx, currentMatch.readId);
    }
}
// ----------------------------------------------------------------------------
// Function filterLoadReads()
// ----------------------------------------------------------------------------
template <typename TReadSeqs>
inline bool isCandidate(TReadSeqs & seq, std::vector<bool>  & filter)
{

//    for (unsigned i = 0; i < length(seq) - 20 + 1; ++i)
//    {
//        hashNext(kmer_shape, begin(seq) + i);
//        if (countOccurrences(indices[bin_no], kmer_shape) > 0)
//            ++common_count;
//        if(common_count >= threshold)
//        {
//            candidate_bins.set(bin_no);
//            break;
//        }
//    }

    return true;
}
// ----------------------------------------------------------------------------
// Function filterLoadReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void filterLoadReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper)
{
    //replace with actual filters
    clear(me.reads.seqs);
    clear(me.reads.names);
    std::vector<bool> filter;
    uint32_t numReads = length(mainMapper.reads.names);
    for (uint32_t i = 0; i< numReads; ++i)
    {
        if (isCandidate(mainMapper.reads.seqs[i], filter))
        {
            appendValue(me.reads.seqs, mainMapper.reads.seqs[i]);
            appendValue(me.reads.names, mainMapper.reads.names[i]);
        }
    }
}

// ----------------------------------------------------------------------------
// Function mapReads()
// ----------------------------------------------------------------------------
template <typename TSpec, typename TConfig, typename TMainConfig>
inline void mapReads(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig>  & mainMapper, uint32_t const & contigOffset)
{
    swap(me.reads, mainMapper.reads);
    _mapReadsImpl(me, mainMapper, me.reads.seqs, contigOffset);
    swap(me.reads, mainMapper.reads);
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
    copyMatches(mainMapper, me, contigOffset);
    appendStats(mainMapper, me);
}

template <typename TSpec, typename TConfig, typename TMainConfig>
inline void runMapper(Mapper<TSpec, TConfig> & me, Mapper<TSpec, TMainConfig> & mainMapper, uint32_t const & contigOffset)
{
    loadContigs(me);
    loadContigsIndex(me);
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
// Function alignMatches()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TContigSeqs>
inline void alignMatches(Mapper<TSpec, TConfig> & me, TContigSeqs & contigSeqs)
{
    typedef MapperTraits<TSpec, TConfig>            TTraits;
    typedef MatchesAligner<LinearGaps, TTraits>     TLinearAligner;
    typedef MatchesAligner<AffineGaps , TTraits>    TAffineAligner;

    start(me.timer);
    setHost(me.cigarSet, me.cigars);
    typename TTraits::TCigarLimits cigarLimits;

    if (me.options.rabema)
        TLinearAligner aligner(me.cigarSet, cigarLimits, me.primaryMatches, contigSeqs, me.reads.seqs, me.options);
    else
        TAffineAligner aligner(me.cigarSet, cigarLimits, me.primaryMatches, contigSeqs, me.reads.seqs, me.options);

    stop(me.timer);
    me.stats.alignMatches += getValue(me.timer);

    if (me.options.verbose > 1)
        std::cerr << "Alignment time:\t\t\t" << me.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function loadAllContigs()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void loadAllContigs(Mapper<TSpec, TConfig> & mainMapper, DisOptions & disOptions)
{
    typedef typename MapperTraits<TSpec, TConfig>::TContigs          TContigs;

    start(mainMapper.timer);
    try
    {
        for (uint32_t i=0; i < disOptions.NUM_OF_BINS; ++i)
        {
            TContigs tmpContigs;
            CharString fileName = disOptions.superContigsIndicesFile;
            append(fileName, std::to_string(i));
            if (!open(tmpContigs, toCString(fileName), OPEN_RDONLY))
                throw RuntimeError("Error while opening reference file.");
            append(mainMapper.contigs.seqs, tmpContigs.seqs);
            append(mainMapper.contigs.names, tmpContigs.names);
        }
    }
    catch (BadAlloc const & /* e */)
    {
        throw RuntimeError("Insufficient memory to load the reference.");
    }
    stop(mainMapper.timer);
    mainMapper.stats.loadContigs += getValue(mainMapper.timer);

    if (mainMapper.options.verbose > 1)
        std::cerr << "Loading reference:\t\t\t" << mainMapper.timer << std::endl;
}

// ----------------------------------------------------------------------------
// Function runDisMapper()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig>
inline void runDisMapper(Mapper<TSpec, TConfig> & mainMapper, DisOptions & disOptions)
{

    Timer<double> timer;

    start(timer);
    configureThreads(mainMapper);

    loadAllContigs(mainMapper, disOptions);

    // Open output file and write header.
    openOutputFile(mainMapper);
    openReads(mainMapper);

    // Process reads in blocks.
    // load reads here

    // classify reads here
    // create mappers and run them on subsets
    // Process reads in blocks.
    while (true)
    {
        if (mainMapper.options.verbose > 1) printRuler(std::cerr);
        loadReads(mainMapper);
        if (empty(mainMapper.reads.seqs)) break;
        initReadsContext(mainMapper, mainMapper.reads.seqs);
        setHost(mainMapper.cigarSet, mainMapper.cigars);
        for (uint32_t i=0; i < disOptions.NUM_OF_BINS; ++i)
        {
            Options options = mainMapper.options;
            set_current_index_file(options, disOptions, i);
            if (!openContigsLimits(options))
                throw RuntimeError("Error while opening reference file.");
            configureMapper<TSpec, TConfig>(options, mainMapper, disOptions.contigOffsets[i]);
        }

        aggregateMatches(mainMapper, mainMapper.reads.seqs);
        rankMatches(mainMapper, mainMapper.reads.seqs);
        if (mainMapper.options.verifyMatches)
            verifyMatches(mainMapper);
        alignMatches(mainMapper);
        writeMatches(mainMapper);
        clearMatches(mainMapper);
        clearAlignments(mainMapper);
        clearReads(mainMapper);
    }
    closeReads(mainMapper);
    closeOutputFile(mainMapper);
    stop(timer);
    if (mainMapper.options.verbose > 0)
        printStats(mainMapper, timer);

 }

// ----------------------------------------------------------------------------
// Function spawnDisMapper()
// ----------------------------------------------------------------------------
template <typename TContigsSize, typename TContigsLen, typename TContigsSum,
          typename TThreading, typename TSequencing, typename TSeedsDistance>
inline void spawnDisMapper(DisOptions & disOptions,
                           TThreading const & /* tag */,
                           TSequencing const & /* tag */,
                           TSeedsDistance const & /* tag */)
{
    disOptions.outputFile = disOptions.superOutputFile;
    typedef ReadMapperConfig<TThreading, TSequencing, TSeedsDistance, TContigsSize, TContigsLen, TContigsSum>  TMainConfig;
    Mapper<void, TMainConfig> disMapper(disOptions);

    runDisMapper(disMapper, disOptions);
}

#endif  // #ifndef APP_YARA_MAPPER_H_
