#include <iostream>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <cstdlib>
#include <stdint.h>
#include <seqan/arg_parse.h>

#include "bamhash_checksum_common.h"

struct Options {
  char version[1024];
  std::vector<seqan::CharString> BAMfiles;
  std::vector<seqan::CharString> Fastq1;
  std::vector<seqan::CharString> Fastq2;
};

int HashBAM(seqan::CharString fname, uint64_t *sum, uint64_t *count) {

  // Open BGZF Stream for reading.
  seqan::Stream<seqan::Bgzf> inStream;
  if (!open(inStream, toCString(fname), "r")) {
    std::cerr << "ERROR: Could not open " << fname << " for reading.\n";
    return 1;
  }

  // Setup name store, cache, and BAM I/O context.
  typedef seqan::StringSet<seqan::CharString> TNameStore;
  typedef seqan::NameStoreCache<TNameStore>   TNameStoreCache;
  typedef seqan::BamIOContext<TNameStore>     TBamIOContext;
  TNameStore      nameStore;
  TNameStoreCache nameStoreCache(nameStore);
  TBamIOContext   context(nameStore, nameStoreCache);

  // Read header.
  seqan::BamHeader header;
  if (readRecord(header, context, inStream, seqan::Bam()) != 0) {
    std::cerr << "ERROR: Could not read header from BAM file " << fname << "\n";
    return 1;
  }
  seqan::clear(header);

  // Define:
  seqan::BamAlignmentRecord record;
  seqan::CharString string2hash;

  // Read record
  while (!atEnd(inStream)) {
    if (readRecord(record, context, inStream, seqan::Bam()) != 0) {
      std::cerr << "ERROR: Could not read record from BAM File " << fname << "\n";
      return 1;
    }
    // Check if flag: reverse complement and change record accordingly
    if (hasFlagRC(record)) {
      reverseComplement(record.seq);
      reverse(record.qual);
    }
    // Check if flag: supplementary and exclude those
    if (!hasFlagSupplementary(record) && !hasFlagSecondary(record)) {
      *count +=1;
      // Construct one string from record
      seqan::append(string2hash, record.qName);
      if(hasFlagLast(record)) {
        seqan::append(string2hash, "/2");
      } else {
        seqan::append(string2hash, "/1");
      }

      seqan::append(string2hash, record.seq);
      seqan::append(string2hash, record.qual);
      seqan::clear(record);

      // Get MD5 hash
      hash_t hex = str2md5(toCString(string2hash), length(string2hash));
      seqan::clear(string2hash);

      hexSum(hex, *sum);
    }
  }
  return 0;
}

int HashFastq(seqan::CharString fname, const char *suffix, uint64_t *sum, uint64_t *cnt) {
  // Define:
  seqan::StringSet<seqan::CharString> idSub;
  seqan::CharString string2hash;
  seqan::CharString id1;
  seqan::CharString seq1;
  seqan::CharString qual1;
  hash_t hex;

  // Open GZStream
  seqan::Stream<seqan::GZFile> gzStream1;

  if (!open(gzStream1, toCString(fname), "r")) {
    std::cerr << "ERROR: Could not open the file: " << fname << " for reading.\n";
    return 1;
  }

  //Setup RecordReader for reading FASTQ file from gzip-compressed file
  seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<> > reader1(gzStream1);

  // Read record
  while (!atEnd(reader1)) {
    if (readRecord(id1, seq1, qual1, reader1, seqan::Fastq()) != 0) {
      if (atEnd(reader1)) {
        std::cerr << "WARNING: Could not continue reading " << fname <<  " at line: " << (*cnt)+1 << "! Is the file corrupt?\n";
        return 1;
      }
      std::cerr << "ERROR: Could not read from " << fname << "\n";
      return 1;
    }

    *cnt +=1;

    // If include id, then cut id on first whitespace
    if (seqan::endsWith(id1,suffix)) {
      seqan::strSplit(idSub, id1, '/', false, 1);
    } else {
      seqan::strSplit(idSub, id1, ' ', false, 1);
    }

    seqan::append(string2hash, idSub[0]);
    seqan::append(string2hash, suffix);
    seqan::append(string2hash, seq1);
    seqan::append(string2hash, qual1);

    // Get MD5 hash
    hex = str2md5(toCString(string2hash), length(string2hash));

    hexSum(hex, *sum);

    seqan::clear(string2hash);
    seqan::clear(idSub);
  }

  return 0;
}

seqan::ArgumentParser::ParseResult
parseCommandLine(Options& options, int argc, char const **argv) {
  unsigned int i;
  seqan::CharString s;
  const char *description = "This program will compute the checksum of one or "
                            "more SAM/BAM files as well as the fastq files used to create it/them. Both "
                            "single and paired-end datasets are supported and can be mixed together. The "
                            "checksum and number of reads/alignments for each file are written to the "
                            "console. In addition, the checksum and number of entries of the BAM files are "
                            "compared to the fastq files and an error message produced if there's a "
                            "mismatch. The return value is 1 if there's an error or mismatch or 0 otherwise.";

  // Setup ArgumentParser.
  seqan::ArgumentParser parser("bamhash_checksum_all");

  setShortDescription(parser, "Compute the checksum of a series of BAM and fastq files.");
  setVersion(parser, BAMHASH_VERSION);
  setDate(parser, "Mar 2015");

  addUsageLine(parser, "[\\fIOPTIONS\\fP] \\fI<-b in1.bam>\\fP \\fI<-1 sample1_1.fastq>\\fP \\fI[-2 sample1_2.fastq]\\fP");
  addDescription(parser, description);

  addOption(parser, seqan::ArgParseOption("b", "bam", "Input SAM/BAM file. Specify this more than once to process multiple files.", seqan::ArgParseArgument::INPUTFILE, "BAM file", 1));
  setRequired(parser, "bam");
  addOption(parser, seqan::ArgParseOption("1", "fastq1", "As above, but for fastq files.", seqan::ArgParseArgument::INPUTFILE, "fastq file", 1));
  setRequired(parser, "fastq1");
  addOption(parser, seqan::ArgParseOption("2", "fastq2", "As above, but the file(s) containing the second mate in each pair.", seqan::ArgParseArgument::INPUTFILE, "fastq file", 1));

  setValidValues(parser, "bam", "bam sam");
  setValidValues(parser, "fastq1", "fq fastq fq.gz fastq.gz");
  setValidValues(parser, "fastq2", "fq fastq fq.gz fastq.gz");

  addSection(parser, "Options");

  // Parse command line.
  seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res;
  }

  for(i=0; i<getOptionValueCount(parser, "bam"); i++) {
    getOptionValue(s, parser, "bam", i);
    options.BAMfiles.push_back(s);
  }
  for(i=0; i<getOptionValueCount(parser, "fastq1"); i++) {
    getOptionValue(s, parser, "fastq1", i);
    options.Fastq1.push_back(s);
  }
  for(i=0; i<getOptionValueCount(parser, "fastq2"); i++) {
    getOptionValue(s, parser, "fastq2", i);
    options.Fastq2.push_back(s);
  }

  return seqan::ArgumentParser::PARSE_OK;
}


int main(int argc, char const **argv) {

  Options info; // Define structure variable
  seqan::ArgumentParser::ParseResult res = parseCommandLine(info, argc, argv); // Parse the command line.
  unsigned int i, err = 0;
  uint64_t globalSum = 0, fileSum = 0;
  uint64_t globalCnt = 0, fileCnt = 0;

  if (res != seqan::ArgumentParser::PARSE_OK) {
    return res == seqan::ArgumentParser::PARSE_ERROR;
  }

  //BAM files
  for(i=0; i<info.BAMfiles.size(); i++) {
    if(HashBAM(info.BAMfiles[i], &fileSum, &fileCnt)) {
      std::cerr << "Error encountered while processing " << info.BAMfiles[i] << "!\n";
      err = 1;
    }
    std::cout << std::hex << fileSum << '\t' << std::dec << fileCnt << '\t' << info.BAMfiles[i] << '\n';
    globalCnt += fileCnt;
    globalSum += fileSum;
    fileCnt = fileSum = 0;
  }

  //fastq1
  for(i=0; i<info.Fastq1.size(); i++) {
    if(HashFastq(info.Fastq1[i], "/1", &fileSum, &fileCnt)) {
      std::cerr << "Error encountered while processing " << info.Fastq1[i] << "!\n";
      err = 1;
    }
    std::cout << std::hex << fileSum << '\t' << std::dec << fileCnt << '\t' << info.Fastq1[i] << '\n';
    if(fileCnt > globalCnt) {
      std::cerr << "Error: more reads in the fastq files than alignments in the BAM file(s).\n";
      err = 1;
    }
    globalSum -= fileSum;
    globalCnt -= fileCnt;
    fileCnt = fileSum = 0;
  }

  //fastq2
  for(i=0; i<info.Fastq2.size(); i++) {
    if(HashFastq(info.Fastq2[i], "/2", &fileSum, &fileCnt)) {
      std::cerr << "Error encountered while processing " << info.Fastq2[i] << "!\n";
      err = 1;
    }
    std::cout << std::hex << fileSum << '\t' << std::dec << fileCnt << '\t' << info.Fastq2[i] << '\n';
    if(fileCnt > globalCnt) {
      std::cerr << "Error: more reads in the fastq files than alignments in the BAM file(s).\n";
      err = 1;
    }
    globalSum -= fileSum;
    globalCnt -= fileCnt;
    fileCnt = fileSum = 0;
  }

  if(err) {
    std::cerr << "Mismatch in the number of reads between the BAM and fastq files!\n";
    return 1;
  } else if(globalSum) {
    std::cerr << "The BAM and fastq files do not match!\n";
    return 1;
  } else {
    std::cerr << "The BAM and fastq files match.\n";
    return 0;
  }
}
