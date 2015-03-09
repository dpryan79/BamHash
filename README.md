BamHash
=======

Hash BAM and FASTQ files to verify data integrity

For each pair of reads in a BAM or FASTQ file we compute a hash value
composed of the readname, whether it is first or last in pair, sequence and quality value.
All the hash values are summed up so the result is independent of the ordering within the files.
The result can be compared to verify that the pair of FASTQ files contain the same read 
information as the aligned BAM file.

Manuscript
==========

In preperation.

Compiling
=========

The only external dependency is on OpenSSL for the MD5 implementation.

Usage
=====

bamhash_checksum_fastq
----------------------

bamhash_checksum_bam
--------------------

bamhash_checksum_all
--------------------

`bamhash_checksum_all` allows the calculation of the checksum of multiple BAM and fastq files, followed by a comparison of these summed values. This facilitates comparing BAM files and their associated fastq files in a single command. Both paired-end and single-end datasets are supported. Multiple BAM and fastq files can be specified. The simplest example usage would be as follows:

    bamhash_checksum_all -b Sample1.bam -1 Sample1_1.fastq -2 Sample1_2.fastq

This command would be appropriate given a paired-end dataset, where `Sample1_1.fastq` contains read #1 from each pair and `Sample1_2.fastq` contains read #2 from each pair. The output could be something similar to:

    4b28c08b9aef33f5	35434414	Sample1.bam
    962044272c27e290	17717207	Sample1_1.fastq
    B5087C646EC75165	17717207	Sample1_2.fastq

"The BAM and fastq files match." would then be printed to standard error. The exit value from `bamhash_chechsum_all` would be 0 in this example. If there's an error or a mismatch in the checksums or number of reads, the return value would be 1. For a single-end dataset, the following command would work:

    bamhash_checksum_all -b Sample2.bam -1 Sample2.fastq

Multiple BAM/fastq files can also be specified, for example:

    bamhash_checksum_all -b Sample1.bam -b Sample2.bam -1 Sample1_1.fastq -1 Sample2.fastq -2 Sample1_2.fastq

Again, `bamhash_checksum_all` will exit with value of 0 if the summed checksum of the fastq files and the number of reads in them matches the BAM files. You can easily capture this value in bash with the `$?` variable.
