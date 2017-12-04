'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-member
# pylint: disable=no-name-in-module
# pylint: disable=ungrouped-imports
from collections import defaultdict
import os
from os.path import splitext
import subprocess
import tempfile

from Bio import Seq, SeqIO, SeqRecord
from pyfaidx import FastaVariant
from pysam import Samfile
import pysam

from synbiochem.utils import io_utils, seq_utils


def parse(reads_filename):
    '''Parse reads file.'''
    return get_reads(reads_filename)


def get_reads(reads_filename, min_length=0):
    '''Gets reads.'''
    reads = []

    if os.path.isdir(reads_filename):
        for dirpath, _, filenames in os.walk(os.path.abspath(reads_filename)):
            for filename in filenames:
                filename = os.path.join(dirpath, filename)
                _get_reads(filename, min_length, reads)
    else:
        _get_reads(reads_filename, min_length, reads)

    return reads


def bin_seqs(barcodes, sequences, evalue=1, ignore_undefined=True):
    '''Bin sequences according to barcodes.'''
    barcode_seqs = defaultdict(dict)

    if barcodes:
        results = seq_utils.do_blast(barcodes, sequences, evalue=evalue,
                                     word_size=4)
        for result in results:
            barcode = 'undefined'

            for alignment in result.alignments:
                barcode = barcodes[alignment.hit_def]

            if barcode != 'undefined' or not ignore_undefined:
                barcode_seqs[barcode][result.query] = sequences[result.query]

    else:
        barcode_seqs['undefined'] = {seq_id: seq
                                     for seq_id, seq in sequences.iteritems()}

    return barcode_seqs


def align(templ_filename, reads, out='align.sam', gap_open=6):
    '''Aligns sequences in barcoded bins.'''
    # Index template:
    subprocess.call(['bwa', 'index', templ_filename])

    # Align and sort:
    sort(mem(reads, templ_filename, gap_open=gap_open), out)


def mem(reads, templ_filename, readtype='pacbio', gap_open=6):
    '''Runs BWA MEM.'''
    reads_file = tempfile.NamedTemporaryFile(delete=False)
    out_file = tempfile.NamedTemporaryFile(delete=False)

    SeqIO.write(reads, reads_file.name, 'fasta')

    with open(out_file.name, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         '-O', str(gap_open),
                         templ_filename, reads_file.name],
                        stdout=out)

    return out_file.name


def get_consensus(sam_filename, templ_filename):
    '''Convert files.'''
    # samtools mpileup -uf ref.fa aln.bam | bcftools view -cg - | vcfutils.pl
    # vcf2fq > cns.fq
    bam_filename = sam_filename + '.bam'
    vcf_filename = sam_filename + '.vcf'
    vcf_call_filename = sam_filename + '_call.vcf'
    fasta_filename = sam_filename + '.fas'
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)

    with open(vcf_filename, 'w') as vcf_file:
        vcf_file.write(pysam.mpileup('-uf', templ_filename, bam_filename))

    # samtools mpileup -uf 3958.fasta CATCCTAGTTGGTACTGCAATACT_3958.sam.bam |
    # bcftools call -c | perl /Applications/bcftools/vcfutils.pl vcf2fq -d 2 >
    # cns.fq

    subprocess.call(['bcftools', 'consensus',
                     '-f', templ_filename,
                     '-o', fasta_filename,
                     vcf_filename])

    consensus = FastaVariant(templ_filename, vcf_filename, sample='consensus',
                             het=True, hom=True)
    print consensus


def sort(in_filename, out_filename):
    '''Custom sorts SAM file.'''
    sam_file = Samfile(in_filename, 'r')
    out_file = Samfile(out_filename, 'wh',
                       template=sam_file,
                       header=sam_file.header)

    for read in sorted([read for read in sam_file],
                       key=lambda x: (-x.query_length,
                                      x.reference_start)):
        out_file.write(read)

    out_file.close()

    return out_filename


def reject_indels(sam_filename, templ_seq, out_filename=None):
    '''Rejects indels.'''
    out_filename = io_utils.get_filename(out_filename)

    sam_file = Samfile(sam_filename, 'r')
    out_file = Samfile(out_filename, 'wh',
                       template=sam_file,
                       header=sam_file.header)

    for read in sam_file:
        if read.cigarstring and str(len(templ_seq)) + 'M' in read.cigarstring:
            out_file.write(read)

    out_file.close()

    return out_filename


def replace_indels(sam_filename, templ_seq, out_filename=None):
    '''Replace indels, replacing them with wildtype.'''
    out_filename = io_utils.get_filename(out_filename)

    with open(out_filename, 'w') as fle:
        SeqIO.write(_replace_indels(sam_filename, templ_seq), fle, 'fasta')

    return out_filename


def _get_reads(filename, min_length, reads):
    '''Gets reads.'''
    _, ext = splitext(filename)

    try:
        with open(filename, 'rU') as fle:
            reads.extend([record
                          for record in SeqIO.parse(fle, ext[1:])
                          if len(record.seq) > min_length])
    except (IOError, ValueError), err:
        print err


def _replace_indels(sam_filename, templ_seq):
    '''Replace indels, replacing them with wildtype.'''
    sam_file = Samfile(sam_filename, 'r')

    for read in sam_file:
        # Perform mapping of nucl indices to remove spurious indels:
        seq = ''.join([read.seq[pair[0]]
                       if pair[0]
                       else templ_seq[pair[1]]
                       for pair in read.aligned_pairs
                       if pair[1] is not None])

        if seq:
            yield SeqRecord.SeqRecord(Seq.Seq(seq), read.qname, '', '')
