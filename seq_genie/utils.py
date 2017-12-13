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
from pysam import Samfile
import pysam

from synbiochem.utils import io_utils


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


def bin_seqs(barcodes, sequences, ignore_undefined=True, search_len=256):
    '''Bin sequences according to barcodes.'''
    barcode_seqs = defaultdict(list)

    max_barcode_len = max([len(barcode) for barcode in barcodes])

    if barcodes:
        for barcode in barcodes.values():
            for sequence in sequences:
                if barcode in sequence.seq[:max_barcode_len + search_len]:
                    if not ignore_undefined or barcode != 'undefined':
                        barcode_seqs[barcode].append(sequence)

    else:
        barcode_seqs['undefined'].extend(sequences)

    return barcode_seqs


def index(filename):
    '''Index file.'''
    subprocess.call(['bwa', 'index', filename])


def align(templ_filename, reads_filename, out='align.sam', gap_open=6):
    '''Aligns sequences in barcoded bins.'''
    # Align and sort:
    sort(mem(reads_filename, templ_filename, gap_open=gap_open), out)


def mem(reads_filename, templ_filename, readtype='ont2d', gap_open=6):
    '''Runs BWA MEM.'''
    out_file = tempfile.NamedTemporaryFile(delete=False)

    with open(out_file.name, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         '-O', str(gap_open),
                         templ_filename, reads_filename],
                        stdout=out)

    return out_file.name


def get_consensus(sam_filename, templ_filename,
                  for_primer, rev_primer):
    '''Convert files.'''
    bam_filename = sam_filename + '.bam'
    fastq_filename = sam_filename + '.fastq'
    fasta_filename = sam_filename + '.fasta'
    pysam.view(sam_filename, '-o', bam_filename, catch_stdout=False)
    pysam.sort('-o', bam_filename, bam_filename)

    with open(fastq_filename, 'w') as fastq:
        proc1 = subprocess.Popen(['samtools', 'mpileup', '-uf',
                                  templ_filename, bam_filename],
                                 stdout=subprocess.PIPE)

        proc2 = subprocess.Popen(['bcftools', 'call', '-c', '--ploidy', '1'],
                                 stdin=proc1.stdout, stdout=subprocess.PIPE)
        proc1.stdout.close()
        proc3 = subprocess.Popen(['perl', '/Applications/bcftools/vcfutils.pl',
                                  'vcf2fq', '-d', '12'],
                                 stdin=proc2.stdout, stdout=subprocess.PIPE)
        proc2.stdout.close()
        stdout, _ = proc3.communicate()
        fastq.write(stdout)

    _convert_consensus(fastq_filename, fasta_filename, for_primer, rev_primer)

    return fasta_filename


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


def _convert_consensus(fastq_filename, fasta_filename, for_primer, rev_primer):
    '''Convert consensus sequence to fasta.'''
    with open(fastq_filename, 'r') as fastq:
        seq_record = SeqIO.read(fastq, 'fastq')

    seq = _pcr(str(seq_record.seq.upper()), for_primer, rev_primer)

    with open(fasta_filename, 'w') as fasta:
        SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq),
                                        id=seq_record.id),
                    fasta, 'fasta')


def _pcr(seq, forward_primer, reverse_primer):
    '''Apply in silico PCR.'''
    for_primer_pos = seq.find(forward_primer.upper())

    rev_primer_pos = \
        seq.find(str(Seq.Seq(reverse_primer).reverse_complement().upper()))

    print str(for_primer_pos) + ' ' + str(rev_primer_pos)

    if for_primer_pos > -1 and rev_primer_pos > -1:
        seq = seq[for_primer_pos:] + \
            seq[:rev_primer_pos + len(reverse_primer)]
    elif for_primer_pos > -1:
        seq = seq[for_primer_pos:]
    elif rev_primer_pos > -1:
        seq = seq[:rev_primer_pos + len(reverse_primer)]

    return seq


def get_mismatches(sam_filename, templ_seq):
    '''Find mismatches in samfile relative to template sequence.'''
    matches = 0
    mismatches = []

    sam_file = Samfile(sam_filename, 'r')

    for read in sam_file:
        for pair in read.aligned_pairs:
            if pair[1] is not None:
                if read.seq[pair[0]] != templ_seq[pair[1]]:
                    mismatches.append(templ_seq[pair[1]] + str(pair[1] + 1) +
                                      read.seq[pair[0]])
                else:
                    matches += 1

    return matches, mismatches


# templ = get_reads('3955.fasta')[0]
# print get_mismatches('cons.sam', templ.seq)
