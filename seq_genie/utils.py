'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
# pylint: disable=too-many-arguments
from collections import defaultdict
from itertools import izip_longest
import os
from os.path import splitext
import subprocess
import tempfile

from Bio import Seq, SeqIO, SeqRecord
from pysam import Samfile, VariantFile
from synbiochem.utils import io_utils, thread_utils

import numpy as np


def get_reads(reads_filename, min_length=0):
    '''Gets reads.'''
    reads = []
    total_reads = 0

    if os.path.isdir(reads_filename):
        for dirpath, _, filenames in os.walk(os.path.abspath(reads_filename)):
            for filename in filenames:
                filename = os.path.join(dirpath, filename)
                total_reads += _get_reads(filename, min_length, reads)
    else:
        total_reads += _get_reads(reads_filename, min_length, reads)

    return reads, total_reads


def bin_seqs(barcodes, sequences, window_size, search_len=32,
             num_threads=0, batch_size=32):
    '''Bin sequences according to barcodes.'''
    barcode_seqs = defaultdict(list)

    num_seqs = len(sequences)

    max_barcode_len = max([len(barcode)
                           for pair in barcodes
                           for barcode in pair])

    if barcodes:
        barcodes = _format_barcodes(barcodes)

        if num_threads:
            thread_pool = thread_utils.ThreadPool(num_threads)

            for idx, seqs in enumerate(_get_batch(sequences, batch_size)):
                _report_barcodes(idx, num_seqs, batch_size, barcode_seqs)

                thread_pool.add_task(_bin_seqs, seqs, max_barcode_len,
                                     search_len, barcodes, barcode_seqs,
                                     window_size)

            thread_pool.wait_completion()
        else:
            for idx, seqs in enumerate(_get_batch(sequences, batch_size)):
                _report_barcodes(idx, num_seqs, batch_size, barcode_seqs)

                _bin_seqs(seqs, max_barcode_len, search_len, barcodes,
                          barcode_seqs, window_size)
    else:
        barcode_seqs['undefined'].extend(sequences)

    return barcode_seqs


def index(filename):
    '''Index file.'''
    subprocess.call(['bwa', 'index', filename])


def mem(templ_filename, reads_filename, out_filename=None,
        readtype='ont2d', gap_open=6):
    '''Runs BWA MEM.'''
    out_file = io_utils.get_filename(out_filename)

    with open(out_file, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         '-O', str(gap_open),
                         templ_filename, reads_filename],
                        stdout=out)

    return out_file


def get_vcf(bam_filename, templ_filename, pcr_offset=0):
    '''Generates a vcf file.'''
    vcf_filename = \
        tempfile.NamedTemporaryFile('w', suffix='.vcf', delete=False).name \
        if pcr_offset else os.path.join(os.path.dirname(bam_filename),
                                        'variants.vcf')

    prc = subprocess.Popen(['samtools',
                            'mpileup',
                            '-uvf',
                            templ_filename,
                            '-t', 'DP',
                            '-o', vcf_filename,
                            bam_filename])

    prc.communicate()

    if pcr_offset:
        vcf_out_filename = os.path.join(os.path.dirname(bam_filename),
                                        'variants.vcf')
        vcf_in = VariantFile(vcf_filename)
        vcf_out = VariantFile(vcf_out_filename, 'w', header=vcf_in.header)

        for rec in vcf_in.fetch():
            rec.pos = rec.pos + pcr_offset
            print rec
            vcf_out.write(rec)

        vcf_out.close()
        return vcf_out_filename

    return vcf_filename


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


def pcr(seq, forward_primer, reverse_primer):
    '''Apply in silico PCR.'''
    for_primer_pos = seq.find(forward_primer.upper())

    rev_primer_pos = \
        seq.find(str(Seq.Seq(reverse_primer).reverse_complement().upper()))

    if for_primer_pos > -1 and rev_primer_pos > -1:
        seq = seq[for_primer_pos:] + \
            seq[:rev_primer_pos + len(reverse_primer)]
    elif for_primer_pos > -1:
        seq = seq[for_primer_pos:]
    elif rev_primer_pos > -1:
        seq = seq[:rev_primer_pos + len(reverse_primer)]

    return seq, for_primer_pos


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


def get_dir(parent_dir, barcodes, ice_id=None):
    '''Get directory from barcodes.'''
    dir_name = os.path.join(parent_dir, '_'.join(barcodes))

    if ice_id:
        dir_name = os.path.join(dir_name, ice_id)

    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

    return dir_name


def _get_reads(filename, min_length, reads):
    '''Gets reads.'''
    _, ext = splitext(filename)

    try:
        with open(filename, 'rU') as fle:
            all_reads = [record for record in SeqIO.parse(fle, ext[1:])]

            reads.extend([record for record in all_reads
                          if len(record.seq) > min_length])

            print 'Reading: %s' % filename
            return len(all_reads)
    except (IOError, ValueError), err:
        print err
        return 0


def _report_barcodes(idx, num_seqs, batch_size, barcode_seqs):
    '''Report barcodes.'''
    seq_lens = [len(seq) for seq in barcode_seqs.values()]

    if seq_lens:
        vals = ((idx + 1) * batch_size,
                num_seqs,
                sum(seq_lens),
                min(seq_lens),
                max(seq_lens),
                np.mean(seq_lens),
                np.median(seq_lens))

        s = 'Seqs: %d/%d\tMatched: %d\tRange: %d-%d\tMean: %.1f\tMedian: %d'
        print s % vals
    else:
        print 'Seqs: %d/%d' % ((idx + 1) * batch_size,  num_seqs)


def _format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    for_brcds = {}
    rev_brcds = {}

    for pair in barcodes:
        for_brcds[pair] = \
            [pair[0], str(Seq.Seq(pair[1]).reverse_complement())]
        rev_brcds[pair] = \
            [pair[1], str(Seq.Seq(pair[0]).reverse_complement())]

    return for_brcds, rev_brcds


def _bin_seqs(seqs, max_barcode_len, search_len, barcodes, barcode_seqs,
              window_size):
    '''Bin a batch of sequences.'''
    for seq in seqs:
        for pairs in barcodes:
            if _check_seq(seq, max_barcode_len, search_len, pairs, barcode_seqs,
                          window_size):
                break


def _check_seq(seq, max_barcode_len, search_len, pairs, barcode_seqs,
               window_size):
    seq_start = seq.seq[:max_barcode_len + search_len]
    seq_end = seq.seq[-(max_barcode_len + search_len):]
    selected_barcodes = [None, None]

    # Check all barcodes:
    for orig, bc_pair in pairs.iteritems():
        _check_pair(orig, bc_pair, seq_start, seq_end,
                    selected_barcodes, window_size)

        if all(selected_barcodes):
            barcode_seqs[tuple(selected_barcodes)].append(seq)
            return True

    return False


def _check_pair(orig, pair, seq_start, seq_end,
                selected_barcodes, window_size):
    '''Check sliding window scores.'''
    if not selected_barcodes[0]:
        for substr in [pair[0][i:i + window_size]
                       for i in xrange(len(pair[0]) - window_size + 1)]:
            if substr in seq_start:
                selected_barcodes[0] = orig[0]
                break

    if not selected_barcodes[1]:
        for substr in [pair[1][i:i + window_size]
                       for i in xrange(len(pair[1]) - window_size + 1)]:
            if substr in seq_end:
                selected_barcodes[1] = orig[1]
                break


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


def _get_batch(iterable, batch_size=64, fillvalue=None):
    '''Get a batch.'''
    args = [iter(iterable)] * batch_size
    return izip_longest(*args, fillvalue=fillvalue)
