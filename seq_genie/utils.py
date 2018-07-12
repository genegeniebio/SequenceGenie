'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
# pylint: disable=too-many-arguments
from collections import defaultdict
import os
from os.path import splitext
import subprocess
import tempfile

from Bio import Seq, SeqIO, SeqRecord
from fuzzywuzzy.fuzz import partial_ratio
from pysam import Samfile, VariantFile
from synbiochem.utils import io_utils, thread_utils


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


def bin_seqs(barcodes, sequences, search_len=32,
             num_threads=0, score_threshold=100):
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

            for idx, seq in enumerate(sequences):
                if idx % 100 == 0:
                    print 'Binning seq: %d/%d' % (idx, num_seqs)

                thread_pool.add_task(_bin_seq, seq, max_barcode_len,
                                     search_len, score_threshold, barcodes,
                                     barcode_seqs)

            thread_pool.wait_completion()
        else:
            for idx, seq in enumerate(sequences):
                if idx % 100 == 0:
                    print 'Binning seq: %d' % idx

                _bin_seq(seq, max_barcode_len, search_len, score_threshold,
                         barcodes, barcode_seqs)
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
        if pcr_offset else bam_filename + '.vcf'

    prc = subprocess.Popen(['samtools',
                            'mpileup',
                            '-uvf',
                            templ_filename,
                            '-t', 'DP',
                            '-o', vcf_filename,
                            bam_filename])

    prc.communicate()

    if pcr_offset:
        vcf_in = VariantFile(vcf_filename)
        vcf_out = VariantFile(bam_filename + '.vcf', 'w', header=vcf_in.header)

        for rec in vcf_in.fetch():
            rec.pos = rec.pos + pcr_offset
            vcf_out.write(rec)

        vcf_out.close()

    return bam_filename + '.vcf'


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


def _format_barcodes(barcodes):
    '''Format barcodes to reduce number of get_rev_complement calls.'''
    form_brcds = {}

    for pair in barcodes:
        forward = [pair[0], str(Seq.Seq(pair[1]).reverse_complement())]
        reverse = [pair[1], str(Seq.Seq(pair[0]).reverse_complement())]
        form_brcds[pair] = [forward, reverse]

    return form_brcds


def _bin_seq(seq, max_barcode_len, search_len, score_thresh, barcodes,
             barcode_seqs):
    '''Bin an individual sequence.'''
    seq_start = seq.seq[:max_barcode_len + search_len]
    seq_end = seq.seq[-(max_barcode_len + search_len):]

    selected_barcodes = \
        _bin_seq_strict(seq_start, seq_end, barcodes) \
        if score_thresh == 100 \
        else _bin_seq_fuzzy(seq_start, seq_end, barcodes, score_thresh)

    if selected_barcodes:
        barcode_seqs[selected_barcodes].append(seq)


def _bin_seq_strict(seq_start, seq_end, barcodes):
    '''Bin seq strict.'''
    for orig, formatted in barcodes.iteritems():
        if formatted[0][0] in seq_start and formatted[0][1] in seq_end:
            return orig

        if formatted[1][0] in seq_start and formatted[1][1] in seq_end:
            return orig

    return None


def _bin_seq_fuzzy(seq_start, seq_end, barcodes, score_threshold,
                   exhaustive=False):
    '''Bin seq fuzzy.'''
    best_barcodes = _bin_seq_strict(seq_start, seq_end, barcodes)

    if not best_barcodes:
        max_scores = score_threshold, score_threshold

        for orig, formatted in barcodes.iteritems():
            for direction in formatted:
                returned_barcodes, max_scores = \
                    _check_scores(orig, direction, seq_start, seq_end,
                                  max_scores)

                if returned_barcodes:
                    best_barcodes = returned_barcodes

                    if not exhaustive:
                        break

    return best_barcodes


def _check_scores(orig, formatted, seq_start, seq_end, max_scores):
    '''Check barcode scores.'''
    barcodes = None

    scores = partial_ratio(formatted[0], seq_start), \
        partial_ratio(formatted[1], seq_end)

    if scores[0] > max_scores[0] and scores[1] > max_scores[1]:
        barcodes = orig
        max_scores = scores

    return barcodes, max_scores


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
