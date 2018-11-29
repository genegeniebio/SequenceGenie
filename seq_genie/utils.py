'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
# pylint: disable=superfluous-parens
import os
from os.path import splitext
import subprocess
import tempfile

from Bio import Seq, SeqIO, SeqRecord
from pysam import Samfile, VariantFile
from synbiochem.utils import io_utils


def get_reads(reads_filename, min_length=0, max_reads=float('inf')):
    '''Gets reads.'''
    reads = []
    total_reads = 0

    if os.path.isdir(reads_filename):
        for dirpath, _, filenames in os.walk(os.path.abspath(reads_filename)):
            for filename in filenames:
                filename = os.path.join(dirpath, filename)
                total_reads += _get_reads(filename, min_length, reads)

                if total_reads > max_reads:
                    break
    else:
        total_reads += _get_reads(reads_filename, min_length, reads)

    return reads, total_reads


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
            print(rec)
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

            print('Reading: %s' % filename)
            return len(all_reads)
    except (IOError, ValueError) as err:
        print(err)
        return 0


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
