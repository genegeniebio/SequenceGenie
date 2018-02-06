'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=no-name-in-module
import os
import sys

from Bio import SeqIO, SeqRecord
import pysam
from synbiochem.utils.mut_utils import get_mutations
from synbiochem.utils.seq_utils import INV_NUCL_CODES

from seq_genie import protein, utils


def _analyse(sam_filename, templ_filename, mut_templ_filename):
    '''Analyse.'''
    out_dir_name = os.path.dirname(sam_filename)
    nucl_filename = os.path.join(out_dir_name, 'nucl.fasta')
    aa_filename = os.path.join(out_dir_name, 'aa.fasta')
    res_filename = os.path.join(out_dir_name, 'out.tsv')

    # Get sequences:
    templ_seq = list(utils.get_reads(templ_filename))[0].seq
    mut_templ_seq = list(utils.get_reads(mut_templ_filename))[0].seq

    # Strip N and spurious mutations:
    nucl_seqs = _strip(sam_filename, templ_seq, mut_templ_seq)

    with open(nucl_filename, 'w') as nucl_file:
        SeqIO.write(nucl_seqs, nucl_file, 'fasta')

    # Translate nucl seqs:
    _trans_nucl(aa_filename, res_filename, templ_seq, nucl_seqs)

    protein.align(templ_filename, nucl_filename)


def _strip(sam_filename, templ_seq, mut_templ_seq):
    '''Replace N and spurious mutations with wildtype.'''
    nucl_seqs = []

    for read in pysam.AlignmentFile(sam_filename, 'r'):
        # Perform mapping of nucl indices to remove N:
        if read.aligned_pairs:
            templ_idx = zip(*read.aligned_pairs)[1]
            prefix = templ_seq[:min([val for val in templ_idx if val])]
            suffix = templ_seq[max(templ_idx) + 1:]
            seq = ''.join([read.seq[pair[0]]
                           if _valid(pair, read, templ_seq, mut_templ_seq)
                           else templ_seq[pair[1]]
                           for pair in read.aligned_pairs
                           if pair[1] is not None])

            if seq:
                nucl_seqs.append(SeqRecord.SeqRecord(prefix + seq + suffix,
                                                     read.qname, '', ''))

    return nucl_seqs


def _valid(pair, read, templ_seq, mut_templ_seq):
    '''Determines whether nucleotide pair is valid.'''
    return pair[0] and \
        read.seq[pair[0]] != 'N' and \
        (read.seq[pair[0]] in INV_NUCL_CODES[mut_templ_seq[pair[1]]] or
         read.seq[pair[0]] == templ_seq[pair[1]])


def _trans_nucl(aa_filename, res_filename, templ_seq, nucl_seqs):
    '''Translate and analyse nucleotide sequences.'''
    with open(aa_filename, 'w') as aa_file, \
            open(res_filename, 'w') as res_file:
        aa_seqs = []
        aa_templ_seq = templ_seq.translate()

        res_file.write('\t'.join(['id', 'mutations', 'seq']) + '\n')

        for nucl_seq in nucl_seqs:
            aa_seq = SeqRecord.SeqRecord(nucl_seq.seq.translate(),
                                         id=nucl_seq.id,
                                         name=nucl_seq.name,
                                         description=nucl_seq.description)
            aa_seqs.append(aa_seq)

            res = [aa_seq.id,
                   get_mutations(str(aa_templ_seq), str(aa_seq.seq)),
                   aa_seq.seq]

            res_file.write('\t'.join(str(val) for val in res) + '\n')

        SeqIO.write(aa_seqs, aa_file, 'fasta')


def main(args):
    '''main method.'''
    # Align:
    sam_filename = protein.align(args[0], args[2])
    _analyse(sam_filename, args[0], args[1])


if __name__ == '__main__':
    main(sys.argv[1:])
