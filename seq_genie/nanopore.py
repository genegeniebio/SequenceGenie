'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
import os
import subprocess
import sys
import tempfile
from Bio import Seq, SeqIO, SeqRecord
from pysam import AlignmentFile

from synbiochem.utils import ice_utils, seq_utils
import matplotlib.pyplot as plt


def convert(directories, min_length=0, fasta_filename='seq.fasta',
            fasta_dir='pass'):
    '''Converts fast5 files to fasta.'''
    with open(fasta_filename, 'w') as fasta_file:
        for directory in directories:
            for dirpath, _, filenames in os.walk(os.path.abspath(directory)):
                for filename in filenames:
                    filename = os.path.join(dirpath, filename)

                    if filename[-6:] == '.fast5' and fasta_dir in filename:
                        subprocess.call('poretools fasta ' + filename,
                                        stdout=fasta_file,
                                        shell=True)
                    elif filename[-6:] == '.fastq' and fasta_dir in filename:
                        with open(filename, 'rU') as fle:
                            records = [record
                                       for record in SeqIO.parse(fle, 'fastq')
                                       if len(record.seq) > min_length]

                            SeqIO.write(records, fasta_file, 'fasta')

    return fasta_filename


def align(templ_filename, seqs_filename, out='align.sam', fmt='fasta'):
    '''Aligns sequences in barcoded bins.'''
    # Index template:
    subprocess.call(['bwa', 'index', templ_filename])

    # Read sequences:
    with open(seqs_filename, 'rU') as fle:
        seqs = {record.id: str(record.seq)
                for record in SeqIO.parse(fle, fmt)}

    # Align and sort:
    _sort(_mem(seqs, templ_filename), out)


def _mem(seqs, templ_filename, readtype='pacbio'):
    '''Runs BWA MEM.'''
    out_file = tempfile.NamedTemporaryFile(delete=False)
    seq_file = tempfile.NamedTemporaryFile(delete=False)

    records = [SeqRecord.SeqRecord(Seq.Seq(seq), seq_id, '', '')
               for seq_id, seq in seqs.iteritems()]

    SeqIO.write(records, seq_file.name, 'fasta')

    with open(out_file.name, 'w') as out:
        subprocess.call(['bwa', 'mem',
                         '-x', readtype,
                         templ_filename, seq_file.name], stdout=out)

    return out_file.name


def _sort(in_filename, out_filename):
    '''Custom sorts SAM file.'''
    sam_file = AlignmentFile(in_filename, 'r')
    out_file = AlignmentFile(out_filename, 'wh',
                             template=sam_file,
                             header=sam_file.header)

    for read in sorted([read for read in sam_file],
                       key=lambda x: (-x.query_length,
                                      x.reference_start)):
        out_file.write(read)

    out_file.close()

    return out_filename


def blast(ice_client, seq):
    '''Does BLAST search against ICE.'''
    results = ice_client.do_blast(seq)

    if results['results']:
        res = results['results'][0]
        print '\t'.join(str(val)
                        for val in [res['entryInfo']['partId'], res['eValue']])


def analyse(seqs_filename):
    '''Analyse reads.'''
    rds = [[seq_id, len(seq), seq]
           for seq_id, seq in seq_utils.read_fasta(seqs_filename).iteritems()]

    # for read in reads:
    #    print '\t'.join(str(val) for val in rds)

    plt.hist([read[1] for read in rds], edgecolor='black')
    plt.title('Read lengths')
    plt.xlabel('Read length')
    plt.ylabel('Frequency')
    plt.show()


def main(args):
    '''main method.'''
    seqs_filename = convert([args[0]], int(args[1]))

    ice_client = ice_utils.ICEClient(args[2], args[3], args[4])

    for ice_id in args[5:]:
        id_seqs = {ice_id: ice_client.get_ice_entry(ice_id).get_seq()}
        templ_filename = ice_id + '.fasta'
        seq_utils.write_fasta(id_seqs, templ_filename)
        align(templ_filename, seqs_filename, out=ice_id + '.sam')


if __name__ == '__main__':
    main(sys.argv[1:])
