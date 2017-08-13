'''
DNA++ (c) DNA++ 2017

All rights reserved.

@author: neilswainston
'''
import os
import subprocess
import sys
from synbiochem.utils import ice_utils, seq_utils
import matplotlib.pyplot as plt


def convert(directories, fasta_filename='seq.fasta', fasta_dir='pass'):
    '''Converts fast5 files to fasta.'''
    fasta_files = []

    for directory in directories:
        for dirpath, _, filenames in os.walk(os.path.abspath(directory)):
            if os.path.split(dirpath)[1] == fasta_dir:
                fasta_file = open(os.path.join(dirpath, fasta_filename), 'w')
                fasta_files.append(fasta_file)

            for filename in filenames:
                filename = os.path.join(dirpath, filename)

                if filename[-6:] == '.fast5' and fasta_dir in filename:
                    subprocess.call('poretools fasta ' + filename,
                                    stdout=fasta_file,
                                    shell=True)

    for fasta_file in fasta_files:
        fasta_file.close()

    return [fasta_file.name for fasta_file in fasta_files]


def blast(ice_client, seq):
    '''Does BLAST search against ICE.'''
    results = ice_client.do_blast(seq)

    if results['results']:
        res = results['results'][0]
        print '\t'.join(str(val)
                        for val in [res['entryInfo']['partId'], res['eValue']])


def main(args):
    '''main method.'''
    reads = []

    for filename in convert(args[3:]):
        for seq_id, seq in seq_utils.read_fasta(filename).iteritems():
            reads.append([seq_id, len(seq), seq])

    for read in reads:
        print '\t'.join(str(val) for val in read)

    plt.hist([read[1] for read in reads], edgecolor='black')
    plt.title('Read lengths')
    plt.xlabel('Read length')
    plt.ylabel('Frequency')
    plt.show()

    ice_client = ice_utils.ICEClient(args[0], args[1], args[2])

    for read in reads:
        blast(ice_client, read[2])


if __name__ == '__main__':
    main(sys.argv[1:])
