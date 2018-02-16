'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
import sys

from seq_genie import pathway


def main(args):
    '''main method.'''
    if args[0] == 'pathway':
        pathway.score_alignments(args[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
