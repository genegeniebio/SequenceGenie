'''
SequenceGenie (c) GeneGenie Bioinformatics Ltd. 2017

All rights reserved.

@author: neilswainston
'''
import sys
import seq_genie


def main(args):
    '''main method.'''
    if args[0] == 'pathway':
        seq_genie.pathway.main(args[1:])


if __name__ == '__main__':
    main(sys.argv[1:])
