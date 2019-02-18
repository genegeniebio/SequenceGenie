'''
SequenceGenie (c) University of Manchester 2018

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import os
from threading import Thread

import numpy as np
import pandas as pd


class ResultsThread(Thread):
    '''Thread-safe class to write results to DataFrames.'''

    def __init__(self, columns, barcodes_df, queue):
        self.__dfs = {}

        self.__dfs['summary'] = barcodes_df
        self.__dfs['summary'].set_index(['forward', 'reverse'], inplace=True)

        # Initialise specific Dataframes:
        for name in ['mutations', 'indels', 'deletions', 'identity', 'depths']:
            self.__dfs[name] = self.__init_df(columns)

        self.__queue = queue
        self.__closed = False
        Thread.__init__(self)

    def run(self):
        '''Run.'''
        while True:
            task = self.__queue.get()

            if task is None:
                break

            self.__update_df(task)
            self.__queue.task_done()

    def write(self, dir_name):
        '''Write result output files.'''
        self.__update_summary()

        for name, df in self.__dfs.items():
            print(dir_name, name)
            df.to_csv(os.path.join(dir_name, name + '.csv'))

    def close(self):
        '''Close.'''
        self.__queue.put(None)

    def __init_df(self, columns):
        '''Initialise a results dataframe.'''
        df = pd.concat([pd.DataFrame(columns=columns),
                        self.__dfs['summary'].copy()],
                       sort=False)
        df.index = self.__dfs['summary'].index
        return df

    def __update_df(self, task):
        '''Update DataFrame.'''
        name, val, col_id, row_ids = task
        self.__dfs[name][col_id].loc[row_ids[0], row_ids[1]] = val

    def __update_summary(self):
        '''Update summary.'''
        self.__dfs['identity'].fillna(0, inplace=True)

        numerical_df = self.__dfs['identity'].drop(
            ['plate_idx'], axis=1).select_dtypes(include=[np.float])

        self.__dfs['summary']['matched_ice_id'] = numerical_df.idxmax(axis=1)
        self.__dfs['summary']['identity'] = numerical_df.max(axis=1)

        for name, df in self.__dfs.items():
            if name != 'summary':
                self.__dfs['summary'][name] = \
                    df.lookup(df.index,
                              self.__dfs['summary']['matched_ice_id'])

        # Remove spurious unidentified entries:
        self.__dfs['summary'] = \
            self.__dfs['summary'][self.__dfs['summary']['identity'] != 0]

        # Sort:
        self.__dfs['summary'].sort_values('identity', ascending=False,
                                          inplace=True)
