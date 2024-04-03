#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday April 2 10:00:09 2024

@author: patrick.woods
"""

#v1  

def features_from_gff(gff, features_column, desired_feature):
    '''
    

    Parameters
    ----------
    gff : String
        The file name corresponding to the gff file to be filtered.
    features_column : Integer
        The column index containing the feature values to extrct.
    desired_feature : String
        A string containing the name of the genetic feature to be extracted from the specified column index such as 'gene' or 'mRNA'.

    Returns
    -------
    file_sub : A pandas DataFrame type object
        A pandas DataFrame type object containing only rows with the specified desired_feature string.

    '''
    import pandas as pd
    import numpy as np
    
    file = pd.read_table(gff, header= None)
    filt1 = file[features_column] == desired_feature
    file_sub = file[filt1]
    return file_sub

#testing out the function
features_from_gff('transcripts.fasta.transdecoder.genomeCentric.gff3', features_column = 2, desired_feature = 'gene')
