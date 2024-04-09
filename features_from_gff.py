#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tuesday April 9 11:43:02 2024

@author: patrick.woods
"""

#v2

def features_from_gff(gff, features_column, desired_feature, contig_column=0, desired_contig='all'):
    '''
    

    Parameters
    ----------
    gff : String
        The file name corresponding to the gff file to be filtered.
    features_column : Integer
        The column index containing the feature values to extrct.
    desired_feature : String
        A string containing the name of the genetic feature to be extracted from the specified column index such as 'gene' or 'mRNA'.
    contig_column : Integer
        The column index containing the contig values to extract. Defaults to 0.
    desired_contig : String
        Defaults to 'all'. If a contig is specified, the function will return a dataframe of specified features from the specified contig.

    Returns
    -------
    file_sub : A pandas DataFrame type object
        A pandas DataFrame type object containing only rows with the specified desired_feature string.

    '''
    import pandas as pd
    import numpy as np
    
    file = pd.read_table(gff, header= None)
    
    feature_filt = file[features_column] == desired_feature
    file_sub_feat = file[feature_filt]
    
    if desired_contig == 'all':
        return file_sub_feat
        
    else:
        contig_filt = file_sub_feat[contig_column] == desired_contig
        file_sub_contig = file_sub_feat[contig_filt]
    
        return file_sub_contig


# Testing out function
features_from_gff('../transcripts.fasta.transdecoder.genomeCentric.gff3', features_column = 2, desired_feature = 'gene', desired_contig='Scaffold_1533')








