#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday April 19 04:28:34 2024

@author: patrick.woods
"""

#v3

def features_from_gff(gff, features_column, desired_feature, contig_column=0, desired_contig='all',start_col=3,stop_col=4, start = None, stop = None):
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
    start_col : Integer
        Defaults to 4 since this is typically the column index of the start column.
    start : Integer
        Defaults to None. Used to define the starting base pair position you would like to filter for.
    stop_col : Integer
        Defaults to 5 since this is typically the column index of the stop column.
    stop : Integer
        Defaults to None. Used to define the stopping base pair position you would like to filter for.
    Returns
    -------
    file_sub_feat : A pandas DataFrame type object.
        A pandas DataFrame type object containing only rows with the specified desired_feature string.
    OR
    
    file_feat_contig_filt : A pandas DataFrame type object.
        A pandas DataFrame type object containing only rows with the specified desired_feature string on the specified contig.
    OR
    
    file_feat_contig_start_stop_filt : A pandas DataFrame type object.
        A pandas DataFrame type object containing only rows with the specified desired_feature string on the specified contig within the start and stop base pair interval.
    OR
    
    file_feat_contig_start_filt : A pandas DataFrame type object.
        A pandas DataFrame type object containing only rows with the specified desired_feature string on the specified contig after the specified start base pair position.
    OR
    
    file_feat_contig_stop_filt : A pandas DataFrame type object.
        A pandas DataFrame type object containing only rows with the specified desired_feature string on the specified contig before the specified stop base pair position.
    '''
    
    import pandas as pd
    
    file = pd.read_table(gff, header= None)
    ### filtering for a feature across the whole genome ###    
    if desired_contig == 'all' and start == None and stop == None:
        
        print('Filtering for every '+str(desired_feature)+' in the gff file.')
        
        feature_filt = file[features_column] == desired_feature
        file_sub_feat = file[feature_filt]
        
        return file_sub_feat
    

    ### filtering for a feature on a specified contig ###
    elif desired_contig != 'all' and start == None and stop == None:
        
        print('Filtering for every ' + str(desired_feature) + ' on contig: ' + str(desired_contig))
        
        feature_filt = file[features_column] == desired_feature
        file_feat_filt = file[feature_filt]
        
        contig_filt = file_feat_filt[contig_column] == desired_contig
        file_feat_contig_filt = file_feat_filt[contig_filt]
    
        return file_feat_contig_filt
    
    
    ### filtering for a feature on a specified contig within a specific start and stop base pair interval ###
    elif desired_contig != 'all' and start != None and stop != None:
        
        print('Filtering for ' + str(desired_feature) + ' on contig: ' + str(desired_contig) + ' using '+str(start)+' and ' +str(stop)  +' as the start and stop positions, respectively.')
        
        feature_filt = file[features_column] == desired_feature
        file_feat_filt = file[feature_filt]
        
        contig_filt = file_feat_filt[contig_column] == desired_contig
        file_feat_contig_filt = file_feat_filt[contig_filt]
    
        start_filt = file_feat_contig_filt[start_col] >= start
        file_feat_contig_start_filt = file_feat_contig_filt[start_filt]
        
        stop_filt = file_feat_contig_start_filt[stop_col] <= stop
        file_feat_contig_start_stop_filt = file_feat_contig_start_filt[stop_filt]
        
        return file_feat_contig_start_stop_filt
    
    ### filtering for a feature on a specified contig with only a starting position specified ###
    elif desired_contig != 'all' and start != None and stop == None:
        
        print('Filtering for every' + str(desired_feature) + ' on contig: ' + str(desired_contig) + ' using '+str(start)+' as a start base pair cutoff position.')
        
        feature_filt = file[features_column] == desired_feature
        file_feat_filt = file[feature_filt]
        
        contig_filt = file_feat_filt[contig_column] == desired_contig
        file_feat_contig_filt = file_feat_filt[contig_filt]
    
        start_filt = file_feat_contig_filt[start_col] >= start
        file_feat_contig_start_filt = file_feat_contig_filt[start_filt]
        
        return file_feat_contig_start_filt
    
    ### filtering for a feature on a specified contig with only a stop position specified ###
    elif desired_contig != 'all' and start == None and stop != None:
        
        print('Filtering for every' + str(desired_feature) + ' on contig: ' + str(desired_contig) + ' using '+str(stop)+' as a stop base pair cutoff position.')
        
        feature_filt = file[features_column] == desired_feature
        file_feat_filt = file[feature_filt]
        
        contig_filt = file_feat_filt[contig_column] == desired_contig
        file_feat_contig_filt = file_feat_filt[contig_filt]
        
        stop_filt = file_feat_contig_filt[stop_col] <= stop
        file_feat_contig_stop_filt = file_feat_contig_filt[stop_filt]
        
        return file_feat_contig_stop_filt
    
### Example code for testing out the function's various parameters ###        
gf_standard=  features_from_gff_v3('transcripts.fasta.transdecoder.genomeCentric.gff3',features_column=2, desired_feature = 'gene')
gf_standard

gf_contig = features_from_gff_v3('transcripts.fasta.transdecoder.genomeCentric.gff3',features_column=2, desired_feature = 'gene', desired_contig='Scaffold_1531')
gf_contig

gf_contig_start = features_from_gff_v3('transcripts.fasta.transdecoder.genomeCentric.gff3',features_column=2, desired_feature = 'gene', desired_contig='Scaffold_1531', start = 5000000)
gf_contig_start

gf_contig_stop = features_from_gff_v3('transcripts.fasta.transdecoder.genomeCentric.gff3',features_column=2, desired_feature = 'gene', desired_contig='Scaffold_1531', stop = 5000000) #####not working
gf_contig_stop

gf_contig_start_stop = features_from_gff_v3('transcripts.fasta.transdecoder.genomeCentric.gff3',features_column=2, desired_feature = 'gene', desired_contig='Scaffold_1531',start=1, stop = 6000000)
gf_contig_start_stop

