"""Returns a set of elements sampled from the second bed file that matches the size distribution of the first bed file"""

import argparse
import pathlib
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

def argparser():
    parser = argparse.ArgumentParser('Matches the size distribtuon of the second to the first bed file.')
    parser.add_argument('-a','--file1',help='Path to the first bedfile.')
    parser.add_argument('-b','--file2',help='Path to the second bedfile.')
    return parser.parse_args()

def ecdf(data: pd.Series):
    #create sorted array of unique elements
    uniques = np.sort(data.unique())
    #compute envenly spaced elements as x elements for the ecdf
    ecdf_x = np.linspace(start=min(uniques),stop=max(uniques),num=60)
    #get size of the raw values
    size = data.shape[0]
    #compute y values for the ecdf
    ecdf_y = []
    for x in ecdf_x:
        #count raw values below or equal to x
        cum_amount = data[data <= x].shape[0]
        data = data[data >= x]
        #save the cummulative amount of values as y values
        ecdf_y.append(cum_amount)
    return ecdf_x.tolist(),ecdf_y



def main():
    #read cli
    args = argparser()

    # read bed files as tables
    file1_df = pd.read_csv(pathlib.Path(args.file1),sep='\t',names=['chr','start','end','source'],usecols=[0,1,2,3])
    file2_df = pd.read_csv(pathlib.Path(args.file2),sep='\t',names=['chr','start','end','source'],usecols=[0,1,2,3])

    file1_df['size'] = file1_df['end'] - file1_df['start']
    file2_df['size'] = file2_df['end'] - file2_df['start']

    #calcuate the log10 of every size to account for the large tail of the size distributions
    file1_df['size'] = file1_df['size'].apply(np.log10)
    file2_df['size'] = file2_df['size'].apply(np.log10)

    file1_orig = file1_df.copy()
    file2_orig = file2_df.copy()

    # Define source dict to parse data sources into ints.
    # This is needed to apply weight according to the data origin in the sampling process.
    source_dict = {'Eichler':4,'GnomAD':3,'DGV':2,'Biobank':1,'Decipher':5}
    rev_source_dict = {4:'Eichler',3:'GnomAD',2:'DGV',1:'Biobank',5:'Decipher'}

    file1_df['source'] = [source_dict[source] for source in file1_df['source']]
    file2_df['source'] = [source_dict[source] for source in file2_df['source']]

    #compute the ecdf for file1 and return the x (bin values) and y values (number of samples in each bin)
    ecdf_x,ecdf_y = ecdf(file1_df['size'])

    matched_elements_patho = []
    matched_elements_non_patho = []
    variant_counter = 0
    #bin the elements of the second file into the same bins
    for (x,y) in zip(ecdf_x,ecdf_y):
        file2_matches = file2_df[file2_df['size'] <= x]
        file1_matches = file1_df[file1_df['size'] <= x]
        file2_df = file2_df[file2_df['size'] > x]
        file1_df = file1_df[file1_df['size'] > x]
        #if the number of matches is less or equal than in file1 take all
        if file2_matches.shape[0] <= y:
            matched_elements_patho.append(file1_matches.sample(file2_matches.shape[0],random_state=42, weights='source'))
            matched_elements_non_patho.append(file2_matches)
            #variant_counter += file2_matches.shape[0]
        else:
            matched_elements_patho.append(file1_matches)
            matched_elements_non_patho.append(file2_matches.sample(y,random_state=42,weights='source'))

            #if variant_counter + file2_matches.shape[0] >= file1_df.shape[0]:
            #    matched_elements.append(file2_matches.sample(file1_df.shape[0]-variant_counter,random_state=42))
            #    variant_counter = file1_df.shape[0]
            #else:
            #    matched_elements.append(file2_matches)
            #    variant_counter += file2_matches.shape[0]

    size_matched_df_patho = pd.concat(matched_elements_patho)
    size_matched_df_non_patho = pd.concat(matched_elements_non_patho)

    # reverse the source mapping
    size_matched_df_patho['source'] = [rev_source_dict[source] for source in size_matched_df_patho['source']]
    size_matched_df_non_patho['source'] = [rev_source_dict[source] for source in size_matched_df_non_patho['source']]

    # plot the size distribution of pathogenic/non-apthogenic variants according to the ecdf
    plt.figure()
    plt.hist(x='size',density=False,histtype='step',bins=ecdf_x,data=file1_orig,label='Pathogenic')
    plt.hist(x='size',density=False,histtype='step',bins=ecdf_x,data=file2_orig,label='Non-Pathogenic')
    #plt.hist(x='size',density=True,histtype='step',data=file2_orig,label='Pathogenic')
    plt.ylabel('Amount of Variants')
    plt.xlabel('Size log10(bp)')
    plt.legend()
    plt.savefig('Pathogenic_and_non_pathogenic_size_dist.png')

    # plot the size distribution of non-apthogenic variants for each data source
    plt.figure()
    plt.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=file2_orig[file2_orig['source'] == 'GnomAD'],label='GnomAD', alpha=0.8)
    plt.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=file2_orig[file2_orig['source'] == 'Eichler'],label='Eichler', alpha=0.8)
    plt.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=file2_orig[file2_orig['source'] == 'Biobank'],label='Biobank', alpha=0.25)
    plt.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=file2_orig[file2_orig['source'] == 'DGV'],label='DGV', alpha=0.25)
    plt.hist(x='size',density=False,histtype='stepfilled',bins=ecdf_x,data=file1_orig,label='Decipher', alpha=0.8)
    plt.ylabel('Number of Variants')
    plt.xlabel('Size log10(bp)')
    plt.legend()
    plt.savefig('Pathogenic_and_non_pathogenic_size_dist_by_source.png')

    size_matched_df_patho.drop('size',axis='columns',inplace=True)
    size_matched_df_patho.to_csv('file_1_size_matched.bed',header=False,index=False,sep='\t')

    size_matched_df_non_patho.drop('size',axis='columns',inplace=True)
    size_matched_df_non_patho.to_csv('file_2_size_matched.bed',header=False,index=False,sep='\t')





















if __name__ == "__main__":
    main()
