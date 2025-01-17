import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import sys
import os
# regular expression for searching infer_exp file
import re

# capturing the variables given by the bash file
acc_path = sys.argv[1]

virus_name = sys.argv[2]


# Function to read accessions from a file and return them as a list
def read_acc_from_file(file_path):
    acc = []
    with open(file_path, 'r') as file:
        for line in file:
            # Strip the newline character and convert the line to an integer
            accession = str(line.strip())
            acc.append(accession)
    return acc

# getting accession list
acc_list = read_acc_from_file(acc_path)

print(acc_list)


from cigar import Cigar

# Cigar function to calculate length of a read from cigar string
def cigar_len(x):
    return Cigar(x).reference_length()


for i in acc_list:

    mismatches = pd.read_csv("work/mismatch/" + i + "_filtered.tsv", sep="\t")
    mismatches = mismatches.drop(columns=['ref', 'base', 'avg_basequality','avg_pos_as_fraction'])
    mismatches = mismatches.set_index("chrom")

    original_depth = pd.read_csv("work/plotting/original_depth/" + i + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
    
    # so I need to create full_range_position which is a temporary 1D dataframe, with a range from the lowest position to the
    # highest position from original_depth. For use later on
    full_range_position = pd.DataFrame({'position': range(original_depth['Position'].min(), original_depth['Position'].max() + 1)})
    
    # I am merging mismatches with full_range_positions, so that if there are any positions from mismatches that have no values
    # I can capture those with this merge function. Ensuring that each position has some value in the mismatch file. 
    mismatches_complete = pd.merge(full_range_position, mismatches, on='position', how='left').fillna(0)
    mismatches_complete = mismatches_complete.drop(columns=['depth','count'])
    
    mismatches_complete.to_csv("work/plotting/mismatches/" + i +".tsv", sep="\t", header=False, index=False)

    # checking if the library is paired-end
    if (os.path.isfile('work/trimmed/'+ i +'_1_val_1.fq') and os.path.isfile('work/trimmed/'+ i +'_2_val_2.fq')) or \
       (os.path.isfile('work/trimmed/'+ i +'_1_val_1.fq.gz') and os.path.isfile('work/trimmed/'+ i +'_2_val_2.fq.gz')):

        # full_range_start and full_range_end serve the same function here as well. I will merge them later on so that every position
        # has some value in the tsv file. It helps with matplotlib otherwise, the plotting would be janky
        original_depth = pd.read_csv("work/plotting/original_depth/" + i + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        full_range_start = pd.DataFrame({'Start': range(original_depth['Position'].min(), original_depth['Position'].max() + 1)})
        full_range_end = pd.DataFrame({'End': range(original_depth['Position'].min(), original_depth['Position'].max() + 1)})



        positive_read1 = pd.read_csv("work/read_start_end/" + i + "_sense_read1.tsv", sep="\t", names=["Start", "CIGAR"])
        
        # applying the cigar read length function
        positive_read1['Length'] = positive_read1['CIGAR'].apply(cigar_len)
        positive_read1['End'] = positive_read1['Start'] + positive_read1['Length']
        
        positive_read1 = positive_read1.drop('Length', axis=1)
        positive_read1 = positive_read1.drop('CIGAR', axis=1)

        positive_read2 = pd.read_csv("work/read_start_end/" + i + "_sense_read2.tsv", sep="\t", names=["Start", "CIGAR"])
        
        # applying the cigar read length function
        positive_read2['Length'] = positive_read2['CIGAR'].apply(cigar_len)
        positive_read2['End'] = positive_read2['Start'] + positive_read2['Length']
        
        positive_read2 = positive_read2.drop('Length', axis=1)
        positive_read2 = positive_read2.drop('CIGAR', axis=1)


        
        negative_read1 = pd.read_csv("work/read_start_end/" +i + "_antisense_read1.tsv", sep="\t", names=["Start", "CIGAR"])
        
        negative_read1['Length'] = negative_read1['CIGAR'].apply(cigar_len)
        negative_read1['End'] = negative_read1['Start'] + negative_read1['Length']
        
        negative_read1 = negative_read1.drop('Length', axis=1)
        negative_read1 = negative_read1.drop('CIGAR', axis=1)

        negative_read2 = pd.read_csv("work/read_start_end/" +i + "_antisense_read2.tsv", sep="\t", names=["Start", "CIGAR"])
        
        negative_read2['Length'] = negative_read2['CIGAR'].apply(cigar_len)
        negative_read2['End'] = negative_read2['Start'] + negative_read2['Length']
        
        negative_read2 = negative_read2.drop('Length', axis=1)
        negative_read2 = negative_read2.drop('CIGAR', axis=1)



        # I need to count the positions for the 5' and 3' peak plots
        counts_pos_1 = positive_read1['Start'].value_counts().reset_index()
        counts_neg_1 = negative_read1['Start'].value_counts().reset_index()

        counts_pos_2 = positive_read2['Start'].value_counts().reset_index()
        counts_neg_2 = negative_read2['Start'].value_counts().reset_index()
        
        # Rename the columns for clarity
        counts_pos_1.columns = ['Start', 'Count']
        counts_neg_1.columns = ['Start', 'Count']
        counts_pos_2.columns = ['Start', 'Count']
        counts_neg_2.columns = ['Start', 'Count']
        
        counts_pos_1 = counts_pos_1.sort_values(by=['Start'])
        counts_neg_1 = counts_neg_1.sort_values(by=['Start'])
        counts_pos_2 = counts_pos_2.sort_values(by=['Start'])
        counts_neg_2 = counts_neg_2.sort_values(by=['Start'])

        # Merging again to ensure each position has a value, 0 if it was 'na' before.
        counts_pos_1_complete = pd.merge(full_range_start, counts_pos_1, on='Start', how='left').fillna(0)
        counts_neg_1_complete = pd.merge(full_range_start, counts_neg_1, on='Start', how='left').fillna(0)
        counts_pos_2_complete = pd.merge(full_range_start, counts_pos_2, on='Start', how='left').fillna(0)
        counts_neg_2_complete = pd.merge(full_range_start, counts_neg_2, on='Start', how='left').fillna(0)
        
        counts_neg_1_complete.to_csv("work/plotting/paired/negative_start/" + i + "_1.tsv", sep="\t", header=False, index=False)
        counts_pos_1_complete.to_csv("work/plotting/paired/positive_start/" + i + "_1.tsv", sep="\t", header=False, index=False)
        counts_neg_2_complete.to_csv("work/plotting/paired/negative_start/" + i + "_2.tsv", sep="\t", header=False, index=False)
        counts_pos_2_complete.to_csv("work/plotting/paired/positive_start/" + i + "_2.tsv", sep="\t", header=False, index=False)



        # I need to count the positions for the 5' and 3' peak plots
        counts_pos_1 = positive_read1['End'].value_counts().reset_index()
        counts_neg_1 = negative_read1['End'].value_counts().reset_index()

        counts_pos_2 = positive_read2['End'].value_counts().reset_index()
        counts_neg_2 = negative_read2['End'].value_counts().reset_index()
        
        # Rename the columns for clarity
        counts_pos_1.columns = ['End', 'Count']
        counts_neg_1.columns = ['End', 'Count']
        counts_pos_2.columns = ['End', 'Count']
        counts_neg_2.columns = ['End', 'Count']
        
        counts_pos_1 = counts_pos_1.sort_values(by=['End'])
        counts_neg_1 = counts_neg_1.sort_values(by=['End'])
        counts_pos_2 = counts_pos_2.sort_values(by=['End'])
        counts_neg_2 = counts_neg_2.sort_values(by=['End'])

        # Merging again to ensure each position has a value, 0 if it was 'na' before.
        counts_pos_1_complete = pd.merge(full_range_end, counts_pos_1, on='End', how='left').fillna(0)
        counts_neg_1_complete = pd.merge(full_range_end, counts_neg_1, on='End', how='left').fillna(0)
        counts_pos_2_complete = pd.merge(full_range_end, counts_pos_2, on='End', how='left').fillna(0)
        counts_neg_2_complete = pd.merge(full_range_end, counts_neg_2, on='End', how='left').fillna(0)
        
        counts_neg_1_complete.to_csv("work/plotting/paired/negative_end/" + i + "_1.tsv", sep="\t", header=False, index=False)
        counts_pos_1_complete.to_csv("work/plotting/paired/positive_end/" + i + "_1.tsv", sep="\t", header=False, index=False)
        counts_neg_2_complete.to_csv("work/plotting/paired/negative_end/" + i + "_2.tsv", sep="\t", header=False, index=False)
        counts_pos_2_complete.to_csv("work/plotting/paired/positive_end/" + i + "_2.tsv", sep="\t", header=False, index=False)

    
    # if the library is single-end then do this
    # Commands are the same as above, just taking into account that they are single, and not paired-end
    elif os.path.isfile('work/trimmed/'+ i +'_trimmed.fq') or os.path.isfile('work/trimmed/'+ i +'_trimmed.fq.gz'):
    
        original_depth = pd.read_csv("work/plotting/original_depth/" + i + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        full_range_start = pd.DataFrame({'Start': range(original_depth['Position'].min(), original_depth['Position'].max() + 1)})
        full_range_end = pd.DataFrame({'End': range(original_depth['Position'].min(), original_depth['Position'].max() + 1)})



        positive = pd.read_csv("work/read_start_end/" + i + "_sense.tsv", sep="\t", names=["Start", "CIGAR"])
        
        positive['Length'] = positive['CIGAR'].apply(cigar_len)
        positive['End'] = positive['Start'] + positive['Length']
        
        positive = positive.drop('Length', axis=1)
        positive = positive.drop('CIGAR', axis=1)


        
        negative = pd.read_csv("work/read_start_end/" +i + "_antisense.tsv", sep="\t", names=["Start", "CIGAR"])
        
        negative['Length'] = negative['CIGAR'].apply(cigar_len)
        negative['End'] = negative['Start'] + negative['Length']
        
        negative = negative.drop('Length', axis=1)
        negative = negative.drop('CIGAR', axis=1)



        
        counts_pos = positive['Start'].value_counts().reset_index()
        counts_neg = negative['Start'].value_counts().reset_index()
        
        # Rename the columns for clarity
        counts_pos.columns = ['Start', 'Count']
        counts_neg.columns = ['Start', 'Count']
        
        counts_pos = counts_pos.sort_values(by=['Start'])
        counts_neg = counts_neg.sort_values(by=['Start'])

        counts_pos_complete = pd.merge(full_range_start, counts_pos, on='Start', how='left').fillna(0)
        counts_neg_complete = pd.merge(full_range_start, counts_neg, on='Start', how='left').fillna(0)
        
        counts_neg_complete.to_csv("work/plotting/single/negative_start/" + i + ".tsv", sep="\t", header=False, index=False)
        counts_pos_complete.to_csv("work/plotting/single/positive_start/" + i + ".tsv", sep="\t", header=False, index=False)


        
        counts_pos = positive['End'].value_counts().reset_index()
        counts_neg = negative['End'].value_counts().reset_index()
        
        # Rename the columns for clarity
        counts_pos.columns = ['End', 'Count']
        counts_neg.columns = ['End', 'Count']

        counts_pos = counts_pos.sort_values(by=['End'])
        counts_neg = counts_neg.sort_values(by=['End'])

        counts_pos_complete = pd.merge(full_range_end, counts_pos, on='End', how='left').fillna(0)
        counts_neg_complete = pd.merge(full_range_end, counts_neg, on='End', how='left').fillna(0)

        counts_neg_complete.to_csv("work/plotting/single/negative_end/" + i + ".tsv", sep="\t", header=False, index=False)
        counts_pos_complete.to_csv("work/plotting/single/positive_end/" + i + ".tsv", sep="\t", header=False, index=False)

    # If the trimmed file was not found in the directory, then print this.
    else:
        print("This accession was neither single ended nor paired. Or this accession was not trimmed.")


# My model to classify the change point analysis using the mean differnce approach
def classifier_model(coverage, top3peaks, chopped, buffer, detection_zone):

    # creating a copy to avoid pointer errors
    # mean_diff is my main array where I will store the mean difference values.
    mean_diff = coverage.copy()
    
    # Chopping off "chopped" number of bases from the ends
    coverage = coverage[chopped:-chopped]    

    # in mean_diff, make the position from 0 to chopped+buffer = 0 as they are not counting now
    for i in range(chopped+buffer):
        mean_diff[i][0]=0

    # in mean_diff, make the end of the genome which is chopped and buffered = 0
    for i in range(len(mean_diff)-(chopped+buffer),len(mean_diff)):
        mean_diff[i][0]=0
    # for i in mean_diff:
    #     print(i)
    
    # range of buffer, len(coverage)-buffer is everything except the chopped portions
    for index in range(buffer,len(coverage)-buffer):
        left_mean = coverage[:index].mean()
        right_mean = coverage[index:].mean()
        diff = abs(right_mean - left_mean)
        # print(diff)
        mean_difference = diff
        # have to add chopped, or else it will be the wrong position
        mean_diff[index+chopped][0]=mean_difference
    
    # getting index of largest mean difference, which corresponds to position
    result = mean_diff.argmax()
    # plt.plot(mean_diff)
    # plt.show()
        
    # this is the "zone" in which we are trying to find if the peaks lie
    if (top3peaks[0][0] - detection_zone) <= result <= (top3peaks[0][0] + detection_zone) or \
        (top3peaks[1][0] - detection_zone) <= result <= (top3peaks[1][0] + detection_zone) or \
        (top3peaks[2][0] - detection_zone) <= result <= (top3peaks[2][0] + detection_zone):
        verdict = True
    else:
        verdict = False

    return result, verdict


# fine tuned variables, which can be edited
chopped = 100
buffer = 225
detection_zone = 100

for acc in acc_list:
    
    if (os.path.isfile('work/trimmed/'+ acc +'_1_val_1.fq') and os.path.isfile('work/trimmed/'+ acc +'_2_val_2.fq')) or \
       (os.path.isfile('work/trimmed/'+ acc +'_1_val_1.fq.gz') and os.path.isfile('work/trimmed/'+ acc +'_2_val_2.fq.gz')):

        original_depth = pd.read_csv("work/plotting/original_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        mismatch = pd.read_csv("work/plotting/mismatches/" + acc + ".tsv", sep="\t", names=["Position", "Count"])

        positive_depth = pd.read_csv("work/plotting/paired/positive_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        negative_depth = pd.read_csv("work/plotting/paired/negative_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])


        positive_start_1 = pd.read_csv("work/plotting/paired/positive_start/" + acc + "_1.tsv", sep="\t", names=["Position", "Count"])
        negative_start_1 = pd.read_csv("work/plotting/paired/negative_start/" + acc + "_1.tsv", sep="\t", names=["Position", "Count"])
        positive_start_2 = pd.read_csv("work/plotting/paired/positive_start/" + acc + "_2.tsv", sep="\t", names=["Position", "Count"])
        negative_start_2 = pd.read_csv("work/plotting/paired/negative_start/" + acc + "_2.tsv", sep="\t", names=["Position", "Count"])

        positive_end_1 = pd.read_csv("work/plotting/paired/positive_end/" + acc + "_1.tsv", sep="\t", names=["Position", "Count"])
        negative_end_1 = pd.read_csv("work/plotting/paired/negative_end/" + acc + "_1.tsv", sep="\t", names=["Position", "Count"])
        positive_end_2 = pd.read_csv("work/plotting/paired/positive_end/" + acc + "_2.tsv", sep="\t", names=["Position", "Count"])
        negative_end_2 = pd.read_csv("work/plotting/paired/negative_end/" + acc + "_2.tsv", sep="\t", names=["Position", "Count"])
        

        #Getting read Orientation from Infer Experiment

        # Open the file in read mode
        with open("work/infer_exp/" + acc + ".infer_experiment.txt", 'r') as file:
            # Iterate through each line in the file
            for line in file:
                # Check if the line contains the first pattern
                match1 = re.search(r'Fraction of reads explained by "1\+\+,1--,2\+-\,2-\+": (\d\.\d+)', line)
                if match1:
                    # Extract the number after the first pattern
                    number1 = match1.group(1)
                    number1 = float(number1)
                
                # Check if the line contains the second pattern
                match2 = re.search(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--": (\d\.\d+)', line)
                if match2:
                    # Extract the number after the second pattern
                    number2 = match2.group(1)
                    number2 = float(number2)

        tolerance = number1 * 0.1
        lower_bound = number1-tolerance
        upper_bound = number1+tolerance

        # similar thing here, if the file is unstranded and within 10% of each other than default to "read_1"
        if number2 >= lower_bound and number2 <= upper_bound:
            read_orientation="read_1"
        elif number1 > number2:
            read_orientation = "read_1"
        elif number1 < number2:
            read_orientation = "read_2"
        
        #5' start peak detection
        if read_orientation=='read_1':
            positive_start2 = pd.read_csv("work/plotting/paired/positive_start/" + acc + "_1.tsv", sep="\t", names=["Position", "Count"])
        elif read_orientation=='read_2':
            positive_start2 = pd.read_csv("work/plotting/paired/positive_start/" + acc + "_2.tsv", sep="\t", names=["Position", "Count"])
        positive_start2 = positive_start2.sort_values(by=['Count'],ascending=False)
        top3peaks = positive_start2.head(3)
        top3peaks = top3peaks.to_numpy()


        #Read coverage Change Point Detection
        depth = pd.read_csv("work/plotting/paired/positive_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        depth = depth.drop(columns=["Virus"])
        depth.set_index("Position", inplace=True)

        coverage = depth['Count'].values.reshape(-1,1)

        result, verdict = classifier_model(coverage, top3peaks, chopped, buffer, detection_zone)

        # Create a list of the dataframes for easy iteration
        dataframes = [original_depth, mismatch, positive_depth, negative_depth,
                      positive_start_1, positive_start_2, negative_start_1, negative_start_2,
                      positive_end_1, positive_end_2, negative_end_1, negative_end_2]
        titles = [virus_name + " Coverage", "Percent Reads with a Mismatch vs "+virus_name,
                  acc + " Sense Reads Coverage",acc + " Antisense Reads Coverage",
                  acc + " Sense Read 1 5\' End", acc + " Sense Read 2 5\' End", acc + " Antisense Read 1 3\' End", acc + " Antisense Read 2 3\' End",
                  acc + " Sense Read 1 3\' End", acc + " Sense Read 2 3\' End", acc + " Antisense Read 1 5\' End", acc + " Antisense Read 2 5\' End"]
        colours = ["magma","black",None, "autumn",None,None,"autumn", "autumn",None,None,"autumn", "autumn"]



        # Create subplots: 4 rows and 2 columns, but only using the first slot for a single graph
        fig = plt.figure()
        axs = []
        gs = gridspec.GridSpec(4, 4)
        axs.append(plt.subplot(gs[0, 0]))
        axs.append(plt.subplot(gs[0, 1]))        
        axs.append(plt.subplot(gs[1, 0]))
        axs.append(plt.subplot(gs[1, 1]))
        axs.append(plt.subplot(gs[2, 0]))
        axs.append(plt.subplot(gs[2, 1]))
        axs.append(plt.subplot(gs[2, 2]))
        axs.append(plt.subplot(gs[2, 3]))
        axs.append(plt.subplot(gs[3, 0]))
        axs.append(plt.subplot(gs[3, 1]))
        axs.append(plt.subplot(gs[3, 2]))
        axs.append(plt.subplot(gs[3, 3]))
        
        # Adding the new axis that spans multiple columns and rows
        axs.append(plt.subplot(gs[0:2, 2:4]))  # This spans rows 0 to 2 and columns 2 to 4

        # Plot each DataFrame on its corresponding subplot
        for i, df in enumerate(dataframes):
            ax = axs[i]
            if i==1:
                df.plot.scatter("Position", "Count",ax=ax, c=colours[i], figsize=(36, 18), s=7)
            else:
                if read_orientation=='read_1':
                    if i==4:
                        ax.plot(top3peaks[0][0],top3peaks[0][1], "ob")
                        ax.plot(top3peaks[1][0],top3peaks[1][1], "ob")
                        ax.plot(top3peaks[2][0],top3peaks[2][1], "ob")
                        ax.legend(['Peaks'])
                elif read_orientation=='read_2':
                    if i==5:
                        ax.plot(top3peaks[0][0],top3peaks[0][1], "ob")
                        ax.plot(top3peaks[1][0],top3peaks[1][1], "ob")
                        ax.plot(top3peaks[2][0],top3peaks[2][1], "ob")
                        ax.legend(['Peaks'])
                if i==2:
                    ax.vlines(result,ymin=0,ymax=df['Count'].max(), color='black', linestyles='dashed', linewidth = 4)

                # if i==0 or i==2 or i==3:
                #     # predicted 5′ end of RNA “stem-loop 2” (SL2)
                #     ax.vlines(10478,ymin=0,ymax=df['Count'].max(), color='blue', linestyles='dashed', linewidth = 1)

                #     # predicted 5′ end of RNA “stem-loop 1” (SL1)
                #     ax.vlines(10394,ymin=0,ymax=df['Count'].max(), color='blue', linestyles='dashed', linewidth = 1)
                    

                #     #position from Andrew's paper: https://www.biorxiv.org/content/10.1101/112904v1.full
                    
                df.plot("Position", "Count", ax=ax, colormap=colours[i], figsize=(36, 18))
            ax.set_title(titles[i])


        axs[12].annotate("Prediction: " + str(verdict), (0.26,0.6), fontsize=35, bbox=dict(boxstyle="round", fc="0.8"))
        axs[12].annotate("sfRNA predicted start: " + str(result), (0.16,0.35), fontsize=30, bbox=dict(boxstyle="round", fc="0.8"))

        # Adjust layout for better spacing
        plt.tight_layout()

        # plt.show()
        fig.savefig(virus_name + '_results/' + acc + '.png')
        plt.close(fig)

        # saving the result to a file in the results folder
        with open(virus_name + '_results/' + 'results.txt',"a+") as f:
            f.write(f"{acc}\t{str(verdict)}\t{str(result)}\n")


    elif os.path.isfile('work/trimmed/'+ acc +'_trimmed.fq') or os.path.isfile('work/trimmed/'+ acc +'_trimmed.fq.gz'):

        original_depth = pd.read_csv("work/plotting/original_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        mismatch = pd.read_csv("work/plotting/mismatches/" + acc + ".tsv", sep="\t", names=["Position", "Count"])
        positive_depth = pd.read_csv("work/plotting/single/positive_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        negative_depth = pd.read_csv("work/plotting/single/negative_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        positive_start = pd.read_csv("work/plotting/single/positive_start/" + acc + ".tsv", sep="\t", names=["Position", "Count"])
        negative_start = pd.read_csv("work/plotting/single/negative_start/" + acc + ".tsv", sep="\t", names=["Position", "Count"])
        positive_end = pd.read_csv("work/plotting/single/positive_end/" + acc + ".tsv", sep="\t", names=["Position", "Count"])
        negative_end = pd.read_csv("work/plotting/single/negative_end/" + acc + ".tsv", sep="\t", names=["Position", "Count"])

        #5' start peak detection
        positive_start2 = pd.read_csv("work/plotting/single/positive_start/" + acc + ".tsv", sep="\t", names=["Position", "Count"])
        positive_start2 = positive_start2.sort_values(by=['Count'],ascending=False)
        top3peaks = positive_start2.head(3)
        top3peaks = top3peaks.to_numpy()


        #Read coverage Change Point Detection
        depth = pd.read_csv("work/plotting/single/positive_depth/" + acc + ".tsv", sep="\t", names=["Virus", "Position", "Count"])
        depth = depth.drop(columns=["Virus"])
        depth.set_index("Position", inplace=True)
        
        coverage = depth['Count'].values.reshape(-1,1)

        result, verdict = classifier_model(coverage, top3peaks, chopped, buffer, detection_zone)
        
        # Create a list of the dataframes for easy iteration
        dataframes = [original_depth, mismatch, positive_depth, negative_depth, positive_start, negative_start, positive_end, negative_end]
        titles = [virus_name + " Coverage","Percent Reads with a Mismatch vs "+virus_name, acc + " Sense Reads Coverage",acc + " Antisense Reads Coverage",
                acc + " Sense Read 5\' End",acc + " Antisense Read 3\' End", acc + " Sense Read 3\' End", acc + " Antisense Read 5\' End"]
        colours = ["magma","black",None, "autumn",None,"autumn",None, "autumn"]
        
        
        
        # Create subplots: 4 rows and 2 columns, but only using the first slot for a single graph
        fig = plt.figure()
        axs = []
        gs = gridspec.GridSpec(4, 3)
        axs.append(plt.subplot(gs[0, 0]))
        axs.append(plt.subplot(gs[0, 1]))        
        axs.append(plt.subplot(gs[1, 0]))
        axs.append(plt.subplot(gs[1, 1]))
        axs.append(plt.subplot(gs[2, 0]))
        axs.append(plt.subplot(gs[2, 1]))
        axs.append(plt.subplot(gs[3, 0]))
        axs.append(plt.subplot(gs[3, 1]))

        axs.append(plt.subplot(gs[0:4, 2:3]))  # This spans rows 0 to 2 and columns 2 to 4

        
        # Plot each DataFrame on its corresponding subplot
        for i, df in enumerate(dataframes):
            ax = axs[i]
            if i==1:
                df.plot.scatter("Position", "Count",ax=ax, c=colours[i], figsize=(27, 18), s=7)
            else:
                if i==4:
                    ax.plot(top3peaks[0][0],top3peaks[0][1], "ob")
                    ax.plot(top3peaks[1][0],top3peaks[1][1], "ob")
                    ax.plot(top3peaks[2][0],top3peaks[2][1], "ob")
                    ax.legend(['Peaks'])
                if i==2:
                    ax.vlines(result,ymin=0,ymax=df['Count'].max(), color='black', linestyles='dashed', linewidth = 4)

                # if i==0 or i==2 or i==3:
                #     # predicted 5′ end of RNA “stem-loop 2” (SL2)
                #     ax.vlines(10478,ymin=0,ymax=df['Count'].max(), color='blue', linestyles='dashed', linewidth = 1)

                #     # predicted 5′ end of RNA “stem-loop 1” (SL1)
                #     ax.vlines(10394,ymin=0,ymax=df['Count'].max(), color='blue', linestyles='dashed', linewidth = 1)
                    

                #     #position from Andrew's paper: https://www.biorxiv.org/content/10.1101/112904v1.full
                    
                df.plot("Position", "Count", ax=ax, colormap=colours[i], figsize=(27, 18))
            ax.set_title(titles[i])

        
        axs[8].annotate("Prediction: " + str(verdict), (0.26,0.6), fontsize=35, bbox=dict(boxstyle="round", fc="0.8"))
        axs[8].annotate("sfRNA predicted start: " + str(result), (0.16,0.35), fontsize=30, bbox=dict(boxstyle="round", fc="0.8"))
        
        # Adjust layout for better spacing
        plt.tight_layout()
        
        # plt.show()
        fig.savefig(virus_name + '_results/' + acc + '.png')
        plt.close(fig)

        # saving the result to a file in the results folder
        with open(virus_name + '_results/' + 'results.txt',"a+") as f:
            f.write(f"{acc}\t{str(verdict)}\t{str(result)}\n")
    
    else:
        print("This accession was neither single ended nor paired. Or this accession was not trimmed. So nothing to plot")
