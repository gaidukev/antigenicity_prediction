
import pandas as pd
import numpy as np
import math
import json

# list of relevant aa indeces as per the paper
aaindex_list = [pd.read_excel("AZAE970101.xlsx", index_col="Unnamed: 0"), pd.read_excel("BENS940104.xlsx", index_col="Unnamed: 0"), pd.read_excel("BETM990101.xlsx", index_col="Unnamed: 0"), pd.read_excel("BONM030103.xlsx", index_col="Unnamed: 0"), pd.read_excel("BONM030104.xlsx", index_col="Unnamed: 0"),
                pd.read_excel("KOLA920101.xlsx", index_col="Unnamed: 0"), pd.read_excel("LUTR910108.xlsx", index_col="Unnamed: 0"), pd.read_excel("MUET010101.xlsx", index_col="Unnamed: 0"), pd.read_excel("TANS760101.xlsx", index_col="Unnamed: 0"), pd.read_excel("ZHAC000106.xlsx", index_col="Unnamed: 0")]

# load positions that have been pre-determined using the logistic regression
pos = []
with open("positions.json", "r") as file:
    pos = json.load(file)

filtered_df = ""


def encode(shape, serum_seq, antigen_seq):
    global aaindex_list
    # shape in 
    # (number of amino acids, number of aaindex)

    global pos
    answer = np.zeros((len(pos), 10))
    
    i = 0
    for position in pos:
        if serum_seq[position] != antigen_seq[position]:
            try:
                for j in range(len(aaindex_list)):
                    df = aaindex_list[j]
                    letter_val = df.loc[serum_seq[position], antigen_seq[position]]
                    if math.isnan(letter_val):
                        letter_val = -(df.loc[antigen_seq[position], serum_seq[position]])
                    #print(answer.shape, i, j)
                    answer[i][j] = letter_val
            except:
                pass
        i = i + 1

            
            

    return answer


sequence_dict = {}
def filter_df(df):
    global filtered_df
    global sequence_dict

    df_new = df
    for i in df.index:
        if i not in sequence_dict.keys():
            df_new = df_new.drop(i, axis=0)

    filtered_df = df_new

def construct_x(serum, sequences, distances, len_sequence):
    # takes in the serum name and all antigen sequences 
    # also: the csv containing all antigenic distances, 
    # start position, length of sequence, 
    # important aaindex for use by aaindex_encode, 
    # and the list of all important positions 
    global aaindex_list
    global pos
    global sequence_dict
    num_pos = len(pos)

    df = pd.read_csv(distances, index_col = 0)
    global filtered_df


    # SHAPE: (num_antigens, 566, 10)
    # where 566 is the standardized length of all the aligned sequences
    # and 10 is the # of aaindeces being examined


    sep = "/"
    with open(sequences, "r") as seq_file:
        


        if len(sequence_dict.keys())==0:
            sequence_dict = json.load(seq_file)

        if type(filtered_df) == str:
            filter_df(df)
        df = filtered_df

            
        x = np.zeros((df.shape[0], num_pos, len(aaindex_list)))




        serum_sequence = sequence_dict[serum]


        i = 0
        keys_list = list(sequence_dict.keys())
        keys_list.sort()
        for key in keys_list:
            if key in df.index.values:
                try: 
                    substring = np.zeros((len_sequence, 10))
                    substring = np.array(encode(substring.shape, serum_sequence, sequence_dict[key]))
                    x[i] = substring
                    i += 1

                except ValueError as e:
                    print(e)
                    
    return x


def scale(vector):
    min = vector[0]
    max = vector[0]
    for val in vector:
        if val < min:
            min = val
        if val > max:
            max = val
        
    for i in range(len(vector)):
        vector[i] = (vector[i] - min) / (max - min)
    return vector



def construct_y(serum, distances):
    # y values are the distances between each antigen-serum pair
    global filtered_df
    df = pd.read_csv(distances)

    if type(filtered_df) == str:
        filter_df(df)
    df = filtered_df
    y = np.zeros(df.shape[0], )
    for antigen_index in range(y.shape[0]):
        #y_val = df.iloc[antigen_index].loc[serum]
        # if y_val > 4:
        #     y[antigen_index] = 1
        # else:
        #     y[antigen_index] = 0
        y[antigen_index] = df.iloc[antigen_index].loc[serum]
    return y
