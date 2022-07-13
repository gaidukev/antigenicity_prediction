# this file encodes antigen-serum pairs into binary x- and y- values
# it is called from inside reg_2.py (file containing logistic regressor)
# so there's no need to call this file on its own

import pandas as pd
import numpy as np

def compare_string(len_sequence, str1, str2):
    # compares two strings (an antigen and a serum)
    # and encodes all differences with 1s
    try:
        assert(len(str1) == len(str2) and len(str1) == len_sequence)
    except AssertionError:
        # suppressing assertion error is probably a bad idea?
        return np.zeros(len_sequence)
    answer = np.zeros(len_sequence)
    for i in range(len_sequence):
        # if dash in that spot, also encode it as a 1
        if str1[i] == "-" or str2[i] == '-' or str1[i] != str2[i]:
            answer[i] = 1
    return answer
        

sequence_dict = {}
# sequence dict contains the names of all the strains 
# and their genetic sequences
# we keep it global so that it's only created ones between different function calls & saves
# computation time with having to parse sequences multiple times

def construct_x(serum, sequences, distances, pos, len_sequence = 1):
    '''
    serum = string of format "BI/15793/68"
    sequences, distances = .csv files
    -> distances: sera are columns, antigens are rows
    Will compare all the antigens against a specific serum & create respective features
    '''
    df = pd.read_csv(distances, index_col = 0)
    x = np.zeros((df.shape[0], len_sequence))
    global sequence_dict

    with open(sequences, "r") as seq_file:
        content = "".join(seq_file.read()).split(">")
        serum_sequence = ""

        if len(sequence_dict.keys())==0:
            # there are no entries in sequence_dict yet
            # ie. the sequences have yet to be parsed.
            for line in content:
                entries = line.split("\n")
                try:
                    name = entries[0].split("/")[0]
                    sequence = "".join(entries[1:])
                    if name not in sequence_dict.keys():
                        # store a sequence under its name
                        sequence_dict[name] = sequence
                except Exception as e:
                    print(e)

        # only grab relevant section of the serum's sequence
        serum_sequence = sequence_dict[serum][pos:pos+len_sequence]


        i = 0
        keys_list = list(sequence_dict.keys())
        keys_list.sort()
        # sorting guarantees same order each time
        for key in keys_list:
            # check if this item is an antigen
            if key in df.index.values:
                try: 
                    substring = np.array(compare_string(len_sequence, sequence_dict[key][pos:pos+len_sequence], serum_sequence))
                    x[i] = substring
                    i += 1

                except ValueError as e:
                    print(e)
        
    return x


def construct_y(serum, distances):
    # 0 if distance <4, 1 otherwise
    # will create a binary vector mapping all antigens to 0 or 1 
    # (depending on how large the genetic distance between the antigen
    # and this specific serum)
    df = pd.read_csv(distances)
    y = np.zeros(df.shape[0], )
    for antigen_index in range(y.shape[0]):

        if df.iloc[antigen_index].loc[serum] >= 4:
            y[antigen_index] = 1
    return y


