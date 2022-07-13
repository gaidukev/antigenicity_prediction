from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import mutual_info_classif
import numpy as np
import pandas as pd
import binary_vectors
from sklearn.feature_selection import SelectKBest
from sklearn import metrics

import random
random.seed()
# seed the random for when we're randomly selecting rows to add / exclude

# note: i'm not sure if any of this is correct

def run_regression(len_sequence, start_pos, distances, alignment, verbose = True):
    if verbose:
        # this is just for debugging purposes
        print("START POS: ", start_pos, "LEN SEQUENCE: ", len_sequence)

    # open the file containing all distances and extract the sera
    dist = pd.read_csv(distances)
    sera = list(dist.columns)

        
    serum = sera[1]
    # construct binary features using the first serum x all antigens
    x_vals = binary_vectors.construct_x(serum, alignment, distances, start_pos, len_sequence)
    y_vals = binary_vectors.construct_y(serum, distances)

    for i in range(2, len(sera)):
        # construct binary features using all other sera x all antigens
        serum = sera[i]
        x_return = binary_vectors.construct_x(serum, alignment, distances, start_pos, len_sequence)
        y_return = binary_vectors.construct_y(serum, distances)

        for val in range(len(y_return)):

            if y_return[val] == 1:
                chance = random.random()
                # only append 1-entries with a certain probability (to make data less skewed & increase the mutal info
                # scores) (is this allowed?)
                if chance > 0.6:
                    x_vals = np.vstack((x_vals, x_return[val])) 
                    y_vals = np.hstack((y_vals, y_return[val]))
            else:
                # always append the 0-entries because there's not that many of them
                x_vals = np.vstack((x_vals, x_return[val])) 
                y_vals = np.hstack((y_vals, y_return[val]))




    x_vals = np.array(x_vals, np.float32)
    y_vals = np.array(y_vals, np.uint8)
    p = np.random.permutation(len(x_vals))
    # create permutation of the x- and y-values (together)

    x_vals, y_vals = x_vals[p], y_vals[p]

    # fit the actual regression
    clf = LogisticRegression().fit(x_vals, y_vals)
    print("SCORE: ", clf.score(x_vals, y_vals))
    y_pred = clf.predict(x_vals)
    metric = metrics.mutual_info_score(y_vals, y_pred)
    # this gives the mutual info
    print(metric)
    return metric








