# PERFORMS FEATURE SELECTION FOR THE CNN
# run this file to generate a .json file of important positions

import find_consurf_intervals
intervals = find_consurf_intervals.find_intervals("new_data_consurf_scores.txt")
import reg_2
import json

distances = "3d_distances_augmented.csv"
alignment = "[augmented] muscle alignment.fa"
important_positions = []
# list that will be dumped into json later
# contains positions in form
# [index_1, index_2, etc]

for interval in intervals:
    # intervals are received in form
    # [(start_pos, interval_length), (start_pos2, interval_length2), ...]
    seq_length = interval[1] - interval[0]

    mutual_info = reg_2.run_regression(seq_length, interval[0], distances, alignment, False)
    if mutual_info > 0: # 1.0e-4:
        # could probably modify thershold of importance (originally 1e-4)
        for i in range(interval[0], interval[1]):
            important_positions.append(i)

with open("positions(3).json", "w") as outfile:
    json.dump(important_positions, outfile)
# store results in positions.json


