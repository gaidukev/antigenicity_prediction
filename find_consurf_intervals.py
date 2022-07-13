# find intervals where consurf scores < 0.99 given a consurf output txt file
# called from imp_pos_to_json.py so no need to call this file on its own
# usually consurf scores are taken from consurf_scores.txt


def find_intervals(filename):


    with open(filename) as consurf_file:
        intervals  = []
        i = 0
        for line in consurf_file:
            try:
                consurf_score =  float(line.split()[2])
                if consurf_score < 0.99:
                    if len(intervals) == 0:
                        intervals.append([i, i+1])
                    else:

                        if intervals[-1:][0][1] == i:
                            intervals[-1:][0][1] = i+1
                        else:
                            intervals.append([i, i+1])
                i+=1
                    


            except Exception as e:
                print("find_consurf_intervals.py error: ", e)
        return intervals

