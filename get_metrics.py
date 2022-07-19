# > 4 = positive
# < 4 = negative

import math
evaled = False
true_pos = 0
true_neg = 0
false_pos = 0
false_neg = 0
four = 4

def set_four(new_four):
    global four
    four = new_four

# evaluate the number of false/true positives/negatives
def evaluate(y_true, y_pred):
    global true_pos, true_neg, false_pos, false_neg, evaled
    if not evaled:
        for i in range(len(y_true)):
            y_t = y_true[i]
            y_p = y_pred[i]
            if y_t > four and y_p > four:
                true_pos += 1
            elif y_t > four and y_p < four:
                false_neg += 1
            elif y_t < four and y_p > four:
                false_pos += 1
            elif y_t < four and y_p < four:
                true_neg += 1
        evaled = True

# find accuracy metrics as per the paper (formulas are also from the paper)
def get_accuracy(y_true, y_pred):
    global true_pos, true_neg, false_pos, false_neg, evaled
    evaluate(y_true, y_pred)

    num = true_pos + true_neg
    denom = true_pos + true_neg + false_pos + false_neg
    return num / denom
    

def get_sensitivity(y_true, y_pred):
    global true_pos, true_neg, false_pos, false_neg, evaled
    evaluate(y_true, y_pred)

    try:
        return true_pos / (true_pos + false_neg)
    except:
        return 0

def get_specificity(y_true, y_pred):
    global true_pos, true_neg, false_pos, false_neg, evaled
    evaluate(y_true, y_pred)

    print("FALSE POS: ", false_pos, "FALSE NEG: ", false_neg, "TRUE POS: ", true_pos, "TRUE NEG: ", true_neg)
    return true_neg / (true_neg + false_pos)

def get_mcc(y_true, y_pred):
    global true_pos, true_neg, false_pos, false_neg, evaled
    evaluate(y_true, y_pred)

    num = true_pos * true_neg - false_pos * false_neg
    denom = (true_pos + false_neg) * (true_pos + false_neg) * (true_neg + false_pos) * (true_neg + false_neg)
    denom = math.sqrt(denom)

    print("FALSE POS: ", false_pos, "FALSE NEG: ", false_neg, "TRUE POS: ", true_pos, "TRUE NEG: ", true_neg)
    try:
        return num / denom
    except:
        return 0