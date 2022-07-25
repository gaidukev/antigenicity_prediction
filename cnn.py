import tensorflow as tf

from tensorflow.keras import layers, models

import aaindex_vectors
import pandas as pd
import numpy as np
import get_metrics
from sklearn import metrics
import time
from math import floor
from pickle import load, dump
from winsound import Beep as beep



noise = {}
noise_rec = {}
#noise_file = open("noise", "wb")
#dump(noise, noise_file)
#noise_file.close()
def find_noise(anti, ser):
    global noise, noise_rec

    if ser in noise_rec:
        noise_rec[ser].append(anti)
    else:
        noise_rec[ser] = [anti]
        
    if anti+"/A" in noise:
        noise[anti+"/A"] += 1
    else:
        noise[anti+"/A"] = 1
    if ser+"/S" in noise:
        noise[ser+"/S"] += 1
    else:
        noise[ser+"/S"] = 1
    

start_time = time.time()
distances = "3d_distances_augmented.csv" # this contains euclidean distances
alignment = "[augmented] muscle alignment.fa" # .fa alignment file

dist = pd.read_csv(distances)
sera = list(dist.columns)


# construct features for first serum x all antigens
serum = sera[1]
start_pos = 0
len_sequence = 566 # standardized sequence length
og_x_vals = aaindex_vectors.construct_x(serum, alignment, distances, len_sequence)
#print(x_vals.any() != 0)
og_y_vals = aaindex_vectors.construct_y(serum, distances)




#print(sera)
for i in range(2, len(sera)):
    # construct features for rest of sera x all antigens
    serum = sera[i]
    new_x_vals = aaindex_vectors.construct_x(serum, alignment, distances, len_sequence)
    if type(new_x_vals) == bool:
        print(serum)
        continue

    # another way to construct features -- using vstack (instead of concat) (they may be the same?)
    # x_vals = np.vstack((x_vals, new_x_vals))
    # y_vals = np.hstack((y_vals, aaindex_vectors.construct_y(serum, distances)))
    og_x_vals = np.concatenate((og_x_vals, new_x_vals))
    og_y_vals = np.concatenate((og_y_vals, aaindex_vectors.construct_y(serum, distances)))

x_vals, y_vals = og_x_vals, og_y_vals
# get a random permutation of x and y values (together)
p = np.random.permutation(len(x_vals))
x_vals, y_vals = x_vals[p], y_vals[p]


# amount of data to use as training data (picked arbitrarily)
# make it < 2800
train_amt = 2000
num_anti = 113

# split into testing and training data
x_train, x_test = x_vals[:train_amt, :, :], x_vals[train_amt:, :, :]


y_train, y_test = y_vals[:train_amt], y_vals[train_amt:]


print(y_test[0], p[0+train_amt], og_y_vals[p[0+train_amt]], "serium:", chr(ord('@')+floor(p[0+train_amt]/num_anti)+2), "anit:", (p[0+train_amt] % num_anti)+2)
print(dist.iat[(p[0+train_amt] % num_anti), floor(p[0+train_amt]/num_anti)+1])
print(sera[floor(p[0+train_amt]/num_anti)+1])
print(dist.iat[(p[0+train_amt] % num_anti), 0])

for i in range(len(y_train)):
    y_train[i] = y_train[i] * .9
# print(x_vals.shape, y_vals.shape, x_train.shape, x_test.shape, y_train.shape, y_test.shape)
# print((x_vals.any()) != 0)

# build model
model = models.Sequential()

model.add(layers.Conv1D(193, (96), strides=3, activation='relu', input_shape = (x_vals.shape[1], x_vals.shape[2]), padding='same'))
model.add(layers.MaxPooling1D(pool_size=(x_vals.shape[2]), padding='same'))
model.add(layers.Dropout(0.1))

model.add(layers.Conv1D(212, (5), activation='relu', strides=2, padding='same'))
model.add(layers.MaxPooling1D(pool_size=(x_vals.shape[2]), padding='same'))
model.add(layers.Dropout(0.165))

model.add(layers.Conv1D(109, (5), activation='relu', strides=3, padding='same'))
model.add(layers.MaxPooling1D(pool_size=(x_vals.shape[2]), padding='same'))
model.add(layers.Dropout(0.1))

model.add(layers.Flatten())

model.add(layers.Dense(256, activation="sigmoid"))#))
model.add(layers.Dropout(0.241))

model.add(layers.Dense(1))


model.compile(optimizer='adam', loss=tf.keras.losses.MeanSquaredLogarithmicError(), 
            metrics=["MeanSquaredLogarithmicError"])

model.summary()

# fit model with 15 epochs
history = model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=100, verbose=0)


ypred = model.predict(x_test)

ypred2 = np.zeros(ypred.shape)
ytrue = np.zeros(ypred.shape)
fless = 0
ntless = 0
ptless = 0
for val in range(len(ytrue)):
    if ypred[val] > 4:
        ypred2[val] = 1
    if y_test[val] > 4:
        ytrue[val] = 1

        
    if y_test[val] > 4 and ypred[val] < 4:
        diff = y_test[val] - ypred[val]
        row =(p[val+train_amt] % num_anti)
        col = floor(p[val+train_amt]/num_anti)+1
        #print("FALSE NEG: ", y_test[val], ypred[val], diff, end=" ")
        #print(sera[col], dist.iat[row, 0], y_test[val] == dist.iat[row, col])
        find_noise(dist.iat[row, 0], sera[col])
        if diff < 2:
            fless += 1
    elif y_test[val] < 4 and ypred[val] > 4:
        diff = ypred[val] - y_test[val]
        row =(p[val+train_amt] % num_anti)
        col = floor(p[val+train_amt]/num_anti)+1
        #print("FALSE POS: ", y_test[val], ypred[val], diff, end=" ")
        #print(sera[col], dist.iat[row, 0], y_test[val] == dist.iat[row, col])
        find_noise(dist.iat[row, 0], sera[col])
        if diff < 2:
            fless += 1
    elif y_test[val] < 4 and ypred[val] < 4:
        diff = abs(ypred[val] - y_test[val])
        #print("TRUE NEG: ", y_test[val], ypred[val], diff)
        if diff < 2:
            ntless += 1
    elif y_test[val] > 4 and ypred[val] > 4:
        if val < 100:
            diff = abs(ypred[val] - y_test[val])
            print("TRUE POS: ", y_test[val], ypred[val], diff)
        diff = abs(ypred[val] - y_test[val])
        if diff < 2:
            ptless += 1
try:

    print("f1 score: ", metrics.f1_score(ytrue, ypred2))
    print("MCC: ", metrics.matthews_corrcoef(ytrue, ypred2))
except:
    pass

print("ACCURACY: ", get_metrics.get_accuracy(y_test, ypred))
print("SENSITIVITY: ", get_metrics.get_sensitivity(y_test, ypred))
print("SPECIFICITY: ", get_metrics.get_specificity(y_test, ypred))
print("\n\n\nTRAINING SET SIZE: ", len(y_train))
print("TESTING SET SIZE: ", len(y_test))
print("FLESS: ", fless, "nTLESS: ", ntless, "pTLESS: ", ptless)
print(time.time() - start_time)
noise_file = open("noise", "rb")
record = load(noise_file)
noise_file.close()
print(noise)
print(noise_rec)
for k,v in noise.items():
    if v >= 4:
        print(k+":", v)
        if k in record:
            record[k] += v
        else:
            record[k] = v
c = {}
for g in noise_rec.values():
    for v in g:
        if v in c:
                c[v] += 1
        else:
                c[v] = 1
for k,v in c.items():
    if v > 1:
        print(k+":", v)
          
noise_file = open("noise", "wb")
dump(record, noise_file)
noise_file.close()
beep(2500, 1000)
# print("MCC: ", get_metrics.get_mcc(y_test, ypred))

