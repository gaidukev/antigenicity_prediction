
import tensorflow as tf

from tensorflow.keras import layers, models

import aaindex_vectors
import pandas as pd
import numpy as np
import get_metrics
from sklearn import metrics

distances = "3d_distances_augmented.csv" # this contains euclidean distances
alignment = "[augmented] muscle alignment.fa" # .fa alignment file

dist = pd.read_csv(distances)
sera = list(dist.columns)


# construct features for first serum x all antigens
serum = sera[1]
start_pos = 0
len_sequence = 566 # standardized sequence length
x_vals = aaindex_vectors.construct_x(serum, alignment, distances, len_sequence)
#print(x_vals.any() != 0)
y_vals = aaindex_vectors.construct_y(serum, distances)




#print(sera)
for i in range(2, len(sera)):
    # construct features for rest of sera x all antigens
    serum = sera[i]
    new_x_vals = aaindex_vectors.construct_x(serum, alignment, distances, len_sequence)

    # another way to construct features -- using vstack (instead of concat) (they may be the same?)
    # x_vals = np.vstack((x_vals, new_x_vals))
    # y_vals = np.hstack((y_vals, aaindex_vectors.construct_y(serum, distances)))
    x_vals = np.concatenate((x_vals, new_x_vals))
    y_vals = np.concatenate((y_vals, aaindex_vectors.construct_y(serum, distances)))


# get a random permutation of x and y values (together)
p = np.random.permutation(len(x_vals))
x_vals, y_vals = x_vals[p], y_vals[p]


# amount of data to use as training data (picked arbitrarily)
# make it < 2800
train_amt = 2500

# split into testing and training data
x_train, x_test = x_vals[:train_amt, :, :], x_vals[train_amt:, :, :]


y_train, y_test = y_vals[:train_amt], y_vals[train_amt:]
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
history = model.fit(x_train, y_train, validation_data=(x_test, y_test), epochs=50)


ypred = model.predict(x_test)

ypred2 = np.zeros(ypred.shape)
ytrue = np.zeros(ypred.shape)
for val in range(len(ytrue)):
    if ypred[val] > 4:
        ypred2[val] = 1
    if y_test[val] > 4:
        ytrue[val] = 1

try:

    print("f1 score: ", metrics.f1_score(ytrue, ypred2))
    print("MCC: ", metrics.matthews_corrcoef(ytrue, ypred2))
except:
    pass

print("ACCURACY: ", get_metrics.get_accuracy(y_test, ypred))
print("SENSITIVITY: ", get_metrics.get_sensitivity(y_test, ypred))
print("SPECIFICITY: ", get_metrics.get_specificity(y_test, ypred))
# print("MCC: ", get_metrics.get_mcc(y_test, ypred))

