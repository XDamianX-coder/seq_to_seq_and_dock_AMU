#!/usr/bin/env python
# coding: utf-8

# General Imports
import os
import pandas as pd
import numpy as np
from rdkit import Chem
import time
import sys
from pathlib import Path
SELFIES_coder_path = Path("../SELFIES_coder")
sys.path.append(SELFIES_coder_path.as_posix())
import SELFIES_coder as SELFIES_CODER
import selfies as sf

SELFIES_CODER.TEST()

# Inital time
start = time.time()


# Know where user actually is
directory_path = os.getcwd()
print("My current directory is : " + directory_path)
folder_name = os.path.basename(directory_path)
print("My directory name is : " + folder_name)

#load file with smiles codes to train and test (in this case one column is enough)
smifile = str(directory_path)+str('/data/zinc20_selected_to_create_model_processed.parquet')
#zinc20_FK_DC_BB_JA_HE_GA_KG_IC_CB_HJ_processed.parquet
data = pd.read_parquet(smifile)
#data = data

from sklearn.model_selection import train_test_split
smiles_train, smiles_test = train_test_split(data["SMILES_canonical"], test_size=0.1, train_size=0.9, random_state=42)
print(smiles_train.shape)
print(smiles_test.shape)


print(data[:10])


data = SELFIES_CODER.get_encoded_SELFIES(data['SMILES_canonical'].to_list())

print(data[2])


# characters that are used in given SMILES dataset along with initial and stopping characters
charset = set("".join(list(data[0]))+"!E")
char_to_int = dict((c,i) for i,c in enumerate(charset))
int_to_char = dict((i,c) for i,c in enumerate(charset))
embed = data[3] + 15 #20
print(str(charset))
print(len(charset), embed)

shap = smiles_train.shape[0]
shap1 = smiles_test.shape[0]


import json
json = json.dumps(data[2])
f = open("SELFIES_to_mol_seq.json","w")
f.write(json)
f.close()

import json 
json = json.dumps(data[1])
f = open("mol_seq_to_SELFIES.json","w")
f.write(json)
f.close()

import json

json = json.dumps(char_to_int)

# open file for writing, "w" 
f = open("mol_seq_to_int.json","w")
f.write(json)
f.close()

import json
json1 = json.dumps(int_to_char)

# open file for writing, "w" 
f = open("int_to_mol_seq.json","w")
f.write(json1)
f.close()


#vectorization of molecular sequence
def vectorize(mol_seqs, shap):
        one_hot =  np.zeros((shap, embed , len(charset)),dtype=np.int8)
        for i,mol_seq in enumerate(mol_seqs):
            #encode the startchar
            one_hot[i,0,char_to_int["!"]] = 1
            #encode the rest of the chars
            for j,c in enumerate(mol_seq):
                one_hot[i,j+1,char_to_int[c]] = 1
            #Encode endchar
            one_hot[i,len(mol_seq)+1:,char_to_int["E"]] = 1
        #Return two, one for input and the other for output
        return one_hot[:,0:-1,:], one_hot[:,1:,:]
mol_seq_train_ = data[0][:int(shap)]
mol_seq_test_ = data[0][int(shap):int(shap+shap1)]
X_train, Y_train = vectorize(mol_seq_train_, shap)
X_test,Y_test = vectorize(mol_seq_test_, shap1)
print(smiles_train.iloc[0])


#going back from vectorized form to redable string
string_test = "".join([int_to_char[idx] for idx in np.argmax(X_train[0,:,:], axis=1)])
print(string_test)

decoded = SELFIES_CODER.convert_back_to_SEFLIES(string_test, data[1])
print(decoded)

decoded_two = SELFIES_CODER.convert_back_to_SEFLIES(mol_seq_train_[0], data[1])
print(decoded_two)

print("Correct encoding-decoding: "+str(decoded == decoded_two))

mol_seq_train_[0]
#Import Keras objects
#LSTM cells special kinds of neural network units that are designed to keep an internal state for longer iterations
from keras.models import Model
from keras.layers import Input
from keras.layers import LSTM
from keras.layers import Dense
from keras.layers import Concatenate
from keras import regularizers
input_shape = X_train.shape[1:]
output_dim = Y_train.shape[-1]
latent_dim = 128
lstm_dim = 128



#encoder-decoder architecture
unroll = False
encoder_inputs = Input(shape=input_shape)
encoder = LSTM(lstm_dim, return_state=True,
                unroll=unroll)
encoder_outputs, state_h, state_c = encoder(encoder_inputs)
states = Concatenate(axis=-1)([state_h, state_c])
neck = Dense(latent_dim, activation="relu")
neck_outputs = neck(states)



decode_h = Dense(lstm_dim, activation="relu")
decode_c = Dense(lstm_dim, activation="relu")
state_h_decoded =  decode_h(neck_outputs)
state_c_decoded =  decode_c(neck_outputs)
encoder_states = [state_h_decoded, state_c_decoded]
decoder_inputs = Input(shape=input_shape)
decoder_lstm = LSTM(lstm_dim,
                    return_sequences=True,
                    unroll=unroll
                   )
decoder_outputs = decoder_lstm(decoder_inputs, initial_state=encoder_states)
decoder_dense = Dense(output_dim, activation='softmax')
decoder_outputs = decoder_dense(decoder_outputs)
#Define the model, that inputs the training vector for two places, and predicts one character ahead of the input
model = Model([encoder_inputs, decoder_inputs], decoder_outputs)
print(model.summary())


from keras.callbacks import History, ReduceLROnPlateau
h = History()
rlr = ReduceLROnPlateau(monitor='val_loss', factor=0.5,patience=15, min_lr=0.000001, verbose=1, min_delta=1e-5) #epsilon=min_delta



from tensorflow.keras.optimizers import RMSprop, Adam
opt=Adam(learning_rate=0.005) #Default 0.001
model.compile(optimizer=opt, loss='categorical_crossentropy')



start_MODEL = time.time()
model.fit([X_train,X_train],Y_train, epochs=200, batch_size=256, shuffle=True, callbacks=[h, rlr], validation_data=[[X_test,X_test],Y_test]) #100
end_MODEL = time.time()
print("The time of execution of above program is :", end_MODEL-start_MODEL)


import pickle
f = open("Blog_history_ZINC_SELFIES.pickle","wb")
pickle.dump(h.history, f)


smiles_to_latent_model = Model(encoder_inputs, neck_outputs)

smiles_to_latent_model.save("mol_seq2lat_ZINC.h5")


latent_input = Input(shape=(latent_dim,))
#reuse_layers
state_h_decoded_2 =  decode_h(latent_input)
state_c_decoded_2 =  decode_c(latent_input)
latent_to_states_model = Model(latent_input, [state_h_decoded_2, state_c_decoded_2])
latent_to_states_model.save("lat2state_ZINC.h5")


#Last one is special, we need to change it to stateful, and change the input shape
inf_decoder_inputs = Input(batch_shape=(1, 1, input_shape[1]))
inf_decoder_lstm = LSTM(lstm_dim,
                    return_sequences=True,
                    unroll=unroll,
                    stateful=True
                   )
inf_decoder_outputs = inf_decoder_lstm(inf_decoder_inputs)
inf_decoder_dense = Dense(output_dim, activation='softmax')
inf_decoder_outputs = inf_decoder_dense(inf_decoder_outputs)
sample_model = Model(inf_decoder_inputs, inf_decoder_outputs)


#Transfer Weights
for i in range(1,3):
    sample_model.layers[i].set_weights(model.layers[i+6].get_weights())
sample_model.save("samplemodel_ZINC.h5")



print(sample_model.summary())


end = time.time()


print("The time of execution of above program is :", end-start)

