import selfies as sf
import pandas as pd

def TEST():
    print('working')


def SMILES_to_SELFIES(SMI):
    
    try:
        mol_SELFIES = sf.encoder(SMI)
        return mol_SELFIES
    except:
        
        val = "Invalid SMILES..."+' : '+str(SMI)
        print(val)
        
def SELFIES_to_SMILES(SELFIES):
    try:
    
        smi = sf.decoder(SELFIES)
        return smi
    except:
        val = "Invalid SELFIES..."+' : '+str(SELFIES)
        print(val)




def drop_possible_none_from_list(list_):
    new_list = []
    try:
        for ele in list_:
            if ele != None:
                new_list.append(ele)
    except:
        print('There was some error...')
    return new_list
            

def SELFIES_length(SELFIE_mol_list): #should be a list of SELFIES
    length_selfies = []
    try:
        for SELFIE_ in SELFIE_mol_list:
            length_selfies.append(SELFIE_.count('['))
    except:
        print('Something went wrong, check source code...')
    return max(length_selfies)


def SELFIES_charset(SELFIES_list): #returns charset of SELFIES elements
    #charset = []
    try:
        element = sf.get_alphabet_from_selfies(SELFIES_list)
        #charset.append(element)
    except:
        print('Error, check the code...')
    return list(element)


def prerare_arbitral_one_hot(): #preparing arbitral one hots to encode SELFIES alphabet
    file = pd.read_csv("../SELFIES_coder/UniCode_char.csv")
    Uni_code_char_file = file['character']
    Uni_code_char_file = pd.DataFrame(Uni_code_char_file, columns=['character'])
    Uni_code_char_file = Uni_code_char_file.dropna(axis='rows')

    starchar = '!' #depends on model used
    endchar = 'E' #depends on model used
    unknown_ch = '?'
    to_be_removed = list(starchar+endchar+unknown_ch+'\\'+'\xad'+'¨'+'"')

    Uni_code_char_list = Uni_code_char_file['character'].to_list()

    for element in to_be_removed:

        Uni_code_char_list.remove(element)
    return Uni_code_char_list


def get_accurate_num_of_one_hot(list_with_all_char, charset): #allows to get accurate length of one hots related with SELFIES alphabet
    single_characters = list_with_all_char[0:len(charset)]

    return single_characters


def make_dictionary(single_characters, SELFIE_val): #lists
    charSELFIES_to_SINGLE_CHAR = {}
    for key in single_characters:
        for value in SELFIE_val:
            charSELFIES_to_SINGLE_CHAR[key] = value
            SELFIE_val.remove(value)
            break
    return charSELFIES_to_SINGLE_CHAR #dictionary of one_hot : SELFIES_char

def change_values_to_keys(dictionary):
    val_key = dict(zip(dictionary.values(), dictionary.keys()))
    return val_key #dictionary of SELFIES_char : one_hot



#Prepare useful arrays of SELFIES characters
def make_array(SELFIES_input):
    arrays = []
    for i in range(len(SELFIES_input)):
        array = SELFIES_input[i].split(']')
        modified_array = []
        for element in range(len(array)):
            try:
                modified_array.append(array[element]+']')
            except:
                print("There was some error...")
        modified_array.pop() #last element contains only ] bracket so it should be removed
        arrays.append(modified_array)
    return arrays


#translate SELFIES characters into single characters given above
#allows to get encoded SELFIES molecules
def translate_SELFIES_array_into_one_hot(array_, SINGLE_CHAR_SELFIE):
    
    whole_list = []
    for i in range(len(array_)):
        new_translated_list = []
        for SELF in array_[i]:
            if SELF in list(SINGLE_CHAR_SELFIE.keys()):
                new_alphabet = SINGLE_CHAR_SELFIE[SELF]
                new_translated_list.append(new_alphabet)
        new_string_as_encoded_SELFIES = ''.join(new_translated_list)
        whole_list.append(new_string_as_encoded_SELFIES)
    return whole_list

#allows to decode encoded string
def convert_back_to_SEFLIES(translated, SELFIE_CHAR_SINGLE):
    new_translated_list_SELF = []
    for char in translated:
        if char in list(SELFIE_CHAR_SINGLE.keys()):
            new_alphabet_ = SELFIE_CHAR_SINGLE[char]
            new_translated_list_SELF.append(new_alphabet_)
            new_string_as_encoded_SELFIES = ''.join(new_translated_list_SELF)
    return new_string_as_encoded_SELFIES




#this function allows to obtain translated SELFIES and compare the results of decoding - 
#if initial data is equal to deciphered data everything is correct
def SELFIES_to_one_hot_and_verify(SELFIES_list_, SINGLE_CHAR_SELFIE, SELFIE_CHAR_SINGLE):
    arrays_test = make_array(SELFIES_list_)
    translation = translate_SELFIES_array_into_one_hot(arrays_test, SINGLE_CHAR_SELFIE)
    decoding_m = []
    for i in range(len(SELFIES_list_)):
        
        decoding_m_ = convert_back_to_SEFLIES(translation[i], SELFIE_CHAR_SINGLE)
        decoding_m.append(decoding_m_)
    if decoding_m == SELFIES_list_:
        print("All works correct, encoding leads to the same data during decoding...")
    else:
        print("Something went wrong...")
    
    return translation


## Whole procedure at once
def get_encoded_SELFIES(SMILES):

    SELFIES = []
    for SMI in SMILES:
        try:
            new_data = SMILES_to_SELFIES(SMI)
            SELFIES.append(new_data)
        except:
            print('Invalid SMILES')
    new_list_ = drop_possible_none_from_list(SELFIES)

    charset_ = SELFIES_charset(new_list_)

    possible_one_hots = prerare_arbitral_one_hot()

    single_hot_enc = get_accurate_num_of_one_hot(possible_one_hots, charset_)

    dictionary_n = make_dictionary(single_hot_enc, charset_)

    dictionary_val_keys = change_values_to_keys(dictionary_n)

    verification = SELFIES_to_one_hot_and_verify(new_list_, dictionary_val_keys, dictionary_n)


    max_length = SELFIES_length(new_list_)
    return verification, dictionary_n, dictionary_val_keys, max_length
