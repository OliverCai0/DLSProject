from tdc.multi_pred import DTI
import math
import pandas as pd
import os

data = DTI(name = 'DAVIS').get_data()

# 'we have selected only proteins with a length between 264 and 1400 residues, which corresponds to 95.7% of the information present in the dataset'
ORIGINAL_NUMBER_OF_ROWS = data.shape[0]
data = data[data['Target'].map(lambda x: len(x) <= 1400 and len(x) >= 264) & data['Drug'].map(lambda x: len(x) >= 38 and len(x) <= 72)]
print('Percentage of information after filtering for SMILE and protein length', data.shape[0] / ORIGINAL_NUMBER_OF_ROWS)

# normalize the binding affinities
data['Y'] = data['Y'].map(lambda x: -math.log10(x))

# encode SMILE characters
smile_char_dictionary = {"#": 1, "]": 2, "1": 3, "3": 4, "-": 5, "2": 6, "(": 7, "[": 8, "S": 9, "5": 10, "I": 11, "=": 12, "O": 13, "C": 14, "6": 15, "N": 16, "F": 17, "4": 18, "n": 19, ")": 20, "o": 21, "H": 22, "c": 23, "s": 24, "Br": 25, "Cl": 26}

# check that there are exactly 26 categories, including the empty_char
assert len(smile_char_dictionary) == 26

# pad all SMILEs less than 72
data['Drug'] = data['Drug'].map(lambda x: x + '0' * (72 - len(x)) if len(x) < 72 else x)

# load subwords
subwords = pd.read_csv(os.path.join(os.getcwd(), "data", 'subword_units_map_uniprot'))



