
import pandas as pd

from sklearn.model_selection import train_test_split

def newFeat(data):
    '''
    Function to calculate ZZ and AC features of the nanopore and then clean the dataset by dropping data 
    corresponding to similar nanopores 
    '''
    
    ZZ_feature = data['mirror_30_clock_vertical'] + data['mirror_90_clock_vertical'] + data['mirror_150_clock_vertical']
    AC_feature = data['mirror_vertical'] + data['mirror_60_clock_vertical'] + data['mirror_120_clock_vertical']

    data['ZZ_features'] = ZZ_feature
    data['AC_features'] = AC_feature
    
    data = data.drop(columns=['pore_id'])
    data = data.drop(columns=['mirror_30_clock_vertical', 'mirror_90_clock_vertical', 'mirror_150_clock_vertical', 'mirror_vertical', 'mirror_60_clock_vertical', 'mirror_120_clock_vertical'])
        
    data = data.drop_duplicates(subset=['unique_id']).reset_index(drop=True)
    
    data = data.drop(columns=['unique_id', 'time'])
    
    return data

# Load all the CSV files 
pore4 = pd.read_csv('input/physicalFeatures/Physical_features-4.csv').reset_index(drop=True)
pore5 = pd.read_csv('input/physicalFeatures/Physical_features-5.csv').reset_index(drop=True)
pore6 = pd.read_csv('input/physicalFeatures/Physical_features-6.csv').reset_index(drop=True)
pore7 = pd.read_csv('input/physicalFeatures/Physical_features-7.csv').reset_index(drop=True)
pore8 = pd.read_csv('input/physicalFeatures/Physical_features-8.csv').reset_index(drop=True)
pore9 = pd.read_csv('input/physicalFeatures/Physical_features-9.csv').reset_index(drop=True)
pore10 = pd.read_csv('input/physicalFeatures/Physical_features-10.csv').reset_index(drop=True)
pore11 = pd.read_csv('input/physicalFeatures/Physical_features-11.csv').reset_index(drop=True)
pore12 = pd.read_csv('input/physicalFeatures/Physical_features-12.csv').reset_index(drop=True)
pore13 = pd.read_csv('input/physicalFeatures/Physical_features-13.csv').reset_index(drop=True)
pore14 = pd.read_csv('input/physicalFeatures/Physical_features-14.csv').reset_index(drop=True)
pore15 = pd.read_csv('input/physicalFeatures/Physical_features-15.csv').reset_index(drop=True)
pore16 = pd.read_csv('input/physicalFeatures/Physical_features-16.csv').reset_index(drop=True)
pore17 = pd.read_csv('input/physicalFeatures/Physical_features-17.csv').reset_index(drop=True)
pore18 = pd.read_csv('input/physicalFeatures/Physical_features-18.csv').reset_index(drop=True)
pore19 = pd.read_csv('input/physicalFeatures/Physical_features-19.csv').reset_index(drop=True)
pore20 = pd.read_csv('input/physicalFeatures/Physical_features-20.csv').reset_index(drop=True)
pore21 = pd.read_csv('input/physicalFeatures/Physical_features-21.csv').reset_index(drop=True)
pore22 = pd.read_csv('input/physicalFeatures/Physical_features-22.csv').reset_index(drop=True)

# Calculate the ZZ and AC features for each nanopore and clean the data
pore4 = newFeat(pore4)
pore5 = newFeat(pore5)
pore6 = newFeat(pore6)
pore7 = newFeat(pore7)
pore8 = newFeat(pore8)
pore9 = newFeat(pore9)
pore10 = newFeat(pore10)
pore11 = newFeat(pore11)
pore12 = newFeat(pore12)
pore13 = newFeat(pore13)
pore14 = newFeat(pore14)
pore15 = newFeat(pore15)
pore16 = newFeat(pore16)
pore17 = newFeat(pore17)
pore18 = newFeat(pore18)
pore19 = newFeat(pore19)
pore20 = newFeat(pore20)
pore21 = newFeat(pore21)
pore22 = newFeat(pore22)

# Concatenate all the above dataframes into a single dataframe
allData = pd.concat([
    pore4, pore5, pore6, pore7, pore8, pore9, pore10, pore11, pore12, pore13, 
    pore14, pore15, pore16, pore17, pore18, pore19, pore20, pore21, pore22
]).reset_index(drop=True)

# Split the complete dataset into training and test sets
data, test = train_test_split(allData, shuffle=True, random_state=42)
data, test = data.reset_index(drop=True), test.reset_index(drop=True)

# Save them in the "input" directory
data.to_csv('input/train.csv', index=False)
test.to_csv('input/test.csv', index=False)

