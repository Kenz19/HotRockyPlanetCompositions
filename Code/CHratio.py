import pandas as pd
import re
import numpy as np

def read_dat_file(file_path): # accounts for both numerical and text data
    '''
    Imports the Abundaces.dat file and converts it into a pandas data frame
    provided that headers have been inserted into the first three columns that 
    do not have headers.
    
    file_path: abundances.dat file path
    '''
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = []
        line_counter = 0

        # Extract header names from the 5th line
        header_line = lines[4]
        headers = re.split(r'\s+', header_line.strip())

        for line in lines:
            line_counter += 1

            if line_counter <= 5:  # Skip the first 4 lines
                continue

            # Split the line using the regular expression pattern (one or more spaces)
            columns = re.split(r'\s+', line.strip())

            # Separate text and numerical columns based on their type
            text_columns = []
            numerical_columns = []
            for column in columns:
                if column is not None and ',' in column:  # Check if the column contains numerical data (comma-separated)
                    numerical_columns.extend(column.split(','))
                else:
                    text_columns.append(column)

            data.append(text_columns + numerical_columns)

        # Create a pandas DataFrame
        df = pd.DataFrame(data, columns=headers)

        return df

def CHedit(desired_ratio, numerator, denominator, abundance_file_path, input_file_path):
    
    # reading in the abundances file
    abundances = read_dat_file(abundance_file_path)
    
    # Converting solar abundances to numeric datatype
    abundances['Solar'] = pd.to_numeric(abundances['Solar'], errors = 'coerce')

    # Grabbing solar abundances
    solar = abundances['Solar'].tolist()
    
    # element symbols in abundance file, use to find the numerator and denominator
    elements = list(abundances['symbol'])
    
    # Establishing initial carbon and hydrogen parameters
    C = solar[elements.index(numerator)] # abundance of numerator element
    xz = solar[5]
    print(xz)
    H = solar[elements.index(denominator)]
    print(H)
    print('Old C/H ratio = ' + str(xz/H))
    print('Old Ca/H ratio = '  + str(C/H))
    
    
    # Establishing current and wanted C/H ratio
    current_CHratio = C/H
    new_CHratio = desired_ratio
    #new_CHratio = float(input('Desired CH ratio: '))
    
    # Calculating factor used to acheive the new CH ratio. C is multipled by this
    x = (1-C)
    factor_numerator = new_CHratio*x
    factor_demoninator = current_CHratio+(new_CHratio*C)
    factor = factor_numerator/factor_demoninator
    
    # Other elements are multiplied by this reducing factor to keep abundance sum 
    # equal to one
    reducing_factor = 1-(C*(factor+1))
    
    # Editing each abundance based on factor. Increasing C and reducing else
    for i in range(len(solar)):
        if i == elements.index(numerator): # 5 refers to carbon
            solar[i] = solar[i]*factor
            
        else: solar[i] = solar[i]*reducing_factor
    
    # needed for converting to astronomical units
    new_H = solar[elements.index(denominator)]
    new_C = solar[5]
    #print(new_C, new_H)
    new_Ca = solar[19]
    print('New Ca/H ratio = ' + str(new_Ca/new_H))
    print('New C/H ratio = ' + str(new_C/new_H))
    
    #print(new_H)
    
    placehold = []
    
    # converting to solar units
    for i in range(len(solar)):
        placehold.append(12 + np.log10(solar[i]/new_H))
        
    # establishing dataframe to be converted to file
    symbol_list = abundances['symbol'].tolist()
    abundance_data = pd.DataFrame(list(zip(symbol_list, placehold, solar)), columns =['Symbol', 'value (astronomical)', 'value (abundance)'])
    
    # now want to compare this dataframe to the input elements in the input file, removing any elements not present from 
    # abundances file
    
    # read second line of the custom.in file
    
    text = open(input_file_path, 'r')
    throwaway = text.readline()
    input_elements = text.readline() # second line contains elements in simulation
    l = input_elements.strip().split(' ')
    #print(l)
    
    abundance_data = abundance_data[abundance_data['Symbol'].isin(l)]
    #print(abundance_data)
    #abundance_data.to_csv('CustomComposition.dat', sep = ' ', header = False, index = False)
    return abundance_data

#CHedit(1, 'Ca', 'H', 'Abundances.dat')
