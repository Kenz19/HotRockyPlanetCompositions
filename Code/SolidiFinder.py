"""
This file will run the entire command line that I have been working with and
make the process more efficient

Firstly the Abundances.dat file will be edited to a desired C/H ratio.

GGchem will then be ran with the new abundances file.

A plot of the mineral abundances against temperature is then produced. T50 is 
also calculated in this step and appended to an external text file

hopefully will loop through a bunch of C/H and get a good T5O vs C/H plot
"""

###############################################################################
import Functions as func
import CHratio as rat
import numpy as np
import shutil
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
###############################################################################

# getting CH value inputs from user
LowestCH = float(input('Lowest C/H limit: '))
HighestCH = float(input('Highest C/H limit: '))
step = float(input('Step between each limit: '))
numerator = input('Numerator element: ')
denominator = input('Denominator element: ')
element = input('Enter element of interest (Symbol, case sensitive): ')

# creating range of CH values
ratios = np.arange(LowestCH, HighestCH + step, step).tolist()
#print(ratios)

# creating file path for storing abundance files for later use if needed
if not os.path.exists('Abundance_files_Ca'):
    os.mkdir('Abundance_files_Ca')

# creating file path for storing Static_conc files for later use if needed
if not os.path.exists('Static_conc_files'):
    os.mkdir('Static_conc_files')

# edits each abundance file and sets it as the file for use in GGchem input file (default.in)
for i in range(len(ratios)):
    abundances = rat.CHedit(ratios[i], numerator, denominator, r'data/Abundances.dat', r'input/default.in')
    abundances.to_csv('Abundance_files_Ca/Abun_' + str(ratios[i]) + '.in', 
                      sep = ' ', header = False, index = False)
    
    # telling GGchem we want a custom abundance file (setting option to 0)
    #func.replace_line('input/default.in', 9, '0\n') # line num should remain the same each time, implement this correctly at later date
    func.replace_line('input/default.in', 10, 'Abundance_files_Ca/Abun_' + str(ratios[i]) + '.in\n')
        
    # running GGchem with new abundance file
    os.system('./ggchem input/default.in')
    
    #shutil.copy('/home/kg543/GGchem/Static_Conc.dat','/home/kg543/GGchem/Static_Conc_files')#/Static_conc_CH' + str(ratios[i]))
    #os.rename('/home/kg543/GGchem/Static_Conc_file/Static_Conc.dat', '/homekg543/GGchem/Static_conc_files/static_conc_CH' + str(ratios[i])) 
    #with open('Static_Conc.dat', 'r') as file:
     #   text = file.read()
      #  with open('c:\\home\\kg543\\GGchem\\Static_conc_files\\Static_Conc' + str(ratios[i])) as file:
       #     file.write(text)
    	
    # running plotting code and getting T50 & CH info	
    new_df = func.condensate_data('Static_Conc.dat')
    
    abundances = open('/home/kg543/GGchem/Abundance_files_Ca/Abun_' + str(ratios[i]) + '.in', 'r')

    # obtaining hydrogen & carbon abundances
    hydrogen_info = abundances.readline().split(' ')
    H = float(hydrogen_info[2]) # change this to be appropriate to what ever element you need

    # throwaway helium and lithium
    He = abundances.readline()
    Li = abundances.readline()
 
    # carbon abundance
    carbon_info = abundances.readline()#.split(' ')
    C_2 = float(carbon_info[2]) # change this to be appropriate to the ratio you need
 
    N = abundances.readline()
    O = abundances.readline()
    F = abundances.readline()
    Na = abundances.readline()
    Mg = abundances.readline()
    Al = abundances.readline()
    Si = abundances.readline()
    P = abundances.readline() 

    # SULPHUR INFORMATION
    sulphur_info = abundances.readline()#.split(' ')
    #C = float(sulphur_info[2])
    Cl = abundances.readline()
    K = abundances.readline() 
    calcium_info = abundances.readline().split(' ')
    C = float(calcium_info[2])

# Obtaining solid data & plotting #############################################
    alldata = new_df 
    alldata = alldata.drop('el', axis=1) # dropping ion column
    keywords = list(alldata.columns) # mineral names

# Correcting the naming errors in GGchem 
    minerals = func.fixGGchem(alldata)

#cfind minerals containing input element
    minerals = func.elfinder(minerals, element)
    solids = func.solidfinder(minerals) # isolating solids

# file information
    file   = 'Static_Conc.dat' # output file from GGchem

    data   = open(file) # opening file
    dummy  = data.readline() # blank line
    dimens = data.readline() # file info
    dimens = np.array(dimens.split()) # splitting file information up
    #print(dimens)
    #sys.exit()

    NELEM  = int(dimens[0]) # total number of elements input
    NMOLE  = int(dimens[1]) # total number of moles present in composition
    NDUST  = int(dimens[2]) # total count of dust (solid) particles
    NPOINT = int(dimens[3]) # number of incramental steps between temp range
    header = data.readline() # column headers
    #print(header)
    data.close()
    dat = np.loadtxt(file,skiprows=3) # all data
    #print(dat)
    keyword = np.array(header.split()) # headers
    #print(keyword)
    #print(keyword)
    NPOINT = len(dat[0:])  



    data = pd.DataFrame(dat)
    data.columns = keyword

    f_column = data["Tg"] # obtaining temperature column from static_conc

    solids = pd.concat([solids,f_column], axis = 1)
    headers = list(solids.columns)

# getting T50 #################################################################

# solids abundance

    solids_edit = solids.drop('Tg', axis = 1)
    solids_edit = 10**solids_edit
    title_solids = list(solids_edit.columns)

    for i in range(len(title_solids)):
        counts = func.elcounter_single(title_solids[i], element)
        solids_edit[title_solids[i]] = solids_edit[title_solids[i]]*counts

    solid_abundance = solids_edit.sum(axis = 1)
    minerals = 10**minerals

    title_minerals = list(minerals.columns)



# accounting for multiple abundancesS

    for i in range(len(title_minerals)):
        count = func.elcounter_single(title_minerals[i], element)
        minerals[title_minerals[i]] = minerals[title_minerals[i]]*count
        total_abundance = minerals.sum(axis = 1)    

# ratio
    solid_percentage = list(solid_abundance/total_abundance)
    solid_percentage = pd.DataFrame(solid_percentage, columns=['solid%'])
    temp = data["Tg"]
    solid_percentage = pd.concat([solid_percentage,temp], axis = 1)

#find row with closest value to 101 in points column
    df_closest = solid_percentage.iloc[(solid_percentage['solid%']-0.5).abs().argsort()[:1]] # T50 found with interpolation

#view results
    T50 = func.interpolate(0.5, solid_percentage, 'solid%', 'Tg')
    #print(T50[0])

# # Pulling C/H from edited abundances.dat ######################################



# # Reading in abundances.dat file

# abundances = func.read_dat_file('Abundances.dat')

# #print(abundances)



# # Converting to workable datatype (float)

# abundances['Solar'] = pd.to_numeric(abundances['Solar'], errors = 'coerce')



# # Pulling solar abundances out

# solar = abundances['Solar'].tolist()

# #print(sum(solar))



# # C and H values

# C = solar[5]

# H = solar[0]



#datatypes = abundances.dtypes



# Putting output in text file #################################################



# saving T50 and C/H ratio to data file

    datafile = open('T50.txt', 'a')
    datafile.write('\n' + str(T50[0]) + ',' + str(C/H) + ','  + str(C_2/H))
    datafile.close()

    solids.plot(x="Tg", y = headers[0:-1])
    plt.ylim(-17,-4)
    #plt.xlim(100, 2000)
    plt.xlabel(r'$T\ \mathrm{[K]}$')
    plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$')

# Save abundances temp plot with appropriate title
    plt.title(element + ' Solids' + ',' + numerator + '/' + denominator + '= ' + str(round(C/H,5)) + ', T50 = ' + str(round(T50[0], 5)) + 'K')

# #plt.savefig('Plots/graph.png')  
    plt.savefig('Plots/' + numerator + '/'+ denominator + 'ratio = ' + str(round(C/H,5)) + '.png')




