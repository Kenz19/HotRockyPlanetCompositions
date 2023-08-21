# imports #####################################################################
import numpy as np
import pandas as pd
import re
import sys
from chempy.chemistry import Substance
from chempy.util.parsing import formula_to_composition
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator





# functions ###################################################################
def fixGGchem(df):
    '''
    The GGchem syntax has some 2 character elements as 2 capital letters. 
    This is not the case for all of them. This function fixes this.
    
    This works for the default input of elements but likely will need to be 
    changed if more elements added to config.
    
    df: dataframe where column headers represent minerals.
    '''
    
    # elements GGchem capitalised in default config
    correct_elements = ['Na', "Mg", "Al", "Si", "Cl", "Ca", "Ti", "Mn", "Fe", 
                        "Cr", "Li", "He", 'Zr'] 
    incorrect_elements = ['NA', 'MG', 'AL', 'SI', 'CL', 'CA', 'TI', 'MN', 'FE',
                          'CR', 'LI', 'HE', 'ZR'] 
    
    # replacing incorrect part in column name
    for i in range(len(incorrect_elements)):
        df.columns = df.columns.str.replace(incorrect_elements[i], 
                                            correct_elements[i]) 
        
    return df

###############################################################################

def elfinder(df, element): 
    '''
    Finds all minerals containing a defined element within a pandas dataframe 
    and seperates them into a new pandas dataframe.
    
    df: Pandas dataframe where column names are minerals
    
    element: Name of element of interest in string format e.g for carbon 
    input: 'C'
    '''
    
    # dictionary: key = atomic number (A), value pairs = no. with that A
    ele = Substance.from_formula(element).composition 
    # isolating atomic number of element
    A = list(ele.keys()) 
    
    df = fixGGchem(df) # correcting capitalisations
    
    new_df = pd.DataFrame() # empty dataframe to contain filtered minerals
    
    # finding composition of each mineral, if mineral contains input element 
    # it is sorted into a new table
    for i in range(len(df.columns)):
        #print(df.columns[i])
        minA = list(formula_to_composition(df.columns[i], prefixes = 'n',
                                      suffixes=('(s)', '(l)', '(g)', '(aq)', 
                                                '(cis)', '(trans)', '[l]', 
                                                '[s]')).keys())
        #print(minA)
        
        for j in range(len(minA)):
            if minA[j] == A[0]:
                # will be more efficient with pd.concat (future change)
                new_df[df.columns[i]] = pd.Series(df[df.columns[i]])
                #print(new_df)
                next
                
            else: continue
    
    # defragmentation of the new dataframe, can be removed when concat 
    # implemented above
    finalframe = new_df.copy() 
    
    return finalframe

###############################################################################

def solidfinder(df):
    '''
    Finds solids/ condensates from GGchem when given a dataframe containing 
    only minerals as the column headers
    
    df: df containing information about minerals from a GGchem output. 
    Column headers must be mineral syntax names
    ''' 
    new_df = df.loc[:, df.columns.str.contains('[s]')]
    
    return new_df

###############################################################################

def elcounter_multiple(df, element): 
    '''
    Counts how many times an element is found in a collection of minerals 
    placed as the column headers of a pandas dataframe.
    
    df: Pandas dataframe where column names are minerals
    
    element: Name of element of interest in string format e.g for carbon 
    input: 'C'
    '''
    
    # dictionary: key = atomic number (A), value pairs = no. with that A
    ele = Substance.from_formula(element).composition 
    
    A = list(ele.keys()) # isolating atomic number of element
    
    df = fixGGchem(df) # correcting capitalisations
    
    count = 0
    
    # finding composition of each mineral, if mineral contains input 
    # element it is sorted into a new table
    for i in range(len(df.columns)):
        
        # atomic numbers of elements within mineral
        minA = list(formula_to_composition(df.columns[i], prefixes = 'n',
                                      suffixes=('(s)', '(l)', '(g)', '(aq)', 
                                                '(cis)', '(trans)', '[l]', 
                                                '[s]')).keys())
        
        # number of each element contained within mineral
        minB = list(formula_to_composition(df.columns[i], prefixes = 'n',
                                      suffixes=('(s)', '(l)', '(g)', '(aq)', 
                                                '(cis)', '(trans)', '[l]', 
                                                '[s]')).values())
        
        for j in range(len(minA)):
            if minA[j] == A[0]:
                # adding on number of element in mineral
                count = count + minB[j] 
                
            else: continue
    
    return count

###############################################################################

def elcounter_single(mineral, element): 
    '''
    Counts how many times an element is found in a single mineral.
    
    mineral: Mineral letter syntax, string e.g 'FeS' Iron Sulphide
    
    element: Name of element of interest in string format e.g for carbon 
    input: 'C'
    '''
    
    # dictionary: key = atomic number (A), value pairs = no. with that A
    ele = Substance.from_formula(element).composition 
    
    A = list(ele.keys()) # isolating atomic number of element
    
    count = 0
    
    # atomic numbers of elements within mineral
    minA = list(formula_to_composition(mineral, prefixes = 'n',
                                       suffixes=('(s)', '(l)', '(g)', '(aq)', 
                                                 '(cis)', '(trans)', '[l]', 
                                                 '[s]')).keys())
        
        # number of each element contained within mineral
    minB = list(formula_to_composition(mineral, prefixes = 'n',
                                      suffixes=('(s)', '(l)', '(g)', '(aq)', 
                                                '(cis)', '(trans)', '[l]', 
                                                '[s]')).values())
    
    for j in range(len(minA)):
            if minA[j] == A[0]:
                count = count + minB[j] # adding on number of element in mineral
                
            else: continue

    return count

###############################################################################

def read_dat_file(file_path): # accounts for both numerical and text data
    '''
    Taking in the GGchem Abundances.dat file and converting it into a pandas
    dataframe
    
    file_path = Abundances.dat
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

            # Split the line using the regular expression pattern 
            # (one or more spaces)
            columns = re.split(r'\s+', line.strip())

            # Separate text and numerical columns based on their type
            text_columns = []
            numerical_columns = []
            for column in columns:
                # Check if the column contains numerical data (comma-separated)
                if column is not None and ',' in column:  
                    numerical_columns.extend(column.split(','))
                else:
                    text_columns.append(column)

            data.append(text_columns + numerical_columns)

        # Create a pandas DataFrame
        df = pd.DataFrame(data, columns=headers)

        return df
    
###############################################################################

def replace_line(file_name, line_num, text):
    '''
    Reads desired line in a file and replaces it with any desired text or data
    '''    
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

###############################################################################

def interpolate(xval, df, xcol, ycol):
    '''
    Compute xval as the linear interpolation of xval where df is a dataframe 
    and df.x are the x coordinates, and df.y are the y coordinates. df.x is 
    expected to be sorted.
    
    '''
    return np.interp([xval], df[xcol], df[ycol])

###############################################################################

def CHeditor(abundance_file_path, new_abundance_file_name, input_file_path, CHratio):  
    '''
    abundance_file_path: File path for GGchem Abundances.dat, string
    
    new_abundance_file_name: Name of new abundance file, string
    
    input_file_path: File path to GGchem sim input file, string
    
    CHratio: desired CHratio in new file, float
    '''
    abundances = read_dat_file(abundance_file_path)
    #print(abundances)

    # Converting solar abundances to numeric datatype
    abundances['Solar'] = pd.to_numeric(abundances['Solar'], errors = 'coerce')

    # Grabbing solar abundances
    solar = abundances['Solar'].tolist()

    # Establishing initial carbon and hydrogen parameters
    C = solar[19] # !!! make this general (search for C & H) C = 5, S = 15, Ca = 19
    H = solar[0]

    # Establishing current and wanted C/H ratio
    current_CHratio = C/H
    new_CHratio = CHratio

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
        if i == 19: # 5 refers to carbon, 15 to sulphur, 19 to calcium
            solar[i] = solar[i]*factor
        
        else: solar[i] = solar[i]*reducing_factor

    # needed for converting to astronomical units
    new_H = solar[0]
    print(new_H)
    
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
    throwaway = text.readline() # ignore
    input_elements = text.readline() # second line contains elements in simulation
    l = input_elements.strip().split(' ') 
    #print(l)
    
    # matching elements in input file
    abundance_data = abundance_data[abundance_data['Symbol'].isin(l)]
    
    #abundance_data.to_csv(new_abundance_file_name, sep = ' ', header = False, index = False)
    
    return abundance_data

###############################################################################

def condensate_data(datafile):
    '''
    An adaptation of the GGchem plotting code used to pull out the data 
    aquired from simulation on each mineral.
    
    dataframe: file path for static_conc.dat
    
    output is a dataframe containing all processed information
    '''
    file   = datafile # output file from GGchem
    
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
    #print(NPOINT)
    
    ############### ESTABLISHING PARAMETERS FROM DATA FILE ####################
    
    bar   = 1.E+6                    # 1 bar in dyn/cm2 (setting units)
    Tg    = dat[:,0]                 # T [K]
    nHtot = dat[:,1]                  # n<H> [cm-3]
    lognH = np.log10(nHtot)          
    press = dat[:,2]                 # p [dyn/cm2] pressure
    Tmin  = np.min(Tg)
    Tmax  = np.max(Tg)
    
    #if (Tmax>4*Tmin): Tmax=4*Tmin
    #if (Tmin<Tmax/3): Tmin=Tmax/3
    #Tmax  = 2800
    #Tmin  = 2000
    
    Narg  = len(sys.argv) # not 100% sure what this is referring too
    #print(Narg)
    if (Narg>1): Tmin=float(sys.argv[1])
    if (Narg>2): Tmax=float(sys.argv[2])
    iii   = np.where((Tg>Tmin) & (Tg<Tmax))[0]
    pmin  = np.min(press[iii])/bar # converting pressure to appropriate units
    pmax  = np.max(press[iii])/bar #      ''                ''
    pmin  = pmin*0.9
    pmax  = pmax*1.1
    if (pmax>pmin*5): 
      pmin = pmin/2.0
      pmax = pmax*2.0
    nHmin = np.min(nHtot[iii])
    nHmax = np.max(nHtot[iii])
    nHmin = nHmin*0.9
    nHmax = nHmax*1.1 
    if (nHmax>nHmin*5): 
      nHmin = nHmin/2.0
      nHmax = nHmax*2.0
    sep = 20
    if (Tmax-Tmin<400): sep=10
    if (Tmax-Tmin>1000): sep=50
    if (Tmax-Tmin>2000): sep=100
    #Tmin  = Tmin*0.95
    #Tmax  = Tmax*1.1
    
    ############## Solid particle densities ###############
    
    ntot  = 0.0*nHtot
    for i in range(3,4+NELEM+NMOLE): # electrons, all atoms, ions and cations
      ntot = ntot + 10**dat[:,i]
    lntot = np.log10(ntot)
    
    solids = []
    smean = []
    ymax = -100.0
    for i in range(4+NELEM+NMOLE,4+NELEM+NMOLE+NDUST,1):
      solid = keyword[i]
      solids.append(solid[1:])
      smean.append(np.mean(dat[iii,i])) 
      ind = np.where(keyword == 'n'+solid[1:])[0]
      if (np.size(ind) == 0): continue
      ind = ind[0]
      yy = dat[:,ind]               # log10 nsolid/n<H>
      ymax = np.max([ymax,np.max(yy[iii])])
      ymin = -99
    indices = np.argsort(smean)
    if (ymax>-99):
      count = 0
      for isolid in reversed(indices):
        solid = solids[isolid]
        ind = np.where(keyword == 'n'+solid)[0]
        if (np.size(ind) == 0): continue
        ind = ind[0]
        yy = dat[:,ind]               # log10 nsolid/n<H>
        ymax = np.max([ymax,np.max(yy[iii])])
        if (np.max(yy[iii])>-99): next#print(solid,ind,np.max(yy[iii]))
        if (np.max(yy[iii])>ymin):
          count = count + 1
      #plt.ylim(ymax-8,ymax+0.3)
      #minorLocator = MultipleLocator(sep)
      #ax.xaxis.set_minor_locator(minorLocator)
      #minorLocator = MultipleLocator(1.0)
      #ax.yaxis.set_minor_locator(minorLocator)
      sz = np.min([9,1+120.0/count])
      col = 1
      if (count>20): 
        sz = np.min([9,1+200.0/count])
        col = 2
      if (count>40): 
        sz = np.min([9,1+250.0/count])
        col = 3
      #leg = plt.legend(loc='best',fontsize=sz,ncol=col,fancybox=True)
      #leg.get_frame().set_alpha(0.7)
    
      for iT in range(0,NPOINT):
        iact = 0
        outp = ' '
        for i in range(4+NELEM+NMOLE+NDUST,4+NELEM+NMOLE+2*NDUST,1):
          #print keyword[i],dat[iT,i]
          if (dat[iT,i]>-200): 
            iact=iact+1
            outp=outp+' '+keyword[i][1:]
        #print(Tg[iT],iact,outp)
    #print(solids)
    ################ WHERE ARE ELEMENTS PRESENT? ###################
    
    ellist = ['H','C','O','N','SI','S','NA','CL','CA','TI','K','AL','MG','FE','LI','F','P','NI','MN','CR','ZN','ZR','RB','CU','B','BR','V','SR','W','el'] # code name
    allist = [' ',' ',' ',' ','Si',' ','Na','Cl','Ca','Ti',' ','Al','Mg','Fe','Li',' ',' ','Ni','Mn','Cr','Zn','Zr','Rb','Cu',' ','Br',' ','Sr',' ','+'] # Chemical name
    exlist = [' He ',' Cl CL Ca CA Cr CR Co Cu CU ',' ',' Ne NE Na NA Ni NI ',' ',' Si SI Sr SR ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' Fe FE ',' ',' ',' ',' ',' ',' ',' ',' ',' Br BR ',' ',' ',' ',' ',' ']
    titels = ['hydrogen','carbon','oxygen','nitrogen','silicon','sulphur','sodium','chlorine','calcium','titanium','potassium','aluminum','magnesium','iron','lithium','fluorine','phosphorus','nickel','manganese','chromium','zinc','zirconium','rubidium','copper','boron','bromine','vanadium','strontium','tungston','charge carriers']
    limits = [4,5,4,6,6,5,6,4,7,8,6,6,6,6,7,6,6,6,6,6,6,6,6,6,6,6,6,6,6,5]   
    condensates = indices
    new_df = pd.DataFrame()
    for i in range(0,30):
      #fig,ax = plt.subplots()
      el = ellist[i]
      al = allist[i]
      ex = exlist[i]
      limit = limits[i]
      titel = titels[i]
      #print(titel+" ...")
      nmax = np.float64(-100)
      nmin = np.float64(0)
      mollist = []
      abulist = []
      maxy = 0.0*dat[:,0]
      for mol in range(3,4+NELEM+NMOLE,1):
        molname = keyword[mol]
        #print(keyword)
        #print(molname)
        ind = str.find(molname,el)
        if (ind < 0): 
          ind = str.find(molname,al)
        if (ind < 0 and el=='el'): 
          ind = str.find(molname,'-')
        if (ind >= 0):
          next1 = molname[ind:ind+2]
          next2 = molname[ind-1:ind+1]
          #print keyword[mol],next1,str.find(ex,next1),len(next1)
          if (len(next1)==1 or str.find(ex,next1)==-1 or molname=='SIS'):
            if (next2!='MN' and next2!='ZN'):
              yy = dat[:,mol]                # log10 nmol [cm-3]
              yy = yy - lntot                # log10 nmol/ntot
              nmax = np.max([nmax,np.max(yy[iii])])
              maxy = maxy + 10**yy
              if (molname=='el'): nmin = np.min([nmin,np.min(yy[iii])])
              mollist.append(mol)   
              abulist.append(np.mean(yy))
      for isolid in condensates:
        solid = solids[isolid]
        #print(solid)
        isol = np.where(keyword == 'n'+solid)[0]
        #print(isol)
        if (np.size(isol) == 0): continue
        isol = isol[0]
        search = el
        if (len(el)==2): search=al
        ind = str.find(solid,search)
        found = 0
        while (ind>=0):
          #print solid,ind
          if (len(search)==2): found=1
          if (found==0):  
            if (ind==len(solid)-1): found=2
          if (found==0):  
            next1 = solid[ind+1]
            #print solid,search,next1
            if (next1.isupper() or next1.isdigit() or next1=='['): found=3
          if (found>0): break  
          ind = solid.find(search,ind+1)
          if (ind<0): break
          #print 'try again with rest ',ind,len(solid),solid[ind:]
        if (found>0):
          yy = dat[:,isol]               # log10 nsolid/n<H>
          yy = yy + lognH - lntot        # log10 nsolid/ntot
          nmax = np.max([nmax,np.max(yy[iii])])
          maxy = maxy + 10**yy
          #print found,isol,keyword[isol],np.max(yy[iii])
          mollist.append(isol)   
          abulist.append(np.max(yy[iii]))
      if (nmax==-100): continue
      count = 0
      indices = np.argsort(abulist)
      maxy = np.log10(maxy)
      nmin = np.min([nmin,np.min(maxy[iii])-limit,nmax-12])
      for ind in reversed(indices):
        mol = mollist[ind]
        abu = abulist[ind]
        molname = keyword[mol]
        if (mol<=4+NELEM+NMOLE):
          yy = dat[:,mol]              # log10 nmol [cm-3]
          yy = yy - lntot              # log10 nmol/ntot
        else:
          yy = dat[:,mol]              # log10 ncond/n<H>
          yy = yy + lognH - lntot      # log10 ncond/ntot 
          molname = molname[1:]
          if (str.find(molname,'[l]')<0):
            #print(Tg[ind], yy[ind], molname)
            molname = molname+'[s]'
            #print(molname)
        #print mol,molname,abu,np.max(yy[iii])
        if (np.max(yy[iii]-maxy[iii])>-limit or molname=='el'):      
        #for i in range(len(solid)):
        #if (molname == solid[i]):
                new_df[molname] = pd.Series(yy)
                new_df = new_df.copy()
                #print(new_df)
                #plt.plot(Tg,yy,label=molname)
                count = count + 1
      #plt.title(titel,fontsize=20)
      #plt.xlabel(r'$T\ \mathrm{[K]}$')
      #plt.ylabel(r'$\mathrm{log}_{10}\ n_\mathrm{mol}/n_\mathrm{tot}$')
      #plt.xlim(Tmin,Tmax)
      #plt.ylim(nmin,nmax+1)
      #minorLocator = MultipleLocator(sep)
      #ax.xaxis.set_minor_locator(minorLocator)
      #minorLocator = MultipleLocator(1.0)
      if (nmax-nmin>50): minorLocator = MultipleLocator(2.0)
      if (nmax-nmin>100): minorLocator = MultipleLocator(5.0)
      if (nmax-nmin>200): minorLocator = MultipleLocator(10.0)
      sz = np.min([9,1+120.0/count])
      col = 1
      if (count>20): 
        sz = np.min([9,1+200.0/count])
        col = 2
      if (count>40): 
        sz = np.min([9,1+250.0/count])
        col = 3
        
    print(new_df)
    return new_df

