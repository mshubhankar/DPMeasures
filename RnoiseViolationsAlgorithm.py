import random
import pandas as pd
import numpy
import string
import datetime
from datetime import date
from collections import defaultdict
import math

def col_in_constraints(constraintSet,df):
    allColomns = []
    for col in df.columns:
        if col in constraintSet:
            allColomns.append(col)
    return allColomns

def harmonic_sum(n,beta):
    sum = 0.0
    for i in range(1,n+1):
        sum += 1.0/(i ** beta)
    return sum

def calculate_all_probs(df,colomnsInConstraints,beta):
    """
    calculate_all_probs - helper function in order to calculate all columns propabilities and store them into a dictionary

    Parameters
    ----------
    df : dataframe
    colomnsInConstraints : list of string
                            a list of all clomumns who are part of a constraint

    Returns
    -------
    a dictionary in which the keys are the columns and the values are lists of propabilities for each unique value

    """
    all_probs = defaultdict(dict)
    i = 1
    for col in colomnsInConstraints:
        temp_sum = harmonic_sum(len(df[col].unique()),beta)
        for cell in df[col].unique():
            value_prob = (1.0/(i ** beta))/float(temp_sum)
            all_probs[col][cell] = value_prob
            i += 1
        i = 1
    return all_probs


def randomize_value(df,data):
    """
    randomize_value - the function that randomize a new value based on its type

    Parameters
    ----------
    df : dataframe
    data : str/int/float/date
        the data which should be randomnly changed into a new value
        
    Returns
    -------
    Generates a new value to be assign into the database
    
    """
    if type(data) is str or isinstance(data, numpy.str_):
        val = data + random.choice(string.ascii_letters)
        
    elif type(data) is numpy.int64 :
        data_string = str(data)
        rnd_digit = random.randint(1,len(data_string))
        digit = int(data_string[rnd_digit-1])
        coin = random.randint(1, 2)
        if coin == 1 :
            if digit == 9 : digit = 0
            else : digit += 1
            data_string = data_string[:rnd_digit-1] + str(digit)+ data_string[rnd_digit-1 + 1:]
        if coin == 2 :
            if digit == 0 : digit = 9
            else : digit -= 1  
            data_string = data_string[:rnd_digit-1] + str(digit)+ data_string[rnd_digit-1 + 1:]
        val = int(data_string)
            
    elif type(data) is float or type(data) is numpy.float64 or type(data) is numpy.float_ or type(data) is numpy.float32:
        data_string = str(data)
        rnd_digit = 0
        while data_string[rnd_digit-1]=='.' or data_string[rnd_digit-1]=='-' or rnd_digit == 0 :
            rnd_digit = random.randint(1,len(data_string))
        digit = int(data_string[rnd_digit-1])
        coin = random.randint(1, 2)
        if coin == 1 :
            if digit == 9 : digit = 0
            else : digit += 1
            data_string = data_string[:rnd_digit-1] + str(digit)+ data_string[rnd_digit-1 + 1:]
        if coin == 2 :
            if digit == 0 : digit = 9
            else : digit -= 1  
            data_string = data_string[:rnd_digit-1] + str(digit)+ data_string[rnd_digit-1 + 1:]
        val = float(data_string)
        
    elif type(data) is date:
        new_day, new_month, new_year = [data.day,date.month,date.year]
        coin = random.randint(1, 3)
        if coin == 1 : 
            if data.month in [1,3,5,7,8,10,12]:
                new_day = random.randint(1,31)
            elif data.month in [4,6,9,11]:
                new_day = random.randint(1,30)
            elif data.month == 2 and data.year % 4 == 0 and data.year % 100 != 0 :
                new_day = random.randint(1,29)
            else : random.randint(1,28)
        if coin == 2 :
            new_month = random.randint(1,12)
        if coin == 3 :
            new_year = random.randint(1921,datetime.datetime.now().year)
        val = datetime.datetime(new_year, new_month, new_day)
    else:
        print("Error: data type not supported. Data: ", data)
        import pdb; pdb.set_trace()
    return val

def replace_value(df,data_col,all_probs):
    # pick random value from all_probs[col] with probability all_probs[col].values()
    
    return numpy.random.choice(list(all_probs[data_col].keys()), 1, list(all_probs[data_col].values()))[0]
    

def flip(p):
    return 1 if random.random() < p else 2


def rand_vio_algorithm(df,colomnsInConstraints,all_probs,typo_prob):
    """
    rand_vio_algorithm - the function chooses a random cell in the database and randomly chooses whether to
                         replace the value with a different value from the column or randomize a new value.

    Parameters
    ----------
    df : dataframe
    colomnsInConstraints : list of string
                            a list of all clomumns who are part of a constraint
    all_probs : dictionary of arrays
                dictionary in which the keys are the columns of the database and the values are
                list of probabilities for each value
        
    Returns
    -------
    Generates a new value to be assign into the database
    
    """
    rand_cell_row = random.choice(df.index)
    rand_cell_col = random.choice(colomnsInConstraints)                  
    rand_cell_data = df.at[rand_cell_row, rand_cell_col]
    
    while(pd.isnull(df.at[rand_cell_row, rand_cell_col])) :
        rand_cell_row = random.choice(df.index)
        rand_cell_col = random.choice(colomnsInConstraints)                  
        rand_cell_data = df.at[rand_cell_row, rand_cell_col]
    
    coin = flip(typo_prob)
    if coin == 1:
        new_val = randomize_value(df,rand_cell_data)
    if coin == 2:
        new_val = replace_value(df,rand_cell_col,all_probs)
    df.at[rand_cell_row,rand_cell_col] = new_val