import networkx as nx
import collections 
import numpy as np
from scipy.stats import cauchy
from sklearn.isotonic import IsotonicRegression  
from sklearn.linear_model import LinearRegression 
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import comb
from collections import defaultdict
import time
import re
import pandasql as psql
import gurobipy as gp
import pandas as pd
from gurobipy import GRB
from collections import Counter

# utility functions
def get_sorted_deg_seq(G): # return ascending degree sequence
    deg_seq = sorted([d for n, d in G.degree()], reverse=False) 
    return deg_seq

def get_deg_his(G, max_deg): # degree histogram from graph G, capped at max_deg
    deg_seq = get_sorted_deg_seq(G)
    degree_count = collections.Counter(deg_seq)
    deg_his = np.zeros(max_deg+1)
    for deg in degree_count:
        deg_his[deg] = degree_count[deg]    
    return deg_his

def deg_seq_to_deg_his(deg_seq, max_deg):
    #assume deg sequence could be non-integer and be bigger than maxDegree
    deg_his = np.zeros(max_deg + 1)
    for deg in deg_seq:
        deg = int(round(deg))
        if deg <= max_deg:
            deg_his[deg]= deg_his[deg]+1
    return deg_his

def his_to_pdf(his):
    area = sum(his)
    if area == 0:
        print("AREA 0!!!!")
        print(his)
    max_deg = len(his)
    pdf = np.zeros(max_deg)
    for i in range(max_deg):
        pdf[i] = his[i]/area
    return pdf
    
def pdf_to_cdf(pdf): # convert probability density function to cumulative distribution
    cdf = np.zeros(len(pdf))
    cdf[0] = pdf[0]
    for i in range(1,len(pdf)):
         cdf[i] = cdf[i-1] + pdf[i]            
    return cdf

def cdf_to_pdf(cdf): # convert cumulative distribution to probability density 
    pdf = np.zeros(len(cdf))
    pdf[0] = cdf[0]
    for i in range(1,len(pdf)):
         pdf[i] = cdf[i] - cdf[i-1]            
    return pdf
    

def dif_deg_his_L1(his1,his2): # compute L1 error between 2 histogram
    #assume the same length
    return sum(abs(his1 - his2))

def dif_deg_his_L2(his1,his2): # compute L2 error between 2 histogram
    return sum(np.square(his1-his2))



def plot_his(trueHis,noisyHis):
    plt.plot(trueHis,'-g', label='trueHis')
    plt.plot(noisyHis,'--r', label='noisyHis')
    plt.legend();
    plt.xscale('log')

    
def plot_cum(trueHis,noisyHis):
    plt.plot(pdfToCdf(trueHis), '3b', label='trueCum')
    plt.plot(pdfToCdf(noisyHis), '2y', label='noisyCum')
    plt.legend();
    plt.xscale('log')

#DP basic functions
def add_laplace(true_counts, sens, epsilon): # Add Laplace noise to true distribution, given sensitivity and epsilon
    scale = 1.0* sens/epsilon
    noisy = true_counts + np.random.laplace(0.0, scale, len(true_counts))
    return noisy

# Smooth out the Laplace noise added to cdf using isotonic regression
def post_process_cdf(noisy_cdf, total_count): 
    ir = IsotonicRegression(y_min=0, y_max=total_count, increasing=True)
    cdf= ir.fit_transform(X=range(len(noisy_cdf)),y=noisy_cdf)   
    return cdf

def post_process_pdf(noisy_pdf, nodes_num):
    cdf = pdf_to_cdf(noisy_pdf)
    cdf = post_process_cdf(cdf, nodes_num)
    pdf = cdf_to_pdf(cdf)
    return pdf


def extend_his(his,max_degree): # extend truncated histogram to have {max_degree+1} points
    #his has a shortern length 
    if (max_degree+1) > len(his):
        his_extended = np.zeros(max_degree + 1)
        his_extended[0:len(his)] = his
        return his_extended
    else:
        return his

    
def sample_prob_list(prob_list):
    normalized = prob_list/sum(prob_list)
    return np.random.choice(len(prob_list),1,p=normalized)[0]

    # normalized = prob_list/sum(prob_list)
    # r = np.random.uniform(0,1,1)
    # s = 0 
    # for i in range(len(prob_list)):
    #     s += normalized[i]
    #     if s >= r:
    #         return i
    # return len(prob_list)-1


#graph transformation/clean up for subgraph counting aglo (e.g. ladder function) 
#this remap the node id, such that node id starts from 0 and increments to the total number of nodes 
def translate(datafile, newDatafile):

    nodeMap = dict()
    
    fin = open(datafile, "r")
    fout = open(newDatafile, "w")
    for ij in fin:
        i,j = ij.split()
        #i = int(i)
        #j = int(j)
        if i not in nodeMap:
            nodeMap[i] = len(nodeMap)
        if j not in nodeMap:
            nodeMap[j] = len(nodeMap)
        
        i_map = nodeMap[i]
        j_map = nodeMap[j]
        if i_map < j_map:
            fout.write(str(i_map)+" "+str(j_map)+"\n")
        else:
            fout.write(str(j_map)+" "+str(i_map)+"\n")

    fout.close()
    fin.close() 




def build_dynamic_queries(constraintSets,df):
    """
    build_dynamic_queries - generates dynamic queries based on the given constraints.
    This function will generate two queries:
    1. unionOfAllTuples - returns the ids of the tuples participating in a violation of the constraints.
    2. unionOfAllPairs - returns pairs (i1,i2) of ids of tuples that jointly violate the constraints.
    
    Parameters
    ----------
    constraintSets : set of strings
        each string represents a constraint from the dcs file
    df : dataframe
        the database frame
        
    Returns
    -------
    list of three string values:
        unionOfAllTuples, unionOfAllPairs are the generated queries
        allColumns is a string consisting of all column names seperated by ','
    """ 
    
    allColumns = ' '.join([str(elem) for elem in df.columns.values.tolist()]).replace(' ',',')

    #Additional conditions for the queries, in order to ignore missing values in the database
    count = 1
    columnsT1 = ""
    for col in df.columns: 
        columnsT1 += "t1."+col
        if count!=len(df.columns) :
            columnsT1+= ' IS NOT NULL AND '
        count+=1
    columnsT1+=" IS NOT NULL "
    columnsT2 = columnsT1.replace('t1','t2')

    count = 0
    for con in constraintSets: 
        if count == 0:
            unionOfAllPairs = " SELECT t1.rowid as t1ctid ,t2.rowid as t2ctid FROM df t1,df t2 WHERE "
            unionOfAllTuples = " SELECT * FROM df t1,df t2 WHERE " 
        else : 
            unionOfAllPairs += " UNION SELECT t1.rowid as t1ctid ,t2.rowid as t2ctid FROM df t1,df t2 WHERE "
            unionOfAllTuples += " UNION SELECT * FROM df t1,df t2 WHERE "
            
        rep = {" ": "_", "&": " and ","not(":"",")":""} 
        rep = dict((re.escape(k), v) for k, v in rep.items()) 
        pattern = re.compile("|".join(rep.keys()))
        con1 = pattern.sub(lambda m: rep[re.escape(m.group(0))], con)     

        # in case the constraint refers to a single tuple
        if "t2" not in con: 
            con1 = re.sub(r'(t1.*?)t1', r'\1t2', con1, 1)
            unionOfAllPairs += con1 +" and t1.ROWID==t2.ROWID and ("+columnsT1+")"
            unionOfAllTuples += con1 +" and t1.ROWID==t2.ROWID and ("+columnsT1+")"
        else:
            unionOfAllPairs += con1 +" and t1.ROWID!=t2.ROWID and ("+columnsT1+" and "+columnsT2+")"
            unionOfAllTuples += con1 +" and t1.ROWID!=t2.ROWID and ("+columnsT1+" and "+columnsT2+")"
        count+=1
        
    return unionOfAllTuples,unionOfAllPairs,allColumns

def constraints_check(df,constraintSets, allColumns, unionOfAllTuples, unionOfAllPairs):
    """
    constraints_check - runs the dynamic queries that have been generated on the database.
    This function will run two queries:
    1. unionOfAllTuples - returns the ids of the tuples participating in a violation of the constraints.
    2. unionOfAllPairs - returns pairs (i1,i2) of ids of tuples that jointly violate the constraints.
    
    Parameters
    ----------
    constraintSets : set of strings
        each string represents a constraint from the dcs file
    unionOfAllTuples : string
    unionOfAllPairs : string    
    df : dataframe
        the database frame
        
    Returns
    -------
    list of two strings and two double variables:
        sdfcWithRep, sdfcNoRep are the results of the unionOfAllPairs and unionOfAllTuples queries, respectively.
        end1-start, end2-start2 are the running times the queries.
    """     
    
    # finds the pairs of tuples that jointly violate a constraint
    start = time.time()    
    violatingPairs =  psql.sqldf("SELECT DISTINCT * FROM (SELECT CASE WHEN t1ctid <= t2ctid THEN t1ctid ELSE t2ctid END AS id1,CASE WHEN t1ctid <= t2ctid THEN t2ctid ELSE t1ctid END AS id2 FROM ("+unionOfAllPairs+")AS A)AS B")
    end1 = time.time()
    
    # finds the tuples that participate in a violation
    # start2 = time.time()
    # violatingTuples = set()
    # for pair in violatingPairs.values:
    #     for item in pair:
    #         violatingTuples.add(item)
    # end2 = time.time()
    
    # return violatingPairs, violatingTuples, end1-start, end2-start2
    return violatingPairs, end1-start


def I_lin_R(rows_violations):
    """
    I_lin_R: computes the measure I^lin_R that is the linear relaxation of the ILP used for computing
    the measure I_R.
    
    - There is a variable x for every tuple in the database such that 0<=x<=1.
    - The constraints are of the form x + y >= 1 where x and y represent two tuples that jointly vioalte a constraint.
    - The objective function is to minimize the sum of all x's.
    
    Parameters
    ----------
    uniquePairsDf : dataframe
        the result of the query that finds all pairs of tuples that jointly violate a constraint.
        
    Returns
    -------
    list of two int variables:
        database_measurer.objVal is the result of the LP.
        end2 - start is the running time of the function.
    """ 
    
    start = time.time()
    # rows_violations = uniquePairsDf.values
    varsDict2 = {}
    database_measurer = gp.Model('Minimal deletions of tuples relaxed')
    database_measurer.setParam('OutputFlag', 0)  # do not show any comments on the screen 
    
    # variables
    for i in rows_violations :
        varsDict2[i[0]] = database_measurer.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x")
        varsDict2[i[1]] = database_measurer.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS, name="x")
    
    # constraints
    for i in rows_violations :
        database_measurer.addConstr(varsDict2[i[0]]+varsDict2[i[1]]>=1, name='con')
    vars= []
    for i in varsDict2:
        vars.append(varsDict2[i])
    
    # objective function
    database_measurer.setObjective(sum(vars), GRB.MINIMIZE)
    
    opt = database_measurer.optimize()
    end2 = time.time()
    return database_measurer.objVal , end2 -start

def check_tuple_violation(df, constraints, vios_epsilon):
    
    allColumns = {}
    orig_columns = df.columns
    for col in df.columns: 
        allColumns[col] = col.replace(' ','_')
    df = df.rename(columns=allColumns)

    distinct_count_dict = {}
    for constraint in constraints:
        allConstraints = build_dynamic_queries([constraint],df)
        violating_pairs, time1 = constraints_check(df,constraints, allColumns, allConstraints[0], allConstraints[1])
        
        distinct_tuples = violating_pairs['id1']._append(violating_pairs['id2']).unique()
        if len(violating_pairs) != 0:
            print("Violating pairs: ", len(violating_pairs))
            
        for tuple_id in distinct_tuples:
            if tuple_id in distinct_count_dict:
                distinct_count_dict[tuple_id] += 1
            else:
                distinct_count_dict[tuple_id] = 1

    if len(distinct_count_dict) == 0:
        return 0
    
    counts = Counter(distinct_count_dict.values()).values()
    values = Counter(distinct_count_dict.values()).keys()
    noisy_counts = np.array(list(counts)) + np.random.laplace(0.0, 1/vios_epsilon, len(counts))
    noisy_counts = noisy_counts.clip(min=0)
    argmax = np.argmax(noisy_counts)
    most_violations_per_tuple = list(values)[argmax]
    return most_violations_per_tuple

    



def get_upperbound(dc_dir, data_file, opt_theta_eps):
    # Calculates the upper bound of the number of violations of the constraints according to Algo 3

    constraints_raw = open(dc_dir +'/dcs.txt', 'r')
    constraints = [line.strip() for line in constraints_raw.readlines()]
    constraints = [x.replace(' ', '_') for x in constraints] #in case the columns names include spaces

    dataset_df = pd.read_csv(data_file)

    constraint_elem = []
    for c in constraints:
        c = c[4:-1].split('&')
        c = [re.split('(!=|>=|<=|>|<|=)', i) for i in c]
        constraint_elem.append(c)
    
    fd_counts = 0 # count of FDs
    all_fd = []
    for constraint_tuples in constraint_elem:
        fd_check = True
        for i in constraint_tuples:
            if i[1] == '=' or i[1] == '!=': # only consider equality constraints for FDs
                pass
            else:
                fd_check = False
                break
        if fd_check:
            fd_counts += 1
            all_fd.append(constraint_tuples)

    epsilon_per_fd = opt_theta_eps/( fd_counts ) # epsilon for each FD

    theta_estimate = 0
    all_equality_attributes = {}
    for constraint_tuples in all_fd:
        fd_equality_attributes = []
        for i in constraint_tuples:
            if i[1] == '=':
                eq_attr = i[0].split('.')[1]
                fd_equality_attributes.append(eq_attr)
                all_equality_attributes[eq_attr] = ''

        #find frequency of the joint distribution of the equality attributes
        freq = dataset_df.groupby(fd_equality_attributes).size().reset_index(name='count')
        
        freq['count'] = freq['count'].apply(lambda x: 0 if x < 0 else x) #make all the counts positive
        freq['count'] += np.random.laplace(0.0, 1/epsilon_per_fd) # add laplace noise
        
        theta_estimate += freq['count'].max() # find max for each FD and add to theta_estimate

    return int(theta_estimate)
    