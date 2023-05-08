##This script is written for automation of H-bond analysis--Krishnendu##

import numpy as np
import os, subprocess
import glob
import pandas as pd
import argparse
import shutil
from functools import reduce


parser = argparse.ArgumentParser(description='Create a H-bond occupancy table for similar systems')
#Required Parameter
parser.add_argument("-dirs",help="The directory names \n if there are 2 systems inside 2 separate directories  system1 & system2 then system should be given as input.",required=True)
#Parameters for VMD H-bond calculation
parser.add_argument("-sel1",default="protein", type=str, required=False,help="<atom selection> (default: protein)")
parser.add_argument("-sel2",default="none", type=str, required=False,help="<atom selection> (default: none)")
parser.add_argument("-dist",default=3.0, type=float, required=False,help="<cutoff distance between donor and acceptor> (default: 3.0)")
parser.add_argument("-ang",default=20, type=float, required=False,help="<angle cutoff> (default: 20)")
parser.add_argument("-polar",default="no", type=str, required=False,help="consider only polar atoms (N, O, S, F)? default: no")
parser.add_argument("-DA",default="both", type=str, required=False, help="sel1 is the donor (D), acceptor (A), or donor and acceptor (both) Only valid when used with two selections, default: both")
parser.add_argument("-i",default=1, type=int, required=False,help="Initial frame for Hbond calculation using VMD")
parser.add_argument("-f",default=-1, type=int, required=False,help="Final frame for Hbond calculation using VMD")
parser.add_argument("-step",default=1, type=int, required=False,help="Steps for Hbond calculation using VMD")
parser.add_argument("-wait",default=2000, type=int, required=False,help="Waiting time for H-bond calculation using VMD")
#Parmeters for creating occupancy table
parser.add_argument("-start",default=1, type=int, help="starting residue",required=False)
parser.add_argument("-stop",default=1000, type=int, help="stop residue",required=False)
parser.add_argument("-occupancy",default=0, type=float, required=False,help="minimum occupancy")
parser.add_argument("-o",default="H-bond-occupancy", type=str, required=False,help="output-filename")
#parameters for sorting the table 
parser.add_argument("-difference",default=15, type=float, required=False,help="If difference of hydrogen bond highest and lowest occupancy value is less than difference then this pair will not be considered")
parser.add_argument("-ival", default=0, type=int, required=False, help="Initial state for computing difference in H-bond occupancy")
parser.add_argument("-fval", default=-1, type=int, required=False, help="Final state for computing difference in H-bond occupancy")


args = parser.parse_args()

dirs = args.dirs
sel1 =args.sel1
sel2 =args.sel2
dist =args.dist
ang  =args.ang
polar=args.polar
DA   =args.DA
initial_frame =args.i
final_frame =args.f
step =args.step
wait =args.wait
start = args.start
stop = args.stop
ival = args.ival
fval = args.fval
difference = args.difference
occupancy = args.occupancy
output =args.o




######################################print parameters used if needed######################
def print_parameters_used():
    print(dirs)
    print(sel1)
    print(sel2)
    print(dist)
    print(ang)
    print(polar)
    print(DA)
    print(initial_frame)
    print(final_frame)
    print(step)
    print(wait)
    print(start)
    print(stop)
    print(difference)
    print(occupancy)
    print(output)
    return None
######################################print parameters used if needed######################



def find_important_H_bond(file,occupancy=occupancy,start=start,stop=stop):
    """
    This function finds important hydrogen bonds in a given file based on occupancy, start and stop
    positions.
    
    :param file: The name of the file containing the data to be analyzed
    :param occupancy: The occupancy threshold for selecting hydrogen bonds. Only hydrogen bonds with
    occupancy greater than or equal to this value will be considered
    :param start: The starting residue number for the search
    :param stop: The stop parameter is the upper limit of the range of residue IDs to consider for
    hydrogen bond analysis. Any hydrogen bond involving a residue with an ID greater than or equal to
    stop will be excluded from the analysis
    :return: two lists: 'o' which contains the occupancy values of the hydrogen bonds that meet the
    specified criteria, and 'ids' which contains the IDs of the atoms involved in those hydrogen bonds.
    """
    f=open(file,"r")
    p=f.readlines()
    ids=[]
    o=[]
    for i in range (len(p)):
        if (
            i > 1
            and (
                (
                    int(p[i].split()[0][3:-5]) > start
                    and int(p[i].split()[0][3:-5]) < stop
                )
                or (
                    int(p[i].split()[1][3:-5]) > start
                    and int(p[i].split()[1][3:-5]) < stop
                )
            )
            and float(p[i].split()[2][:-1]) >= occupancy
        ):
            d=int(p[i].split()[0][3:-5])
            a=int(p[i].split()[1][3:-5])
            if d > a:
                ids.append(p[i].split()[1][:-5]+" "+p[i].split()[0][:-5])
            if a >= d:
                ids.append(p[i].split()[0][:-5]+" "+p[i].split()[1][:-5])
            o.append(float(p[i].split()[2][:-1]))
    return o,ids

def get_unique_occupancy(pairs,array):
    """
    The function calculates the total occupancy of unique pairs in an array.
    
    :param pairs: It is a list of unique values that represent the different types of occupancy that are
    being considered
    :param array: The `array` parameter is a 2D numpy array with two rows. The first row contains
    numerical values representing the occupancy of different rooms, and the second row contains strings
    representing the type of room
    :return: The function `get_unique_occupancy` returns an array `occ` which contains the sum of the
    values in `array[0]` for each unique element in `pairs` that is found in `array[1]`.
    """
    occ=np.zeros(len(pairs))
    for i in range (len(pairs)):
        for j in range (len(array[1])):
            if pairs[i]==array[1][j]:  
                occ[i]=occ[i]+array[0][j]
    return occ

def process_dataframe(DataFrame, difference):
    """
    The function removes rows from a pandas DataFrame where the difference between the minimum and
    maximum values in each row is less than a specified threshold.
    
    :param DataFrame: A pandas DataFrame containing numerical data
    :param difference: The minimum difference between the maximum and minimum values in a row of the
    DataFrame that is required for the row to be kept in the new DataFrame. If the difference is less
    than this value, the row will be removed from the new DataFrame
    :return: a new data frame with rows removed where the difference between the minimum and maximum
    values in the row is less than the specified difference.
    """
    minimum= DataFrame.min(numeric_only=True, axis=1)
    maximum= DataFrame.max(numeric_only=True, axis=1)
    diff=maximum-minimum
    remove_index=np.where(diff< difference)[0].tolist()
    return DataFrame.drop(remove_index, axis=0)

#############################Check if all necessary softwares are present#####################
if shutil.which('vmd') is None:
    "print VMD not present in your system please install VMD"
#############################Check if all necessary softwares are present#####################

#######################################Actual Code##########################################
VMD=shutil.which('vmd')
Pymol=shutil.which('pymol')
TCL_Path=VMD[:-7]+"lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl"
print(VMD)
F=glob.glob(dirs+"*")
F.sort()
XTCS=glob.glob(dirs+"*/*.xtc")
XTCS.sort()
PDBfiles=glob.glob(dirs+"*/*.pdb")
PDBfiles.sort()




print("#################################Performing H-Bond Calculation using VMD#################################")
for i in range (len(F)):
    print(f"system{str(i + 1)}")
    with open(F[i]+'/input_script.vmd', 'w') as f:
        str1="mol addfile "+XTCS[i]+" first "+str(initial_frame)+" last "+str(final_frame)+" step "+str(step)+" waitfor "+str(wait)+"\n"
        if sel2!="none":
            str2="source "+TCL_Path+" \n set sel1 [atomselect top "+sel1+" ] \n"
            str3 = (
                f"hbonds -sel1 $sel1 -sel2 $sel2 -dist {str(dist)} -ang {str(ang)} -polar "
                + polar
                + " -DA "
                + DA
                + " -writefile yes -plot no -outfile "
                + F[i]
                + "/hbonds.dat -type pair -detailout "
                + F[i]
                + "/hbonds-details.dat \n"
            )
        else:
            str2="source " +TCL_Path+" \n set sel1 [atomselect top "+sel1+" ] \n set sel2 [atomselect top "+sel2+" ] \n"
            str3 = (
                f"hbonds -sel1 $sel1 -dist {str(dist)} -ang {str(ang)} -polar "
                + polar
                + " -DA "
                + DA
                + " -writefile yes -plot no -outfile "
                + F[i]
                + "/hbonds.dat -type pair -detailout "
                + F[i]
                + "/hbonds-details.dat \n"
            )
        str4="exit"
        f.write(str1+str2+str3+str4)
        f.close()
    subprocess.run(VMD+" "+PDBfiles[i]+" -dispdev none -e "+F[i]+"/input_script.vmd",shell=True)

Files=[F[i]+"/hbonds-details.dat" for i in range (len(F))]
print("################################# H-Bond Calculation is over #################################")

Pairs=[]
unique_pairs=[]
Systems=[]

for i in range (len(F)):
    Systems.append([])
    a=find_important_H_bond(Files[i])
    Pairs.append(a)
    unique_pairs=unique_pairs+a[1]
unique_pairs=np.unique(unique_pairs)

for i in range (len(Pairs)):
    Systems[i]=get_unique_occupancy(unique_pairs,Pairs[i])

Resid1=[]
Resid2=[]
for i in range (len(unique_pairs)):
    Resid1.append(unique_pairs[i].split()[0])
    Resid2.append(unique_pairs[i].split()[1])

Resid1=pd.Series(Resid1)
Resid2=pd.Series(Resid2)
for i in range (len(Systems)):
    Systems[i]=pd.Series(np.around(Systems[i]))

colname = [f"System{str(i + 1)}" for i in range (len(Systems))]
rowname=[Resid1,Resid2]
rowname.extend(iter(Systems))
rowname=np.array(rowname, )
rows=rowname.T

df = pd.DataFrame()
df["#Resid1"]=rows[:, 0]
df["Resid2"]=rows[:, 1]
for i in range (len(colname)):
    df[f"System{str(i + 1)}"] = np.array(rows[:, 2+i], dtype=float)

df["Difference"] = (
    df[f"System{str(F[fval][-1])}"] - df[f"System{str(F[ival][-1])}"]
)




df_new=process_dataframe(df, difference)
df_final=df_new.sort_values("Difference")
df_final.to_csv(output+".csv",index=False,sep="\t")

subprocess.run(Pymol+" "+PDBfiles[0],shell=True)
exit()
#



