##This script is written for automation of H-bond analysis--Krishnendu##

import numpy as np
import os, subprocess
import glob
import pandas as pd
import argparse
import shutil



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
occupancy = args.occupancy
output =args.o


def find_important_H_bond(file,occupancy=occupancy,start=start,stop=stop):
    f=open(file,"r")
    p=f.readlines()
    ids=[]
    o=[]
    for i in range (len(p)):
        if i >1: 
            if ((int(p[i].split()[0][3:-5])>start and int(p[i].split()[0][3:-5])<stop) or (int(p[i].split()[1][3:-5])>start and int(p[i].split()[1][3:-5])<stop)) and float(p[i].split()[2][:-1])>=occupancy:
                d=int(p[i].split()[0][3:-5])
                a=int(p[i].split()[1][3:-5])
                if d > a:
                    ids.append(p[i].split()[1][:-5]+" "+p[i].split()[0][:-5])
                if a >= d:
                    ids.append(p[i].split()[0][:-5]+" "+p[i].split()[1][:-5])
                o.append(float(p[i].split()[2][:-1]))
    return o,ids

def get_unique_occupancy(pairs,array):
    occ=np.zeros(len(pairs))
    for i in range (len(pairs)):
        for j in range (len(array[1])):
            if pairs[i]==array[1][j]:  
                occ[i]=occ[i]+array[0][j]
    return occ


############################Test if required software is present #####################################
if shutil.which('vmd')==None:
    "print VMD not present in your system please install VMD"

########################################################Actual Code##########################################################################
VMD=shutil.which('vmd')
TCL_Path=VMD[:-7]+"lib/vmd/plugins/noarch/tcl/hbonds1.2/hbonds.tcl"

F=glob.glob(dirs+"*")
F.sort()
XTCS=glob.glob(dirs+"*/*.xtc")
XTCS.sort()
PDBfiles=glob.glob(dirs+"*/*.pdb")
PDBfiles.sort()


print("#################################Performing H-Bond Calculation using VMD#################################")
for i in range (len(F)):
    print("system"+str(i+1))
    with open(F[i]+'/input_script.vmd', 'w') as f:
        if sel2!="none":
            str1="mol addfile "+XTCS[i]+" first "+str(initial_frame)+" last "+str(final_frame)+" step "+str(step)+" waitfor 2000 \n"
            str2="source "+TCL_Path+" \n set sel1 [atomselect top "+sel1+" ] \n"
            str3="hbonds -sel1 $sel1 -sel2 $sel2 -dist "+str(dist)+" -ang "+str(ang)+" -polar "+polar+" -DA "+DA+" -writefile yes -plot no -outfile "+F[i]+"/hbonds.dat -type pair -detailout "+F[i]+"/hbonds-details.dat \n"
            str4="exit"
        else:
            str1="mol addfile "+XTCS[i]+" first "+str(initial_frame)+" last "+str(final_frame)+" step "+str(step)+" waitfor 2000 \n"
            str2="source " +TCL_Path+" \n set sel1 [atomselect top "+sel1+" ] \n set sel2 [atomselect top "+sel2+" ] \n"
            str3="hbonds -sel1 $sel1 -dist "+str(dist)+" -ang "+str(ang)+" -polar "+polar+" -DA "+DA+" -writefile yes -plot no -outfile "+F[i]+"/hbonds.dat -type pair -detailout "+F[i]+"/hbonds-details.dat \n"
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
    Systems[i]=pd.Series(Systems[i])

colname=["#Resid1","Resid2"]
for i in range (len(Systems)):
    colname.append("System"+str(i+1))
rowname=[Resid1,Resid2]
for j in range (len(Systems)):
    rowname.append(Systems[j])
rowname=np.array(rowname)
rows=rowname.T

df = pd.DataFrame(rows, columns=colname)
print(df)

df.to_csv(output+".csv",index=False,sep="\t")

exit()
