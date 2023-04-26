from pymol.cgo import *
from pymol import cmd
import numpy  as np
import pandas as pd


def draw_hbond(filename, cutoff, color_positive="red", color_negative="blue",gap=0,width=3):
    f=pd.read_csv(filename, sep="\t")
    resid1=f["#Resid1"].apply(lambda x:int(x[3:]))
    resid2=f["Resid2"].apply(lambda x:int(x[3:]))
    diff=np.array(f["Difference"],dtype=float)
    for i in range (len(resid1)):
        if diff[i] >cutoff:
            cmd.select("s1", "resid "+str(resid1[i])+" and name CA")
            cmd.select("s2", "resid "+str(resid2[i])+" and name CA")
            cmd.distance("d1", "s1", "s2")
            cmd.show("spheres", "s1 or s2")
            cmd.hide("labels", "d1")
            cmd.set("dash_gap", str(gap))
            cmd.color(color_positive, "d1")
            cmd.set("dash_width", str(width))
        elif diff[i] <-cutoff:
            print(i)
            cmd.select("g1", "resid "+str(resid1[i])+" and name CA")
            cmd.select("g2", "resid "+str(resid2[i])+" and name CA")
            cmd.distance("d2", "g1", "g2")
            cmd.show("spheres", "g1 or g2")
            cmd.hide("labels", "d2")
            cmd.set("dash_gap", str(gap))
            cmd.color(color_negative, "d2")
            cmd.set("dash_width", str(width))
cmd.extend("draw_hbond", draw_hbond)

#draw_hbond()
