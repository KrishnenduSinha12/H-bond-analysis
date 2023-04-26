# This code gives gives a comparison table of hydrogen bond occupancy of similar systems (eg Apo and Ligand Bound proteins, proteins with different mutations).


# Installation:
  copy the H-bond-analysis.py file to your directory and use it   

# How to use the code:
1) create folders corresponding to the systems for eg you have 3systems Apo Bound and intermediate, creates folders with name system1 system2 and system3 
2) Inside each system folder there should be one pdbfile and one xtcfile.
3) To run the analysis script use python H-bond-analysis.py -dirs dirname (here system) 
4) for details use python H-bond-analysis --help
5) It will write a csv file containing all the information 
6) use the vis.py script to visualize the H-bond-occupancy difference between 2 systems. for that open a pdb using pymol. load the vis.py script as run vis.py. Then call the draw_hbond function as draw_hbond("H-bond-occupancy.csv", cutoff=20, color_positive="red", color_negative="blue")
