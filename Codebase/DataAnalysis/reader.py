import uproot
import numpy as np
import ROOT
file = uproot.open("results/ATLAS/WZ_graphs_onshell.root") 
gr = file["Exp_0"].all_members
print(gr)
file = open("results/text/compLimt.txt", 'w')

print(gr["fX"])
for i in zip(gr["fX"], gr["fY"]):
    print(i)
# print(len(file["Exp_0"].values()[0]))
# print(len(file["Exp_0"].values()[1]))
# gr = file["fID_gr"].all_members
# print(len(gr["fX"]))
# print(len(gr["fY"]))

# print(gr)
# print(file["Exp_0"].values()[0])


# ['CLs_gr;1', 'CLsexp_gr;1', 'clsu1s_gr;1', 'clsu2s_gr;1', 'clsd1s_gr;1', 'clsd2s_gr;1', 'upperLimit_gr;1', 
#  'expectedUpperLimit_gr;1', 'fIDList;1', 'fID_gr;1', 'forbiddenFunction$mirror;3', 'forbiddenFunction$mirror;2', 
#  'forbiddenFunction$mirror;1', 'SubGraphs;1', 'SubGraphs/CLs_Contour_0;1', 'SubGraphs/CLsexp_Contour_0;1', 
#  'SubGraphs/clsu1s_Contour_0;1', 'SubGraphs/clsu2s_Contour_0;1', 'SubGraphs/clsd1s_Contour_0;1', 
#  'SubGraphs/clsd2s_Contour_0;1', 'SubGraphs/CLs_Contour_0_Up;1', 'SubGraphs/CLsexp_Contour_0_Up;1', 
#  'SubGraphs/clsu1s_Contour_0_Up;1', 'SubGraphs/clsu2s_Contour_0_Up;1', 'SubGraphs/clsd1s_Contour_0_Up;1', 
#  'SubGraphs/clsd2s_Contour_0_Up;1', 'SubGraphs/CLs_Contour_0_Down;1', 'SubGraphs/CLsexp_Contour_0_Down;1', 
#  'SubGraphs/clsu1s_Contour_0_Down;1', 'SubGraphs/clsu2s_Contour_0_Down;1', 'SubGraphs/clsd1s_Contour_0_Down;1', 
#  'SubGraphs/clsd2s_Contour_0_Down;1', 'Band_1s_0;1', 'Band_2s_0;1', 'Obs_0;1', 'Exp_0;1', 'FinalCurves;1', 'Obs_0_Up;1', 
#  'Exp_0_Up;1', 'Obs_0_Down;1', 'Exp_0_Down;1']
exit()
file = open("results/text/compLimt.txt", 'w')

PNN = open("results/text/PNNPCA_FS_MLMGridSig.txt", 'r')
NN = open("results/text/NN_FS_MLMGridSig.txt", 'r')
MaxOut = open("results/text/MaxOutPCA_FS_MLMGridSig.txt", 'r')

PNN.readline()
NN.readline()
MaxOut.readline()


file.write("m1    m2    PNNnbkg    PNNnsig    PNNnexpsig    NNnbkg    NNnsig    NNnexpsig    MaxOutnbkg    MaxOutnsig    MaxOutnexpsig\n")

for line_PNN_F, line_NN_F, line_Maxout_F in zip(PNN.readlines(), NN.readlines(), MaxOut.readlines()):
    line_PNN = [line.replace("\n", "") for line in line_PNN_F.split("    ")]
    line_NN = [line.replace("\n", "") for line in line_NN_F.split("    ")]
    line_MaxOut = [line.replace("\n", "") for line in line_Maxout_F.split("    ")]
    
    addLine = f"{line_PNN[0]}    {line_PNN[1]}    {line_PNN[2]}    {line_PNN[3]}    {line_PNN[4]}    {line_NN[2]}    {line_NN[3]}    {line_NN[4]}    {line_MaxOut[2]}    {line_MaxOut[3]}    {line_MaxOut[4]}\n"
    file.write(addLine)
file.close()
    
    
