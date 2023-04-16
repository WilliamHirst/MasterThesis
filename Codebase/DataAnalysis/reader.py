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
    
    
