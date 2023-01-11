import matplotlib.pyplot as plt
from Plot_stuff.ROCM import *



def HM(model, X, Y, W, columns):
    columns_s = columns[Y.to_numpy() == 1] 
    unique_c = columns_s.unique()


    print(unique_c)
    for c in unique_c:
        index_i = (columns == c).to_numpy()
        print(index_i)
        X_i = X[index_i]
        Y_i = Y[index_i]
        W_i = W[index_i]
        auc = plotRoc(  Y_i, 
                        model.predict(X_i, batch_size=8192), 
                        W_i,
                        "",
                        plot = False,
                        return_score = True)
        

        
        
