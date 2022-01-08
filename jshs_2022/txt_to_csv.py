import os
import pandas as pd

def txt_to_csv(dirname):
    l = os.listdir(dirname)
    for i in range(len(l)):
        if l[i] != ".DS_Store":
            try:
                txt_to_csv(dirname + "/" + l[i])
            except:
                if l[i][-3:] == "txt":
                    temp = pd.read_csv(dirname + "/" + l[i], sep="\t")
                    temp.to_csv(dirname + "/" + l[i][:-3] + "csv")
                    os.remove(dirname + "/" + l[i])

txt_to_csv("/Users/jaisharma/Documents/jshs_2022")
