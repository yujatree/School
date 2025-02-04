import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt(f'./data/15050000')
data = (data.flatten())[0:10000].reshape(100,100)
data = (data > 0.5).astype(int)
print(data)

#def cluster_size(data):
#    clusterIDs = -1 * np.ones(data.shape).astype(int)
#    idx = 0
#
#    # initialize
#    for i in range(data.shape[0]):
#        for j in range(data.shape[1]):
#            if (data[i,j] == 1):
#                clusterIDs[i,j] = idx
#                idx +=1
#
#    while(True):
#        history = clusterIDs.copy()
#
#        for i in range(data.shape[0]):
#            for j in range(data.shape[1]):
#                if (clusterIDs[i,j] != -1):
#
#                    if (i == 0):
#                        if (clusterIDs[data.shape[0]-1,j] != -1): clusterIDs[i,j], clusterIDs[data.shape[0]-1,j] = min(clusterIDs[i,j], clusterIDs[data.shape[0]-1,j]), min(clusterIDs[i,j], clusterIDs[data.shape[0]-1,j])
#                        if (clusterIDs[i+1,j] != -1): clusterIDs[i,j], clusterIDs[i+1,j] = min(clusterIDs[i,j], clusterIDs[i+1,j]), min(clusterIDs[i,j], clusterIDs[i+1,j])
#                    elif (i == data.shape[0]-1):
#                        if (clusterIDs[i-1,j] != -1): clusterIDs[i,j], clusterIDs[i-1,j] = min(clusterIDs[i,j], clusterIDs[i-1,j]), min(clusterIDs[i,j], clusterIDs[i-1,j])
#                        if (clusterIDs[0,j] != -1): clusterIDs[i,j], clusterIDs[0,j] = min(clusterIDs[i,j], clusterIDs[0,j]), min(clusterIDs[i,j], clusterIDs[0,j])
#                    else:
#                        if (clusterIDs[i-1,j] != -1): clusterIDs[i,j], clusterIDs[i-1,j] = min(clusterIDs[i,j], clusterIDs[i-1,j]), min(clusterIDs[i,j], clusterIDs[i-1,j])
#                        if (clusterIDs[i-1,j] != -1): clusterIDs[i,j], clusterIDs[i-1,j] = min(clusterIDs[i,j], clusterIDs[i-1,j]), min(clusterIDs[i,j], clusterIDs[i-1,j])
#
#                    if (j == 0):
#                        if (clusterIDs[i,data.shape[1]-1] != -1): clusterIDs[i,j], clusterIDs[i,data.shape[1]-1] = min(clusterIDs[i,j], clusterIDs[i,data.shape[1]-1]), min(clusterIDs[i,j], clusterIDs[i,data.shape[1]-1])
#                        if (clusterIDs[i,j+1] != -1): clusterIDs[i,j], clusterIDs[i,j+1] = min(clusterIDs[i,j], clusterIDs[i,j+1]), min(clusterIDs[i,j], clusterIDs[i,j+1])
#                    elif (j == data.shape[1]-1):
#                        if (clusterIDs[i, j-1] != -1): clusterIDs[i,j], clusterIDs[i,j-1] = min(clusterIDs[i,j], clusterIDs[i,j-1]), min(clusterIDs[i,j], clusterIDs[i,j-1])
#                        if (clusterIDs[i,0] != -1): clusterIDs[i,j], clusterIDs[i,0] = min(clusterIDs[i,j], clusterIDs[i,0]), min(clusterIDs[i,j], clusterIDs[i,0])
#                    else:
#                        if (clusterIDs[i, j-1] != -1): clusterIDs[i,j], clusterIDs[i,j-1] = min(clusterIDs[i,j], clusterIDs[i,j-1]), min(clusterIDs[i,j], clusterIDs[i,j-1])
#                        if (clusterIDs[i,j+1] != -1): clusterIDs[i,j], clusterIDs[i,j+1] = min(clusterIDs[i,j], clusterIDs[i,j+1]), min(clusterIDs[i,j], clusterIDs[i,j+1])
#
#        if (np.all(clusterIDs == history)): break
#
#    for key, id in enumerate(np.unique(clusterIDs)):
#        if (id != -1):
#            clusterIDs[clusterIDs == id] = key
#
#    print(clusterIDs)
#    print(np.unique(clusterIDs, return_counts = True))
#
#cluster_size(data)
