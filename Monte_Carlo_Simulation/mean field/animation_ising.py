import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

itemlist = sorted(item for item in os.listdir('./one_data')
                        if not item.startswith('.'))

fig, ax = plt.subplots(figsize = (10,10))

def animateSnapshot(i):
    ax.clear()

    data = np.genfromtxt(f'./one_data/{i:08d}', delimiter = '')
    reshaped_data = (data.flatten())[:10000].reshape(100,100)
    
    image = ax.imshow(reshaped_data,
                      cmap = 'binary')
    title = ax.set_title(f"{i:08d}")

    return image, title

animation = FuncAnimation(fig, animateSnapshot,
                          repeat = True,
                          ##### EDIT HERE #####
                          frames = range(0, 200000, 1000))
                          ##### EDIT HERE #####
animation.save("./one_data/result.gif", dpi = 100, writer = PillowWriter(fps = 7))

print(" Created an Animation ! >:D ")
