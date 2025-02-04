import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import os

itemlist = sorted(item for item in os.listdir('./animation_correct')
                        if not item.startswith('.'))

fig, ax = plt.subplots(figsize = (10,10))

print(np.genfromtxt(f'./animation_correct/00000000')

def animateSnapshot(i):
    ax.clear()

    data = np.genfromtxt(f'./animation_correct/{i:08d}')
    reshaped_data = (data.flatten())[0:10000].reshape(100,100)

    image = ax.imshow(reshaped_data, delimiter = '',
                      cmap = 'binary')
    title = ax.set_title(f"{i:08d}")

    return image, title

animation = FuncAnimation(fig, animateSnapshot,
                          repeat = True,
                          ##### EDIT HERE #####
                          frames = range(0,10000, 1000))
                          ##### EDIT HERE #####
animation.save("./animation_correct/result.gif", dpi = 30, writer = PillowWriter(fps = 30))

print(" Created an Animation ! >:D ")

