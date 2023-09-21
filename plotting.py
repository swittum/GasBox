import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

print('Entering Python Script')

fig, ax = plt.subplots(figsize=(5, 5))

with h5py.File('./pos_animation.h5', 'r') as file:
    animation = np.array(file['pos_animation'])

    def animate(i):
        circles = [plt.Circle(particle, radius=.01) for particle in animation[i]]
        ax.clear()
        for circle in circles:
            ax.add_patch(circle)
        return circles

    ani = FuncAnimation(fig, animate, frames=4000, interval=.01, blit=True)
    writer = FFMpegWriter(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    ani.save("animation.mp4", writer=writer)

    

    