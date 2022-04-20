import matplotlib.pyplot as plt
import numpy as np

def plot_circle(ax, x0, y0, r, **kwargs):
    th = np.linspace(0, 2*np.pi, 100)
    y = np.sin(th)*r + y0
    x = np.cos(th)*r + x0
    ax.plot(x, y, **kwargs)
