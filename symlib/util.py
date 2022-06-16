import matplotlib.pyplot as plt
import numpy as np

def plot_circle(ax, x0, y0, r, **kwargs):
    th = np.linspace(0, 2*np.pi, 100)
    y = np.sin(th)*r + y0
    x = np.cos(th)*r + x0
    ax.plot(x, y, **kwargs)

def set_units_x(x, h, scale, param):
    dx = np.zeros(x.shape)
    for dim in range(3): dx[:,dim] = x[:,dim] - h["x"][dim]
    dx *= scale*1e3/param["h100"]
    return dx

def set_units_v(v, h, scale, param):
    v *= np.sqrt(scale)
    dv = np.zeros(v.shape)
    for dim in range(3): dv[:,dim] = v[:,dim] - h["v"][dim]
    return dv

def set_units_param(scale, param):
    return param["mp"]/param["h100"], param["eps"]*scale/param["h100"]

def set_units_halos(h, scale, param):
    h = np.copy(h)    
    for dim in range(3):
        x0 = np.copy(h["x"][0,:,dim])
        v0 = np.copy(h["v"][0,:,dim])

        for hi in range(len(h)):
            h["x"][hi,:,dim] -= x0
            h["v"][hi,:,dim] -= v0
            if hi == 0:
                h["x_core"][hi,:,dim] = 0
                h["v_core"][hi,:,dim] = 0
            else:
                h["x_core"][hi,:,dim] -= x0
                h["v_core"][hi,:,dim] -= v0
                
            h["x"][hi,:,dim] *= 1e3*scale/param["h100"]
            h["x_core"][hi,:,dim] *= 1e3*scale/param["h100"]
            
            if dim == 0:
                h["rvir"][hi,:] *= 1e3*scale/param["h100"]
                h["mvir"][hi,:] *= 1/param["h100"]

    
    for hi in range(len(h)):
        invalid = h[hi]["rvir"] < 0
        h[hi, invalid]["rvir"] = -1
        h[hi, invalid]["x"] = -1
        h[hi, invalid]["v"] = -1
        h[hi, invalid]["x_core"] = -1
        h[hi, invalid]["v_core"] = -1
        
    return h
