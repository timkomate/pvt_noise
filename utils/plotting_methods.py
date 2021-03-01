import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def plot_synthetic(t,ccf,distance,model, xlim = None):
    fig, axs = plt.subplots(2)
    fig.suptitle("Station separation: {:.2f} km".format(distance))
    axs[0].plot(t,ccf)
    if xlim != None:
        axs[0].set_xlim(xlim)
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("Amplitude")
    axs[1].plot(model["freq"], model["phase_vel"])
    axs[1].set_xlabel("Frequency [Hz]")
    axs[1].set_ylabel("Phase velocity [km/s]")
    plt.show()

def plot_results(freqs,c_branches0,distance,title,
            model = None,c_branches1 = None,c_branches2 = None, freq_zeros = None, wlengths = np.array([]),
            fmin = None, fmax = None):
    legend = []
    if c_branches2 is None:
        plt.plot(freqs, c_branches0.transpose(), c = "red")
        legend.append("Measured phase velocites")
    else:
        plt.plot(freqs, c_branches0.transpose(), c = "red")
        plt.plot(freqs, c_branches1.transpose(), c = "gray")
        plt.plot(freq_zeros, c_branches2.transpose(), c = "green")
        legend.append("Initial phase velocities")
        legend.append("Corrected phase velocities")
    if model is not None:
        plt.plot(model["freq"], model["phase_vel"], 
            linestyle = "dashed")
    if wlengths.size:
        plt.vlines(
            x = wlengths*4/distance,
            ymin = 0,
            ymax = 8, 
            colors = "black", 
            linestyles = "dashed",
        )
        for i in wlengths:
            legend.append("{} wavelength".format(i))
    #plt.legend(legend)
    plt.ylim([0,8])
    plt.xlim([fmin,fmax])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase velocity [km/s]")
    plt.title(title)
    plt.show()

def plot_pv_period(ccf,pv1, distance, pv2= None, pv3 = None, title = "", bg_model = None,
            fmin = None, fmax = None, wlengths = np.array([]),
            s1 = None, s2 = None):
    fig = plt.figure(figsize = [16,9], dpi = 100)

    if(s1 != None and s2 != None):
        ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, fig = fig)
        ax2 = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan = 2, fig = fig)
        ax3 = plt.subplot2grid((3,3), (1, 2), rowspan=2, fig = fig)

        lons = [s1[0], s2[0]]
        lats = [s1[1], s2[1]]

        map = Basemap(llcrnrlon=8.,llcrnrlat=37.5,urcrnrlon=32.,urcrnrlat=54.5,
                resolution='i', projection='merc', ax = ax3)
        map.drawmapboundary(fill_color='gray')
        map.fillcontinents(color='white')
        map.drawcountries()
        map.drawcoastlines()
        x,y = map(lons,lats)
        map.plot(x, y, color='black')
        map.plot(x, y, marker='v',color='green')
    else:
        ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, fig = fig)
        ax2 = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan = 3, fig = fig)
        ax3 = None
    
    #fig.set_size_inches(10, 15)
    ax1.plot(ccf[0], ccf[1], 'black')
    ax2.plot(1/pv1[0], pv1[1].transpose(), 'gray')
    if(pv2 != None):
        ax2.plot(1/pv2[0], pv2[1].transpose(), 'red')
    if(pv3 != None):
        ax2.plot(pv3[0], pv3[1].transpose(), "green")
    if(type(bg_model) != type(None)):
        ax2.plot(1/bg_model["freq"], bg_model["phase_vel"], "blue",
            linestyle = "dashed")
    if wlengths.size:
        ax2.vlines(
            x = distance/(wlengths*4),
            ymin = 0,
            ymax = 8, 
            colors = "black", 
            linestyles = "dashed",
        )
    
    ax1.set_xlim([0,np.max(ccf[0])])
    ax1.set_title(title)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Amplitude")
    
    ax2.set_ylim([0,8])
    ax2.set_xlim([1/fmax,1/fmin])
    ax2.set_xlabel("Period [s]")
    ax2.set_ylabel("Phase velocity [km/s]")
    return fig

def plot_pv(ccf,pv1, distance, pv2= None, pv3 = None, title = "", bg_model = None,
            fmin = None, fmax = None, wlengths = np.array([]),
            s1 = None, s2 = None):
    fig = plt.figure(figsize = [16,9], dpi = 100)

    if(s1 != None and s2 != None):
        ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, fig = fig)
        ax2 = plt.subplot2grid((3,3), (1,0), colspan=2, rowspan = 2, fig = fig)
        ax3 = plt.subplot2grid((3,3), (1, 2), rowspan=2, fig = fig)

        lons = [s1[0], s2[0]]
        lats = [s1[1], s2[1]]

        map = Basemap(llcrnrlon=8.,llcrnrlat=37.5,urcrnrlon=32.,urcrnrlat=54.5,
                resolution='i', projection='merc', ax = ax3)
        map.drawmapboundary(fill_color='gray')
        map.fillcontinents(color='white')
        map.drawcountries()
        map.drawcoastlines()
        x,y = map(lons,lats)
        map.plot(x, y, color='black')
        map.plot(x, y, marker='v',color='green')
    else:
        ax1 = plt.subplot2grid((3,3), (0,0), colspan=3, fig = fig)
        ax2 = plt.subplot2grid((3,3), (1,0), colspan=3, rowspan = 3, fig = fig)
        ax3 = None
    
    #fig.set_size_inches(10, 15)
    ax1.plot(ccf[0], ccf[1], 'black')
    ax2.plot(pv1[0], pv1[1].transpose(), 'gray')
    if(pv2 != None):
        ax2.plot(pv2[0], pv2[1].transpose(), 'red')
    if(pv3 != None):
        ax2.plot(pv3[0], pv3[1].transpose(), "green")
    if(type(bg_model) != type(None)):
        ax2.plot(bg_model["freq"], bg_model["phase_vel"], "blue",
            linestyle = "dashed")
    if wlengths.size:
        ax2.vlines(
            x = wlengths*4/distance,
            ymin = 0,
            ymax = 8, 
            colors = "black", 
            linestyles = "dashed",
        )
    
    ax1.set_xlim([0,np.max(ccf[0])])
    ax1.set_title(title)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Amplitude")
    
    ax2.set_ylim([0,8])
    ax2.set_xlim([fmin,fmax])
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Phase velocity [km/s]")
    return fig


def plot_pv_diff(freqs,c_branches,distance,title,
            model = None,c_branches2 = None, wlengths = np.array([]),
            fmin = None, fmax = None, branch = None):
    if branch is None:
        branch = c_branches.shape[0] // 2
    plt.plot(freqs, c_branches[branch] - c_branches2[branch], c = "red")
    plt.show()