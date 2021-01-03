import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_synthetic(t,ccf,distance,model):
    fig, axs = plt.subplots(2)
    fig.suptitle("Station separation: {:.2f} km".format(distance))
    axs[0].plot(t,ccf)
    axs[0].set_xlim([-1*distance,10*distance])
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

def plot_pv_diff(freqs,c_branches,distance,title,
            model = None,c_branches2 = None, wlengths = np.array([]),
            fmin = None, fmax = None, branch = None):
    if branch is None:
        branch = c_branches.shape[0] // 2
    plt.plot(freqs, c_branches[branch] - c_branches2[branch], c = "red")
    plt.show()

    