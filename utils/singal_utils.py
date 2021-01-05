import numpy as np
import matplotlib.pyplot as plt
import scipy.signal

def calculate_simmetric_part(t,ccf, plot = False):
    ccf_acausal = ccf[t <= 0]
    ccf_causal = ccf[t >= 0]
    ccf_simmetric = 1./2 * (ccf_causal + np.flip(ccf_acausal))
    t_simmetric = t[t >= 0]
    ccf_simmetric = scipy.signal.detrend(
        data = ccf_simmetric,
        type="linear"
    )
    ccf_simmetric = ccf_simmetric - np.mean(ccf_simmetric)
    if plot:
        fig, axs = plt.subplots(4)
        fig.suptitle("Original – causal – acausal – simmetric")
        axs[0].plot(t,ccf)
        axs[1].plot(t_simmetric,ccf_causal)
        axs[2].plot(-1*t_simmetric, ccf_acausal)
        axs[3].plot(t_simmetric, ccf_simmetric)
        plt.show()
    return ccf_simmetric,t_simmetric

def compute_taper(count, width):
    result = np.ones(count)
    offset = int(min(count, int(width)))
    result[0:offset] = np.linspace(
        start = -1,
        stop = 1,
        num = offset
    )
    return (np.sin(result * np.pi / 2) + 1) / 2

def compute_window(freq, value, count, index, time_step):
    omega = 2 * np.pi * freq
    alpha = np.power(value,2) * omega * time_step
    result = np.arange(
        start = 0,
        stop = count
    )
    result = result - index
    result = np.power(result * time_step, 2)
    result = result * np.power(omega, 2) / (4 * alpha)
    return np.exp(-result)

def compute_gfilter(f,t,time_step, value, plot = False):
    omega = 2 * np.pi * f
    alpha = np.power(value,2) * omega * time_step
    amp = omega / (2 * np.sqrt(np.pi * alpha))
    yfil = amp * np.exp(-omega*omega / (4 * alpha) * np.power(t,2))
    yfil = yfil * np.cos(omega * t)
    if plot:
        plt.plot(t, yfil)
        plt.show()
    return yfil

def get_zero_crossings(x,y):
    guess = np.where(np.diff(np.signbit(y)))[0]
    zeros = np.zeros(guess.shape)
    i = 0
    for z in guess:
        xi = [x[z], x[z+1]]
        yi = [y[z], y[z+1]]
        a, b = np.polyfit(xi, yi, 1)
        zero = -b/a
        zeros[i] = zero
        i += 1
    return zeros