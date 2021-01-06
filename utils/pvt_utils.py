import scipy.special
import scipy.io
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
import time
import pyfftw
import utils.plotting_methods
import utils.singal_utils

def pick_closest_curve(model,c_branches,freqs,fmin):
    model_frequencies = model["freq"].to_numpy()
    fmodel_index = np.argwhere(model_frequencies >= fmin)[0]
    fmeasured_index = np.argmin(freqs >= fmin)
    pv_model = model["phase_vel"].to_numpy()[fmodel_index]
    pv_measured = c_branches[:,fmeasured_index]
    difference = np.abs(pv_model - pv_measured)
    difference[np.isnan(difference)] = np.inf
    branch_idex = np.argmin(difference)
    return branch_idex

def local_search(c_original, p_measured, freqs, cdiff, distance, max_freq = None):
    c = np.copy(c_original)
    if max_freq is None:
        max_freq = np.inf
    branches_num = c.shape[0]
    freqs_num = c.shape[1]
    #looping through the branches:
    for i in np.arange(branches_num):
        #looping through the frequencies
        for j in np.arange(freqs_num):
            if freqs[j] > max_freq: 
                break
            if (np.isnan(c[i,j])):
                continue
            j_plus = j+1 if j < branches_num else branches_num
            j_minus = j-1 if j < 0 else 0
            c_up = c[i,j_plus]
            c_down = c[i,j_minus]
            #print(c_up,c_down)
            max_diff = np.min([c_up,c_down]) * 0.5
            c[i,j] = calculate_velocity(p_measured[j],freqs[j],distance,c[i,j],cdiff,max_diff)
    return c

def get_direction(phase,frequency,distance,c_guess,cdiff):
    c = np.array([c_guess - cdiff, c_guess+cdiff])
    H = scipy.special.hankel2(0,(frequency*2*np.pi / c * distance))
    Hp = -np.angle(H)
    differences = np.abs(Hp-phase)
    derivative = (differences[0] - differences[-1]) / (2*cdiff)
    return np.sign(derivative)

def calculate_velocity(phase, frequency, distance, c_guess, cdiff,max_diff):
    if np.isnan(max_diff):
        max_iteration = 1000
    else:
        max_iteration = max_diff/cdiff
    #print(max_iteration)
    i = 0
    d = np.inf
    c0 = c_guess
    c1 = c_guess
    c_original = c_guess
    if np.abs(phase) > np.pi:
        phase += np.pi * 2
    direction = get_direction(phase, frequency,distance,c0,cdiff)
    while True:
        H = scipy.special.hankel2(0,(frequency*2*np.pi / c1 * distance))
        Hp = -np.angle(H)
            
        difference = np.abs(Hp - phase)       
        if d < difference:
            return c0
        d = difference
        c0 = c1
        c1 = c1 + direction * cdiff
        i += 1
        #print(max_iteration,i,cdiff)
        if i > max_iteration:
            return c_original
    return c0

def measure_phase(freqs, t,ccf,dt,distance,gamma,gammaw):
    p = np.zeros(
        shape = freqs.shape
    )
    a = np.zeros(
        shape = freqs.shape
    )
    i = 0
    #df = 1./((ccf.size) * dt)
    for f in freqs:
        pyfftw.interfaces.cache.enable()
        #Constructing gauss filter
        yfil = utils.singal_utils.compute_gfilter(f,t,dt,gamma)
        #Convolve the CCF with the gaussian filter
        ccf_h = scipy.signal.convolve(
            in1 = ccf,
            in2 = yfil,
            mode = "full",
            method = "fft"
        )
        max_index = np.argmax(np.abs(ccf_h))
        
        #Setting up the weighting function
        weights = utils.singal_utils.compute_window(
            freq = f, 
            value = gammaw,
            count = ccf_h.size, 
            index = max_index, 
            time_step = dt
        )
        ccf_w = ccf_h * weights
        ccf_w = ccf_w[ccf.size//2 - 1:-1*ccf.size//2]
        weights = weights[ccf.size//2:-1*ccf.size//2 + 1]
        
        if False:
            #fig = plt.figure(1)
            plt.clf()
            plt.plot(t,ccf / np.max(np.abs(ccf)))
            plt.plot(t,weights / np.max(np.abs(weights)))
            plt.plot(t,ccf_w / np.max(np.abs(ccf_w)))
            plt.xlim([-1*distance,distance])
            plt.show(block = False)
            plt.pause(0.05)
        
        #measure the phase
        ccf_w = ccf_w
        spectra = pyfftw.interfaces.numpy_fft.fft(ccf_w)
        fft_freq = np.fft.fftfreq(ccf_w.size,dt)
        #freq_index = int(np.round(f/df))
        freq_index = np.argmin(np.abs(fft_freq - f))
        p[i] = np.angle(spectra[(freq_index) - 1])  + np.pi
        a[i] = np.real(spectra[(freq_index) - 1])
        #print("{} {} {}".format(f,p[i],freq_index - 1))
        #p[i] = np.angle(spectra[(2*freq_index)])  + np.pi
        i += 1
    return [p, a]

def high_freq_approx(phase,freqs,distance,branch_num,min_vel,max_vel):
    p = -np.unwrap(phase)
    p_branches = np.tile(p , reps = (branch_num,1))
    c_branches = np.tile(2*np.pi*distance*freqs , reps = (branch_num,1))
    branches = np.arange(
        start = -(branch_num-1),
        stop = (branch_num+1),
        step = 2) * np.pi
    p_branches = (p_branches.transpose() - branches).transpose()
    p_branches = p_branches + (np.pi / 4)
    c_branches = c_branches / p_branches
    c_branches[c_branches < min_vel] = np.nan
    c_branches[c_branches > max_vel] = np.nan
    return c_branches

def real_part_crossings(amplitudes,freqs,distance,branch_num,min_vel,max_vel,plot = False):
    freq_zeros = utils.singal_utils.get_zero_crossings(freqs,amplitudes)
    jn_zeros = scipy.special.jn_zeros(0,freq_zeros.size+branch_num)

    crossings = jn_zeros[0:branch_num]

    jn_zeros = np.insert(jn_zeros, 0, np.flip(-crossings))
    
    c_zeros = np.zeros((branch_num,freq_zeros.size))

    adjust=0

    for i in np.arange(branch_num):
        
        c_zeros[i,:] = (2*np.pi*freq_zeros*distance)/jn_zeros[2*i-adjust:freq_zeros.size + 2*i-adjust]
        adjust = 1
    
    c_zeros[c_zeros < min_vel] = np.nan
    c_zeros[c_zeros > max_vel] = np.nan

    
    if plot:
        plt.plot(freqs, amplitudes)
        plt.show()

    return [freq_zeros, c_zeros]



