import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import obspy.signal.tf_misfit
import scipy.io
import scipy.signal
import scipy.fftpack
import multiprocessing
import pyfftw
from timeit import default_timer as timer

import utils.plotting_methods
import utils.pvt_driver
import utils.singal_utils
import utils.synthetic_utils
import utils.parameter_init 
import utils.io_methods
import utils.pvt_utils as pvt
import glob

if __name__ == "__main__":
    param = utils.parameter_init.Config("config.cfg")
    
    dt = param.dt
    distance = param.distance

    fmax = 1./param.min_period
    fmin = 1./param.max_period
    cdiff = param.cdiff
    branch_num = param.branch_num 
    freq_hankel = 1./param.h_period
    min_vel = param.min_vel
    max_vel = param.max_vel

    gamma_a,gamma_b = np.polyfit(param.gamma[0:2],param.gamma[2:],1)
    gammaw_a,gammaw_b = np.polyfit(param.gammaw[0:2],param.gammaw[2:],1)

    gamma = gamma_a*distance + gamma_b
    gammaw = gammaw_a*distance + gammaw_b

    wlength_a,wlength_b = np.polyfit(param.wlength[0:2],param.wlength[2:],1)
    wlength = wlength_a*distance + wlength_b
    wlength = int(wlength/dt)
    """model = pd.read_csv(
        filepath_or_buffer = param.background_model,
        delimiter = " ",
        header = None,
        comment = "#",
        names=["mode", "period", "phase_vel", "group_vel"]
    )
    model["freq"] = 1./model["period"]"""
    
    
    model = pd.read_csv(
        filepath_or_buffer = param.background_model,
        delimiter = " ",
        header = None,
        comment = "#",
        names=["freq","phase_vel"]
    )
    t,ccf,freqs,phase,amp = utils.synthetic_utils.calculate_synthetic_ccf(model, distance, param)
    
    ccf = utils.singal_utils.downweight_ends(ccf,wlength)
    #ccf = pvt.add_noise(ccf,1)
    utils.plotting_methods.plot_synthetic(t,ccf,distance,model)
    obspy.signal.tf_misfit.plot_tfr(
        st=ccf[(t>0) & (t < distance/1.5)],
        w0= 6,
        dt = 0.2,
        fmin = 1/200,
        fmax=1,
        mode="power"
    )

    [p,a] = pvt.measure_phase(
        freqs = freqs,
        t = t,
        ccf = ccf,
        dt = dt,
        distance= distance,
        gamma = gamma,
        gammaw= gammaw
    )

    c_branches0 = pvt.high_freq_approx(
        phase = p,
        freqs = freqs,
        distance= distance,
        branch_num=branch_num,
        min_vel=min_vel,
        max_vel=max_vel
    )

    c_branches1 = pvt.local_search(
        c_original = c_branches0,
        p_measured = -p,
        freqs = freqs,
        cdiff = cdiff,
        distance= distance,
        max_freq = freq_hankel
    )

    freq_zeros, c_branches2 = pvt.real_part_crossings(
        amplitudes = a,
        freqs=freqs,
        distance=distance,
        branch_num=branch_num,
        min_vel=min_vel,
        max_vel=max_vel
    )
    
    center_branch1 = pvt.pick_closest_curve(
        model = model,
        c_branches= c_branches1,
        freqs = freqs,
        fmin = fmin
    )
    
    center_branch2 = pvt.pick_closest_curve(
        model = model,
        c_branches= c_branches2,
        freqs = freq_zeros,
        fmin = fmin
    )
    
    if param.plot:
        utils.plotting_methods.plot_results(
            freqs=freqs,
            c_branches0 = c_branches0,
            title = "station separation: {:.2f} km".format(distance),
            distance = distance,
            model= model,
            c_branches1=c_branches1,
            c_branches2=c_branches2,
            freq_zeros=freq_zeros,
            wlengths=np.array([1,2,3]),
            fmin = fmin,
            fmax = fmax
        )
    
    mask1 = np.arange(
        start= center_branch1-param.branch_to_save,
        stop = center_branch1+param.branch_to_save+1,
        step = 1
    )

    mask2 = np.arange(
        start= center_branch2-param.branch_to_save,
        stop = center_branch2+param.branch_to_save+1,
        step = 1
    )

    c0, c1, c2 = c_branches0[mask1,:], c_branches1[mask1,:],c_branches2[mask2,:]
    
    utils.io_methods.save_results_ascii(
        freqs=freqs,
        pv = c1,
        filename = "ascii_data/hankel{:.0f}.xy".format(distance)
    )

    utils.io_methods.save_results_ascii(
        freqs=freqs,
        pv = c0,
        filename = "ascii_data/high_freq{:.0f}.xy".format(distance)
    )

    utils.io_methods.save_results_ascii(
        freqs=freq_zeros,
        pv = c2,
        filename = "ascii_data/zero_crossings{:.0f}.xy".format(distance)
    )

    output1 = open("ascii_data/real_part{:.0f}.xy".format(distance),"w")
    for i in np.arange(freqs.size):
        output1.write("{} {}\n".format(freqs[i],a[i]))
    output1.close()
    
    output2 = open("ascii_data/crossings{:.0f}.xy".format(distance), "w")
    for i in freq_zeros:
        output2.write("{} 0\n".format(i))
    output2.close()

    output3 = open("ascii_data/bensen{:.0f}.xy".format(distance), "w")
    bensen_criterion= (3*np.nanmean(model["phase_vel"]))/distance # in frequency
    output3.write("{} {}\n".format(bensen_criterion,min_vel))
    output3.write("{} {}".format(bensen_criterion,max_vel))
    output3.close()