import plotting_methods
import obspy.signal.tf_misfit
import utils.pvt_utils as pvt
import utils.parameter_init
import utils.io_methods
import utils.synthetic_utils
import utils.singal_utils
import numpy as np
import pandas as pd
from timeit import default_timer as timer

def calculate(t,ccf,freqs,dt,distance,model,config):
    fmax = 1./config.min_period
    fmin = 1./config.max_period
    cdiff = config.cdiff
    branch_num = config.branch_num 
    freq_hankel = 1./config.h_period
    min_vel = config.min_vel
    max_vel = config.max_vel

    gamma_a,gamma_b = np.polyfit(config.gamma[0:2],config.gamma[2:],1)
    gammaw_a,gammaw_b = np.polyfit(config.gammaw[0:2],config.gammaw[2:],1)

    gamma = gamma_a*distance + gamma_b
    gammaw = gammaw_a*distance + gammaw_b

    #print(gamma, gammaw)
    
    [p,a] = pvt.measure_phase(
        freqs = freqs,
        t = t,
        ccf = ccf,
        dt = dt,
        distance= distance,
        gamma = gamma,
        gammaw= gammaw
    )

    c_branches0 = np.zeros(p.shape)
    c_branches1 = np.zeros(p.shape)
    c_branches2 = np.zeros(p.shape)

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
    
    if config.plot:
        plotting_methods.plot_results(
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
        start= center_branch1-config.branch_to_save,
        stop = center_branch1+config.branch_to_save+1,
        step = 1
    )

    mask2 = np.arange(
        start= center_branch2-config.branch_to_save,
        stop = center_branch2+config.branch_to_save+1,
        step = 1
    )
    
    return c_branches0[mask1,:], c_branches1[mask1,:],c_branches2[mask2,:], freq_zeros



def run(path):
    param = utils.parameter_init.Config("config.cfg")
    model = pd.read_csv(
        filepath_or_buffer = param.background_model,
        delimiter = " ",
        header = None,
        comment = "#"
    )
    model.columns = ["freq", "phase_vel"]
    t,ccf,distance,dt,nstack,n1,s1,n2,s2 = utils.io_methods.read_measured_data(path)
    ccf, t = utils.singal_utils.calculate_simmetric_part(t,ccf)
    df = 1./((ccf.size) * dt)
    
    freqs = np.arange(
        start = 1./param.max_period,
        stop =  1./param.min_period,
        step= df
    )
    taper = utils.singal_utils.compute_taper(
        count = ccf.size,
        width = 1./dt * param.taper_length
    )
    ccf = ccf * taper

    #Padding with zeros. Single-sided CCF -> Double-sided CCF
    t,ccf = utils.synthetic_utils.pad_zeros(t,ccf)

    obspy.signal.tf_misfit.plot_tfr(
        st=ccf[(t>0) & (t < distance/1.5)],
        w0= 6,
        dt = 0.2,
        fmin = 1/200,
        fmax=1,
        mode="power"
    )

    c_branches = calculate(t,ccf,freqs,dt,distance,model,param)
    #print(c_branches.shape,timer() - start)
    """pvt.save_pv(
        c_branches= c_branches,
        distance= distance,
        model = model,
        freqs = freqs,
        gamma=8,
        gammaw=21,
        lat1 = 

    )"""