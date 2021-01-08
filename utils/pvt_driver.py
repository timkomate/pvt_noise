import utils.plotting_methods
import obspy.signal.tf_misfit
import utils.pvt_utils as pvt
import utils.parameter_init
import utils.io_methods
import utils.synthetic_utils
import utils.singal_utils
import utils.model_utils as mu
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
    try:
        param = utils.parameter_init.Config("config.cfg")

        dt = param.dt
        distance = param.distance

        gamma_a,gamma_b = np.polyfit(param.gamma[0:2],param.gamma[2:],1)
        gammaw_a,gammaw_b = np.polyfit(param.gammaw[0:2],param.gammaw[2:],1)

        gamma = gamma_a*distance + gamma_b
        gammaw = gammaw_a*distance + gammaw_b

        wlength_a,wlength_b = np.polyfit(param.wlength[0:2],param.wlength[2:],1)
        wlength = wlength_a*distance + wlength_b
        wlength = int(wlength/dt)

        model = pd.read_csv(
            filepath_or_buffer = param.background_model,
            delimiter = " ",
            header = None,
            comment = "#"
        )
        model.columns = ["freq", "phase_vel"]
        [t,ccf,distance,dt,nstack,n1,s1,n2,s2] = utils.io_methods.read_measured_data(path)
        ccf, t = utils.singal_utils.calculate_simmetric_part(t,ccf)
        df = 1./((ccf.size) * dt)
        
        freqs = np.arange(
            start = 1./param.max_period,
            stop =  1./param.min_period,
            step= df
        )
        
        ccf = utils.singal_utils.downweight_ends(
            data = ccf,
            wlength = wlength
        )
        
        #Padding with zeros. Single-sided CCF -> Double-sided CCF
        t,ccf = utils.synthetic_utils.pad_zeros(t,ccf)

        if param.plot:
            obspy.signal.tf_misfit.plot_tfr(
                st=ccf[(t>0) & (t < distance/1.5)],
                w0= 6,
                dt = 0.2,
                fmin = 1/200,
                fmax=1,
                mode="power"
            )

        c_branches0, c_branches1, c_branches2, freq_zeros = calculate(t,ccf,freqs,dt,distance,model,param)
        data = utils.io_methods.read_json("stationsall.json")
        
        lon1 = data[n1][s1]["longitude"]
        lat1 = data[n1][s1]["latitude"]

        lon2 = data[n2][s2]["longitude"]
        lat2 = data[n2][s2]["latitude"]

        pattern = "eigen_lon{}_lat{}_S.eigen"
        points = mu.find_points(lon1,lat1, lon2,lat2)
        bg_model = mu.average_model(points, param.model_path, pattern)
        model_name = "{}/{}-{}/{}-{}.xy".format(param.save_path,s1,s2,s1,s2)
        
        #print(c_branches.shape,timer() - start)
        fname = "pv_{}_{}_{}_{}_{:.2f}km_{}.mat".format(n1,s1,n2,s2,distance,nstack)
        utils.io_methods.save_pv_format(
            save_path = param.save_path,
            filename = fname,
            c_branches= c_branches1,
            distance= distance,
            model = model,
            freqs = freqs,
            gamma=gamma,
            gammaw=gammaw,
            lat1 = lat1,
            lat2= lat2,
            lon1 = lon1,
            lon2 = lon2,
            nstack = nstack,
            n1 = n1,
            s1 = s1,
            n2 = n2,
            s2 = s2,
            c_zeros = c_branches2,
            f_zeros = freq_zeros
        )
        utils.io_methods.save_bg_model(model_name,bg_model)
    except:
        return