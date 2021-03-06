import utils.plotting_methods
import obspy.signal.tf_misfit
import utils.pvt_utils as pvt
import utils.parameter_init
import utils.io_methods
import utils.synthetic_utils
import utils.singal_utils
from utils.setup_logger import logger
import utils.model_utils as mu
import utils.pvt_exceptions
import numpy as np
import pandas as pd
from timeit import default_timer as timer
import os

import matplotlib.pyplot as plt

def calculate(t,ccf,freqs,dt,distance,model,config):
    debug = False
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

    """center_branch1 = pvt.pick_closest_curve(
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
    
    mask1 = np.arange(
        start= center_branch1-config.branch_to_save,
        stop = center_branch1+config.branch_to_save+1,
        step = 1
    )

    mask2 = np.arange(
        start= center_branch2-config.branch_to_save,
        stop = center_branch2+config.branch_to_save+1,
        step = 1
    )"""

    #c_branches1[c_branches1 < min_vel] = np.nan
    #c_branches1[c_branches1 > max_vel] = np.nan
    if debug:
        fig1 = utils.plotting_methods.plot_pv(
            ccf = ([1],[1]),
            fmin=1./config.max_period,
            fmax=1./config.min_period,
            pv1 = (freqs, c_branches1),
            pv2 = (freq_zeros,c_branches2),
            distance=distance,
            #wlengths= np.array([1,2,3]),
            #title="{}.{}-{}.{}: {:.2f}km".format(n1,s1,n2,s2,distance),
            #bg_model=bg_model
        )
        plt.show()
    
    return c_branches0, c_branches1,c_branches2, freq_zeros



def run(path):
    start = timer()
    try:
        print("working on: {}".format(path))
        param = utils.parameter_init.Config("config.cfg")

        gamma_a,gamma_b = np.polyfit(param.gamma[0:2],param.gamma[2:],1)
        gammaw_a,gammaw_b = np.polyfit(param.gammaw[0:2],param.gammaw[2:],1)

        wlength_a,wlength_b = np.polyfit(param.wlength[0:2],param.wlength[2:],1)
        
        data = utils.io_methods.read_json("stationsall.json")

        [t,ccf,distance,dt,nstack,n1,s1,n2,s2] = utils.io_methods.read_measured_data(path)

        if(distance < param.min_distance):
            raise utils.pvt_exceptions.StationsTooClose(n1,s1,n2,s2)

        lon1 = data[n1][s1]["longitude"]
        lat1 = data[n1][s1]["latitude"]
        lon2 = data[n2][s2]["longitude"]
        lat2 = data[n2][s2]["latitude"]

        pattern = "eigen_lon{}_lat{}_S.eigen"
        points = mu.find_points(lon1,lat1, lon2,lat2)
        bg_model = mu.average_model(points, param.model_path, pattern)
        
        gamma = gamma_a*distance + gamma_b
        gammaw = gammaw_a*distance + gammaw_b
        wlength = wlength_a*distance + wlength_b
        wlength = int(wlength/dt)

        

        fname = "pv_{}_{}_{}_{}_{:.2f}km_{}.mat".format(n1,s1,n2,s2,distance,nstack)
        folder = "{}/{}-{}/".format(param.save_path,s1,s2)
        full_name = "{}/{}".format(folder,fname)
        
        if(not param.overwrite and os.path.isfile(full_name)):
            raise utils.pvt_exceptions.DataExcist(full_name)
                
        ccf, t = utils.singal_utils.calculate_simmetric_part(t,ccf,param.ccf_part)
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

        _, c_branches1, c_branches2, freq_zeros = calculate(t,ccf,freqs,dt,distance,bg_model,param)

        if(param.plot):
            fig1 = utils.plotting_methods.plot_pv(
                pv1 = (freqs, c_branches1),
                pv3= (freq_zeros,c_branches2),
                distance=distance,
                title="{}.{}-{}.{}: {:.2f}km".format(n1,s1,n2,s2,distance)
            )
            plt.figure.savefig("./figg")
        
        model_name = "{}/{}-{}/{}-{}.xy".format(param.save_path,s1,s2,s1,s2)
        
        if not os.path.exists(folder):
            os.makedirs(folder)
        
        utils.io_methods.save_pv_format(
            filename = full_name,
            c_branches= c_branches1,
            distance= distance,
            model = bg_model,
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
        logger.info("{}::{:.2f}::{}::{}".format(path.split("/")[-1],distance, timer() - start, 1))
    except utils.pvt_exceptions.DataExcist as e:
        print(e)
        logger.info("{}::{:.2f}::{}::{}".format(path.split("/")[-1],distance, timer() - start, 2))
        return
    except utils.pvt_exceptions.StationsTooClose as e:
        print(e)
        logger.info("{}::{:.2f}::{}::{}".format(path.split("/")[-1],distance, timer() - start, 3))
        return
    except IndexError as e:
        print(e)
        logger.info("{}::{:.2f}::{}::{}".format(path.split("/")[-1],distance, timer() - start, 4))
        return
    except:
        logger.info("{}::{:.2f}::{}::{}".format(path.split("/")[-1],distance, timer() - start, 5))
        print("Unknown error at:{}".format(path))
        return