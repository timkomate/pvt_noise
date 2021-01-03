import plotting_methods
import parameter_init as param
import obspy.signal.tf_misfit
import pvt_noise_utils as pvt
import numpy as np
import pandas as pd
from timeit import default_timer as timer

def run(path):
    model = pd.read_csv(
        filepath_or_buffer = param.background_model,
        delimiter = " ",
        header = None,
        comment = "#"
    )
    model.columns = ["freq", "phase_vel"]
    t,ccf,distance,dt,nstack,n1,s1,n2,s2 = pvt.read_measured_data(path)
    ccf, t = pvt.calculate_simmetric_part(t,ccf)
    df = 1./((ccf.size) * dt)
    
    freqs = np.arange(
        start = 1./param.max_period,
        stop =  1./param.min_period,
        step= df
    )
    taper = pvt.compute_taper(
        count = ccf.size,
        width = 1./dt * param.taper_length
    )
    ccf = ccf * taper

    #Padding with zeros. Single-sided CCF -> Double-sided CCF
    t,ccf = pvt.pad_zeros(t,ccf)

    obspy.signal.tf_misfit.plot_tfr(
        st=ccf[(t>0) & (t < distance/1.5)],
        w0= 6,
        dt = 0.2,
        fmin = 1/200,
        fmax=1,
        mode="power"
    )

    c_branches = pvt.calculate(t,ccf,freqs,dt,distance,model)
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