import numpy as np
import scipy.signal
from utils.parameter_init import Config
import utils.plotting_methods
import utils.singal_utils
import matplotlib.pyplot as plt

def add_noise(ccf,percentage):
    max = np.max(np.abs(ccf))
    noise = np.random.normal(loc = 0, scale = max, size=(ccf.size)) 
    return ccf + noise * (percentage/100) 

def calculate_synthetic_ccf(model, distance, config):
    fmax = 1./config.min_period
    fmin = 1./config.max_period
    df = config.df
    #distance = config.distance
    dt = config.dt
    taper_length = config.taper_length
    ccf_spectra = synthetic_spectra(model,df,dt,distance)
    plot_synthetic = False
    t,ccf,freqs,dt,ccf_spectra = synthetic_ccf(ccf_spectra, df, dt, plot = False)

    """ taper = compute_taper(
        count = ccf.size,
        width = 1./dt * taper_length
    )
    ccf = ccf * taper """
    ccf = utils.singal_utils.downweight_ends(ccf,taper_length)
    
    #Padding with zeros. Single-sided CCF -> Double-sided CCF
    t,ccf = pad_zeros(t,ccf)

    if plot_synthetic:
        utils.plotting_methods.plot_synthetic(
            t = t,
            ccf=ccf,
            distance=distance,
            model=model
        )

    #Discarding redundant frequencies/periods from the ccf_spectra    
    ccf_spectra = ccf_spectra[(fmin <= freqs) & (freqs <= fmax)]
    freqs = freqs[(fmin <= freqs) & (freqs <= fmax)]
    phase = -np.angle(ccf_spectra)
    real_part = np.real(ccf_spectra)
    return [t,ccf,freqs,phase,real_part]

def compute_taper(count, width):
    result = np.ones(count)
    offset = int(min(count, int(width)))
    result[0:offset] = np.linspace(
        start = -1,
        stop = 1,
        num = offset
    )
    return (np.sin(result * np.pi / 2) + 1) / 2

def pad_zeros(t,ccf):
    l = t > 0
    z = np.zeros(t[l].size)
    ccf_double_sided = np.insert(
        arr = ccf,
        obj = 0,
        values = z
    )
    n_times = -np.flip(t[l])
    t_double_sided = np.insert(
        arr = t,
        obj = 0,
        values = n_times
    )
    return [t_double_sided, ccf_double_sided]

def synthetic_spectra(model, df, dt, distance):
    #pv = np.flip(model["phase_vel"].to_numpy())
    #freqs_pv = np.flip(model["freq"].to_numpy())
    pv = (model["phase_vel"].to_numpy())
    freqs_pv = (model["freq"].to_numpy())

    n = int(np.round(1./ (df*dt)))
    l = int(n//2 + 1)

    freqs = np.linspace(
        start = 1,
        stop = l,
        num = l
    )

    freqs = freqs * df

    pv_int = np.interp(
        x = freqs,
        xp = freqs_pv,
        fp = pv
    )
    
    spectra = scipy.special.hankel2(0,(freqs*2*np.pi / pv_int) * distance)
    spectra = np.pad(
        array = spectra,
        pad_width=(0,l - spectra.size),
        mode = "constant"
    )
    mask = np.isnan(spectra)
    spectra[mask] = 0
    return spectra

def synthetic_ccf(spectra, df, dt, plot = False):
    n = int(np.round(1./ (df*dt)))
    #l = int(n//2 + 1)

    spectra = np.append(
        arr= spectra,
        values = np.conjugate(np.flip(spectra[:-1]))
    )

    t = np.arange(
        start = 0,
        stop = n,
        step = 1
    ) * dt

    freqs = np.fft.fftfreq(t.size,dt)

    ccf = np.fft.ifft(spectra ,n )
    if plot:
        plt.plot(t,ccf)
        plt.show()
        plt.plot(freqs,np.angle(spectra))
        plt.show()
    return [t,np.real(ccf), freqs, dt, spectra]