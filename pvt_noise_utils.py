import scipy.special
import scipy.io
import scipy.signal
import numpy as np
import matplotlib.pyplot as plt
import time
import parameter_init as param
import pyfftw
import plotting_methods

def nextpow2(x):
    return 1<<(x-1).bit_length()

def add_noise(ccf,percentage):
    max = np.max(np.abs(ccf))
    noise = np.random.normal(loc = 0, scale = max, size=(ccf.size)) 
    return ccf + noise * (percentage/100) 

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

def save_synthetic_ccf(t,ccf,dt,distance):
    dd = {
        "Station1": "S1",
        "Station2": "S2",
        "Dist": distance,
        "dtnew": dt,
        "nstack": 1,
        "lagsx1x2": t,
        "cross12": ccf
    }
    
    scipy.io.savemat(
        "./synthetic.mat",
        dd
    )

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
        yfil = compute_gfilter(f,t,dt,gamma)
        #Convolve the CCF with the gaussian filter
        ccf_h = scipy.signal.convolve(
            in1 = ccf,
            in2 = yfil,
            mode = "full",
            method = "fft"
        )
        max_index = np.argmax(np.abs(ccf_h))
        
        #Setting up the weighting function
        weights = compute_window(
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

def calculate_synthetic_ccf(model, distance):
    fmax = 1./param.min_period
    fmin = 1./param.max_period
    df = param.df
    #distance = param.distance
    dt = param.dt
    taper_length = param.taper_length
    ccf_spectra = synthetic_spectra(model,df,dt,distance)
    plot_synthetic = False
    t,ccf,freqs,dt,ccf_spectra = synthetic_ccf(ccf_spectra, df, dt, plot = False)

    taper = compute_taper(
        count = ccf.size,
        width = 1./dt * taper_length
    )
    ccf = ccf * taper
    
    #Padding with zeros. Single-sided CCF -> Double-sided CCF
    t,ccf = pad_zeros(t,ccf)

    if plot_synthetic:
        plotting_methods.plot_synthetic(
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

def calculate(t,ccf,freqs,dt,distance,model):
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

    #print(gamma, gammaw)
    
    [p,a] = measure_phase(
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

    c_branches0 = high_freq_approx(
        phase = p,
        freqs = freqs,
        distance= distance,
        branch_num=branch_num,
        min_vel=min_vel,
        max_vel=max_vel
    )

    c_branches1 = local_search(
        c_original = c_branches0,
        p_measured = -p,
        freqs = freqs,
        cdiff = cdiff,
        distance= distance,
        max_freq = freq_hankel
    )

    freq_zeros, c_branches2 = real_part_crossings(
        amplitudes = a,
        freqs=freqs,
        distance=distance,
        branch_num=branch_num,
        min_vel=min_vel,
        max_vel=max_vel
    )

    center_branch1 = pick_closest_curve(
        model = model,
        c_branches= c_branches1,
        freqs = freqs,
        fmin = fmin
    )
    
    center_branch2 = pick_closest_curve(
        model = model,
        c_branches= c_branches2,
        freqs = freq_zeros,
        fmin = fmin
    )
    
    if param.plot:
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
        start= center_branch1-param.branch_to_save,
        stop = center_branch1+param.branch_to_save+1,
        step = 1
    )

    mask2 = np.arange(
        start= center_branch2-param.branch_to_save,
        stop = center_branch2+param.branch_to_save+1,
        step = 1
    )
    
    return c_branches0[mask1,:], c_branches1[mask1,:],c_branches2[mask2,:], freq_zeros

def real_part_crossings(amplitudes,freqs,distance,branch_num,min_vel,max_vel,plot = False):
    freq_zeros = get_zero_crossings(freqs,amplitudes)
    jn_zeros = scipy.special.jn_zeros(0,freq_zeros.size+branch_num)
    crossings = np.linspace(
        start = 2.40 - branch_num * 3.12,
        stop = 2.40-3.12,
        num= branch_num
    )
    jn_zeros = np.insert(jn_zeros, 0, crossings)
    
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

def get_zero_crossings(x,y):
    guess = np.where(np.diff(np.signbit(y)))[0]
    zeros = np.zeros(guess.shape)
    print(guess.size, zeros.size)
    i = 0
    for z in guess:
        xi = [x[z], x[z+1]]
        yi = [y[z], y[z+1]]
        a, b = np.polyfit(xi, yi, 1)
        zero = -b/a
        zeros[i] = zero
        i += 1
    return zeros

def read_measured_data(path):
    matfile = scipy.io.loadmat(
        file_name = path,
        squeeze_me = True
    )
    distance = matfile["Dist"] / 1000 #m to km
    dt = matfile["dtnew"]
    t = matfile["lagsx1x2"]
    ccf = matfile["cross12"]
    nstack = matfile["nstack"]
    s1 = matfile["Station1"]
    s2 = matfile["Station2"]
    try:
        cutvec = matfile["cutvec"]
        return [t,ccf,distance,dt,nstack,s1,s2,cutvec]
    except:
        n1 = matfile["Network1"]
        n2 = matfile["Network2"]
        return [t,ccf,distance,dt,nstack,n1,s1,n2,s2]
    
def save_pv(c_branches,distance,model,freqs,gamma,gammaw,
            lat1,lat2,lon1,lon2,nstack,filename):
    dd = {
        "crayan": c_branches,
        "Dist": distance,
        "frayan": freqs,
        "gamma1": gamma,
        "gammaw": gammaw,
        "lat1": lat1,
        "lat2": lat2,
        "lon1": lon1,
        "lon2": lon2,
        "nstack": nstack,
        "fs": model["freq"].to_numpy(),
        "pvs": model["phase_vel"].to_numpy()
    }
    
    scipy.io.savemat(
        file_name=filename,
        mdict=dd
    )

def calculate_simmetric_part(t,ccf):
    ccf_acausal = ccf[t <= 0]
    ccf_causal = ccf[t >= 0]
    ccf_simmetric = 1./2 * (ccf_causal + np.flip(ccf_acausal))
    t_simmetric = t[t >= 0]
    ccf_simmetric = scipy.signal.detrend(
        data = ccf_simmetric,
        type="linear"
    )
    ccf_simmetric = ccf_simmetric - np.mean(ccf_simmetric)
    if param.plot:
        fig, axs = plt.subplots(4)
        fig.suptitle("Original – causal – acausal – simmetric")
        axs[0].plot(t,ccf)
        axs[1].plot(t_simmetric,ccf_causal)
        axs[2].plot(-1*t_simmetric, ccf_acausal)
        axs[3].plot(t_simmetric, ccf_simmetric)
        plt.show()
    return ccf_simmetric,t_simmetric

def save_results_ascii(freqs,pv,filename):
    shape = pv.shape
    output = open(filename,"w")
    for i in np.arange(shape[0]):
        for j in np.arange(shape[1]):
            output.write("{} {}\n".format(freqs[j], pv[i,j]))
        output.write(">\n")
    output.close()
