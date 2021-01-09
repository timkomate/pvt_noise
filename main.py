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
    pool = multiprocessing.Pool(processes=4)
    #files = ["./input_data/CCF_GR_MOX_HU_PSZ_ZZ_2377_WN.mat","./input_data/CCF_GR_MOX_HU_SOP_ZZ_722_WN.mat", "./input_data/CCF_GR_MOX_HU_TRPA_ZZ_1006_WN.mat", "./input_data/CCF_HU_PSZ_CR_ZAG_ZZ_441_WN.mat"]
    files = glob.glob("{}*.mat".format(param.input_path))
    print(files)
    pool.map(utils.pvt_driver.run, files)

#TODO