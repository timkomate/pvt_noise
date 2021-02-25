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
    pool = multiprocessing.Pool(processes=param.cpus)
    files = glob.glob("{}*.mat".format(param.input_path))
    pool.map(utils.pvt_driver.run, files)

#TODO
#Remove branches_to_save variable: OK
#