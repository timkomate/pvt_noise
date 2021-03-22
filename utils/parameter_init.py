import configparser

class Config(object):
    def __init__(self,config_name):
        config = configparser.ConfigParser()
        config.read(config_name)
        #Data case
        self.plot = config.getboolean("CASE", "plot")
        self.overwrite = config.getboolean("CASE", "overwrite")
        self.save_plots = config.getboolean("CASE", "save_plots")
        self.fileformat = config.get("CASE", "fileformat")

        #Sythetic
        self.background_model = config.get("SYNTHETIC", "background_model")
        self.dt = config.getfloat("SYNTHETIC", "dt")
        self.df = config.getfloat("SYNTHETIC", "df")
        self.distance = config.getfloat("SYNTHETIC", "distance")
        self.add_noise = config.getboolean("SYNTHETIC", "add_noise")
        self.noise = config.getfloat("SYNTHETIC", "noise")
        #General
        self.ccf_part = config.get("GENERAL", "ccf_part")
        self.cpus = config.getint("GENERAL", "cpus")
        self.max_period = config.getfloat("GENERAL", "max_period")
        self.min_period =  config.getfloat("GENERAL", "min_period")
        self.branch_num = config.getint("GENERAL", "branch_num")
        self.h_period = config.getfloat("GENERAL", "h_period")
        self.min_vel = config.getfloat("GENERAL", "min_vel")
        self.max_vel = config.getfloat("GENERAL", "max_vel")
        self.min_distance = config.getfloat("GENERAL", "min_distance")
        self.max_distance = config.getfloat("GENERAL", "max_distance")
        self.cdiff = config.getfloat("GENERAL", "cdiff")
        self.gamma = list(map(float,config.get("GENERAL", "gamma").split(",")))
        self.gammaw = list(map(float,config.get("GENERAL", "gammaw").split(",")))
        self.wlength = list(map(float,config.get("GENERAL", "wlength").split(",")))
        #PATHS
        self.input_path = config.get("PATHS", "input_path")
        self.save_path = config.get("PATHS", "save_path")
        self.model_path = config.get("PATHS", "model_path")