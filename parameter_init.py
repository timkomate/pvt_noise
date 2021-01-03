import configparser

config = configparser.ConfigParser()
config.read("./config.cfg")

#Data case
real_data_case = config.getboolean("CASE", "real_data_case")
plot = config.getboolean("CASE", "plot")

#Sythetic
background_model = config.get("SYNTHETIC", "background_model")
dt = config.getfloat("SYNTHETIC", "dt")
df = config.getfloat("SYNTHETIC", "df")
distance = config.getfloat("SYNTHETIC", "distance")

#General
max_period = config.getfloat("GENERAL", "max_period")
min_period =  config.getfloat("GENERAL", "min_period")
branch_num = config.getint("GENERAL", "branch_num")
branch_to_save = config.getint("GENERAL", "branch_to_save")
taper_length = config.getfloat("GENERAL", "taper_length")
h_period = config.getfloat("GENERAL", "h_period")
min_vel = config.getfloat("GENERAL", "min_vel")
max_vel = config.getfloat("GENERAL", "max_vel")
cdiff = config.getfloat("GENERAL", "cdiff")
gamma = list(map(float,config.get("GENERAL", "gamma").split(",")))
gammaw = list(map(float,config.get("GENERAL", "gammaw").split(",")))