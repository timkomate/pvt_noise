import logging
from datetime import datetime
import os

now = datetime.now()
dt_string = now.strftime("%Y-%m-%d-%H:%M:%S")

if not os.path.exists("./logs"):
            os.makedirs("./logs")

log_format = "%(asctime)s::%(filename)s::%(message)s"
logging.basicConfig(
    filename="./logs/{}.log".format(dt_string), 
    level='INFO', 
    format=log_format
)

logger = logging.getLogger("pvt_noise")