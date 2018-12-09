# coding=utf-8
import json
import logging
import logging.config
import os

LOG_OUT_DIR = 'log/'


def setup_logging(default_path='config/logging.json', default_level=logging.INFO, env_key='LOG_CFG'):
    """
    Sets up logging configuration
    :param default_path:
    :param default_level:
    :param env_key:
    :return:
    """
    path = default_path
    log_cfg = os.getenv(env_key)
    python_env = os.environ.get('PYTHON_ENV', 'dev')

    if log_cfg:
        path = log_cfg

    if os.path.exists(path) and python_env != 'prod':
        with open(path, 'rt') as f:
            config = json.load(f)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(
            level=default_level,
            format='[%(asctime)s] %(filename)s:%(lineno)s - %(levelname)s %(message)s',
            datefmt='%y-%m-%d %H:%M:%S'
        )


def get_logger(name=None):
    """
    Sets up and returns a properly configured logging object
    :param name:
    :return:
    """
    if name is None:
        name = __name__

    logger = logging.getLogger(name)
    # common.setLevel(logging.DEBUG)
    return logger


class CustomFormatter(logging.Formatter):
    def __init__(self, fmt=None, datefmt="%y-%m-%d %H:%M:%S", use_color=True):
        super().__init__(fmt, datefmt=datefmt)
        self.fmt = fmt
        self.use_color = use_color
        self.colormap = {
            'WARNING': 35,  # Magenta
            'INFO': 32,  # Green
            'DEBUG': 34,  # Blue
            'CRITICAL': 41,  # Background Red
            'ERROR': 31  # Red
        }

    def format(self, record):
        level_name = record.levelname
        level_message = record.msg
        if self.use_color and level_name in self.colormap:
            level_name_color = "\033[1;{}m{}:\033[0m".format(self.colormap.get(level_name), level_name)
            level_message_color = "\033[1;{}m{}\033[0m".format(self.colormap.get(level_name), level_message)
            record.levelname = level_name_color
            record.msg = level_message_color

        return logging.Formatter.format(self, record)
