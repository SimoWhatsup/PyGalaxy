# coding=utf-8
import os
import logging

LOG_OUT_DIR = 'log/'


def init_logger(name):
    """
    Sets up and returns a properly configured logging object
    :param name:
    :return:
    """
    survey_logger = logging.getLogger(name)
    survey_logger.setLevel(logging.DEBUG)
    # Prevents duplicate log entries after reinitialization.
    if not survey_logger.handlers:
        if not os.path.isdir(LOG_OUT_DIR):
            os.makedirs(LOG_OUT_DIR)

        fh = logging.FileHandler(LOG_OUT_DIR + name + '.out')
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        FORMAT = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
        formatter = logging.Formatter(FORMAT, datefmt="%y-%m-%d %H:%M:%S")
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        survey_logger.addHandler(fh)
        survey_logger.addHandler(ch)

    return survey_logger
