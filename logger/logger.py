# coding=utf-8
import logging


def init_logger(name):
    """
    Sets up and returns a properly configured logging object
    :param name:
    :return:
    """
    surveyLogger = logging.getLogger(name)
    surveyLogger.setLevel(logging.DEBUG)
    # Prevents duplicate log entries after reinitialization.
    if not surveyLogger.handlers:
        fh = logging.FileHandler('log/' + name + '.log')
        fh.setLevel(logging.DEBUG)
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        FORMAT = '%(asctime)s - %(name)s - %(levelname)s: %(message)s'
        formatter = logging.Formatter(FORMAT, datefmt="%y-%m-%d %H:%M:%S")
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        surveyLogger.addHandler(fh)
        surveyLogger.addHandler(ch)

    return surveyLogger
