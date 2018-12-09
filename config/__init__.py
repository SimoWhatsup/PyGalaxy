# coding=utf-8
import importlib
import os
import sys

# env is either dev (by default) or other (if set by environment variable)
env = os.environ.get('PYTHON_ENV', 'dev')
conf_name = 'config.{name}'.format(name=env)

try:
    config = importlib.import_module(conf_name)
except ImportError:
    sys.stderr.write("Error: Can't find the file '%s.py' in the directory." % conf_name)
    sys.exit(1)
