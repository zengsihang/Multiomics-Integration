r"""
GLUE (Graph-Linked Unified Embedding)
"""
# from os.path import dirname, basename, isfile, join
# import glob
# modules = glob.glob(join(dirname(__file__), "*.py"))
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and not f.endswith('__init__.py')]

try:
    from importlib.metadata import version
except ModuleNotFoundError:
    from pkg_resources import get_distribution
    version = lambda name: get_distribution(name).version

from . import data, genomics, graph, models, num, plot
from .utils import config, log


name = "scglue"
__version__ = version(name)
