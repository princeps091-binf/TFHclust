# %%
# run from project root (ease of module access)

import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
from matplotlib.cm import viridis
import os
from importlib import reload
import src.detect_cluster
reload(src.detect_cluster)
import src.integrate_replicate
reload(src.integrate_replicate)

import matplotlib.patches as pat  # Patches like pat.Polygon()
from matplotlib.collections import PolyCollection  # Collections of patches
