from pyquil.quil import Program
import pyquil.api as api
from pyquil.gates import *

from grove.pyvqe.vqe import VQE
from scipy.optimize import minimize
import numpy as np

import sys

def ground_state_energy(ansatz):
    print(ansatz)
