"""WCS-L: a well-cost model for EGS campaigns with learning, a spread rate and random trouble time.

Built on a validated re-implementation of the Sandia Well Cost Simplified model (``wcsl.wcs``).
"""
from . import wcs
from .design import DesignOption, evaluate as evaluate_designs
from .calibration import Calibration, calibrate, calibrate_cape17, load_fervo
from .experience import experience, steady_state_ceiling
from .learning import (bm_factor, brett_millheim, fit_brett_millheim, fit_quadratic, fit_wright,
                       learning_rate, wright)
from .model import ACTIVITIES, WCSL, Uncertainty, WCSLParams, WellPlan, activity_hours, percentiles, uniform_campaign
from .npt import NPTModel
from .spread import SpreadRate

__version__ = "0.1.0"
