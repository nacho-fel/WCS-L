"""Non-productive time (trouble): random events whose frequency falls with experience.

    N_n ~ Poisson(lambda_inf + (lambda0 - lambda_inf) exp(-k_lambda (E_n - 1))),   D_j ~ LogNormal(log(median), sigma)

lambda_inf is the long-run trouble rate that experience does not remove (0 in the first calibration; 0.35
per well from the band recalibration: Icelandic wells had trouble in 24 of 77 wells, 31 %).

Replaces the WCS's zero trouble time plus flat contingency with a right-skewed distribution.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class NPTModel:
    lambda0: float = 1.5        # mean events on the first well
    k_lambda: float = 0.3       # decay of the event rate with experience
    median_days: float = 1.5    # median duration of one event
    sigma: float = 0.8          # log-standard deviation of the duration
    lambda_inf: float = 0.0     # long-run event rate per well that never goes away

    def rate(self, E):
        return self.lambda_inf + (self.lambda0 - self.lambda_inf) * np.exp(-self.k_lambda * (np.asarray(E, float) - 1.0))

    def mean_event_days(self) -> float:
        return self.median_days * np.exp(self.sigma ** 2 / 2)

    def expected_days(self, E):
        return self.rate(E) * self.mean_event_days()

    def sample(self, E, rng: np.random.Generator, size: int) -> np.ndarray:
        """Trouble days, shape (size, len(E))."""
        lam = self.rate(E)
        counts = rng.poisson(lam, size=(size, len(lam)))
        total = np.zeros(counts.shape)
        kmax = int(counts.max()) if counts.size else 0
        mu = np.log(self.median_days)
        for j in range(kmax):
            d = rng.lognormal(mu, self.sigma, size=counts.shape)
            total += np.where(counts > j, d, 0.0)
        return total
