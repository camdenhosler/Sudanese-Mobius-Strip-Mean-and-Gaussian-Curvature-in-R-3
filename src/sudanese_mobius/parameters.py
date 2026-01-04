#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class GridParameters:
    resolution: int = 201
    u_min: float = -0.5 * np.pi
    u_max: float =  0.5 * np.pi
    v_min: float = 0.0
    v_max: float = 2.0 * np.pi

    @property
    def u(self):
        return np.linspace(self.u_min, self.u_max, self.resolution)

    @property
    def v(self):
        return np.linspace(self.v_min, self.v_max, self.resolution)