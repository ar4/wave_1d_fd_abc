"""Propagate a 1D wavefield using different absorbing boundary
conditions.
"""
import numpy as np
from wave_1d_fd_abc import pml

class Propagator(object):
    """A finite difference propagator for the 1D wave equation."""
    def __init__(self, model, abc_width, dx, dt=None):
        self.nx = len(model)
        self.dx = np.float32(dx)
        self.abc_width = abc_width
        max_vel = np.max(model)
        if dt:
            self.dt = dt
        else:
            self.dt = 0.6 * self.dx / max_vel
        self.nx_padded = self.nx + 2*self.abc_width
        self.model_padded = np.pad(model,
                                   (self.abc_width, self.abc_width),
                                   'edge')
        self.wavefield = [np.zeros(self.nx_padded, np.float32),
                          np.zeros(self.nx_padded, np.float32)
                         ]
        self.current_wavefield = self.wavefield[0]
        self.previous_wavefield = self.wavefield[1]


class Pml(Propagator):
    """Perfectly Matched Layer."""
    def __init__(self, model, pml_width, dx, dt=None, sigma=None):
        abc_width = pml_width + 4
        super(Pml, self).__init__(model, abc_width, dx, dt)
        self.lpml = [np.zeros(abc_width, np.float32),
                     np.zeros(abc_width, np.float32)
                    ]
        self.current_lpml = self.lpml[0]
        self.previous_lpml = self.lpml[1]

        self.rpml = [np.zeros(abc_width, np.float32),
                     np.zeros(abc_width, np.float32)
                    ]
        self.current_rpml = self.rpml[0]
        self.previous_rpml = self.rpml[1]

        if sigma:
            self.sigma = sigma
        else:
            self.sigma = np.ones(abc_width)

        self.sigma_x = np.gradient(self.sigma)

    def step(self, num_steps, sources=None, sources_x=None):
        """Propagate wavefield."""

        num_sources = sources.shape[0]
        source_len = sources.shape[1]
        pml.pml.step(self.current_wavefield, self.previous_wavefield,
                     self.current_lpml, self.current_rpml,
                     self.previous_lpml, self.previous_rpml,
                     self.sigma, self.sigma_x,
                     self.model_padded, self.dt, self.dx,
                     sources, sources_x, num_steps,
                     self.nx_padded, num_sources, source_len,
                     self.abc_width)

        if num_steps%2 != 0:
            self.current_wavefield, self.previous_wavefield = \
                    self.previous_wavefield, self.current_wavefield
            self.current_lpml, self.previous_lpml = \
                    self.previous_lpml, self.current_lpml
            self.current_rpml, self.previous_rpml = \
                    self.previous_rpml, self.current_rpml

        return self.current_wavefield[self.abc_width: \
                                      self.nx_padded-self.abc_width]
