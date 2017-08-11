"""Propagate a 1D wavefield using different absorbing boundary
conditions.
"""
import numpy as np
from wave_1d_fd_abc import oneway, oneway_plain, twoway_plain
from wave_1d_fd_pml.propagators import Propagator
from wave_1d_fd_pml.propagators import Pml2

class Oneway(Propagator):
    """One-way absorbing boundary condition."""
    def __init__(self, model, dx, dt=None, abc_width=10, method=3):
        super(Oneway, self).__init__(model, dx, dt=dt, abc_width=abc_width)
        self.method = method

    def step(self, num_steps, sources, sources_x):
        """Propagate wavefield."""

        num_sources = sources.shape[0]
        source_len = sources.shape[1]
        oneway.oneway.step(self.current_wavefield, self.previous_wavefield,
                     self.model_padded, self.dt, self.dx,
                     sources, sources_x, num_steps,
                     self.abc_width, self.pad_width, self.method)

        if num_steps%2 != 0:
            self.current_wavefield, self.previous_wavefield = \
                    self.previous_wavefield, self.current_wavefield

        return self.current_wavefield[self.total_pad: \
                                      self.nx_padded-self.total_pad]


class Oneway_plain(Propagator):

    def step(self, num_steps, sources, sources_x):
        """Propagate wavefield."""

        num_sources = sources.shape[0]
        source_len = sources.shape[1]
        oneway_plain.oneway_plain.step(self.current_wavefield, self.previous_wavefield,
                     self.model_padded, self.dt, self.dx,
                     sources, sources_x, num_steps,
                     self.total_pad)

        if num_steps%2 != 0:
            self.current_wavefield, self.previous_wavefield = \
                    self.previous_wavefield, self.current_wavefield

        return self.current_wavefield[self.total_pad: \
                                      self.nx_padded-self.total_pad]


class Twoway_plain(Propagator):

    def step(self, num_steps, sources, sources_x):
        """Propagate wavefield."""

        num_sources = sources.shape[0]
        source_len = sources.shape[1]
        twoway_plain.twoway_plain.step(self.current_wavefield, self.previous_wavefield,
                     self.model_padded, self.dt, self.dx,
                     sources, sources_x, num_steps,
                     self.total_pad)

        if num_steps%2 != 0:
            self.current_wavefield, self.previous_wavefield = \
                    self.previous_wavefield, self.current_wavefield

        return self.current_wavefield[self.total_pad: \
                                      self.nx_padded-self.total_pad]


