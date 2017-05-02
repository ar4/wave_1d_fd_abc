"""Propagate a 1D wavefield using different absorbing boundary
conditions.
"""
import numpy as np
from wave_1d_fd_abc import pml
from numba import jit

class Propagator(object):
    """A finite difference propagator for the 1D wave equation."""
    def __init__(self, model, dx, dt=None, abc_width=20, pad_width=8):
        self.nx = len(model)
        self.dx = np.float32(dx)
        self.abc_width = abc_width
        self.pad_width = pad_width
        self.total_pad = self.abc_width + self.pad_width
        max_vel = np.max(model)
        if dt:
            self.dt = dt
        else:
            self.dt = 0.6 * self.dx / max_vel
        self.nx_padded = self.nx + 2*self.total_pad
        self.model_padded = np.pad(model,
                                   (self.total_pad, self.total_pad),
                                   'edge')
        self.wavefield = [np.zeros(self.nx_padded, np.float32),
                          np.zeros(self.nx_padded, np.float32)
                         ]
        self.current_wavefield = self.wavefield[0]
        self.previous_wavefield = self.wavefield[1]

    def steps(self, num_steps, sources, sources_x, interval=1):
        num_intervals = int(np.floor(num_steps/interval))
        saved_steps = np.zeros([num_intervals, self.nx_padded])
        for i in range(num_intervals):
            self.step(interval, sources[:, i*interval:(i+1)*interval],
                      sources_x)
            saved_steps[i, :] = self.current_wavefield[:]
        return saved_steps


class Pml(Propagator):
    """Perfectly Matched Layer."""
    def __init__(self, model, dx, dt=None, pml_width=20, profile=None):
        super(Pml, self).__init__(model, dx, dt, pml_width)
        self.phi = [np.zeros(self.nx_padded, np.float32),
                    np.zeros(self.nx_padded, np.float32)
                   ]
        self.current_phi = self.phi[0]
        self.previous_phi = self.phi[1]

        if profile is None:
            profile = dx*np.arange(-1, pml_width-1, 1,
                                   dtype=np.float32)
            profile[:1] = 0

        self.sigma = np.zeros(self.nx_padded, np.float32)
        self.sigma[self.total_pad-1:self.pad_width-1:-1] = profile
        self.sigma[-self.total_pad:-self.pad_width] = profile
        self.sigma[:self.pad_width] = self.sigma[self.pad_width]
        self.sigma[-self.pad_width:] = self.sigma[-self.pad_width-1]

        self.sigma_x = np.gradient(self.sigma)

    def step(self, num_steps, sources, sources_x):
        """Propagate wavefield."""

        num_sources = sources.shape[0]
        source_len = sources.shape[1]
        pml.pml.step(self.current_wavefield, self.previous_wavefield,
                     self.current_phi, self.previous_phi,
                     self.sigma, self.sigma_x,
                     self.model_padded, self.dt, self.dx,
                     sources, sources_x, num_steps,
                     self.nx_padded, num_sources, source_len,
                     self.abc_width, self.pad_width)

        if num_steps%2 != 0:
            self.current_wavefield, self.previous_wavefield = \
                    self.previous_wavefield, self.current_wavefield
            self.current_pml, self.previous_pml = \
                    self.previous_pml, self.current_pml

        return self.current_wavefield[self.total_pad: \
                                      self.nx_padded-self.total_pad]


    def py_step(self, num_steps, sources, sources_x):
        """Propagate wavefield."""

        @jit(nopython=True)
        def first_x_deriv(f, i, dx):
            return (f[i+1] - f[i-1])/(2*dx)

        @jit(nopython=True)
        def second_x_deriv(f, i, dx):
            return (-735*f[i-8]+15360*f[i-7]
                    -156800*f[i-6]+1053696*f[i-5]
                    -5350800*f[i-4]+22830080*f[i-3]
                    -94174080*f[i-2]+538137600*f[i-1]
                    -924708642*f[i+0]
                    +538137600*f[i+1]-94174080*f[i+2]
                    +22830080*f[i+3]-5350800*f[i+4]
                    +1053696*f[i+5]-156800*f[i+6]
                    +15360*f[i+7]-735*f[i+8])/(302702400*dx**2)

        @jit(nopython=True)
        def fd_interior(f, fp, model_padded, dt, dx, i):
            f_xx = second_x_deriv(f, i, dx)
            fp[i] = model_padded[i]**2 * dt**2 * f_xx + 2 * f[i] - fp[i]

        @jit(nopython=True)
        def fd_pml(f, fp, phi, phip, sigma, sigma_x, model_padded,
                   dt, dx, i):
            f_xx = second_x_deriv(f, i, dx)
            f_x = first_x_deriv(f, i, dx)
            phi_x = first_x_deriv(phi, i, dx)
            sigma_xx = (sigma[i+1] - 2*sigma[i] + sigma[i-1])/dx**2

            phi_xt = f_xx - sigma_xx*phi[i] - sigma_x[i] * phi_x

            fp[i] = ((model_padded[i]**2 * dt**2 * phi_xt + 2 * f[i])
                     / (1 - sigma_x[i] * dt / 2) - fp[i])
            phip[i] = dt * f_x + phi[i] * (1 - dt * sigma_x[i])

        @jit(nopython=True)
        def numba_step(f, fp, phi, phip, sigma, sigma_x, model_padded,
                       dt, dx, pad_width, abc_width, nx, num_steps):

            lpml_start = pad_width
            interior_start = lpml_start + abc_width
            rpml_start = interior_start + nx
            rpad_start = rpml_start + abc_width

            lpml = range(lpml_start, interior_start)
            interior = range(interior_start, rpml_start)
            rpml = range(rpml_start, rpad_start)

            for step in range(num_steps):

                # left PML
                #for i in lpml:
                #    fd_pml(f, fp, phi, phip, sigma, sigma_x,
                #           model_padded, dt, dx, i)

                # interior (no PML)
                for i in interior:
                    fd_interior(f, fp, model_padded, dt, dx, i)

                # right PML
                for i in rpml:
                    fd_pml(f, fp, phi, phip, sigma, sigma_x,
                           model_padded, dt, dx, i)

                # Add sources
                for i in range(sources.shape[0]):
                    sx = sources_x[i] + abc_width
                    source_amp = sources[i, step]
                    fp[sx] += (model_padded[sx]**2 * dt**2 * source_amp)

                tmp = f
                f = fp
                fp = tmp
                tmp = phi
                phi = phip
                phip = tmp

        numba_step(self.current_wavefield, self.previous_wavefield,
                   self.current_phi, self.previous_phi,
                   self.sigma, self.sigma_x,
                   self.model_padded, self.dt, self.dx,
                   self.pad_width, self.abc_width, self.nx,
                   num_steps)

        return self.current_wavefield[self.total_pad: \
                                      self.nx_padded-self.total_pad]
