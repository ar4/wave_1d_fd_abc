"""Test the propagators."""
import pytest
import numpy as np
from wave_1d_fd_abc.propagators import (Pml)

def ricker(freq, length, dt, peak_time):
    """Return a Ricker wavelet with the specified central frequency."""
    t = np.arange(-peak_time, (length)*dt - peak_time, dt, dtype=np.float32)
    y = ((1.0 - 2.0*(np.pi**2)*(freq**2)*(t**2))
         * np.exp(-(np.pi**2)*(freq**2)*(t**2)))
    return y

def green(x0, x1, dx, dt, t, v, v0, f):
    """Use the 1D Green's function to determine the wavefield at a given
    location and time due to the given source.
    """
    y = np.sum(f[:np.maximum(0, int((t - np.abs(x1-x0)/v)/dt))])*dt*dx*v0/2
    return y

@pytest.fixture
def model_one():
    """Create a model with one reflector, and the expected wavefield."""
    N = 100
    rx = int(N/2)
    model = np.ones(N, dtype=np.float32) * 1500
    model[rx:] = 2500
    max_vel = 2500
    dx = 5
    dt = 0.001
    nsteps = np.ceil(0.27/dt).astype(np.int)
    source = ricker(25, nsteps, dt, 0.05)
    sx = 35
    expected = np.zeros(N)
    # create a new source shifted by the time to the reflector
    time_shift = np.round((rx-sx)*dx / 1500 / dt).astype(np.int)
    shifted_source = np.pad(source, (time_shift, 0), 'constant')
    # reflection and transmission coefficients
    r = (2500 - 1500) / (2500 + 1500)
    t = 1 + r

    # direct wave
    expected[:rx] = np.array([green(x*dx, sx*dx, dx, dt,
                                    (nsteps+1)*dt, 1500, 1500,
                                    source) for x in range(rx)])
    # reflected wave
    expected[:rx] += r*np.array([green(x*dx, (rx-1)*dx, dx, dt,
                                       (nsteps+1)*dt, 1500, 1500,
                                       shifted_source) for x in range(rx)])
    # transmitted wave
    expected[rx:] = t*np.array([green(x*dx, rx*dx, dx, dt,
                                      (nsteps+1)*dt, 2500, 1500,
                                      shifted_source) for x in range(rx, N)])
    return {'model': model, 'dx': dx, 'dt': dt, 'nsteps': nsteps,
            'sources': np.array([source]), 'sx': np.array([sx]),
            'expected': expected}

@pytest.fixture
def model_two():
    """Create a random model and verify small wavefield at large t."""
    N = 100
    np.random.seed(0)
    model = np.random.random(N).astype(np.float32) * 3000 + 1500
    max_vel = 4500
    dx = 5
    dt = 0.6 * dx / max_vel
    nsteps = np.ceil(1.0/dt).astype(np.int)
    num_sources = 10
    sources_x = np.zeros(num_sources, dtype=np.int)
    sources = np.zeros([num_sources, nsteps], dtype=np.float32)
    for sourceIdx in range(num_sources):
        sources_x[sourceIdx] = np.random.randint(N)
        sources[sourceIdx, :] = ricker(25, nsteps, dt, 0.05)
    expected = np.zeros(N, dtype=np.float32)
    return {'model': model, 'dx': dx, 'dt': dt, 'nsteps': nsteps,
            'sources': sources, 'sx': sources_x, 'expected': expected}

@pytest.fixture
def versions():
    """Return a list of implementations."""
    return [Pml]

def test_one_reflector(model_one, versions):
    """Verify that the numeric and analytic wavefields are similar."""

    for v in versions:
        _test_version(v, model_one, atol=8.0)


def test_allclose(model_two, versions):
    """Verify that all implementations produce similar results."""

    for v in versions:
        _test_version(v, model_two, atol=8.0)


def _test_version(version, model, atol):
    """Run the test for one implementation."""
    v = version(model['model'], model['dx'], model['dt'])
    y = v.step(model['nsteps'], model['sources'], model['sx'])
    assert np.allclose(y, model['expected'], atol=atol)
