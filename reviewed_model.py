"""Wilson–Cowan numerical checks added in November 2026.

Equations and default parameters follow the course notebook and the Neuromatch
Wilson–Cowan tutorials. See THIRD_PARTY_NOTICES.md for attribution and licensing.
All reviewed simulations use the same baseline-subtracted sigmoid.
"""
from itertools import product
import numpy as np
from scipy.optimize import root
from scipy.special import expit


DEFAULT = dict(tau_E=1., tau_I=2., a_E=1.2, a_I=1., theta_E=2.8,
               theta_I=4., wEE=9., wEI=4., wIE=13., wII=11.,
               I_ext_E=0., I_ext_I=0.)


def response(x, gain, threshold):
    return expit(gain * (x - threshold)) - expit(-gain * threshold)


def derivative(x, gain, threshold):
    value = expit(gain * (x - threshold))
    return gain * value * (1 - value)


def drift(state, p, extra_e=0., synaptic_noise=0.):
    state = np.asarray(state)
    noise = np.broadcast_to(synaptic_noise, state.shape)
    e, i = state[..., 0], state[..., 1]
    e_drive = p['wEE']*e-p['wEI']*i+p['I_ext_E']+extra_e+noise[..., 0]
    i_drive = p['wIE']*e-p['wII']*i+p['I_ext_I']+noise[..., 1]
    return np.stack(((-e+response(e_drive,p['a_E'],p['theta_E']))/p['tau_E'],
                     (-i+response(i_drive,p['a_I'],p['theta_I']))/p['tau_I']),axis=-1)


def jacobian(state, p):
    e, i = state
    de = derivative(p['wEE']*e-p['wEI']*i+p['I_ext_E'],p['a_E'],p['theta_E'])
    di = derivative(p['wIE']*e-p['wII']*i+p['I_ext_I'],p['a_I'],p['theta_I'])
    return np.array([[(-1+p['wEE']*de)/p['tau_E'], -p['wEI']*de/p['tau_E']],
                     [p['wIE']*di/p['tau_I'], (-1-p['wII']*di)/p['tau_I']]])


def fixed_points(p, grid_size=7, residual_tolerance=1e-8):
    """Multistart root search; accept converged, bounded, low-residual roots."""
    lower = -expit(-np.array([p['a_E']*p['theta_E'],p['a_I']*p['theta_I']]))
    upper = 1 + lower
    found = []
    for initial in product(np.linspace(.02,.96,grid_size),repeat=2):
        solution = root(lambda x:drift(x,p),initial,jac=lambda x:jacobian(x,p))
        candidate = solution.x
        residual = np.linalg.norm(drift(candidate,p),ord=np.inf)
        if not solution.success or not np.isfinite(candidate).all() or residual>residual_tolerance:
            continue
        if np.any(candidate<lower-1e-8) or np.any(candidate>upper+1e-8):
            continue
        if any(np.linalg.norm(candidate-r['position'])<1e-5 for r in found):
            continue
        eig = np.linalg.eigvals(jacobian(candidate,p))
        found.append({'position':candidate, 'residual':float(residual),
                      'eigenvalues':eig, 'stable':bool(np.all(eig.real < -1e-8))})
    return sorted(found,key=lambda x:x['position'][0])


def simulate(p, duration=100., dt=.1, initial=(.05,.05), *, pulse_amplitude=0.,
             pulse_window=(1000.,3000.), noise='none', sigma=0., tau_noise=5.,
             seed=0, clip=False):
    """Euler rates with optional additive noise or an exact-discretized OU input.

    initial may have shape (2,) or (n,2). Sigma is a rate diffusion amplitude for
    additive noise, but an input stationary standard deviation for OU noise;
    identical numerical sigma values do not imply matched perturbations.
    """
    if noise not in ('none','additive','ou'):
        raise ValueError('noise must be none, additive, or ou')
    if duration<=0 or dt<=0 or sigma<0 or tau_noise<=0:
        raise ValueError('duration, dt and tau_noise must be positive; sigma nonnegative')
    initial=np.asarray(initial,dtype=float)
    if initial.shape[-1:]!=(2,) or not np.isfinite(initial).all():
        raise ValueError('initial must be finite with final dimension 2')
    time=np.arange(0,duration,dt)
    rates=np.empty((len(time),)+initial.shape);rates[0]=initial
    rng=np.random.default_rng(seed)
    eta=np.zeros_like(initial)
    decay=np.exp(-dt/tau_noise)
    for k in range(len(time)-1):
        pulse=pulse_amplitude if pulse_window[0]<=time[k]<pulse_window[1] else 0.
        if noise=='ou':
            eta=decay*eta+sigma*np.sqrt(1-decay**2)*rng.standard_normal(initial.shape)
        new=rates[k]+dt*drift(rates[k],p,pulse,eta)
        if noise=='additive':new+=sigma*np.sqrt(dt)*rng.standard_normal(initial.shape)
        rates[k+1]=np.clip(new,0,1) if clip else new
    return time,rates
