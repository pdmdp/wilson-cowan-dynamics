"""Validate the model and reproduce the reviewed memory/phase-space figures."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reviewed_model import DEFAULT, drift, jacobian, fixed_points, simulate


def main():
    p=DEFAULT.copy()
    roots=fixed_points(p)
    assert len(roots)==3 and sum(r['stable'] for r in roots)==2
    point=np.array([.3,.2]);eps=1e-6
    numerical=np.column_stack([(drift(point+eps*np.eye(2)[j],p)-drift(point-eps*np.eye(2)[j],p))/(2*eps) for j in range(2)])
    assert np.allclose(jacobian(point,p),numerical,atol=1e-7)
    _,det=simulate(p)
    for mode in ['additive','ou']:
        _,zero=simulate(p,noise=mode,sigma=0)
        assert np.array_equal(det,zero),mode
    _,noise_a=simulate(p,noise='additive',sigma=.01,seed=8)
    _,noise_b=simulate(p,noise='additive',sigma=.01,seed=8)
    assert np.array_equal(noise_a,noise_b)
    print('Root residuals, Jacobian, zero-noise limits and seed checks passed.',flush=True)

    amplitudes=np.linspace(0,.50,51)
    time,rates=simulate(p,duration=5500,initial=np.full((51,2),.05),pulse_amplitude=amplitudes,clip=True)
    after=rates[time>=4700,:,0].mean(axis=0)
    during=rates[(time>=1100)&(time<3000),:,0].mean(axis=0)
    flags=after>.5
    threshold=float(amplitudes[np.flatnonzero(flags)[0]]) if flags.any() else None
    assert threshold is not None and 0<threshold<.5
    index=int(np.flatnonzero(flags)[0])
    showcase=[index-1,index,index+1]
    dt_checks=[]
    for j in showcase:
        t_half,half=simulate(p,duration=5500,dt=.05,pulse_amplitude=float(amplitudes[j]),clip=True)
        end=float(half[t_half>=4700,0].mean())
        assert bool(end>.5)==bool(flags[j])
        dt_checks.append({'amplitude':float(amplitudes[j]),'post_rate_dt_0_05':end})
    print('Memory sweep threshold:',threshold,flush=True)
    plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans'})
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    axes[0].plot(amplitudes,during,'o-',ms=3,label='During pulse',color='#317c9c')
    axes[0].plot(amplitudes,after,'s-',ms=3,label='After pulse',color='#156b58')
    axes[0].axvline(threshold,color='#8b9296',linestyle=':',label=f'First persistent response: {threshold:.2f}')
    axes[0].set(xlabel='Excitatory input pulse amplitude',ylabel='Mean excitatory rate',title='Input sweep')
    axes[0].legend(fontsize=8)
    for j in showcase:axes[1].plot(time,rates[:,j,0],label=f'Input {amplitudes[j]:.2f}')
    axes[1].axvspan(1000,3000,color='#d6bd70',alpha=.18)
    axes[1].set(xlabel='Time (ms)',ylabel='Excitatory rate',title='Persistent response after stimulus removal')
    axes[1].legend(fontsize=9)
    for ax in axes:ax.spines[['top','right']].set_visible(False)
    out=Path(__file__).resolve().parent
    (out/'assets').mkdir(exist_ok=True);(out/'results').mkdir(exist_ok=True)
    fig.savefig(out/'assets/input-gated-memory.png',dpi=180,facecolor='white')
    plt.close(fig)
    _,trajectories=simulate(p,duration=60,initial=np.array([[.05,.05],[.4,.1],[.7,.5]]))
    fig,ax=plt.subplots(figsize=(6.2,4.5),layout='constrained')
    e,i=np.meshgrid(np.linspace(0,.97,110),np.linspace(0,.8,110))
    field=drift(np.stack([e,i],axis=-1),p)
    ax.contour(e,i,field[...,0],levels=[0],colors=['#317c9c'])
    ax.contour(e,i,field[...,1],levels=[0],colors=['#c46d43'],linestyles='dashed')
    for j in range(3):ax.plot(trajectories[:,j,0],trajectories[:,j,1],lw=1.2,alpha=.75)
    for r in roots:
        ax.scatter(*r['position'],s=65,c='#156b58' if r['stable'] else 'white',edgecolors='#182a34',zorder=5)
    ax.set(xlabel='Excitatory rate',ylabel='Inhibitory rate',title='Validated fixed points and trajectories',xlim=(-.02,.98),ylim=(-.02,.8))
    ax.spines[['top','right']].set_visible(False)
    fig.savefig(out/'assets/phase-plane.png',dpi=180,facecolor='white')
    result={'provenance':'November 2026 numerical review; separate from the coursework outputs.',
      'parameters':p,'checks':['root convergence and residuals','analytic vs finite-difference Jacobian','zero-noise equivalence','seed reproducibility','dt=0.05 memory classification near threshold'],
      'roots':[{'position':r['position'].tolist(),'residual':r['residual'],'stable':r['stable'],'eigenvalues_real':r['eigenvalues'].real.tolist(),'eigenvalues_imag':r['eigenvalues'].imag.tolist()} for r in roots],
      'memory':{'duration_ms':5500,'dt_ms':.1,'pulse_window_ms':[1000,3000],'post_window_ms':[4700,5500],'persistence_threshold':.5,'first_persistent_amplitude':threshold,'amplitudes':amplitudes.tolist(),'during_pulse':during.tolist(),'after_pulse':after.tolist(),'step_size_check':dt_checks}}
    (out/'results/review.json').write_text(json.dumps(result,indent=2)+'\n')
    print('Reviewed figures and results written.',flush=True)


if __name__=='__main__':main()
