import numpy
import math
import astropy.units as u
import astropy.constants as const
import rebound
from galpy.orbit import Orbit
from galpy.util import bovy_coords
from galpy.potential import KeplerPotential
import matplotlib.animation as animation
from IPython.display import HTML
import matplotlib as mpl
import matplotlib.pyplot as plt
from galpy.util import bovy_plot
bovy_plot.bovy_print(axes_labelsize=17.,text_fontsize=12.,
                     xtick_labelsize=15.,ytick_labelsize=15.)
_SAVE_GIFS= False

R0= 8246.7*u.pc
vo= 220.*u.km/u.s
MSgrA= 4.261*10**6*u.Msun
sim= rebound.Simulation()
sim.units= ('AU', 'yr', 'Msun')
# GC
sim.add(m=MSgrA.to_value(u.Msun))
# S2
sim.add(m=0.,
        a=(125.058*u.mas*R0).to_value(u.AU,equivalencies=u.dimensionless_angles()),
        e=0.884649,
        inc=(134.567*u.deg).to_value(u.rad),
        omega=(66.263*u.deg).to_value(u.rad),
        Omega=(228.171*u.deg).to_value(u.rad),
        T=(8.37900*u.yr).to_value(u.yr)-0.35653101)

ogc= Orbit([0.,0.,0.,0.,0.,0.],ro=R0,vo=vo)

pt= sim.particles[1]
o= Orbit([(ogc.ra()*u.deg+(pt.y*u.AU/R0).to(u.deg,equivalencies=u.dimensionless_angles())).to_value(u.deg),
          (ogc.dec()*u.deg+(pt.x*u.AU/R0).to(u.deg,equivalencies=u.dimensionless_angles())).to_value(u.deg),
          (ogc.dist()*u.kpc+pt.z*u.AU).to_value(u.kpc),
          (ogc.pmra()*(u.mas/u.yr)+(pt.vy*u.AU/u.yr/R0).to(u.mas/u.yr,equivalencies=u.dimensionless_angles())).to_value(u.mas/u.yr),
          (ogc.pmdec()*(u.mas/u.yr)+(pt.vx*u.AU/u.yr/R0).to(u.mas/u.yr,equivalencies=u.dimensionless_angles())).to_value(u.mas/u.yr),
          (ogc.vlos()*(u.km/u.s)+pt.vz*u.AU/u.yr).to_value(u.km/u.s)],radec=True,ro=R0,vo=vo)

periods= 2.
times= numpy.linspace(0.,periods*16.0455,100001)*u.yr # 2 periods
kp= KeplerPotential(amp=MSgrA,ro=R0)
o.integrate(times,kp)

mpl.rcParams['figure.figsize'] = (5., 10.)
o.plot(d1='(ra-{})*{}'.format(ogc.ra(),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec(),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
      xlabel=r'$\Delta \mathrm{RA}\,(\mathrm{arcsec})$',
      ylabel=r'$\Delta \mathrm{Dec}\,(\mathrm{arcsec})$')

plt.plot([0.],[0.],'k+',ms=30.)
plt.xlim(0.055,-0.075)
plt.ylim(-0.03,0.21)

from galpy.potential.DissipativeForce import DissipativeForce
from astropy.constants import c
c = c.to_value(u.km/u.s)
class SchwarzschildPrecessionForce(DissipativeForce):
    def __init__(self,amp=1.,fsp=1.,gamma=1.,beta=1.,ro=None,vo=None):
        DissipativeForce.__init__(self,amp=amp,ro=ro,vo=vo,
                                  amp_units='mass')
        self._fsp= fsp
        self._gamma= gamma
        self._beta= beta
    # Following functions implement the vec{r}/r and \vec{v} parts, respectively
    def _force_firstterm(self,r,v):
        return 1./(c/self._vo)**2./r**3.*(2.*(self._gamma+self._beta)*self._amp/r-self._gamma*v**2.)
    def _force_secondterm(self,r,vr):
        return 2.*(1.+self._gamma)/(c/self._vo)**2./r**2.*vr
    # Now compute the three projections of the forve
    def _Rforce(self,R,z,phi=0.,t=0.,v=None):
        r= numpy.sqrt(R**2.+z**2.)
        vr= R/r*v[0]+z/r*v[2]
        vmag= numpy.sqrt(v[0]**2.+v[1]**2.+v[2]**2.)
        return self._fsp*(self._force_firstterm(r,vmag)*R+self._force_secondterm(r,vr)*v[0])
    def _zforce(self,R,z,phi=0.,t=0.,v=None):
        r= numpy.sqrt(R**2.+z**2.)
        vr= R/r*v[0]+z/r*v[2]
        vmag= numpy.sqrt(v[0]**2.+v[1]**2.+v[2]**2.)
        return self._fsp*(self._force_firstterm(r,vmag)*z+self._force_secondterm(r,vr)*v[2])
    def _phiforce(self,R,z,phi=0.,t=0.,v=None):
        r= numpy.sqrt(R**2.+z**2.)
        vr= R/r*v[0]+z/r*v[2]
        return self._fsp*(self._force_secondterm(r,vr)*v[1]*R)

sp= SchwarzschildPrecessionForce(amp=MSgrA,ro=R0,fsp=1.)

periods = 4.
times= numpy.linspace(0.,periods*16.0455,1001)*u.yr
o.integrate(times,kp)
osp= o()
osp.integrate(times,kp+sp)

mpl.rcParams['figure.figsize'] = (5., 10.)
osp.plot(
        d1='(ra-{})*{}'.format(ogc.ra(),
        (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
        d2='(dec-{})*{}'.format(ogc.dec(),
        (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
        xlabel=r'$\Delta \mathrm{RA}\,(\mathrm{arcsec})$',
        ylabel=r'$\Delta \mathrm{Dec}\,(\mathrm{arcsec})$')

plt.plot([0.],[0.],'k+',ms=30.)
plt.xlim(0.075,-0.095)
plt.ylim(-0.03,0.21)
plt.title("Orbite de S2 autour de SgrA*")

plt.plot(o.time()*u.yr,o.vlos(times))
plt.xlabel(r'$t\ (\mathrm{year})$')
plt.ylabel(r'$V_{\mathrm{los}}\,(\mathrm{km\,s}^{-1})$');
print(min(o.vlos(times)), max(o.vlos(times)))

fe = 1.36*10**14 # Fréquence d'émission de S2
vr = 7700.
v = 4000.
fr = (fe*(1-(vr/c)+(vr**2/c**2)-(v**2/(2*c**2))))
print("Frequency emmited by S2:", "{:.2e}".format(fe))
print("Frequency received on Earth: "+"{:.2e}".format(fr))
print("Frequency shift:", "{:.2e}".format(fe-fr))

nperiods= 120.
times= numpy.linspace(0.,nperiods*16.0455,1000001)*u.yr # n périodes
o.integrate(times[times < 16.0455*u.yr],kp)
osp= o()
osp.integrate(times,kp+sp)

def delta_ra(ra):
    return (ra-ogc.ra())*(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())
def delta_dec(dec):
    return (dec-ogc.dec())*(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())

mpl.rcParams['figure.figsize'] = (6., 7.)

fig, ax= plt.subplots()
line2 = plt.plot([0.],[0.],'k+',ms=10.)

plt.xlabel(r'$\Delta \mathrm{RA}\,(\mathrm{arcsec})$')
plt.ylabel(r'$\Delta \mathrm{Dec}\,(\mathrm{arcsec})$')
plt.xlim(0.075,-0.25)
plt.ylim(-0.19,0.26)

nstride= round(len(times)/nperiods)+1

line,= ax.plot(delta_ra(osp.ra(times[:nstride])),delta_dec(osp.dec(times[:nstride])),'-')

orbitValue= ax.annotate(r'$\mathrm{year}=%.0f$' % (2010.),
                 (0.05,0.05),xycoords='axes fraction',
                 horizontalalignment='left',verticalalignment='bottom',size=18.)
redshiftValue= ax.annotate(r'$\mathrm{Frequency\ shift}=%s$' % (0.),
                 (0.05,0.05),xycoords='axes fraction',
                 horizontalalignment='left',verticalalignment='top',size=18.)

plt.legend(frameon=False,fontsize=18.,loc='upper right')

def update(i):
    minRange= int(round((i+0.5)*nstride))
    maxRange= int(round((i+1.5)*nstride))
    line.set_xdata(delta_ra(osp.ra(times[minRange:maxRange])))
    line.set_ydata(delta_dec(osp.dec(times[minRange:maxRange])))
    orbitValue.set_text(r'$\mathrm{Year}=%.0f$' % (round(2010.35653101+times[i*nstride].to_value(u.yr),-1)))
    v = o.vlos(times[i])
    # Changing v for Vlos
    fr = (fe*(1-(vr/c)+(vr**2/c**2)-(v**2/(2*c**2))))
    shift = "{:.6e}".format(fe-fr)
    redshiftValue.set_text(r'$\mathrm{Frequency\ shift}=%s$' % shift)
    progress = (round(i*100/nperiods))
    if progress%10==0:
        print(progress, "%")
    return (line,)

anim = animation.FuncAnimation(fig,update,frames=len(times)//nstride,interval=1,blit=True,repeat=True)
plt.tight_layout()
anim.save('S2-SgrA.gif',writer='ffmpeg',dpi=200)
plt.close()

def distance(ra, dec):
    return abs(delta_dec(dec)/delta_ra(ra))

def fr(v):
    return (fe*(1-(vr/c)+(vr**2/c**2)-(v**2/(2*c**2))))

fig, ax= plt.subplots()

plt.xlabel(r'$\mathrm{Time}\,(\mathrm{years})$')
plt.ylabel(r'$\mathrm{Frequency Received}\,(\mathrm{Hz})$')

line,= ax.plot(times[:nstride],fr(o.vlos(times[:nstride])),'-')
