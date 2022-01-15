from ast import increment_lineno
import numpy
import astropy.units as u
import rebound
from galpy.orbit import Orbit
from galpy.util import coords         # instead of galpy.util.bovy_coords
from galpy.potential import KeplerPotential
from galpy.util import bovy_plot                  #bovy_plot
import matplotlib.animation as animation
from IPython.display import HTML
import matplotlib.pyplot as plt
bovy_plot.bovy_print(axes_labelsize=17., text_fontsize=12., xtick_labelsize=15., ytick_labelsize=15.)

from galpy.potential.DissipativeForce import DissipativeForce
from astropy.constants import c

_SAVE_GIFS = False

# Getting the initial conditions using rebound

# Paramètres de l'orbite de l'étoile S2
R0= 8246.7*u.pc
vo= 220.*u.km/u.s
MSgrA= 4.261*10**6*u.Msun
sim= rebound.Simulation()
sim.units= ('AU', 'yr', 'Msun')
# GC
sim.add(m=MSgrA*(u.Msun))
# S2
sim.add(m=0.,
        a=(125.058*u.mas*R0).to_value(u.AU,equivalencies=u.dimensionless_angles()),
        e=0.884649,
        inc=(134.567*u.deg)*(u.rad),
        omega=(66.263*u.deg)*(u.rad),
        Omega=(228.171*u.deg)*(u.rad),
        T=(8.37900*u.yr)*(u.yr)-0.35653101) # time since 2010's apocenter (calculated: 2010+0.35653101)

ogc = Orbit([0.,0.,0.,0.,0.,0.],ro=R0,vo=vo)

pt = sim.particles[1]
o = Orbit([(ogc.ra()*u.deg+(pt.y*u.AU/R0).to(u.deg,equivalencies=u.dimensionless_angles()))*(u.deg),
          (ogc.dec()*u.deg+(pt.x*u.AU/R0).to(u.deg,equivalencies=u.dimensionless_angles()))*(u.deg),
          (ogc.dist()*u.kpc+pt.z*u.AU)*(u.kpc),
          (ogc.pmra()*(u.mas/u.yr)+(pt.vy*u.AU/u.yr/R0).to(u.mas/u.yr,equivalencies=u.dimensionless_angles()))*(u.mas/u.yr),
          (ogc.pmdec()*(u.mas/u.yr)+(pt.vx*u.AU/u.yr/R0).to(u.mas/u.yr,equivalencies=u.dimensionless_angles()))*(u.mas/u.yr),
          (ogc.vlos()*(u.km/u.s)+pt.vz*u.AU/u.yr)*(u.km/u.s)], radec=True,ro=R0,vo=vo)

times= numpy.linspace(0.,2.*16.0455,100001)*u.yr # 2 periods
kp= KeplerPotential(amp=MSgrA,ro=R0)
o.integrate(times,kp)

figsize=(5,10)
o.plot(d1='(ra-{})*{}'.format(ogc.ra(),(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec(),(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       xlabel="Delta RA (arcsec)",
       ylabel="Delta Dec (arcsec)")
plt.plot([0.],[0.],'k+',ms=30.)
plt.xlim(0.055,-0.075)
plt.ylim(-0.03,0.21)

figsize=(8,5)
o.plot(d1='t*1e9',
        d2='(ra-{})*{}'.format(ogc.ra(),(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
      xlabel="t (depuis 2010)",
      ylabel="Delta RA (arcsec)")

figsize=(8,5)
o.plot(d1='t*1e9',
       d2='(dec-{})*{}'.format(ogc.dec(),(u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
      xlabel="t (depuis 2010)",
      ylabel="Delta Dec (arcsec)")

plt.plot(o.time()*u.yr,o.vlos(times))
plt.xlabel("t (depuis 2010)")
plt.ylabel("V_los (km/s-1)")

c= c*(u.km/u.s)
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


sp= SchwarzschildPrecessionForce(amp=MSgrA,ro=R0,fsp=5.)
times= numpy.linspace(0.,4.*16.0455,1001)*u.yr # 4 periods
o.integrate(times,kp)
osp= o()
osp.integrate(times,kp+sp)
figsize=(8,5)
plt.plot(o.time()*u.yr,o.r(times)*u.AU)
plt.plot(osp.time().to(u.yr),osp.r(times)*(u.AU))
plt.xlabel("t (depuis 2010)")
plt.ylabel("r (AU)")
plt.ylim(0.,2000.)

figsize=(5,10)
o.plot(d1='(ra-{})*{}'.format(ogc.ra()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
      xlabel="Delta RA (arcsec)",
      ylabel="Delta Dec (arcsec)")
osp.plot(d1='(ra-{})*{}'.format(ogc.ra()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
        overplot=True)
plt.plot([0.],[0.],'k+',ms=30.)
plt.xlim(0.075,-0.095)
plt.ylim(-0.03,0.21)


sp= SchwarzschildPrecessionForce(amp=MSgrA,ro=R0,fsp=1.)
times= numpy.linspace(0.,4.*16.0455,1001)*u.yr # 4 periods
o.integrate(times,kp)
osp= o()
osp.integrate(times,kp+sp)
figsize=(5,11)
o.plot(d1='(ra-{})*{}'.format(ogc.ra()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       xlabel="Delta RA (arcsec)",
       ylabel="Delta Dec (arcsec)")
osp.plot(d1='(ra-{})*{}'.format(ogc.ra()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
       d2='(dec-{})*{}'.format(ogc.dec()*(u.deg),
                              (u.deg/u.arcsec).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles())),
        overplot=True)
plt.plot([0.],[0.],'k+',ms=30.)
plt.xlim(0.075,-0.095)
plt.ylim(-0.03,0.21)

figsize=(8,5)
plt.plot(o.time().to(u.yr),osp.vlos(times)-o.vlos(times))
plt.xlabel("t (depuis 2010)")
plt.ylabel("delta V_los (km/s-1)")

figsize=(8,5)
plt.plot(o.time().to(u.yr),(osp.ra(times)-o.ra(times))
     *(u.deg/u.mas).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles()))
plt.xlabel("t (depuis 2010)")
plt.ylabel("delta RA (mas)")

figsize=(8,5)
plt.plot(o.time().to(u.yr),(osp.dec(times)-o.dec(times))
     *(u.deg/u.mas).to(u.dimensionless_unscaled,equivalencies=u.dimensionless_angles()))
plt.xlabel("t (depuis 2010)")
plt.ylabel("delta RA (mas)")

