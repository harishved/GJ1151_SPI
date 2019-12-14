# Sub-Alfvenic flux for dipole field configuration
# All units are cgs unless specified otherwise
#
import matplotlib.pyplot as plt
import numpy as np

alpha = 1		# Interaction strength of obstacle with flow
B0 = 100		# Surface magnetic field in Gauss
Bp0 = 1.0		# Planetary magnetic field
theta_M = np.pi/2 	# Angle between planetary field and stellar field
Omega_star = 2.0*np.pi/(125*24*3600)	# Rotation angular velocity of the star
Rp = 6400e5		# Planetary radius in cm
R_star = 1.3e10		# Radius of the star
M_star_msun = 0.168	# Stellar mass in units of solar mass
nbase = 1e6		# Plasma density at the plabnet (only used to computed the Alfven speed)		# Does not enter the energy / flux determination but for that (so plasma temp does not matter)
Tcor = 1e6		# Base temperature in the corona
vsound = 180e5*(Tcor/2e6)**0.5	# Coronal sound speed
rsonic = 2e11/R_star*M_star_msun*(Tcor/2e6)**-1	# Radius of the sonic point
#
dvec = np.linspace(5,50,1000)*R_star	# Orbital distance of the planet
n0 = nbase*(dvec/R_star)**-2		# Plasma density at the planet
B_r = B0*(dvec/R_star)**-2		# B-field at the planet Eqn B8 Turnpenney 2018
#
tp = np.load("parker_wind_vel_sol_lin.npz")	# Numerically pre computed wind solution 
vsw = vsound*np.interp(dvec/R_star/rsonic,tp['r'],tp['v'])	# Wind velocity
#
plt.plot(dvec/R_star/rsonic,vsw/1e5,'k')
plt.xlabel("Distance / sonic radius"); plt.ylabel("Sound speed [km/s]");
plt.minorticks_on(); 
plt.title("Sonic radius = %.1f Stellar radius"%rsonic)
plt.tight_layout(); plt.savefig("sound_speed.pdf"); plt.close()
#
v_orb = 30e5*(M_star_msun/(dvec/1.5e13))**0.5	# Orbital speed of planet
v_rot = dvec*Omega_star				# Corotation speed
B_phi = B_r * dvec*Omega_star/vsw		# Azimuthal field Eqn B9 Turnpenney 2018
B_tot = (B_r**2+B_phi**2)**0.5			# Total B-field at planet
#
B_ang = np.arctan2(B_phi,B_r)		# Angle the B-field makes with the radial direction
#					# First terms of RHS of eqn B11 Turnpenney 2018
vrel = (v_orb**2+vsw**2)**0.5		# Relative vel;ocioty bewteen wind and obstacle
vrel_angle = np.arctan2(v_orb,vsw)	# Angle bewteen the radial vactor and relative velocity
#					# Secind term of RHS of eqn B11 Turnpenney 2018
theta = np.absolute(B_ang - vrel_angle)	# Eqn B11 of Turnpenney 2018
#
plt.plot(dvec/R_star,B_ang*180/np.pi,'r',label=r"$\theta_B$");
plt.plot(dvec/R_star,vrel_angle*180/np.pi,'b',label=r"$\theta_v$");
plt.plot(dvec/R_star,theta*180/np.pi,'k');
plt.xlabel("Distance / stellar radius"); plt.ylabel("Theta [deg]")
plt.tight_layout(); plt.savefig("theta.pdf"); plt.close()
#
#
ct2 = (np.sin(theta))**2		# Geomatric factor in efficiency 
v_alf = 2.18e11*B_tot*n0**-0.5		# Alfven speed at the planet
M_A = vrel/v_alf			# Alfven mach number
#
Bp = Bp0 * M_star_msun**0.5/(dvec/1.5e13)**1.5 / 365.
Rp_eff = Rp * (Bp/B_tot)**(1./3.) * (3*np.cos(theta_M/2))**0.5	# Effective pkabetary radius (Saur 2013)
Rp_eff[Rp_eff<Rp]=Rp			# Reff cannot be smaller than Rplanet
#
plt.plot(dvec/R_star,Bp)
plt.xlabel("Distance / stellar radius"); 
plt.ylabel("Planetary B-field [Gauss]")
plt.minorticks_on(); 
plt.tight_layout(); 
plt.savefig("Bp.pdf"); 
plt.close()
#
S_mks = 2*np.pi*(Rp_eff/1e2)**2 * (v_alf/1e2) * (alpha*M_A)**2 * (B_tot/1e4)**2/(4*np.pi*1e-7) * ct2
S = S_mks * 1e7	# MKS tO CGS units
#
#
period_marker = np.array([0.5,1,2,3,5,7,10,15,20])
dvec_marker = ((period_marker/365.)*M_star_msun**0.5)**(1./1.5) *1.5e13/R_star
#
#
plt.figure(figsize=(4,5.5))
ax2 = plt.subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
ax1 = plt.subplot2grid((3,1),(0,0),rowspan=1,colspan=1)
#
#
ax2.plot(dvec/R_star,np.log10(S),label="ST-model")
ax2.plot(dvec/R_star,np.log10(S/ct2/M_A*0.25),label="LZ-model")
ax1.plot(dvec/R_star,np.log10(M_A))
#ax2.fill_between(dvec/R_star,np.log10(S),np.log10(S/M_A/ct2*0.25),hatch='x',color="yellow",alpha=0.2)
#
# Compute the beam solid angle of emitter /4pi
tp1 = (np.cos(np.arccos(0.3)-0.3/2) - np.cos(np.arccos(0.3)+0.3/2)) / 2
tp2 = (np.cos(np.arccos(0.7)-0.7/2) - np.cos(np.arccos(0.7)+0.7/2)) / 2
#
#
bsa_min = min(tp1,tp2)
bsa_max = max(tp1,tp2)
#
# Min and max power based on beaming solid angle and 1-10% emission efficieny
pow_min = 2e21/0.01*bsa_min * 2.8e6*B0/(167.-120.)/1e6
pow_max = 2e21/0.1*bsa_max * 2.8e6*B0/(167.-120.)/1e6
#
ax2.fill_between([np.amin(dvec)/R_star,np.amax(dvec)/R_star],[np.log10(pow_min),np.log10(pow_min)],[np.log10(pow_max),np.log10(pow_max)],color="cyan",alpha=0.25)
#
ax11 = ax1.twiny()
#
ax11.set_xticks(dvec_marker)
ax11.set_xticklabels(period_marker)
ax11.set_xlim([np.amin(dvec)/R_star, np.amax(dvec)/R_star])
ax1.set_xlim([np.amin(dvec)/R_star, np.amax(dvec)/R_star])
ax2.set_xlim([np.amin(dvec)/R_star, np.amax(dvec)/R_star])
#
ax11.set_xlabel("Orbital period [days]")
ax1.tick_params(top=False,which='both')
ax11.tick_params(top=False,which='minor')
ax1.set_ylabel(r"${\rm log}_{10}(M_A)$")
#
ax2.set_ylabel(r"${\rm log}_{10}$(Star-ward flux [ergs/s])")
ax2.set_xlabel("Distance / Stellar radius")
ax2.set_ylim([21,25])
ax1.set_ylim([-3,0])
ax1.plot([np.amin(dvec)/R_star,np.amax(dvec)/R_star],[22,22],'k--')
ax1.set_xticklabels([])
#
ax2.legend(loc=3)
ax2.text(x=25,y=24.4,s=r"$B_\ast=100\,{\rm G}$",fontsize=14)
ax2.text(x=25,y=24.0,s=r"$n_{\rm base} = 10^6\,{\rm cm}^{-3}$",fontsize=14)
ax1.text(x=-4,y=0.7,s="Open field",fontsize=16)
plt.tight_layout()
plt.savefig("open.pdf")
plt.close()
