
'''Program to calculate the secular and diffusion timescales of a test
particle, due to perturbations from a single planet. The equations and
explanations are provided in Sec. 5.4 of Pearce & Wyatt (2014), and are
based on e.g. Murray & Dermott (1999) and Termaine (1993).

To use the program, simply change the values in the 'User Inputs' section
just below. You should not have to change anything outside of that 
section.

Note: the secular timescale breaks down if the planet and test particle
have the exact same location!'''

############################### Libraries ###############################
import numpy as np
from scipy import integrate

############################## User Inputs ##############################
# Star mass (Solar masses)
mStar_mSun = 1.

# Planet mass (Jupiter masses)
mPlt_mJup = 1.

# Planet semimajor axis (au)
aPlt_au = 5.2

# Test particle semimajor axis (au)
aTest_au = 6.

#################### Constants (don't change these!) ####################
# Jupiter mass in Solar masses
mJup_mSun = 9.55e-4

############################### Functions ###############################
def GetLaplaceCoefficient(j, s, alpha):
	'''Get the Laplace coefficient b_s^j(alpha), given by Equation 7 in
	Pearce et al. (2021)'''
	
	laplaceCoefficientIntegrand = lambda psi: np.cos(j*psi) / (1.0 - 2*alpha*np.cos(psi) + alpha**2)**s
	
	laplaceCoefficient = 1.0 / np.pi * integrate.quad(laplaceCoefficientIntegrand, 0, 2*np.pi)[0]
	
	return laplaceCoefficient

#------------------------------------------------------------------------
def GetSecularTimescale_yr(mStar_mSun, mPlt_mJup, aPlt_au, aTest_au):
	'''Get the secular timescale for a test particle at semimajor axis 
	aTest, that is being perturbed by a planet of mass mPlt and 
	semimajor axis aPlt. Uses Equation 17 in Pearce & Wyatt (2014) if 
	aPlt <= aTest (i.e. internal perturber), or Equation 6 in 
	Pearce et al. (2021) if aPlt > aTest (external perturber).'''

	# Catch if planet has zero mass
	if mPlt_mJup == 0: return np.inf

	# Get planet mass over star mass
	mPlt_mSun = mPlt_mJup * mJup_mSun
	mPlt_mStar = mPlt_mSun / mStar_mSun
		
	# Alpha is the ratio of aPlt to aTest, defined to be 0 < alpha < 1
	alpha = min(aPlt_au/aTest_au, aTest_au/aPlt_au)
		
	# Planet orbital period
	tPlt_yr = GetOrbitalPeriod_yr(aPlt_au, mStar_mSun)
		
	# Laplace coefficient (integral evaluted numerically)
	b1_32 = GetLaplaceCoefficient(1.0, 1.5, alpha)
		
	# For internal perturber
	if aPlt_au <= aTest_au:
		secularTime_yr = 4.*tPlt_yr / mPlt_mStar * alpha**-2.5 / b1_32
		
	# For external perturber
	else:
		secularTime_yr = 4.*tPlt_yr / mPlt_mStar * alpha**-0.5 / b1_32

	return secularTime_yr

#------------------------------------------------------------------------
def GetDiffusionTimescale_yr(mStar_mSun, mPlt_mJup, aPlt_au, aTest_au):
	'''Get the diffusion timescale for a test particle at semimajor axis 
	aTest, that is being perturbed by a planet of mass mPlt and 
	semimajor axis aPlt. Equation 18 in Pearce & Wyatt (2014).'''

	# Catch if planet has zero mass
	if mPlt_mJup == 0: return np.inf
	
	# Get planet mass over star mass
	mPlt_mSun = mPlt_mJup * mJup_mSun

	# Alpha is the ratio of aPlt to aTest, defined to be 0 < alpha < 1
	alpha = min(float(aPlt_au)/aTest_au, float(aTest_au)/aPlt_au)
	
	# Planet orbital period
	tPlt_yr = GetOrbitalPeriod_yr(aPlt_au, mStar_mSun)
	
	# Get the diffusion timescale
	diffusionTime_yr = 0.01 * tPlt_yr * alpha**0.5 * (mPlt_mSun / mStar_mSun)**-2
	
	return diffusionTime_yr

#------------------------------------------------------------------------
def GetOrbitalPeriod_yr(a_au, mStar_mSun):
	'''Get the orbital period of a body orbiting a star'''

	# Mean motion
	n_radPerYr = abs(4*np.pi**2 * mStar_mSun / a_au**3)**.5

	# Orbital period
	period_yr = 2*np.pi / n_radPerYr

	return period_yr

################################ Program ################################

secularTime_yr = GetSecularTimescale_yr(mStar_mSun, mPlt_mJup, aPlt_au, aTest_au)

print('Secular time: %.2e yr' % secularTime_yr)

diffusionTime_yr = GetDiffusionTimescale_yr(mStar_mSun, mPlt_mJup, aPlt_au, aTest_au)

print('Diffusion time: %.2e yr' % diffusionTime_yr)


