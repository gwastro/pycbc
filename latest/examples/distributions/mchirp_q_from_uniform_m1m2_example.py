import matplotlib.pyplot as plt
from pycbc import distributions
from pycbc import conversions
import numpy as np

# Create chirp mass and mass ratio distribution object that is uniform
# in mass1 and mass2
minmc = 5
maxmc = 60
mc_distribution = distributions.MchirpfromUniformMass1Mass2(mc=(minmc,maxmc))
# generate q in a symmetric range [min, 1/min] to make mass1 and mass2
# symmetric
minq = 1/4
maxq = 1/minq
q_distribution = distributions.QfromUniformMass1Mass2(q=(minq,maxq))

# Take 100000 random variable samples from this chirp mass and mass ratio
# distribution.
n_size = 100000
mc_samples = mc_distribution.rvs(size=n_size)
q_samples = q_distribution.rvs(size=n_size)

# Convert chirp mass and mass ratio to mass1 and mass2
m1 = conversions.mass1_from_mchirp_q(mc_samples['mc'],q_samples['q'])
m2 = conversions.mass2_from_mchirp_q(mc_samples['mc'],q_samples['q'])

# Check the 1D marginalization of mchirp and q is consistent with the 
# expected analytical formula
n_bins = 200
xq = np.linspace(minq,maxq,100)
yq = ((1+xq)/(xq**3))**(2/5)
xmc = np.linspace(minmc,maxmc,100)
ymc = xmc

plt.figure(figsize=(10,10))
# Plot histograms of samples in subplots
plt.subplot(221)
plt.hist2d(mc_samples['mc'], q_samples['q'], bins=n_bins, cmap='Blues')
plt.xlabel('chirp mass')
plt.ylabel('mass ratio')
plt.colorbar(fraction=.05, pad=0.05,label='number of samples')

plt.subplot(222)
plt.hist2d(m1, m2, bins=n_bins, cmap='Blues')
plt.xlabel('mass1')
plt.ylabel('mass2')
plt.colorbar(fraction=.05, pad=0.05,label='number of samples')

plt.subplot(223)
plt.hist(mc_samples['mc'],density=True,bins=100,label='samples')
plt.plot(xmc,ymc*mc_distribution.norm,label='$P(M_c)\propto M_c$')
plt.xlabel('chirp mass')
plt.ylabel('PDF')
plt.legend()

plt.subplot(224)
plt.hist(q_samples['q'],density=True,bins=n_bins,label='samples')
plt.plot(xq,yq*q_distribution.norm,label='$P(q)\propto((1+q)/q^3)^{2/5}$')
plt.xlabel('mass ratio')
plt.ylabel('PDF')
plt.legend()

plt.tight_layout()
plt.show()