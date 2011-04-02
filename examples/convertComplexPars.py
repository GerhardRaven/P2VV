from P2VV import convertComplexPars
from math import sqrt

amplitudes = True

if amplitudes :
  # transversity amplitudes:
  #        Re(0) Im(0)    Re(par)       Im(par)         Re(perp)       Im(perp)       Re(S)          Im(S)
  vals = ((1.,   0.   ), (-6.09615e-01, -1.38858e-01), (-6.35889e-01,  1.46244e-01), (-1.74615e-01,  2.64948e-01))

  covs = ((0.,   0.,       0.,           0.,             0.,           0.,             0.,           0.       ),  # Re(0)
          (0.,   0.,       0.,           0.,             0.,           0.,             0.,           0.       ),  # Im(0)
          (0.,   0.,       5.817e-05,   -3.991e-05,     -1.843e-05,   -2.226e-05,      1.170e-06,   -2.039e-05),  # Re(par)
          (0.,   0.,      -3.991e-05,    2.542e-04,      1.702e-05,    1.227e-04,      6.763e-06,    2.945e-05),  # Im(par)
          (0.,   0.,      -1.843e-05,    1.702e-05,      5.820e-05,    2.611e-05,      1.411e-06,    3.162e-05),  # Re(perp)
          (0.,   0.,      -2.226e-05,    1.227e-04,      2.611e-05,    1.502e-04,      2.743e-06,    4.393e-05),  # Im(perp)
          (0.,   0.,       1.170e-06,    6.763e-06,      1.411e-06,    2.743e-06,      2.061e-05,   -9.360e-06),  # Re(S)
          (0.,   0.,      -2.039e-05,    2.945e-05,      3.162e-05,    4.393e-05,     -9.360e-06,    2.066e-04))  # Im(S)

  # convert amplitudes from cartesian to polar (squared magnitudes)
  # normalize on first three amplitudes
  vals1, covs1 = convertComplexPars(vals, covs, inType='cart', outType='polar',
      normalize = 3, magnSquared = True)

  # convert amplitudes from polar (squared magnitudes) to cartesian
  # normalize on first amplitude
  vals2, covs2 = convertComplexPars(vals1, covs1, inType='polar',
      outType='cart', normalize = 1, magnSquared = True)

else :
  # lambda CPV
  vals = ((1.37452,     3.96158e-01),)
  covs = (( 2.527e-01, -6.181e-01  ),
          (-6.181e-01,  1.824      ))

  # convert lambda from cartesian to polar
  vals1, covs1 = convertComplexPars(vals, covs, inType='cart', outType='polar',
      normalize = 0, magnSquared = False)

  # convert lambda from polar to cartesian
  vals2, covs2 = convertComplexPars(vals1, covs1, inType='polar',
      outType='cart', normalize = 0, magnSquared = False)


# print initial cartesian values and covariance matrix
print 'initial cartesian values and covariance matrix'
for i, val in enumerate(vals) :
  for j, val in enumerate(val) :
    print '{0:9.3g} +- {1:.1g}  '.format(val, sqrt(covs[2*i+j][2*i+j])),
  print
print

for covs_i in covs :
  print '(',
  for cov in covs_i :
    print '{0:9.3g}'.format(cov),
  print ')'
print

# print polar values and covariance matrix
print 'polar values and covariance matrix'
for i, val in enumerate(vals1) :
  for j, val in enumerate(val) :
    print '{0:9.3g} +- {1:.1g}  '.format(val, sqrt(covs1[2*i+j][2*i+j])),
  print
print

for covs_i in covs1 :
  print '(',
  for cov in covs_i :
    print '{0:9.3g}'.format(cov),
  print ')'
print

# print final cartesian values and covariance matrix
print 'final cartesian values and covariance matrix'
for i, val in enumerate(vals2) :
  for j, val in enumerate(val) :
    print '{0:9.3g} +- {1:.1g}  '.format(val, sqrt(covs2[2*i+j][2*i+j])),
  print
print

for covs_i in covs2 :
  print '(',
  for cov in covs_i :
    print '{0:9.3g}'.format(cov),
  print ')'
print

