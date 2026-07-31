# Standard flavour-density-matrix history

These numerical tables were generated from the public FortEPiaNO 1.3 code
using its `ini/explanatory.ini` Standard-Model configuration:

- three active flavours and zero primordial lepton asymmetry;
- normal mass ordering with the mixing parameters recorded in `ini.log`;
- 500 logarithmic time samples;
- 30 selected generalized Gauss--Laguerre momentum nodes below \(y=20\);
- finite-temperature QED corrections through the configured order;
- diagonal collision integrals and the public code's configured
  off-diagonal damping prescription.

The independent run completed with `Total errors: 0` after FortEPiaNO's
7,225 bundled tests had passed.  Its final serialized value is
\(N_{\rm eff}=3.043523644\).  Cosmic reconstructs the same quantity directly
from `Tr rho` and the recorded order-121 quadrature rule.

The tables are generated numerical reference data, not a copy or translation
of FortEPiaNO source code.  FortEPiaNO is described by Bennett et al.,
JCAP 04 (2021) 073, arXiv:2012.02726.  The QKE equations and full collision
operator are discussed by Froustey, Pitrou & Volpe, JCAP 12 (2020) 015,
arXiv:2008.01074.
