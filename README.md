This script calculates second order perturbation energy (E2) using Fock and density matrices in natural atomic orbital (NAO) basis, chemist’s localized
property-optimized orbitals (CLPO)[^1][^2] to Lewis hybrid orbitals (LHO) and LHO to NAO transformation matrices as outputted by JANPA[^3][^4] program. See this paper[^5] by Nikolaienko for the more details on calculation procedure.

> To get needed files the following options should be used when running JANPA:
>
> -doFock -Fock_NAO_File -SDS_NAO_File -CLPO2LHO_File -LHO2NAO_File

The script has been tested with Python 3.8.10 and JANPA 2.02, but should work
with any Python 3.x version

[^1]: Nikolaienko T. Y., Bulavin L. A. Int. J. Quantum Chem., 2019, 119, e25798.
DOI: [10.1002/qua.25798](https://dx.doi.org/10.1002/qua.25798)
[^2]: Nikolaienko T. Y. Phys. Chem. Chem. Phys., 2019, 21, 5285–5294.
DOI: [10.1039/C8CP07276K](https://dx.doi.org/10.1039/C8CP07276K)
[^3]: Nikolaienko T. Y., Bulavin L. A., Hovorun D. M. Comput. Theor. Chem.,
2014, 1050, 15–22. DOI: [10.1016/j.comptc.2014.10.002](https://dx.doi.org/10.1016/j.comptc.2014.10.002)
[^4]: https://janpa.sourceforge.net
[^5]: Yu. Nikolaienko T., Kryachko E. S., Dolgonos G. A. J. Comput. Chem.,
2018, 39, 1090–1102. DOI: [10.1002/jcc.25061](https://dx.doi.org/10.1002/jcc.25061)
***
Please cite this paper if you use this script in your work:

A.S. Gazizov, E.A. Kuznetsova, A.Z. Kamaletdinov, A.V. Smolobochkin, O.A. Lodochnikova, D.P. Gerasimova, A.R. Burilov, M.A. Pudovik, *Org. Chem. Front.*, **2023**

***

