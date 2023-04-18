#!/usr/bin/env python3
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# This script calculates second order perturbation energy (E2) using Fock and
# density matrices in natural atomic orbital (NAO) basis, chemist’s localized
# property-optimized orbitals (CLPO) [1,2] to Lewis hybrid orbitals (LHO) and LHO
# to NAO transformation matrices as outputted by JANPA [3,4] program. See this
# paper [5] by Nikolaienko for the more details on calculation procedure.
#
# To get needed files the following options should be used when running JANPA:
#
# -doFock -Fock_NAO_File -SDS_NAO_File -CLPO2LHO_File -LHO2NAO_File
# 
# The script has been tested with Python 3.8.10 and JANPA 2.02, but should work
# with any Python 3.x version
# ----------------------------------------------------
# [1] Nikolaienko T. Y., Bulavin L. A. Int. J. Quantum Chem., 2019, 119, e25798.
# DOI: 10.1002/qua.25798.
# [2] Nikolaienko T. Y. Phys. Chem. Chem. Phys., 2019, 21, 5285–5294.
# DOI: 10.1039/C8CP07276K.
# [3] Nikolaienko T. Y., Bulavin L. A., Hovorun D. M. Comput. Theor. Chem.,
# 2014, 1050, 15–22. DOI: 10.1016/j.comptc.2014.10.002.
# [4] https://janpa.sourceforge.net
# [5] Yu. Nikolaienko T., Kryachko E. S., Dolgonos G. A. J. Comput. Chem.,
# 2018, 39, 1090–1102. DOI: 10.1002/jcc.25061.
# ----------------------------------------------------
# Please cite this paper if you use this script in your work:
# 
# A. S. Gazizov, A. V. Smolobochkin, T. S. Rizbayeva, S. Z. Vatsadze, A. R. Burilov, O. G. Sinyashin, I. V. Alabugin, J. Org. Chem., 2023
# 
#-----------------------------------------------------------------
#-----------------------------------------------------------------
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
#-----------------------------------------------------------------
#-----------------------------------------------------------------
__version__ = '0.1-alpha'
#-----------------------------------------------------------------
import numpy as np
from tabulate import tabulate
import configargparse
#-----------------------------------------------------------------
def loadMatrices(Fock_nao_f, sds_nao_f, cplo_2_lho_f, lho_2_nao_f):
    # We first load whole file as a matrix
    tmp_M = np.loadtxt(Fock_nao_f, skiprows=3, dtype = 'str')
    # Here is total number of orbitals
    nNOs = len(tmp_M) #260

    # Here we load Fock matrix in NAO basis
    Fock_nao = np.loadtxt(Fock_nao_f, skiprows=3, usecols=range(nNOs))
    # Here we load CLPO to LHO transformation matrix
    cplo_2_lho = np.loadtxt(cplo_2_lho_f, skiprows=3, usecols=range(nNOs))
    # Here we load LHO to NAO transformation matrix
    lho_2_nao = np.loadtxt(lho_2_nao_f, skiprows=3, usecols=range(nNOs))
    # Desity matrix in NAO basis (eigenvalues are orbital occupancies)
    sds_nao = np.loadtxt(sds_nao_f, skiprows=3, usecols=range(nNOs))
    # Here we load orbitals' names
    orb_names = np.loadtxt(cplo_2_lho_f, dtype = 'str', skiprows=3, usecols = nNOs)

    # Here we calculate the NAO -> CPLO transformation matrix...
    nao_2_cplo = np.linalg.inv(lho_2_nao) @ np.linalg.inv(cplo_2_lho) 
    # And apply it to Fock and density matrices
    sds_cplo = nao_2_cplo.T @ sds_nao @ nao_2_cplo
    Fock_cplo = nao_2_cplo.T @ Fock_nao @ nao_2_cplo

    return Fock_cplo, sds_cplo, orb_names
#-----------------------------------------------------------------
def doE2Analysis(Fock_matrix, D_matrix, orbital_names = None, qCT_threshold = 0.01, E_threshold = 0.1):
    lowest_donor_occ = 1.0
    highest_accpt_occ = 1.0

    ii_Name = 'unknown'
    jj_Name = 'unknown'
    qCT = -1
    E2 = -1

    titles = ['Donor orbital','Acceptor orbital', 'Donor occupancy', 'Acceptor occupancy', 'Charge transfer (e)', 'E2 energy (kcal/mol)']
    interaction = []
    result = []
    result.append(titles)

    for i in range(len(Fock_matrix)):
        for j in range(len(Fock_matrix)):
            if not i == j:
                # Orbital energies
                Eii = Fock_matrix[i,i]
                Ejj = Fock_matrix[j,j]
                # Orbital occupancies
                Dii = D_matrix[i,i]
                Djj = D_matrix[j,j]
                # Orbital names, if any
                if not orbital_names is None:
                    ii_Name = orbital_names[i]
                    jj_Name = orbital_names[j]

                interaction.append(ii_Name + ' (' + str(i + 1) + ')')
                interaction.append(jj_Name + ' (' + str(j + 1) + ')')
                interaction.append(round(Dii, 4))
                interaction.append(round(Djj, 4))

                if (Dii > lowest_donor_occ) and (Djj < highest_accpt_occ) and (Djj < Dii):
                    #if (Ejj > Eii):
                    Eij = Fock_matrix[i,j]
                    qCT1 = 2*(Eij/(Eii - Ejj))**2
                    E2 = qCT1*(Ejj - Eii)*627.509
                    Dij = D_matrix[i,j]
                    qCT = Dij*Dij/Dii
                    if qCT > qCT_threshold:
                        interaction.append(round(qCT, 4))
                        interaction.append(round(E2, 2))
                        result.append(interaction)
                interaction = []
    return result
#-----------------------------------------------------------------
#-----------------------------------------------------------------
if __name__ == "__main__":
    parser = configargparse.ArgParser(description = '-----------------------------------------------------------------\nThis script calculates second order perturbation energy (E2) using Fock and\ndensity matrices in natural atomic orbital (NAO) basis, chemist’s localized\nproperty-optimized orbitals (CLPO) [1,2] to Lewis hybrid orbitals (LHO) and LHO\nto NAO transformation matrices as outputted by JANPA [3] program. See this [4]\npaper by Nikolaienko for the more details on calculation procedure.\n\nTo get needed files the following options should be used when running JANPA:\n\n-doFock -Fock_NAO_File -SDS_NAO_File -CLPO2LHO_File -LHO2NAO_File\n\n----------------------------------------------------\n[1] Nikolaienko T. Y., Bulavin L. A. Int. J. Quantum Chem., 2019, 119, e25798.\nDOI: 10.1002/qua.25798.\n[2] Nikolaienko T. Y. Phys. Chem. Chem. Phys., 2019, 21, 5285–5294.\nDOI: 10.1039/C8CP07276K.\n[3] Nikolaienko T. Y., Bulavin L. A., # Hovorun D. M. Comput. Theor. Chem.,\n2014, 1050, 15–22. DOI: 10.1016/j.comptc.2014.10.002.\n[4] Yu. Nikolaienko T., Kryachko E. S., Dolgonos G. A. J. Comput. Chem.,\n2018, 39, 1090–1102. DOI: 10.1002/jcc.25061.\n----------------------------------------------------\nPlease cite this paper if you use this script in your work:\n\nlfahgakfljghsdlkfgh\n\n-----------------------------------------------------------------', formatter_class = configargparse.RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--version', action='version', version = "%(prog)s ("+__version__+")", help = 'print script version and exit')
    parser.add_argument('-c', '--config-file', is_config_file = True, help = 'optional config file')
    parser.add_argument('-F', '--Fock-matrix', required = True, metavar = 'FILE', default = 'Fock_NAO.txt', help = 'path to file with Fock matrix in NAO basis. Default: %(default)s')
    parser.add_argument('-SDS', '--SDS-matrix', required = True, metavar = 'FILE', default = 'SDS_NAO.txt', help = 'path to file with density matrix in NAO basis. Default: %(default)s')
    parser.add_argument('-C2L', '--CLPO2LHO-matrix', required = True, metavar = 'FILE', default = 'CLPO2LHO.txt', help = 'path to file with chemist’s localized property-optimized orbitals (CLPO) to Lewis hybrid orbitals (LHO) transformation matrix. Default: %(default)s')
    parser.add_argument('-L2N', '--LHO2NAO-matrix', required = True, metavar = 'FILE', default = 'LHO2NAO.txt', help = 'path to file with Lewis hybrid orbitals (LHO) to natural atomic orbital (NAO) transformation matrix. Default: %(default)s')
    parser.add_argument('-O', '--output-file', required = True, metavar = 'FILE', help = 'path to output file', default = 'E2_output.txt')
#-----------------------------------------------------------------
    args = parser.parse_args()

    try:
        Fock_CPLO, SDS_CPLO, ORB_NAMES = loadMatrices(args.Fock_matrix, args.SDS_matrix, args.CLPO2LHO_matrix, args.LHO2NAO_matrix)
    except FileNotFoundError as err:
        print('ERROR!!!')
        print(str(err))
        exit(-1)

    print('STARTED calculation\n')

    result = doE2Analysis(Fock_CPLO, SDS_CPLO, orbital_names = ORB_NAMES, qCT_threshold = 0.01, E_threshold = 0.1)
    result_table = tabulate(result, headers = 'firstrow', stralign = 'left', numalign = 'center', floatfmt = ('','',".4f",".4f",".4f",".2f"))
    print(result)

    with open(args.output_file, 'w') as f:
            f.write(result_table)
    print('\nFINISHED')
#-----------------------------------------------------------------