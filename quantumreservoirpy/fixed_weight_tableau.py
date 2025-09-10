import numpy as np
import stim
from bitarray import bitarray
from qiskit.quantum_info import Pauli, PauliList


def fixed_weight_tableau(n_qubits, n_meas, weight, XYZ = False):

    """Generates a tableau in dict form, with n_meas pure Z stabilizers of given weight
    and pure X destabilizes of arbitrary weight
    This can be given as the argument "tableau" to the Stabilizer constructor
    
    If XYZ is False, will yield random Z stabilizers only"""

    if n_meas >= n_qubits:
        raise Exception(
            "n_meas should be strictly less than n_qubits"
        )
    if weight >= n_qubits:
        raise Exception(
            "Your Paulis are overweight :( weight should be strictly less than n_qubits"
        )
    
       
    # bistring generation in lexicographic order: https://stackoverflow.com/a/58072652
    def kbits(n, k):
        limit=1<<n
        val=(1<<k)-1
        while val<limit:
            yield "{0:0{1}b}".format(val,n)
            minbit=val&-val #rightmost 1 bit
            fillbit = (val+minbit)&~val  #rightmost 0 to the left of that bit
            val = val+minbit | (fillbit//(minbit<<1))-1

    Stabz = np.zeros(n_qubits)

    for val in kbits(n_qubits, weight):
        valarray = bitarray(val).tolist()
        Stabz = np.vstack((Stabz, valarray))

    Stabz = np.delete(Stabz, 0, 0)
    Stabz = Stabz.astype(int)
    
    rng = np.random.default_rng()
    rng.shuffle(Stabz)
          
    if XYZ:

        # Currently I'm letting Stim take care of choosing a set of linearly independent stabilizers
        # but I'm leaving this rank function here in case we need it later.

        # def gf2_rank(rows):
        #     """
        #     Find rank of a matrix over GF2.

        #     The rows of the matrix are given as nonnegative integers, thought
        #     of as bit-strings.

        #     This function modifies the input list. Use gf2_rank(rows.copy())
        #     instead of gf2_rank(rows) to avoid modifying rows.
        #     """
        #     rank = 0
        #     while rows:
        #         pivot_row = rows.pop()
        #         if pivot_row:
        #             rank += 1
        #             lsb = pivot_row & -pivot_row
        #             for index, row in enumerate(rows):
        #                 if row & lsb:
        #                     rows[index] = row ^ pivot_row
        #     return rank

        Paulis = ['X', 'Y', 'Z']
        signs = ['+', '-']

        def generate_random_Paulistring():
            PauliArray = rng.choice(Paulis, n_qubits)
            mask = rng.choice(Stabz).astype(bool)
            newPaulistring = rng.choice(signs) + ''.join(np.where(mask, PauliArray, 'I'))
            return newPaulistring

        
        Stabxyz = PauliList(generate_random_Paulistring())
        print(Stabxyz)

        while (Stabxyz.size < 3*n_qubits) & (Stabxyz.size < 2**(weight-1)): # generate enough strings that at least n_qubit of them are linearly independent
            
            commute = False
            already_here = False
            while not commute or already_here:
                newPaulistring = Pauli(generate_random_Paulistring())
                commute = newPaulistring.commutes(Stabxyz).all()
                already_here = Stabxyz.equiv(newPaulistring).any()
        
            Stabxyz = Stabxyz.insert(0, newPaulistring)
            print(Stabxyz)

        Stab = Stabxyz.to_labels()

    else:
        Stab = (3*Stabz).tolist() # In Stim, 0=I, 1=X, 2=Y, 3=Z

    stimStabz = [stim.PauliString(stab) for stab in Stab]

    tableau = stim.Tableau.from_stabilizers(stimStabz, allow_redundant=True, allow_underconstrained=not XYZ)

    # Translate to Qiskit and choose n_meas stabilizers
    stabilizer = [str(tableau.z_output(i)).replace('_','I') for i in range(n_meas)]
    destabilizer = [str(tableau.x_output(i)).replace('_','I') for i in range(n_meas)]
    tableau_dict = {"stabilizer" : stabilizer, "destabilizer" : destabilizer}

    return(tableau_dict)

# print(fixed_weight_tableau(10,9,9,XYZ=True))
