import numpy as np
import stim


def fixed_weight_tableau(n_qubits, n_meas, weight):

    """Generates a tableau in dict form, with n_meas pure Z stabilizers of given weight
    and pure X destabilizes of arbitrary weight
    This can be given as the argument "tableau" to the Stabilizer constructor"""

    if n_meas >= n_qubits:
        raise Exception(
            "n_meas should be strictly less than n_qubits"
        )
    if weight >= n_qubits:
        raise Exception(
            "Your Paulis are overweight :( weight should be strictly less than n_qubits"
        )

    Stabz = np.sum(np.array([np.eye(n_qubits, k=-i) for i in range(weight)]), 0)
    Stabz = np.delete(Stabz, slice(1,weight-1), 0) # delete all underweight bitstrings except the first one

    # Initialize bistring generation in lexicographic order: https://stackoverflow.com/a/58072652
    val = (1<<weight)-1
    valarray = np.fromstring("{0:0{1}b}".format(val,n_qubits), dtype='u1') - ord('0')

    # This loop guarantees a full-rank matrix, necessary to get a valid stabilizer tableau:
    # it iterates over all possible bitstrings of the given weight and adds the bitstring to
    # the list of stabilizer only if it is linearly independent from all the others
    while(np.linalg.matrix_rank(Stabz)<n_qubits):

        currentrank = np.linalg.matrix_rank(Stabz)
        Stabz = np.vstack((Stabz, valarray))

        while(np.linalg.matrix_rank(Stabz)==currentrank):
            minbit=val&-val #rightmost 1 bit
            fillbit = (val+minbit)&~val  #rightmost 0 to the left of that bit
            val = val+minbit | (fillbit//(minbit<<1))-1
            valarray = np.fromstring("{0:0{1}b}".format(val,n_qubits), dtype='u1') - ord('0')
            Stabz[-1] = valarray

    stimStabz = (3*Stabz.astype(int)).tolist() # In Stim, 0=I, 1=X, 2=Y, 3=Z
    stimStabz = [stim.PauliString(stimStabz[i]) for i in range(n_qubits)]

    tableau = stim.Tableau.from_stabilizers(stimStabz)

    # Translate to Qiskit
    stabilizer = [str(tableau.z_output(i)).replace('_','I') for i in range(n_qubits-n_meas,n_qubits)]
    destabilizer = [str(tableau.x_output(i)).replace('_','I') for i in range(n_qubits-n_meas, n_qubits)]
    tableau_dict = {"stabilizer" : stabilizer, "destabilizer" : destabilizer}

    return(tableau_dict)