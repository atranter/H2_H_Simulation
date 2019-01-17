from openfermion import *
from openfermionpsi4 import run_psi4
import os
import sys

class OpenFermionWrapper():
    ''' 
    This class is designed to act as a wrapper for the OpenFermion package.
    This wrapper's focus is to create an interface between a client program
    and the OpenFermion syntax with an increased ease of use. 
    The intended use case for this wrapper is when a client program would
    like to define a molecule and create the hamiltonian, energy eigenstates, 
    or a quantum circuit written in QASM of the time evolution of the
    hamiltonian through the Trotter-Suzuki decomposition method. 
    Given the increased complexity around these use-cases, this wrapper is 
    created to keep in mind that client programs may want to restructure the
    created data. Hence, in contrast to normal programming design, there will
    be no private members of this class.
    
    Typical Progression:
    Set variables using setter functions
    load_molecule
    create_hamiltonians
    perform_transform
    generate_circuit
    
    Author: William A. Simon
    Date: 12/7/2018
    
    TODO:
    support BKSF, BKTree, Parity mappings
    get N lowest energy states
     '''
    def __init__(self):
        self.name = None
        self.geometry = None
        self.basis = None
        self.multiplicity = None
        self.charge = None
        self.description = ''
        self.filename = None

        self.plugin = "psi4"

        self.molecule = None

        self.active_space_start = None
        self.active_space_stop = None

        self.molecular_hamiltonian = None
        self.fermion_hamiltonian = None
        self.n_qubits = None

        self.ground_state_energy = None

        self.qubit_hamiltonian_bk = None
        self.qubit_hamiltonian_jw = None
        self.qubit_hamiltonian_bksf = None
        self.qubit_hamiltonian_bktree = None

        self.jw_circuit = None
        self.bk_circuit = None

    def set_name(self, name):
        self.name = name
    def set_geometry(self, geometry):
        self.geometry = geometry
    def set_basis(self, basis): 
        self.basis = basis
    def set_multiplicity(self, multiplicity): 
        self.multiplicity = multiplicity
    def set_charge(self, charge):
        self.charge = charge
    def set_description(self, description): 
        self.description = description
    def set_active_space_start(self, orbital): 
        self.active_space_start = orbital
    def set_active_space_stop(self, orbital): 
        self.active_space_stop = orbital
    def set_plugin(self, plugin):
        if((plugin != 'pyscf') and (plugin != "pyscf")):
            sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
                           Plugin not recognized\n''')
        self.plugin = plugin

    def mapping_error(self, mapping):
        ''' 
        ARGS:
            mapping - a string
        RETURNS:
            None
        A simple helper function to print and exit with an error related
        to an incorrect mapping request. 
        '''
        sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
             Mapping: {} was not recognized'''.format(mapping))

    def set_parameters(self, geometry, basis, multiplicity, charge):
        ''' 
        ARGS:
            geometry - Geometry of the molecule in the OpenFermion/Psi4 syntax
            basis - a string defining the chosen basis
            multiplicity - integer
            charge - integer
        RETURNS:
            None

        A helper function typically called by this class to set the necessary
        variables defining the desired molecule if they have not been set yet.
        This function calls the simple setter functions for the input variables.
        '''
        if(self.geometry == None):
            self.set_geometry(geometry)
        if(self.basis == None):
            self.set_basis(basis)
        if(self.multiplicity == None):
            self.set_multiplicity(multiplicity)
        if(self.charge == None):
            self.set_charge(charge)

    def check_parameters(self):
        ''' 
        ARGS:
            None

        RETURNS:
            None

        Helper function to make sure the necessary parameters have been
        set. If they have not been, system exits with an error message.

        ''' 
        if((self.geometry == None)     or (self.basis == None) or
           (self.multiplicity == None) or (self.charge == None)):
            sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
                            Must define all molecular parameters\n''')

    def load_molecule(self, geometry=None, basis=None, multiplicity=None, 
       charge=None, forceCalculation=False):
        ''' 
        ARGS:
            Necessary: None
            Optional: geometry - Geometry of the molecule in the 
                                 OpenFermion/Psi4 syntax
                      basis - a string defining the chosen basis
                      multiplicity - integer
                      charge - integer
                      forceCalculation - a boolean (preset to false) that will
                                         force the run_psi4 method when set
                                         to True instead of loading the saved
                                         integral calculations

        RETURNS:
            None

        This function is designed to handle the generation of a molecule 
        primarily using the MolecularData.load() method of OpenFerimon. For
        this method to work properly, integral calculations must be made prior
        through the use of psi4 (and pyscf in later versions). The MolecularData
        object is then stored under the self.molecule parameter of this class.

        '''

        self.set_parameters(geometry, basis, multiplicity, charge)
        self.check_parameters()

        # Generate and populate instance of MolecularData data structure
        if(self.description != ''):
            molecule = MolecularData(self.geometry, self.basis, 
                                     self.multiplicity, self.charge,
                                     self.description)
        else:
            molecule = MolecularData(self.geometry, self.basis, 
                                     self.multiplicity, self.charge)

        # Determine if integrals have been previously generated or not
        if(not(os.path.exists(molecule.filename + '.hdf5')) or 
              (forceCalculation)):

            if(self.plugin == 'psi4'):
                # Run Psi4 calculation protocol
                run_scf, run_mp2, run_cisd, run_ccsd, run_fci = 1, 0, 0, 0, 0

                molecule = run_psi4(molecule, run_scf=run_scf, run_mp2=run_mp2,
                                    run_cisd=run_cisd, run_ccsd=run_ccsd,
                                    run_fci=run_fci)
            else:
                # Run PySCF calculation protocol
            	sys.exit('''\n\n --- OpenFermionWrapper Error ---\n 
                                Only Psi4 compatability is implemented''')
            # Save molecule for future, so that regeneration is not required
            molecule.save()
            molecule.load()
        else:
            molecule.load()

        self.molecule = molecule

    def create_hamiltonians(self):
        ''' 
        ARGS:
            None
        RETURNS: 
            None

        This function is designed to generate the fermion hamiltonian and
        qubit hamiltonians for both of the transformations (Bravyi-Kitaev and
        Jordan-Wigner) allowed in OpenFermion. To reduce logical errors, checks
        are made prior to the transform to ensure the molecule has been
        generated. If not, this function will call the load_molecule method. 
        Additionally, relevant values are stored in the molecular_hamiltonian,
        fermion_hamiltonian, and n_qubits members of this class.

        '''

        # Try to load the molecule if the client has not yet loaded it 
        if(not self.molecule):
            self.load_molecule()

        if((self.active_space_start) and (self.active_space_stop)):
        # Get the Hamiltonian in an active space.
            self.molecular_hamiltonian = self.molecule.get_molecular_hamiltonian(
                    occupied_indices=range(self.active_space_start),
                    active_indices=range(self.active_space_start, 
                                         self.active_space_stop))
        else:
            self.molecular_hamiltonian = self.molecule.get_molecular_hamiltonian()
        
        self.fermion_hamiltonian = get_fermion_operator(self.molecular_hamiltonian)
        
        self.n_qubits = count_qubits(self.fermion_hamiltonian)

    def set_ground_state_energy(self):
        ''' 
        ARGS:
            None
        RETURNS:
            None

        This function uses the fermion hamiltonian to calculate the ground state
        energy of the given molecule using the get_ground_state method of
        OpenFermion. The energy is stored as a floating point number in the 
        ground_state_energy variable of this class. 
        If the fermion hamiltonain is not generated, this function will call
        the create_hamiltonians method prior to calculating the energy. 
        '''
        if(not self.fermion_hamiltonian):
            self.create_hamiltonians()

        sparse_matrix = get_sparse_operator(self.fermion_hamiltonian)
        self.ground_state_energy = get_ground_state(sparse_matrix)[0]

    def perform_transform(self, mapping):
        ''' 
        ARGS:
            mapping - a string determining the desired tranformation method
        RETURNS:
            None

        This function performs the approprate transformation to the 
        FermionOperator representing the current molecule. 
        The currently supported mappings are:
            Bravyi-Kitaev ("BK")
            Jordan-Wigner ("JW")
        The Qubit Hamiltonians are stored under the qubit_hamiltonian_bk or
        qubit_hamiltonian_jw variables. 
        '''

        # Attempt to catch user logic error
        if(not self.fermion_hamiltonian):
            self.create_hamiltonians()

        if(mapping == "BK"):
            self.qubit_hamiltonian_bk = bravyi_kitaev(self.fermion_hamiltonian)
        elif(mapping == "JW"):
            self.qubit_hamiltonian_jw = jordan_wigner(self.fermion_hamiltonian)
        else:
            self.mapping_error(mapping)

    def generate_circuit(self, mapping=None):
        ''' 
        ARGS:
            Optional - mapping - a string determining the desired tranformation 
        RETURNS:
            None

        This function uses the built-in cability in OpenFermion to create
        a quantum circuit in QASM that represents the exponentiation of the 
        hamiltonian over time using the Trotter-Suzuki decomposition method.
        In keeping with the theme, checks are made to attempt to catch user
        logic error and still result in success. Note: If a mapping is
        not defined, recursive calls are made to use all implemented 
        transformations. 
        '''
        if(mapping == 'BK'):
            # Attempt to catch user logic error
            if(not self.qubit_hamiltonian_bk):
                self.perform_transform(mapping)
            self.bk_circuit = pauli_exp_to_qasm([self.qubit_hamiltonian_bk])
        elif(mapping == 'JW'):
            # Attempt to catch user logic error
            if(not self.qubit_hamiltonian_jw):
                self.perform_transform(mapping)
            self.jw_circuit = pauli_exp_to_qasm([self.qubit_hamiltonian_jw])
        else:
            # Attempt to catch user logic error
            self.generate_circuit('BK')
            self.generate_circuit('JW')

