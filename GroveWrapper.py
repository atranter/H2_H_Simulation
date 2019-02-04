from pyquil import Program, get_qc
from pyquil.api import WavefunctionSimulator
import pyquil.api as api
from pyquil.gates import *

from forestopenfermion import qubitop_to_pyquilpauli

from grove.pyvqe.vqe import VQE
from scipy.optimize import minimize
import numpy as np

import sys

from pyquil.paulis import sZ, exponentiate

from OpenFermionWrapper import OpenFermionWrapper

def get_H2_hamil(bondlength):
	bondlength = bondlength[0]
	molecule = OpenFermionWrapper()
	geometry = [("H", (0,0,0)), ("H", (0,0,bondlength))]
	molecule.load_molecule(geometry=geometry, basis="sto-3g", charge=0,
						   multiplicity=1, forceCalculation=True)
	molecule.perform_transform("JW")
	return molecule.qubit_hamiltonian_jw


def get_circuit(qubit_hamil):
	# convert to pauli operator syntax
	pyquil_generator = qubitop_to_pyquilpauli(qubit_hamil)
	
	# create empty circuit
	pyquil_program = Program()
	
	# exponentiate each term in the hamiltonian and append to circuit
	for term in pyquil_generator.terms:
		pyquil_program += exponentiate(term)

	return pyquil_program, pyquil_generator


def callable(bondlength):
	hamil = get_H2_hamil(bondlength)
	circuit, hamil = get_circuit(hamil)
	return circuit


def find_best():
	qvm = api.QVMConnection()

	vqe = VQE(minimizer=minimize, minimizer_kwargs={'method': 'nelder-mead'})
	initial = [float(2)]
	hamil = get_H2_hamil([0.7414])
	circuit, hamil = get_circuit(hamil)

	print(vqe.vqe_run(callable, hamil, initial, samples=5, qvm=qvm))


def ground_state_energy(qubit_hamil):
	# get a connection with the local server
	qvm = api.QVMConnection()

	# circuit is the exponentiation of the hamiltonian, while hamil is the 
	# hamiltonian in the pauli operator syntax
	circuit, hamil= get_circuit(qubit_hamil)

	# create an instance of the vqe
	vqe = VQE(minimizer = minimize, 
			  minimizer_kwargs={'method': 'nelder-mead'})
	hamil = sZ(1)
	return (vqe.expectation(circuit, hamil, None, qvm))

def expectation():
	lengths = np.linspace(0.1, 2, num = 50)
	expectation_values = []
	for length in lengths:
		hamil = get_H2_hamil([float(length)])
		expectation_values.append(ground_state_energy(hamil))
	for value in zip(lengths, expectation_values):
		print(value)

def execute():
	bondlength = float(0.5455)
	hamil = get_H2_hamil([bondlength])
	circuit, hamil = get_circuit(hamil)
	qc = get_qc('9q-square-qvm')
	# circuit = Program(H(0), CNOT(0,1))
	half_circuit = Program()
	for i, gate in enumerate(circuit):
		half_circuit += gate
		if(i == len(circuit)/2):
			break

	wfn = WavefunctionSimulator().wavefunction(circuit)
	print(wfn[0])

	result = qc.run_and_measure(half_circuit, trials=10)
	print(result[0])
	print(result[1])
	print(result[2])
	print(result[3])	
	print(result[4])
	print(result[5])
	print(result[6])
	print(result[7])
	print(result[8])

# find_best()
# expectation()
execute()