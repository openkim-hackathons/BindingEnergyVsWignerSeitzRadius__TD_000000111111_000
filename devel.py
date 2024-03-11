from kim_python_utils.ase import CrystalGenomeTest
from numpy import multiply
from math import pi
from crystal_genome_util.aflow_util import get_stoich_reduced_list_from_prototype
from copy import deepcopy
from ase.build import bulk

class BindingEnergyVsWignerSeitzRadius(CrystalGenomeTest):
    def _calculate(self, structure_index: int, a_min_frac: float, a_max_frac: float, N: int):
        """
        structure_index:
            KIM tests can loop over multiple structures (i.e. crystals, molecules, etc.). 
            This indicates which is being used for the current calculation.

        a_min_frac:
            A fraction which indicates the smallest lattice constant for which the Test
            will attempt to compute the cohesive energy.  The smallest lattice spacing
            will be calculated as a_min = a_min_frac*a_0.  Must be strictly greater than
            zero and strictly less than one.

        a_max_frac:
            A fraction which indicates the largest lattice constant for which the Test
            will attempt to compute the cohesive energy.  The largest lattice spacing
            will be calculated as a_max = a_max_frac*a_0.  Must be strictly greater than
            one.

        N
            The number of lattice spacings sampled
        """
        atoms = self.atoms[structure_index]

        average_wigner_seitz_radius = []
        binding_potential_energy_per_atom = []
        binding_potential_energy_per_formula = []
        a_frac_range = a_max_frac - a_min_frac
        a_frac_step = a_frac_range/(N-1)
        original_cell = deepcopy(atoms.cell) # have to deepcopy! Otherwise it refers to the same object
        natoms = len(atoms)

        for i in range(N):
            a_frac = a_frac_step*i + a_min_frac
            print("evaluating a_frac = " + str(a_frac) + " ...")
            atoms.set_cell(multiply(original_cell,a_frac),scale_atoms = True)
            average_wigner_seitz_radius.append(((atoms.get_volume()/natoms)*(3/4)/pi)**(1/3))            
            pe_per_atom = atoms.get_potential_energy()/natoms
            binding_potential_energy_per_atom.append(pe_per_atom)
            binding_potential_energy_per_formula.append(pe_per_atom*sum(get_stoich_reduced_list_from_prototype(self.prototype_label)))

        self._add_property_instance("binding-energy-relation-crystal")
        self._add_common_crystal_genome_keys_to_current_property_instance(structure_index,write_stress=False,write_temp=False) # last two default to False
        self._add_key_to_current_property_instance("average-wigner-seitz-radius",average_wigner_seitz_radius,"angstrom")
        self._add_key_to_current_property_instance("binding-potential-energy-per-atom",binding_potential_energy_per_atom,"eV")
        self._add_key_to_current_property_instance("binding-potential-energy-per-formula",binding_potential_energy_per_formula,"eV")

        with open("debug"+str(structure_index)+".txt","w") as f:
            for r,u in zip(average_wigner_seitz_radius,binding_potential_energy_per_atom):
                f.write(str(r)+" "+str(u)+"\n")


# This queries for equilibrium structures in this prototype and builds atoms
# test = BindingEnergyVsWignerSeitzRadius(model_name="MEAM_LAMMPS_KoJimLee_2012_FeP__MO_179420363944_002", stoichiometric_species=['Fe','P'], prototype_label='AB_oP8_62_c_c')
                
# Alternatively, for debugging, give it atoms object or a list of atoms objects
atoms = bulk('NaCl','rocksalt',a=4.58)
test = BindingEnergyVsWignerSeitzRadius(model_name="Sim_LAMMPS_EIM_Zhou_2010_BrClCsFIKLiNaRb__SM_259779394709_000", atoms=atoms)
test(a_min_frac=0.75, a_max_frac=1.5, N=51)