
import RotorCalculation
import numpy as np

import multiprocessing


# setup
# type of graph
type_of_graph = 1
# 1 - rotating csys
# 0 - non rotating csys

include_excitation_in_graph = 0
# 1 - yes
# 0 - yes

# definition of system
SDM = np.zeros(4)  # matrix of mechanical data
# radial stiffness [N/m]
SDM[0] = 10000000.0
# rotating damping [N*m/s], not taken into account yet. Will be in future versions.
SDM[1] = 0
# non rotating damping [N*m/s]
SDM[2] = 400
# mass of rotor [kg]
SDM[3] = 5

# analysis definition
# length of simulation time
calculation_time = 0.2
# time step of analysis [s]
dt = 0.000005
# distance between geometrical center and center of mass. Excitations of vibrations.
mass_eccentricity = 0.01

# end of setup

calculation = RotorCalculation.CalculationOfRotor()

if __name__ == '__main__':
    p_1 = multiprocessing.Process(target=calculation.calculate_rotor_sym_for_range_of_frequencies,
                                  args=(45, 100, 50, SDM, calculation_time, dt, mass_eccentricity,
                                        type_of_graph, include_excitation_in_graph))
    p_2 = multiprocessing.Process(target=calculation.calculate_rotor_sym_for_range_of_frequencies,
                                  args=(100, 155, 50, SDM, calculation_time, dt, mass_eccentricity,
                                        type_of_graph, include_excitation_in_graph))
    p_3 = multiprocessing.Process(target=calculation.calculate_rotor_sym_for_range_of_frequencies,
                                  args=(155, 210, 50, SDM, calculation_time, dt, mass_eccentricity,
                                        type_of_graph, include_excitation_in_graph))
    p_4 = multiprocessing.Process(target=calculation.calculate_rotor_sym_for_range_of_frequencies,
                                  args=(210, 265, 50, SDM, calculation_time, dt, mass_eccentricity,
                                        type_of_graph, include_excitation_in_graph))
    p_5 = multiprocessing.Process(target=calculation.calculate_rotor_sym_for_range_of_frequencies,
                                  args=(265, 320, 50, SDM, calculation_time, dt, mass_eccentricity,
                                        type_of_graph, include_excitation_in_graph))

    p_1.start()
    p_2.start()
    p_3.start()
    p_4.start()
    p_5.start()

    p_1.join()
    p_2.join()
    p_3.join()
    p_4.join()
    p_5.join()
