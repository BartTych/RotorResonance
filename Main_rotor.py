
import RotorCalculation
import Data_storage
import numpy as np
import matplotlib.pyplot as plt

import multiprocessing



# definition of system
# setup
#type of graph
type_of_graph = 1
# 1 - rotating csys
# 0 - non rotating csys

include_excitation = 1
#1 - yes
#0 - yes

SDM = np.zeros((4))
# stifness [N/m]
SDM[0] = 10005350.0
# rotating damping [N*m/s]
SDM[1] = 0
# nonrotating damping [N*m/s]
SDM[2] = 400
# mass [kg]
SDM[3] = 5

# analysis definition
calculation_time = 0.2
phase = 0.0
dt = 0.000005
mass_ecentricity = 0.01

# end of setup

calculation = RotorCalculation.CalculationOfRotor()
#calculation.caltulate(45,100, 1000, SDM, calculation_time, dt, mass_ecentricity)

if __name__ == '__main__':
    p_1 = multiprocessing.Process(target=calculation.caltulate, args=(45,100, 50, SDM, calculation_time, dt, mass_ecentricity,type_of_graph,include_excitation))
    p_2 = multiprocessing.Process(target=calculation.caltulate, args=(100,155, 50, SDM, calculation_time, dt, mass_ecentricity,type_of_graph,include_excitation))
    p_3 = multiprocessing.Process(target=calculation.caltulate, args=(155,210, 50, SDM, calculation_time, dt, mass_ecentricity, type_of_graph,include_excitation))
    p_4 = multiprocessing.Process(target=calculation.caltulate, args=(210,265, 50, SDM, calculation_time, dt, mass_ecentricity, type_of_graph,include_excitation))
    p_5 = multiprocessing.Process(target=calculation.caltulate, args=(265,320, 50, SDM, calculation_time, dt, mass_ecentricity, type_of_graph,include_excitation))

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

