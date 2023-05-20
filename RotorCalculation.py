import numpy as np
import Data_storage
from Ploting import Plot
from csys_rot import Ratate_csys

class CalculationOfRotor:

    """
    Calculate transient analysis of Jeffcott rotor:

    Assumptions:
    - analysis is linear
    - rotor is rotating in 3D space around axis Z, plane XY is where center of geometry and mass are moving
    - force created by inertia is applied to center of geometry and position of that center, is affected
    only by rotation of rotor not by deflection of rotor. It is assumed that deflections are small compared
     to unbalance. Similar to assumption in beams where deflection is assumed to small relative to beam span
      and do not change geometry of load.
    - there is only non rotating damping, so damping is applied to geometrical center and is proportional to
    speed of that center in plane xy
    - stiffness is proportional to deflection of center of geometry and is applied at that center in direction
    opposite to deformation
    - so all forces are acting on center of geometry.

    X matrix - used to store current state of the system:
                 excitation                                 response
        row 0 [0][0][0] x_poz_exc  [0][0][1] y_poz_exc   [0][1][0] x_poz_r [0][1][1] y_poz_r
        row 1 [1][0][0] x_vel_exc  [1][0][1] y_vel_exc   [1][1][0] x_vel_r [1][1][1] y_vel_r
        row 2 [2][0][0] x_acc_exc  [2][0][1] y_vel_exc   [2][1][0] x_acc_r [2][1][1] y_acc_r
        row 3 [3][0][0] x_Fsum_exc [3][0][1] y_Fsum_exc  [3][1][0] x_Fsum_r[3][1][1] y_Fsum_r
        row 4 [4][0][0] x_Fa_exc   [4][0][1] y_Fa_exc    [4][1][0] x_Fa_r  [4][1][1] y_Fa_r
        row 5 [5][0][0] x_Fk_exc   [5][0][1] y_Fk_exc    [5][1][0] x_Fk_r  [5][1][1] y_Fk_r
        row 6 [6][0][0] x_Fd_exc   [6][0][1] y_Fd_exc    [6][1][0] x_Fd_r  [6][1][1] y_Fd_r
    """

    def __init__(self, SDM, mass_eccentricity, calculation_time, type_of_graph, include_excitation, dt):

        self.dt = dt
        self.include_excitation = include_excitation
        self.type_of_graph = type_of_graph
        self.calculation_time = calculation_time
        self.mass_eccentricity = mass_eccentricity
        self.SDM = SDM

    def calculate_rotor_sym_for_range_of_frequencies(self, start_hz, end_hz, number_of_simulations):

        for frequency in np.linspace(start_hz, end_hz, number_of_simulations, False):
            exc_x_rot, exc_y_rot, res_x_rot, res_y_rot = self.__calculate_rotor_for_one_frequency(frequency)

            Plot.prepare_plot(frequency, self.include_excitation, exc_x_rot, exc_y_rot, res_x_rot, res_y_rot)

    def __calculate_rotor_for_one_frequency(self, Hz):
        t = 0

        excitation_log = Data_storage.DataStorage(self.SDM, self.dt, 0, True, False, False, False, False, False, False,
                                                  False)
        response_log = Data_storage.DataStorage(self.SDM, self.dt, 1, True, False, False, False, False, False, False,
                                                False)

        # X_matrix matrix - see doc
        X_matrix = np.zeros((7, 2, 2))

        for i in range(int(self.calculation_time / self.dt)):
            self.__calculate_one_cycle_of_jeff_linear_rotor_sym(t,X_matrix, Hz)
            phase = self.__calculate_phase_angle(t, Hz)

            excitation_log.log_data(X_matrix, i, phase, for_every_n=1)
            response_log.log_data(X_matrix, i, phase, for_every_n=1)
            t += self.dt

        if self.type_of_graph:
            exc_x_rot, exc_y_rot = Ratate_csys.rotate_excitation_by_phase_angle_of_loging(excitation_log)
            res_x_rot, res_y_rot = Ratate_csys.rotate_response_by_phase_angle_of_loging(response_log)
        else:
            exc_x_rot, exc_y_rot = excitation_log.get_pos_data()
            res_x_rot, res_y_rot = response_log.get_pos_data()

        return exc_x_rot, exc_y_rot, res_x_rot, res_y_rot

    def __calculate_one_cycle_of_jeff_linear_rotor_sym(self, T, X, Hz):
        self.__cal_new_poz(X, T, Hz)
        self.__cal_new_velocity(X, T,Hz)

        # s - stiffness
        self.__cal_new_s_forces_n(X)
        # d - damping
        self.__cal_new_d_forces(X)

        # a - acceleration
        self.__cal_new_a_forces(X)
        # f - force
        self.__cal_new_f_sum(X)

        self.__cal_new_accelerations(X, T, Hz)

    @staticmethod
    def __calculate_phase_angle(T, Hz):
        return 2 * np.pi * T * Hz


    def __calc_position_of_exc(self, T, Hz):
        return self.mass_eccentricity * np.cos(2 * np.pi * Hz * T),\
               self.mass_eccentricity * np.sin(2 * np.pi * Hz * T)


    def __calc_velocity_of_exc(self, T, Hz):
        return -self.mass_eccentricity * 2 * np.pi * Hz * np.sin(2 * np.pi * Hz * T),\
               self.mass_eccentricity * 2 * np.pi * Hz * np.cos(2 * np.pi * Hz * T)


    def __calc_acceleration_of_exc(self, T, Hz):
        return -self.mass_eccentricity * (2 * np.pi * Hz)**2 * np.cos(2 * np.pi * Hz * T),\
               -self.mass_eccentricity * (2 * np.pi * Hz)**2 * np.sin(2 * np.pi * Hz * T)

    def __cal_new_poz(self, X, T, Hz):
        # calc new position of excitation
        X[0][0][0], X[0][0][1] = self.__calc_position_of_exc(T, Hz)

        # calc new position of response
        X[0][1][0] = X[0][1][0] + X[1][1][0] * self.dt
        X[0][1][1] = X[0][1][1] + X[1][1][1] * self.dt

    def __cal_new_velocity(self, X, T, Hz):
        # calc new velocity of excitation
        X[1][0][0], X[1][0][1] = self.__calc_velocity_of_exc(T, Hz)

        # calc new velocity of response
        X[1][1][0] = X[1][1][0] + X[2][1][0] * self.dt
        X[1][1][1] = X[1][1][1] + X[2][1][1] * self.dt

    def __cal_new_accelerations(self, X, T, Hz):
        # excitation
        X[2][0][0], X[2][0][1] = self.__calc_acceleration_of_exc(T, Hz)

        # response
        X[2][1][0] = X[3][1][0]/self.SDM[3]
        X[2][1][1] = X[3][1][1]/self.SDM[3]


    def __cal_new_s_forces_n(self, X):
        # stiffness forces acting on excitation
        X[5][0][0] = self.SDM[0] * X[0][1][0]
        X[5][0][1] = self.SDM[0] * X[0][1][1]
        # stiffness forces acting on excitation
        X[5][1][0] = -self.SDM[0] * X[0][1][0]
        X[5][1][1] = -self.SDM[0] * X[0][1][1]


    def __cal_new_d_forces(self, X):
        X[6][1][0] = -self.SDM[2] * X[1][1][0]
        X[6][1][1] = -self.SDM[2] * X[1][1][1]

    def __cal_new_a_forces(self, X):
        X[4][1][0] = -self.SDM[3] * X[2][0][0]
        X[4][1][1] = -self.SDM[3] * X[2][0][1]

    @staticmethod
    def __cal_new_f_sum(X):
        # sum for excitation
        X[3][0][0] = X[4][0][0] + X[5][0][0] + X[6][0][0]
        X[3][0][1] = X[4][0][1] + X[5][0][1] + X[6][0][1]

        # sum for response
        X[3][1][0] = X[4][1][0] + X[5][1][0] + X[6][1][0]
        X[3][1][1] = X[4][1][1] + X[5][1][1] + X[6][1][1]
