import numpy as np
import Data_storage
import matplotlib.pyplot as plt
from matplotlib import figure
import gc
import os


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

    def __init__(self):
        pass

    def calculate_rotor_sym_for_range_of_frequencies(self, start_hz, end_hz, number_of_simulations, SDM,
                                                     calculation_time, dt, mass_eccentricity, type_of_graph,
                                                     include_excitation):

        for n in np.linspace(start_hz, end_hz, number_of_simulations, False):
            exc_x_rot, exc_y_rot, res_x_rot, res_y_rot = self.calculate_rotor(n, dt, calculation_time,
                                                                              SDM, mass_eccentricity, type_of_graph)
            # prepare and save plot
            self.prepare_plot(n, include_excitation, exc_x_rot, exc_y_rot, res_x_rot, res_y_rot)

    def calculate_rotor(self, Hz, dt, calculation_time, SDM, mass_eccentricity, rotate_result):
        t = 0
        excitation_log = Data_storage.log_storage_rotor(SDM, dt, 0, True, False, False, False, False, False, False,
                                                        False)
        response_log = Data_storage.log_storage_rotor(SDM, dt, 1, True, False, False, False, False, False, False,
                                                      False)
        X = np.zeros((7, 2, 2))
        for i in range(int(calculation_time / dt)):
            self.calculate_one_cycle_of_jeff_linear_rotor_sym(t, dt, X, SDM, Hz, mass_eccentricity)
            phase = self.calculate_phase_angle(t, Hz)

            excitation_log.log_data(X, i, phase, for_every_n=1)
            response_log.log_data(X, i, phase, for_every_n=1)
            t += dt
        # rotation of calculated data csys
        if rotate_result:
            response_log.rotate_response_by_phase_angle_of_loging()
            excitation_log.rotate_excitation_by_phase_angle_of_loging()

        # get calculated data
        exc_x_rot, exc_y_rot = excitation_log.get_pos_data()
        res_x_rot, res_y_rot = response_log.get_pos_data()

        return exc_x_rot, exc_y_rot, res_x_rot, res_y_rot

    @staticmethod
    def prepare_plot(Hz, include_excitation, exc_x_rot, exc_y_rot, res_x_rot, res_y_rot):
        fig = figure.Figure(figsize=(25, 25))
        ax = fig.subplots(1)
        ax.set_aspect('equal')
        ax.plot(res_x_rot, res_y_rot, label="position of response", color="orangered")
        if include_excitation:
            ax.scatter(exc_x_rot, exc_y_rot, s=15, label="position of excitation")
        ax.set_xlabel('position x [m]', fontsize=20)
        ax.set_ylabel('position y [m]', fontsize=20)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        fig.legend(loc='upper right', fontsize=20)
        file = f'{os.getcwd()}/results'
        if not os.path.exists(file):
            os.mkdir(file)
        fig.savefig(f"{file}/{round(Hz, 2):.2f}_Hz.png")

        ax.clear()
        plt.close('all')
        gc.collect()

        print(f"calculated frequency :{Hz}")

    def calculate_one_cycle_of_jeff_linear_rotor_sym(self, T, dt, X, SDM, Hz, mass_eccentricity):
        self.cal_new_poz(X, T, dt, Hz, mass_eccentricity)
        self.cal_new_velocity(X, T, dt, Hz, mass_eccentricity)

        # s - stiffness
        self.cal_new_s_forces_n(X, SDM)
        # d - damping
        self.cal_new_d_forces(X, SDM)

        # a - acceleration
        self.cal_new_a_forces(X, SDM)
        # f - force
        self.cal_new_f_sum(X)

        self.cal_new_accelerations(X, SDM, T, dt, Hz, mass_eccentricity)

    @staticmethod
    def calculate_phase_angle(T, Hz):
        return 2 * np.pi * T * Hz

    @staticmethod
    def calc_position_of_exc(T, Hz, ecc):
        return ecc * np.cos(2 * np.pi * Hz * T), ecc * np.sin(2 * np.pi * Hz * T)

    @staticmethod
    def calc_velocity_of_exc(T, Hz, ecc):
        return -ecc * 2 * np.pi * Hz * np.sin(2 * np.pi * Hz * T), ecc * 2 * np.pi * Hz * np.cos(2 * np.pi * Hz * T)

    @staticmethod
    def calc_acceleration_of_exc(T, Hz, ecc):
        return -ecc * (2 * np.pi * Hz)**2 * np.cos(2 * np.pi * Hz * T),\
               -ecc * (2 * np.pi * Hz)**2 * np.sin(2 * np.pi * Hz * T)

    def cal_new_poz(self, X, T, dt, Hz, ecc):
        # calc new position of excitation
        X[0][0][0], X[0][0][1] = self.calc_position_of_exc(T, Hz, ecc)

        # calc new position of response
        X[0][1][0] = X[0][1][0] + X[1][1][0] * dt
        X[0][1][1] = X[0][1][1] + X[1][1][1] * dt

    def cal_new_velocity(self, X, T, dt, Hz, ecc):
        # calc new velocity of excitation
        X[1][0][0], X[1][0][1] = self.calc_velocity_of_exc(T, Hz, ecc)

        # calc new velocity of response
        X[1][1][0] = X[1][1][0] + X[2][1][0] * dt
        X[1][1][1] = X[1][1][1] + X[2][1][1] * dt

    def cal_new_accelerations(self, X, SDM, T, dt, Hz, ecc):
        # excitation
        X[2][0][0], X[2][0][1] = self.calc_acceleration_of_exc(T, Hz, ecc)

        # response
        X[2][1][0] = X[3][1][0]/SDM[3]
        X[2][1][1] = X[3][1][1]/SDM[3]

    @staticmethod
    def cal_new_s_forces_n(X, SDM):
        # stiffness forces acting on excitation
        X[5][0][0] = SDM[0] * X[0][1][0]
        X[5][0][1] = SDM[0] * X[0][1][1]
        # stiffness forces acting on excitation
        X[5][1][0] = -SDM[0] * X[0][1][0]
        X[5][1][1] = -SDM[0] * X[0][1][1]

    @staticmethod
    def cal_new_d_forces(X, SDM):
        X[6][1][0] = -SDM[2] * X[1][1][0]
        X[6][1][1] = -SDM[2] * X[1][1][1]

    @staticmethod
    def cal_new_a_forces(X, SDM):
        X[4][1][0] = -SDM[3] * X[2][0][0]
        X[4][1][1] = -SDM[3] * X[2][0][1]

    @staticmethod
    def cal_new_f_sum(X):
        # sum for excitation
        X[3][0][0] = X[4][0][0] + X[5][0][0] + X[6][0][0]
        X[3][0][1] = X[4][0][1] + X[5][0][1] + X[6][0][1]

        # sum for response
        X[3][1][0] = X[4][1][0] + X[5][1][0] + X[6][1][0]
        X[3][1][1] = X[4][1][1] + X[5][1][1] + X[6][1][1]
