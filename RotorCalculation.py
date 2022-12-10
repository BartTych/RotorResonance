import numpy as np
import Data_storage
import matplotlib.pyplot as plt
from matplotlib import figure
import gc
class CalculationOfRotor:

    '''
    Calculate one cycle of transient analysis of Jeffcott rotor
    Analysis is linear

    Assumptions:
    - rotor is rotating in 3D space around axis Z, plane XY is where center of geometry and mass are moving
    - force created by inertia is applied to center of geometry and position of that center, is affected
    only by rotation of rotor not by deflection of rotor. It is assumed that deflections are small compared to unbalance.
    Similar to assumption in beams where deflection is assumed to small relative to beam span and do not change geometry
     of load.
    - there is only non rotating damping, so damping is applied to geometrical center and is proportional to
    speed of that center in plane xy
    - stiffness is proportional to deflection of center of geometry and is applied at that center in direction opposite to
    deformation
    - so all forces are acting on center of geometry.

    X matrix:
                 excitation                                 response
        wiersz 0 [0][0][0] x_poz_exc  [0][0][1] y_poz_exc   [0][1][0] x_poz_c [0][1][1] y_poz_c
        wiersz 1 [1][0][0] x_vel_exc  [1][0][1] y_vel_exc   [1][1][0] x_vel_c [1][1][1] y_vel_c
        wiersz 2 [2][0][0] x_acc_exc  [2][0][1] y_vel_exc   [2][1][0] x_acc_c [2][1][1] y_acc_c
        wiersz 3 [3][0][0] x_Fsum_exc [3][0][1] y_Fsum_exc  [3][1][0] x_Fsum_c[3][1][1] y_Fsum_c
        wiersz 4 [4][0][0] x_Fa_exc   [4][0][1] y_Fa_exc    [4][1][0] x_Fa_c  [4][1][1] y_Fa_c
        wiersz 5 [5][0][0] x_Fk_exc   [5][0][1] y_Fk_exc    [5][1][0] x_Fk_c  [5][1][1] y_Fk_c
        wiersz 6 [6][0][0] x_Fd_exc   [6][0][1] y_Fd_exc    [6][1][0] x_Fd_c  [6][1][1] y_Fd_c


    macierz SDM - for rotor

    n - number of spring

    # sztywnosc
    SDM[0][n] = 3005350.0

    # tlumienie rotating
    SDM[1][n] = 0 - narazie 0 dla prostego modelu
    bo to wprowadza zlozone zachowanie do modelu
    i wymaga rozbudowania modelu obliczeniowego

    # tlumienie nonrotating
    SDM[2][n] = 820.0

    # masa
    SDM[3][n] = 5
    '''

    def __init__(self):
        pass

    def calculate_rotor(self,Hz, dt, calculation_time, SDM, mass_ecentricity, rotete_result):
        T = 0
        excitation_log = Data_storage.log_storage_rotor(SDM, dt, 0, True, False, False, False, False, False, False,
                                                        False)
        response_log = Data_storage.log_storage_rotor(SDM, dt, 1, True, False, False, False, False, False, False,
                                                      False)

        X = np.zeros((7, 2, 2))
        for i in range(int(calculation_time / dt)):
            self.calculate_one_cycle_of_jeff_linear_rotor(T, dt, X, SDM, Hz, mass_ecentricity)
            phase = self.calculate_phase_angle(T, Hz)

            excitation_log.log_data(X, i, phase, for_every_n=1)
            response_log.log_data(X, i, phase, for_every_n=1)
            T += dt

        # exc_x, exc_y = excitation_log.get_pos_data()
        # res_x, res_y = response_log.get_pos_data()

        # rotation of calculated data
        if rotete_result:
            response_log.rotate_response_by_phase_angle_of_loging()
            excitation_log.rotate_excitation_by_phase_angle_of_loging()

        # get calculated data
        exc_x_rot, exc_y_rot = excitation_log.get_pos_data()
        res_x_rot, res_y_rot = response_log.get_pos_data()

        return exc_x_rot, exc_y_rot,res_x_rot, res_y_rot

    def prepare_plot(self, Hz, SDM,exc_x_rot,exc_y_rot, res_x_rot, res_y_rot):
        #fig, ax = plt.subplots(1, figsize=(25, 25))
        fig = figure.Figure(figsize=(25, 25))
        ax = fig.subplots(1)

        ax.set_aspect('equal')
        # ax.plot(res_x_rot, res_y_rot, label="pozycja odpowiedzi",color="magenta")
        ax.plot(res_x_rot, res_y_rot, label="pozycja odpowiedzi", color="orangered")
        #ax.scatter(exc_x_rot,exc_y_rot, s = 15, label="pozycja wymuszenia")
        ax.set_xlabel('pozycja x [m]', fontsize=20)
        ax.set_ylabel('pozycja y [m]', fontsize=20)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        fig.legend(loc='upper right', fontsize=20)
        file = r'/Users/bart/python/Rezonans_p3_11/wyniki/rotor'
        fig.savefig(f"{file}/{round(Hz, 2):.2f}_Hz_{SDM[2]}_tlumienie.png")

        ax.clear()

        plt.figure().clear()
        fig.clear()


        plt.cla()
        plt.clf()
        plt.close('all')

        gc.collect()

        print(f"gotowe {Hz}")

    def caltulate(self,start,end, ilosc, SDM, calculation_time, dt, mass_ecentricity):

        for n in np.linspace(start,end,ilosc,False):
            Hz = n

            # container to calculate sym for given Hz and model args
            exc_x_rot, exc_y_rot, res_x_rot, res_y_rot = self.calculate_rotor(n, dt, calculation_time, SDM, mass_ecentricity, rotete_result=True)

            # prepare and save plot
            self.prepare_plot(n, SDM,exc_x_rot, exc_y_rot, res_x_rot, res_y_rot)
            del(exc_x_rot)
            del(exc_y_rot)
            del(res_x_rot)
            del(res_y_rot)


    def calculate_one_cycle_of_jeff_linear_rotor(self,T, dt, X, SDM, Hz, mass_ecentricity):
        self.cal_new_poz(X, T, dt, Hz, mass_ecentricity)
        self.cal_new_velocity(X, T, dt, Hz, mass_ecentricity)

        #s - stifness
        self.cal_new_s_forces_n(X, SDM)
        #d - damping
        self.cal_new_d_forces(X, SDM)

        # a - acceleration
        self.cal_new_a_forces(X, SDM)
        # f - force
        self.cal_new_f_sum(X)

        self.cal_new_accelerations(X, SDM, T, dt, Hz, mass_ecentricity)




    def calculate_phase_angle(self, T, Hz):
        return 2 * np.pi * T * Hz

    def calc_position_of_exc(self, T, Hz, ecc):
        return ecc * np.cos(2 * np.pi * Hz * T),ecc * np.sin(2 * np.pi * Hz * T)

    def calc_velocity_of_exc(self, T, Hz, ecc):
            return -ecc * 2 * np.pi * Hz * np.sin(2 * np.pi * Hz * T),ecc * 2 * np.pi * Hz * np.cos(2 * np.pi * Hz * T)

    def calc_acceleration_of_exc(self, T, Hz, ecc):
            return -ecc * (2 * np.pi * Hz)**2 * np.cos(2 * np.pi * Hz * T),-ecc * (2 * np.pi * Hz)**2 * np.sin(2 * np.pi * Hz * T)

    def cal_new_poz(self, X, T, dt, Hz, ecc):
        #licze nowa pozycje wymuszenia
        X[0][0][0], X[0][0][1] = self.calc_position_of_exc(T, Hz, ecc)

        #licze nowa pozycje odpowiedzi, x i y
        X[0][1][0] = X[0][1][0] + X[1][1][0] * dt
        X[0][1][1] = X[0][1][1] + X[1][1][1] * dt


    def cal_new_velocity(self, X, T, dt, Hz, ecc):
        # licze nowa predkosc wymuszenia
        X[1][0][0], X[1][0][1] = self.calc_velocity_of_exc(T, Hz, ecc)

        #licze nowa predkosc odpowiedzi, x i y
        X[1][1][0] = X[1][1][0] + X[2][1][0] * dt
        X[1][1][1] = X[1][1][1] + X[2][1][1] * dt

    def cal_new_accelerations(self, X, SDM,T, dt, Hz, ecc):
        # licze nowa predkosc wymuszenia
        X[2][0][0], X[2][0][1] = self.calc_acceleration_of_exc(T, Hz, ecc)

        # licze nowe przyspieszenie odpowiedzi, x i y
        X[2][1][0] = X[3][1][0]/SDM[3]
        X[2][1][1] = X[3][1][1]/SDM[3]

    def cal_new_s_forces_n(self, X, SDM):
        #licze sily sprezystosci dzialajace na wymuszenie
        X[5][0][0] = SDM[0] * X[0][1][0]
        X[5][0][1] = SDM[0] * X[0][1][1]

        X[5][1][0] = -SDM[0] * X[0][1][0]
        X[5][1][1] = -SDM[0] * X[0][1][1]

    def cal_new_d_forces(self, X, SDM):
        # licze sily tlumienia dzialjace na wymuszenie
        # narzie bedzie zero tylko tutaj

        #licze sily tlumienia dzilajace na odpowiedz

        X[6][1][0] = -SDM[2] * X[1][1][0]
        X[6][1][1] = -SDM[2] * X[1][1][1]

    def cal_new_a_forces(self, X, SDM):
        # licze dla wymuszenia
        # podstawa jest przyspieszenie srodka masy wymuszenia
        # wiec jest to wymuszenie
        # a odpowiedz srodka masy to co innego
        # tego mi brakuje
        X[4][1][0] = -SDM[3] * X[2][0][0]
        X[4][1][1] = -SDM[3] * X[2][0][1]

    def cal_new_f_sum(self,X):
        #licze sume dla wymuszenia
        X[3][0][0] = X[4][0][0] + X[5][0][0] + X[6][0][0]
        X[3][0][1] = X[4][0][1] + X[5][0][1] + X[6][0][1]

        #licze sume dla odpowiedzi
        X[3][1][0] = X[4][1][0] + X[5][1][0] + X[6][1][0]
        X[3][1][1] = X[4][1][1] + X[5][1][1] + X[6][1][1]
