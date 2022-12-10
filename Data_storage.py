import matplotlib.pyplot as plt
import numpy
import numpy as np
import bisect

class amplitude_frequency_data:

    amplitude:list[float] # dane amplitudy
    frequency:list[float] # dane czestotliwosci

    min_fr = 0
    max_fr = 0
    excitation_A = 0
    direction: bool

    def __init__(self, excitation_A:float, min_fr:int, max_fr:int):
        self.excitation_A = excitation_A
        self.min_fr = min_fr
        self.max_fr = max_fr


    def add_data(self, amplitude:list, frequency:list):
        self.amplitude = amplitude
        self.frequency = frequency

    def get_data(self):
        return self.amplitude, self.frequency



    def get_frequency_range(self):
        return self.min_fr, self.max_fr

class amplitudeData:
    amplitude: list[float]  # dane amplitudy

    min_fr = 0
    max_fr = 0
    excitation_A = 0
    direction: bool

    def __init__(self, excitation_A: float, min_fr: int, max_fr: int):
        self.excitation_A = excitation_A
        self.min_fr = min_fr
        self.max_fr = max_fr

    def add_data(self, amplitude):
        self.amplitude = amplitude

    def get_data(self):
        return self.amplitude


class FFT_result:
    amplitude: list[float]  # dane amplitudy
    czestosci = []

    min_fr = 0
    max_fr = 0
    excitation_A = 0
    direction: bool

    def __init__(self, excitation_A: float, min_fr: int, max_fr: int):
        self.excitation_A = excitation_A
        self.min_fr = min_fr
        self.max_fr = max_fr

    def add_data(self, amplitude, czestosci):
        self.amplitude = amplitude
        self.czestosci = czestosci

    def get_data(self):
        return self.amplitude, self.czestosci


class log_storage:
    ''' Clasa twozy log danych z macierzy X co n-ty krok
        schemat macierzy X:
                  kol 0       kol 1     ..   kol n
        wiersz 0 poz_exc    poz_node_1  ..   poz_node_n
        wiersz 1 vel_exc    vel_node_1  ..   vel_node_n
        wiersz 2 acc_exc    acc_node_1  ..   acc_node_n
        wiersz 3 F_exc      F_node 1    ..   F_node_n
        wiersz 4 D_exc      D_node 1    ..   D_node_n

        F - force
        D - damping force
     '''

    pos = bool
    vel = bool
    acc = bool

    force = bool
    dumping = bool
    inertia = bool
    force_control_sum = bool
    power = bool

    time = bool

    pos_log = list
    vel_log = list
    acc_log = list

    force_log = list
    dumping_log = list
    inertia_log = list

    force_control_sum_log = list

    SDM = np.zeros((3, 2))
    dt = float

    def __init__(self, SDM, dt, position, velocity, acceleration, force, dumping, inertia, force_control_sum,power):
        self.SDM = SDM
        self.dt = dt
        self.pos = position
        self.acc = acceleration
        self.vel = velocity

        self.force = force
        self.dumping = dumping
        self.inertia = inertia
        self.force_control_sum = force_control_sum
        self.power = power

        self.time = []
        self.cycle = 0

        self.phase_log = []
        self.pos_log = []
        self.vel_log = []
        self.acc_log = []

        self.force_log = []
        self.dumping_log = []
        self.inertia_log = []
        self.force_control_sum_log =[]
        self.moment_power = []

    def log_data(self, X, i, phase,T, for_every_n, node_number):
        if np.mod(i, for_every_n) == 0:

            self.log_phase(phase)
            if self.pos:
                self.log_position(X, node_number)
            if self.vel:
                self.log_velocity(X, node_number)
            if self.acc:
                self.log_acceleration(X, node_number)
            if self.force:
                self.log_force(X, node_number)
            if self.dumping:
                self.log_dumping_force(X, node_number)
            if self.inertia:
                self.log_inertia(X, node_number)
            if self.force_control_sum:
                self.log_force_control_sum(X, node_number)
            if self.power:
                self.log_power(X, node_number)

                self.log_time(T)

            #self.cycle = i


    def clear_all_data(self):
        self.pos_log.clear()
        self.pos_log.clear()
        self.phase_log.clear()
        #self.x_acc_log.clear()
        #self.y_acc_log.clear()

    def log_phase(self, phase):
        self.phase_log.append(phase)

    def log_position(self,X ,node_number):
        self.pos_log.append(X[0][node_number])

    def log_velocity(self, X, node_number):
        self.vel_log.append(X[1][node_number])

    def log_acceleration(self,X, node_number):
        self.acc_log.append(X[2][node_number])

    def log_force(self,X, node_number):
        self.force_log.append(X[3][node_number])

    def log_dumping_force(self,X, node_number):
        self.dumping_log.append(X[4][node_number])

    def log_inertia(self, X, node_number):
        self.inertia_log.append(-X[2][node_number] * self.SDM[2][0])

    def log_force_control_sum(self,X,node_number):
        self.force_control_sum_log.append(X[3][node_number] + X[4][node_number] - X[2][node_number] * self.SDM[2][0])
    def log_power(self, X, node_number):
        self.moment_power.append(X[3][node_number] * X[1][node_number])
    def log_time(self,T):
        self.time.append(T)
    def get_time(self):
        return self.time
    def get_cycle(self):
        return self.cycle

    def get_pos_data(self):
        return self.pos_log
    def get_vel_data(self):
        return self.vel_log
    def get_acc_data(self):
        return self.acc_log
    def get_force_data(self):
        return self.force_log
    def get_damping_data(self):
        return self.dumping_log
    def get_inertia_data(self):
        return self.inertia_log
    def get_force_control_sum(self):
        return self.force_control_sum_log
    def get_moment_power(self):
        return self.moment_power

    def get_restored_X_matrix_part(self):
        X_partial = np.zeros((5, 1))
        if len(self.pos_log) != 0:
            X_partial[0] = self.pos_log[-1]
            X_partial[1] = self.vel_log[-1]
            X_partial[2] = self.acc_log[-1]
            X_partial[3] = self.force_log[-1]
            X_partial[4] = self.dumping_log[-1]

        return X_partial

    def get_last_two_position_logs(self):
        return self.pos_log[-1], self.pos_log[-2]
    def get_last_two_velocity_logs(self):
        return self.vel_log[-1], self.vel_log[-2]

    '''calculate numeric integration of moment power thru whole recorded period'''
    def calculate_power_per_cycle(self):
        integration = 0.0

        for n in self.moment_power:
            integration += n * self.dt
        return integration

    def calculate_normalization_coeff_V_vs_X(self):
        #licze zakres polozenia
        max_x = max(self.pos_log)
        min_x = min(self.pos_log)

        max_v = max(self.vel_log)
        min_v = min(self.vel_log)

        return 1/((max_v - min_v) / (max_x - min_x))

    def calculate_normalization_coeff_V_vs_X_by_frequencu(self, Fr):
        return 1 / (2 * np.pi * Fr)

    def normalize_velocity_by_coefficent(self, coeff):
        self.vel_log = [n * coeff for n in self.vel_log]


    def rotate_all_logged_data_by_phase_angle(self):
        for i,n in enumerate(self.phase_log):
            self.__rotate_response_data_point_by_phase_angle(i, n)

    def get_phase_angle_of_last_log(self):
        return self.phase_log[-1]

    def __rotate_response_data_point_by_phase_angle(self, i, n):
        """ Ta metoda dziala tylko jesli osie sa zeskalowane i odpowiedz jest okręgiem """
        x = self.pos_log[i]
        y = self.vel_log[i]

        # licze promien odpowiedzi
        r = np.hypot(x, y)
        # licze kat fazowy stanu odpowiedzi
        response_phase = np.angle([complex(x,y)])

        # odejmuje kat wymuszenia od kata odpowiedzi
        new_phase = response_phase + n
        # licze nowe x i y
        # przypisuje wartosci w listach danych
        self.pos_log[i] = r * np.cos(new_phase)[0]
        self.vel_log[i] = r * np.sin(new_phase)[0]

    def return_max_absolute_value_stored_of_pos(self):
        m = abs(max(self.pos_log))
        mi = abs(min(self.pos_log))

        return abs(max(m,mi))

    def return_max_absolute_value_stored_of_vel(self):
        m = abs(max(self.vel_log))
        mi = abs(min(self.vel_log))

        return abs(max(m, mi))

    def trim_all_data_by_end_time(self,T):
        # licze index danych odpowiedni dla T
        index = bisect.bisect_left(self.time,T)
        if self.pos:
            self.pos_log = self.pos_log[:index]
        if self.vel:
            self.vel_log = self.vel_log[:index]
        if self.acc:
            self.acc_log = self.acc_log[:index]

        self.phase_log = self.phase_log[:index]
        if self.force:
            self.force_log = self.force_log[:index]
        if self.dumping:
            self.dumping_log = self.dumping_log[:index]
        if self.inertia_log:
            self.inertia_log = self.inertia_log[:index]
        if self.force_control_sum:
            self.force_control_sum_log = self.force_control_sum_log[:index]
        if self.power:
            self.moment_power = self.moment_power[:index]


        #przycinam wszystkie listy z danymi
class log_storage_rotor:

    ''' Clasa twozy log danych z macierzy X co n-ty krok
        schemat macierzy X dla rotora:

         X matrix:
                 wymuszenie                                 odpowiedz
        wiersz 0 [0][0][0] x_poz_exc  [0][0][1] y_poz_exc   [0][1][0] x_poz_c [0][1][1] y_poz_c
        wiersz 1 [1][0][0] x_vel_exc  [1][0][1] y_vel_exc   [1][1][0] x_vel_c [1][1][1] y_vel_c
        wiersz 2 [2][0][0] x_acc_exc  [2][0][1] y_vel_exc   [2][1][0] x_acc_c [2][1][1] y_acc_c
        wiersz 3 [3][0][0] x_Fsum_exc [3][0][1] y_Fsum_exc  [3][1][0] x_Fsum_c[3][1][1] y_Fsum_c
        wiersz 4 [4][0][0] x_Fa_exc   [4][0][1] y_Fa_exc    [4][1][0] x_Fa_c  [4][1][1] y_Fa_c
        wiersz 5 [5][0][0] x_Fk_exc   [5][0][1] y_Fk_exc    [5][1][0] x_Fk_c  [5][1][1] y_Fk_c
        wiersz 6 [6][0][0] x_Fd_exc   [6][0][1] y_Fd_exc    [6][1][0] x_Fd_c  [6][1][1] y_Fd_c

        x skladowa x danej wartosci
        y skladowa y danej wartosci
        F - force
        Fa - sila wywolana przyspieszeniem
        Fk - sila wywolana sztywnoscia
        Fd - dila wywolana tlumieniem
     '''
    log_column = int
    pos = bool
    vel = bool
    acc = bool

    force = bool
    dumping = bool
    inertia = bool
    force_control_sum = bool
    power = bool

    x_pos_log = list
    y_pos_log = list

    x_vel_log = list
    y_vel_log = list

    x_acc_log = list
    y_acc_log = list

    x_force_log = list
    y_force_log = list
    x_inertia_log = list
    y_inertia_log = list
    x_dumping_log = list
    y_dumping_log = list

    x_force_control_sum_log = list
    y_force_control_sum_log = list

    SDM = np.zeros((3, 2))
    dt = float

    def __init__(self, SDM, dt, log_column, position, velocity, acceleration, force, dumping, inertia, force_control_sum,power):
        self.SDM = SDM
        self.dt = dt
        self.log_column = log_column
        self.pos = position
        self.acc = acceleration
        self.vel = velocity

        self.force = force
        self.dumping = dumping
        self.inertia = inertia
        self.force_control_sum = force_control_sum
        self.power = power

        self.phase_angle = []

        self.x_pos_log = []
        self.y_pos_log = []

        self.x_vel_log = []
        self.y_vel_log = []

        self.x_acc_log = []
        self.y_acc_log = []

        self.x_force_log = []
        self.y_force_log = []

        self.x_inertia_log = []
        self.y_inertia_log = []

        self.x_dumping_log = []
        self.y_dumping_log = []

        self.x_force_control_sum_log =[]
        self.y_force_control_sum_log =[]

        self.moment_power = []

    def log_data(self, X, i, phase_angle, for_every_n):
        if np.mod(i, for_every_n) == 0:
            self.log_phase_angle(phase_angle)

            if self.pos:
                self.log_position(X)
            if self.vel:
                self.log_velocity(X)
            if self.acc:
                self.log_acceleration(X)
            if self.force:
                self.log_force(X)
            if self.dumping:
                self.log_dumping_force(X)
            if self.inertia:
                self.log_inertia(X)
            if self.force_control_sum:
                self.log_force_control_sum(X)
            if self.power:
                self.log_power(X)
    def clear_all_data(self):
        self.x_pos_log.clear()
        self.y_pos_log.clear()
        self.phase_angle.clear()
        #self.x_acc_log.clear()
        #self.y_acc_log.clear()


    def scale_pos_log(self, scale):
        self.x_pos_log = self.x_pos_log * scale
        self.y_pos_log = self.y_pos_log * scale

    def rotate_excitation_by_phase_angle_of_loging(self):
        for i, n in enumerate(self.phase_angle):
            self.rotate_excitation_data_point_by_phase_angle(i)

    def rotate_response_by_phase_angle_of_loging(self):
        for i, n in enumerate(self.phase_angle):
            self.rotate_response_data_point_by_phase_angle(i, n)

    def rotate_excitation_data_point_by_phase_angle(self, i):
        x = self.x_pos_log[i]
        y = self.y_pos_log[i]

        r = np.sqrt(x**2 + y**2)

        self.x_pos_log[i] = r
        self.y_pos_log[i] = 0

    def rotate_response_data_point_by_phase_angle(self, i, n):
        x = self.x_pos_log[i]
        y = self.y_pos_log[i]

        # licze promien odpowiedzi
        r = np.sqrt(x**2 + y**2)
        # licze kat fazowy stanu odpowiedzi
        response_phase = np.angle([complex(x,y)])

        # odejmuje kat wymuszenia od kata odpowiedzi
        new_phase = response_phase - n
        # licze nowe x i y
        # przypisuje wartosci w listach danych
        self.x_pos_log[i] = r * np.cos(new_phase)
        self.y_pos_log[i] = r * np.sin(new_phase)


    def log_phase_angle(self, angle):
        self.phase_angle.append(angle)

    def log_position(self, X):
        self.x_pos_log.append(X[0][self.log_column][0])
        self.y_pos_log.append(X[0][self.log_column][1])

    def log_velocity(self, X):
        self.x_vel_log.append(X[1][self.log_column][0])
        self.y_vel_log.append(X[1][self.log_column][1])

    def log_acceleration(self, X):
        self.x_acc_log.append(X[2][self.log_column][0])
        self.y_acc_log.append(X[2][self.log_column][1])

    def log_force(self, X):
        self.x_force_log.append(X[5][self.log_column][0])
        self.y_force_log.append(X[5][self.log_column][1])

    def log_dumping_force(self, X):
        self.x_dumping_log.append(X[6][self.log_column][0])
        self.y_dumping_log.append(X[6][self.log_column][1])

    def log_inertia(self, X):
        self.x_inertia_log.append(X[4][self.log_column][0])
        self.y_inertia_log.append(X[4][self.log_column][1])

    def log_force_control_sum(self, X):
        self.x_force_control_sum_log.append(X[3][self.log_column][0] )
        self.y_force_control_sum_log.append(X[3][self.log_column][1] )

    def log_power(self, X):
        pass
        #self.moment_power.append(X[3][self.log_column] * X[1][self.log_column])

    def get_pos_data(self):
        return [n for n in self.x_pos_log], [n for n in self.y_pos_log]
    def get_vel_data(self):
        return self.x_vel_log, self.y_vel_log
    def get_acc_data(self):
        return self.x_acc_log, self.y_acc_log
    def get_force_data(self):
        return self.x_force_log, self.y_force_log
    def get_damping_data(self):
        return self.x_dumping_log, self.y_dumping_log
    def get_inertia_data(self):
        return self.x_inertia_log, self.y_inertia_log
    def get_force_control_sum(self):
        return self.x_force_control_sum_log, self.y_force_control_sum_log
    def get_moment_power(self):
        return self.moment_power

