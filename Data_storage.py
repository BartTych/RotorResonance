import numpy as np


class LogStorageRotor:
    # operacje na danych pawinny byc gdzie indziej
    # nie w tej klasie
    # ona jest odpowiedzialna za przechowywanie danych i to wszystko
    """
     class for log of data for X matrix made every n-th step
        schematics of X matrix:

                  excitation                                 response
        row 0 [0][0][0] x_poz_exc  [0][0][1] y_poz_exc   [0][1][0] x_poz_r [0][1][1] y_poz_r
        row 1 [1][0][0] x_vel_exc  [1][0][1] y_vel_exc   [1][1][0] x_vel_r [1][1][1] y_vel_r
        row 2 [2][0][0] x_acc_exc  [2][0][1] y_vel_exc   [2][1][0] x_acc_r [2][1][1] y_acc_r
        row 3 [3][0][0] x_Fsum_exc [3][0][1] y_Fsum_exc  [3][1][0] x_Fsum_r[3][1][1] y_Fsum_r
        row 4 [4][0][0] x_Fa_exc   [4][0][1] y_Fa_exc    [4][1][0] x_Fa_r  [4][1][1] y_Fa_r
        row 5 [5][0][0] x_Fk_exc   [5][0][1] y_Fk_exc    [5][1][0] x_Fk_r  [5][1][1] y_Fk_r
        row 6 [6][0][0] x_Fd_exc   [6][0][1] y_Fd_exc    [6][1][0] x_Fd_r  [6][1][1] y_Fd_r
        x - component x of given value
        y - component y of given value
        Fa - force created because of acceleration
        Fk - force created because of stiffness
        Fd - force created because of damping
     """
    # niewiem czy to wogule jest potrzebne
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


    def __init__(self, SDM, dt, log_column, position, velocity, acceleration, force, dumping, inertia,
                 force_control_sum, power):
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

        self.x_force_control_sum_log = []
        self.y_force_control_sum_log = []

        self.moment_power = []

    # log ma mylaca nazwe bo to jest wewnetrzna funkcja
    # to co moge tez zmienic to jak decyduje co ma byc logowane
    # powinna to poprostu lista jako parametr ce mega zwiekszy czytelnocs

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

    def scale_pos_log(self, scale):
        self.x_pos_log = self.x_pos_log * scale
        self.y_pos_log = self.y_pos_log * scale

    # te jest cos co niepowinno byc w tej klasie
    # to sa modyfikacj na denych
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

        r = np.sqrt(x**2 + y**2)
        response_phase = np.angle([complex(x, y)])

        new_phase = response_phase - n
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
        self.x_force_control_sum_log.append(X[3][self.log_column][0])
        self.y_force_control_sum_log.append(X[3][self.log_column][1])

    def log_power(self, X):
        pass
        # self.moment_power.append(X[3][self.log_column] * X[1][self.log_column])

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

    def get_phase_angle(self):
        return self.phase_angle
