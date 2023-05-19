
from Data_storage import DataStorage
import numpy as np

class Ratate_csys:

    @staticmethod
    def rotate_excitation_by_phase_angle_of_loging(dataToRotate: DataStorage):
        X, Y = dataToRotate.get_pos_data()

        for i, angle in enumerate(X):
            Ratate_csys.__rotate_excitation_data_point_by_phase_angle(X, Y, i)
        return X, Y

    @staticmethod
    def rotate_response_by_phase_angle_of_loging(dataToRotate: DataStorage):
        phaseData = dataToRotate.get_phase_angle()
        X,Y = dataToRotate.get_pos_data()

        for i, angle in enumerate(phaseData):
            Ratate_csys.__rotate_response_data_point_by_phase_angle(X,Y, i, angle)
        return X, Y

    @staticmethod
    def __rotate_excitation_data_point_by_phase_angle(X, Y, index):
        x = X[index]
        y = Y[index]

        r = np.sqrt(x ** 2 + y ** 2)

        X[index] = r
        Y[index] = 0

    @staticmethod
    def __rotate_response_data_point_by_phase_angle(X, Y, index, angle):
        x = X[index]
        y = Y[index]

        r = np.sqrt(x ** 2 + y ** 2)
        response_phase = np.angle([complex(x, y)])

        new_phase = response_phase - angle
        X[index] = r * np.cos(new_phase)
        Y[index] = r * np.sin(new_phase)
