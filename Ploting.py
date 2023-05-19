import matplotlib.pyplot as plt
from matplotlib import figure
import gc
import os


class Plot:

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