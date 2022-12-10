import numpy as np
import matplotlib.pyplot as plt

a_x = np.linspace(-4, 0, 3)
a_y = np.linspace(-2, 2, 3)

print(np.stack((a_x, a_y), axis=1))

