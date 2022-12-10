import copy
import os
from decimal import Decimal

import numpy as np

import matplotlib.pyplot as plt
import Data_storage
import RotorCalculation
from tqdm import tqdm

def flat_top(n, L):
    a0 = 0.21557895
    a1 = 0.41663158
    a2 = 0.277263158
    a3 = 0.083578947
    a4 = 0.006947368
    w = a0 - a1 * np.cos((2*np.pi*n)/(L-1)) + a2 * np.cos((4*np.pi*n)/(L-1))\
        -a3 * np.cos((6*np.pi*n)/(L-1)) + a4 * np.cos((8*np.pi*n)/(L-1))
    return w

#x = 5
#y = 5

#import numpy as np
#import matplotlib.pyplot as plt
#%matplotlib inline
def calculate_phase_angle(T, Hz):
    return 2 * np.pi * T * Hz

def cal_pozytion_of_excitation(T, Hz1, amplitude, phi):
    return amplitude * np.sin(2 * np.pi * ( Hz1 * T + phi ))#  + amplitude * np.sin(2*np.pi* T * Hz2)\
          # + amplitude * np.sin(2*np.pi* T * 170)

def cal_velocity_of_excitation(T, Hz1 , amplitude, phi):
    return amplitude * 2 * np.pi * Hz1 * np.cos(2 * np.pi * (T * Hz1 + phi / 360))#

def pole_wektorowe(r_wymuszenia, Hz, T, range):

    # potrzebne jest polozenie wymuszenia - przesuniecie sily sprezystosci
    # i predkosc wymuszenia - przesuniecie sily tlumienia
    exc_x = cal_pozytion_of_excitation(T, Hz, r_wymuszenia, 0)
    exc_v = cal_velocity_of_excitation(T, Hz, r_wymuszenia, 0)

    x, y = np.meshgrid(np.linspace(-range, range, 30),np.linspace(-range, range, 30))

    # tu licze pozycje y bez zeskalowania, potrzebny mi jest do tego wspolczynnik skalowania
    # z poprzedniego kroku
    # podczas liczenia wykozystuje te nowe y1
    # ale na koncu wrzucam wyliczone wektory dla pozycji zeskalowanego x,y


    m = 1.0
    k = 10000000.0
    ni = 4000.0

    # vektor na x
    u = y
    # vektor na y
    v = (-(x-exc_x) * k + -(y-exc_v) * ni) / m

    return x,y,u,v
    #plt.quiver(x,y,u,v)

# co moge zrobic to zrobic pole wektorowe dla danego momentu w czasie czyli tak naprawde w zaleznosci od tego gdzie jest polozenie wymuszenia
# wtedy moge zrobic wzory tego jaka jest pochodna x i wyplotowac to jako wektor na kirunku x
# jaka jest pochodna na y i jaki jest wektor na y
# jest cool bo x i v w pelni opisuja to co sie dzieje z masa
# pole wektorowe bedzie sie zmieniac w czasie bo polozenie wymuszenia tez bedzie sie zmieniac
# pytanie jak najlepiej to wyswietlic. ej to pole bedzie sie zmieniac z ukladzie wspolrzednych ktory nieobraca sie razem z wymuszeniem
# a ajk to bedzie wygladac jak bede obracac ukladem razem z wymuszeniem ?
# jak zrobie pole wektorowe to wyzwaniem bedzie obrocenie go do nowego wirujacego ukladu wspolrzednych
# eee nic trudnego, skladam oba wekory razem i potem tylko obracam je odpowiednio i rozkladam na nowe wektory jednostkowe
# wiec krok jeden to napisanie rownan pola wektorowego z parametrem aktualniego polozenia wymuszenia

# jakie sily beda wchodzic do pola wektorowego
# sztywnosc sprezyny i plumienie od predkosci
# reakcja czyli sila przyspieszenia niebedzie tam wchodzic
# czyli to bedzie pole sil zewnetrznych dzialajacych na mase

# na kierunku x
# zmiane w przemieszczeniu wywoluje predkosc , czyli wektor na x jest poprostu wartoscia predkosci
# na kierunku y
# sila od sztywnosci i tlumienia wywoluje przyspieszenie czyli zmiane predkosc


def cal_vector_field_in_rotated_coordinates(x, y, rot, T, Hz, r_wymuszenia, SDM, nor_coeff):

    u = copy.deepcopy(x)
    v = copy.deepcopy(y)

    #tymczasowo daje T = 0
    exc_x = cal_pozytion_of_excitation(T, Hz, r_wymuszenia, 0)
    exc_v = cal_velocity_of_excitation(T, Hz, r_wymuszenia, 0)

    #m = 1.0
    #k = 10.0
    #ni = 5.0

    k = SDM[0][0]

    # tlumienie
    ni = SDM[1][0]

    # masa
    m = SDM[2][0]


    for index,n in np.ndenumerate(x):

        # tu jest blad bo mam 2 rozne rzeczy wyliczone
        # dla xx jest policzony vektor predkosci
        # a dla yy jest policzony wektor przyspieszenia

        # wiec
        # albo przyspieszenie wszedzie
        # albo predkosc wszedzie
        # opcja predkosc wszedzie

        xx = y[index] / nor_coeff
        #v = (-(x-exc_x) * k + -(y-exc_v) * ni) / m
        yy = (-(x[index] - exc_x) * k + -(y[index] - exc_v) * ni) / m
        # yy = yy * nor_coeff
        # u = y
        # vektor na y
        # v = ( -(x - exc_x) * k + -(y - exc_v) * ni) / m

        r = np.hypot( xx , yy)
        ang = np.angle([complex(xx, yy)])

        u[index] = r * np.cos(ang + rot)[0]
        v[index] = r * np.sin(ang + rot)[0]

    return u, v

def obrucenie_ukladu_wspolrzednych(x, y, rot):
    """ tu prawdopodobnie jest blad i ta metoda obraca prawidlowo tylko macierze ktore sa kwadratowe.
     musza miec taki sam kszatalt """

    x1 = copy.deepcopy(x)
    y1 = copy.deepcopy(y)

    for index, n in np.ndenumerate(x):
        r = np.hypot(x[index], y[index])
        ang = np.angle([complex(x[index], y[index])])
        x1[index] = r * np.cos(ang + rot)[0]
        y1[index] = r * np.sin(ang + rot)[0]

    return x1, y1

def Calculate_rotated_vector_field(Hz, r_wym, T, SDM, nor_coeff, rot_by_change_of_csys, range):

    x, y, u, v = pole_wektorowe(r_wym, Hz, T, range)

    x1 = x

    #print(f"macierz y: {y}" )

    y1 = y# / nor_coeff
    #y =  y1
    #print(f"macierz podzielenie: {np.divide(y,y1)}")
    #fig, ax = plt.subplots(1)
    #ax.scatter(x, y, s = 0.1)

    #wyglada ze to dziala ok
    rot_ange = calculate_phase_angle(T, Hz)
    #rot_ange = 0

    x2, y2 = obrucenie_ukladu_wspolrzednych(x1, y1, -rot_ange)

    #x2 = x1
    #y2 = y1
    # test skalowania po obrocie
    y2 = y2 / nor_coeff


    #ax.scatter(x2, y2, s = 0.1)
    #ax.set_aspect('equal')
    #plt.show()

    u1, v1 = cal_vector_field_in_rotated_coordinates(x2, y2, rot_by_change_of_csys, T, Hz, r_wym, SDM, nor_coeff)

    # u i v sa tymczasowo
    # doceleowo u1 i v1
    return x, y, u1, v1

'''
for i in np.linspace(0, 6, 100):

    Hz = 40
    r_wym = 0.0005

    x, y, u, v = pole_wektorowe(r_wym, Hz, 0)
    x1, y1 = obrucenie_ukladu_wspolrzednych(x, y, i)
    u1, v1 = cal_vector_field_in_rotated_coordinates(x1, y1, i, 0, Hz, r_wym)

    #plt.show()
    fig, ax = plt.subplots(1, figsize = (25, 25))
    # ax.scatter(res_x, res_y, s = 0.1)
    # ax.scatter(exc_x,exc_y, s = 0.1)

    ax.quiver(x, y, u, v, label = "sily wymuszenia", color = 'red')
    ax.quiver(x, y, u1, v1, label = "sily wymuszenia w obruconym ukladzie")


    #ax.scatter(x,y, s = 4, label="punkty bazowe pola wektorowego")
    #ax.scatter(x1,y1, s = 4, label="punkty przesuniete pola wektorowego")
    ax.set_aspect('equal')
    ax.set_xlabel('pozycja x [m]', fontsize = 20)
    ax.set_ylabel('predkosc y [m/s]', fontsize = 20)
    ax.tick_params(axis = 'x', labelsize = 20)
    ax.tick_params(axis = 'y', labelsize = 20)
    fig.legend(loc = 'upper right', fontsize = 20)

    exc_x = cal_pozytion_of_excitation(i, Hz, r_wym, 0)
    exc_v = cal_velocity_of_excitation(i, Hz, r_wym, 0)

    file = r'/Users/bart/python/Rezonans/Wyniki/pole_vectorowe/'
    plt.savefig(f"{file}/x_{round(i, 5):.5f}_v.png")
    plt.close()
'''
'''
for i in range(100):

    x, y, u, v = Calculate_rotated_vector_field(30, 0.005, i / 10)

    fig, ax = plt.subplots(1, figsize=(25, 25))
    ax.quiver(x, y, u, v, label="sily wymuszenia", color='red')
    ax.set_xlabel('pozycja x [m]', fontsize=20)
    ax.set_ylabel('predkosc y [m/s]', fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    fig.legend(loc='upper right', fontsize=20)

    file = r'/Users/bart/python/Rezonans/Wyniki/pole_vectorowe/'
    plt.savefig(f"{file}/x_{round(i, 5):.5f}_v.png")
    plt.close()
'''

