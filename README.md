# RotorResonance

Script for plotting position of Jeffcott rotor geometric center
in transient with excitation at given frequency.

The Jeffcott rotor is a single disk symmetrically mounted
on a uniform elastic shaft, used in rotodynamics as basic model for
analysis of rotating systems.

I made this script to better understand what is happening during
resonance of mechanical system.

Basic assumption is that rotor is excited by unbalance of rotor.
Script is plotting position of rotor center of geometry, in rotating
csys (coordinate system) or non rotating system at frequency 
defined by user.

user is defining properties of analysed system nad type of plot
in main_rotor.py file, in setup section.

Plot in rotating csys:
Plot is rotated in such a way that excitation (unbalance ) is always
positioned on x-axis positive side. So you can see what happens to 
response (rotor center of gravity) from time zero (response at point
x=0,y=0) to defined end time of calculation. Excitation and response
are rotating ccw (counter clock wise). 

Plot in non rotating csys:
Response is plotted form point of observer looking at system 
form non rotating perspective. Response also start at point x=0,y=0
but excitation is no longer at fix point as for rotating csys. 


 