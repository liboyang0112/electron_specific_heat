#!/home/boyang/miniconda3/bin/python

from sympy import *

x = symbols('x')
vx = symbols('vx')
y = symbols('y')
z = symbols('z')
#kT = symbols('kT', positive=True)
skTby2pim = symbols('A', positive=True)
n = symbols('n', positive=True)
m = symbols('m', positive=True)
kT = skTby2pim*skTby2pim*2*pi*m
a = n*sqrt(m/2/pi/kT)
b = m/2/kT
eptdf = m*a*m/2/pi/kT*x*(x**2+y**2+z**2)*exp(-b*((x-vx)**2+y**2+z**2))/2
gammadf = a*x*exp(-b*((x-vx)**2))
ptdf = m*a*x*x*exp(-b*((x-vx)**2))
ept = expand(integrate(integrate(integrate(eptdf,(y,-oo,oo)),(z,-oo,oo)),(x,0,oo)))
gamma = expand(integrate(gammadf,(x,0,oo)))
pt = expand(integrate(ptdf,(x,0,oo)))
print(gamma)
print()
print(ept)
print()
print(pt)
