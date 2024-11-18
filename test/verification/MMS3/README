This case solves all equations in the plasma fluid model
without an spatial dependence to check time accuracy:

dne/dt + dGe/dx = k ne
dni/dt + dGi/dx = k ne
d2phidx2 = e(ne-ni)/eps0
dEe/dt + dGE/dx = S_E + Jheat - inel_term - el_term
Ge=mue ne E-De dne/dx = 0
Gi=mui ni E-Di dni/dx = 0
GE = 5/3 mue Ee E - 5/3 De dEe/dx = 0
S_E = -10.0*Ee
k=-5.0
n0=1e6

Initially

ne=ni=Ee=n0*(1+sin(pi*x/L))

Final solution after time t

ne=ni=n0*(1+sin(pi*x/L)) exp(-5.0*t)
Ee=n0*(1+sin(pi*x/L)) exp(-10.0*t)

