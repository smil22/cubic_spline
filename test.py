import library
import numpy as np
import matplotlib.pyplot as plt

#Dots cloud construction
f = lambda x: 1/(1+x**2)
a,b,n = -5,5,50
xi = np.linspace(a,b,n)
yi = [f(t) for t in xi]

S = np.array([])
for x in xi:
    Sx = library.cubic_spline(xi,yi,x)
    S = np.append(S,Sx)
plt.scatter(xi,yi)
plt.plot(xi,S,'-r')
plt.show()


