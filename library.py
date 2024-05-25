import numpy as np
import scipy.linalg as la


def spline_function(ui,vi,wi,xi,yi,x):
    """This function calculates the value in x of the spline function
     determined for the interval [xi,xi+1]."""
    Sx = yi+ui*(x-xi)+(1/2)*vi*(x-xi)**2+(1/6)*wi*(x-xi)**3
    return Sx

def find_interval(x,xi):
    """This function returns two indices, that of the element of
     start and end of the interval to which x belongs in
     the list ordered in ascending order, xi."""
    n = len(xi)
    for i in range(n-1):
        if xi[i] <= x and x <= xi[i+1]:
            j = i+1
            return i,j
    return 0,0

def cubic_spline(xi,yi,x):
    """This function implements the interpolation method by
     cubic splines for the point cloud (xi,yi)."""
    n = len(xi)
    hi = np.zeros((1,n-1))[0]
    alpha = np.zeros((1,n-2))[0]
    beta = np.zeros((1,n-2))[0]
    bi = np.zeros((1,n-2))[0]
    for i in range(n-1):
        hi[i] = xi[i+1]-xi[i]
        if i >= 1:
            alpha[i-1] = (hi[i-1]+hi[i])/3
            bi[i-1] = ((yi[i+1]-yi[i])/hi[i])-((yi[i]-yi[i-1])/hi[i-1])
            beta[i-1] = hi[i]/6
    betai = beta[0:n-3]
    
    A = np.diag(alpha)+np.diag(betai,-1)+np.diag(betai,1)
    vi = la.solve(A,bi)
    vi = np.append(0,vi)
    vi = np.append(vi,0)
    wi = np.zeros((1,n))[0]
    ui = np.zeros((1,n))[0]
    for i in range(n-1):
        wi[i] = (vi[i+1]-vi[i])/hi[i]
        ui[i] = ((yi[i+1]-yi[i])/hi[i])-(1/2)*vi[i]*hi[i]-(1/6)*wi[i]*hi[i]**2
    i,j = find_interval(x,xi)
    u,v,w = ui[i],vi[i],wi[i]
    Sx = spline_function(u,v,w,xi[i],yi[i],x)
    return Sx