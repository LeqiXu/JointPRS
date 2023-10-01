cimport cython
cimport numpy as np
from libc.math cimport sqrt, pow, cosh, sinh, exp, log
from libc.stdlib cimport rand
import numpy as np

cdef extern from "limits.h":
    np.int64_t RAND_MAX

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double psi(double x, double alpha, double lam):
    return -alpha*(cosh(x)-1.0)-lam*(exp(x)-x-1.0)

@cython.boundscheck(False)
@cython.wraparound(False) 
cdef double dpsi(double x, double alpha, double lam):
    return -alpha*sinh(x)-lam*(exp(x)-1.0)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double g(double x, double sd, double td, double f1, double f2):
    if (x >= -sd) and (x <= td):
        return 1.0
    elif x > td:
        return f1
    elif x < -sd:
        return f2

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double gigrnd(double p, double a, double b):
    cdef double lam = p
    cdef double omega = sqrt(a*b)
    cdef double alpha, x, t, s, eta, zeta, theta, xi, p_, r, td, sd, q, rnd, f1, f2
    cdef int swap = 0

    if lam < 0:
        lam = -lam
        swap = 1

    alpha = sqrt(pow(omega,2)+pow(lam,2))-lam

    # find t
    x = -psi(1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        t = 1.0
    elif x > 2.0:
        if (alpha == 0.0) and (lam == 0.0):
            t = 1.0
        else:
            t = sqrt(2.0 / (alpha + lam))
    elif x < 0.5:
        if (alpha == 0.0) and (lam == 0.0):
            t = 1.0
        else:
            t = log(4.0 / (alpha + 2.0 * lam))

    # find s
    x = -psi(-1.0, alpha, lam)
    if (x >= 0.5) and (x <= 2.0):
        s = 1.0
    elif x > 2.0:
        if (alpha == 0.0) and (lam == 0.0):
            s = 1.0
        else:
            s = sqrt(4.0 / (alpha * cosh(1) + lam))
    elif x < 0.5:
        if (alpha == 0.0) and (lam == 0.0):
            s = 1.0
        elif alpha == 0.0:
            s = 1.0 / lam
        elif lam == 0.0:
            s = log(1.0 + 1.0 / alpha + sqrt(1.0 / pow(alpha, 2) + 2.0 / alpha))
        else:
            s = min(1.0 / lam, log(1.0 + 1.0 / alpha + sqrt(1.0 / pow(alpha, 2) + 2.0 / alpha)))

    eta = -psi(t, alpha, lam)
    zeta = -dpsi(t, alpha, lam)
    theta = -psi(-s, alpha, lam)
    xi = dpsi(-s, alpha, lam)

    p_ = 1.0/xi
    r = 1.0/zeta

    td = t-r*eta
    sd = s-p_*theta
    q = td+sd

    while True:
        U = rand() / float(RAND_MAX)
        V = rand() / float(RAND_MAX)
        W = rand() / float(RAND_MAX)
        if U < q/(p_+q+r):
            rnd = -sd+q*V
        elif U < (q+r)/(p_+q+r):
            rnd = td-r*log(V)
        else:
            rnd = -sd+p_*log(V)

        f1 = exp(-eta-zeta*(rnd-t))
        f2 = exp(-theta+xi*(rnd+s))
        if W*g(rnd, sd, td, f1, f2) <= exp(psi(rnd, alpha, lam)):
            break

    rnd = exp(rnd)*(lam/omega+sqrt(1.0+pow(lam,2)/pow(omega,2)))
    if swap:
        rnd = 1.0/rnd

    rnd = rnd/sqrt(a/b)
    return rnd

@cython.boundscheck(False)
@cython.wraparound(False)
cdef double calculate_xx(np.ndarray[double, ndim=2] m_inv_cov_gg, np.ndarray[double, ndim=1] beta_all_jj_group):
    cdef np.int64_t i,j,k
    cdef np.int64_t size = beta_all_jj_group.shape[0]
    cdef np.ndarray[double, ndim=1] temp = np.zeros(size)
    cdef double result = 0.0
    
    # First calculate: m_inv_cov_gg @ beta_all_jj_group
    for i in range(size):
        for k in range(size):
            temp[i] += m_inv_cov_gg[i, k] * beta_all_jj_group[k]

    # Then calculate: beta_all_jj_group.T @ (m_inv_cov_gg @ beta_all_jj_group)
    for j in range(size):
        result += beta_all_jj_group[j] * temp[j]
    
    return result            

@cython.boundscheck(False)
@cython.wraparound(False)
def psi_update(np.int64_t p_tot, np.ndarray[np.int64_t, ndim=1] map_nonzero_grp, 
               dict uni_nonzero_grp, dict m_inv_cov, np.ndarray[double, ndim=2] beta_all, 
               np.ndarray[double, ndim=1] xx, np.ndarray[double, ndim=1] psi, double a, 
               np.ndarray[np.int64_t, ndim=1] n_grp, np.ndarray[double, ndim=1] delta):
    cdef np.ndarray[np.int64_t, ndim=1] group
    cdef np.int64_t jj, gg
    cdef double temp_psi
    for jj in range(p_tot):
        gg = map_nonzero_grp[jj]
        group = uni_nonzero_grp[gg]
        if len(group) == 1:
            xx[jj] = (beta_all[jj, group] ** 2) * m_inv_cov[gg]
        else:
            xx[jj] = calculate_xx(m_inv_cov[gg], beta_all[jj, group])
        while xx[jj] > 0:
            try:
                temp_psi = gigrnd(a-0.5*n_grp[jj], 2.0*delta[jj], xx[jj])
            except:
                continue
            else:
                break
        psi[jj] = temp_psi
    
    return psi
