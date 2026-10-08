import numpy as np
from numba import njit # optional (on-the-fly compilation)

# Hide compilation warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
import warnings
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)

@njit # optional (on-the-fly compilation)
def orbital_to_cartesian(ell, GM):
    
    """
    Converts orbital elements to cartesian coordinates
    :param ell: orbital elements [a, e, i, ω, Ω,  λ]
    :param GM: gravitational parameter
    :return: cartesian coordinates and speeds ([x, y, z], [vx, vy, vz])
    """

    # Parsing orbital elements to original function parameters
    a = ell[0]
    e = ell[1]
    i = ell[2]
    omega = ell[3]
    OMEGA = ell[4]
    l = ell[5]

    k = e*np.cos(omega+OMEGA)
    h = e*np.sin(omega+OMEGA)
    q = np.sin(i/2)*np.cos(OMEGA)
    p = np.sin(i/2)*np.sin(OMEGA)

    # The original function converted to python (not my code, don't ask me how it works)
    n = np.sqrt(GM / a ** 3)
    fle = l-k*np.sin(l)+h*np.cos(l)
    corf = (l - fle + k * np.sin(fle) - h * np.cos(fle)) / (1 - k * np.cos(fle) - h * np.sin(fle))
    shift = 1e-8

    while abs(corf) > shift:
        fle += corf
        corf = (l - fle + k * np.sin(fle) - h * np.cos(fle)) / (1 - k * np.cos(fle) - h * np.sin(fle))
        shift *= 1.1

    lf = -k * np.sin(fle) + h * np.cos(fle)
    sam1 = -k * np.cos(fle) - h * np.sin(fle)
    asr = 1 / (1 + sam1)

    phi = np.sqrt(1 - k*k - h*h)
    psi = 1/(1+phi)

    x1 = a * (np.cos(fle) - k -psi * h * lf)
    y1 = a * (np.sin(fle) - h + psi * k * lf)

    vx1 = n * asr * a * (-np.sin(fle) - psi * h * sam1)
    vy1 = n * asr * a * (np.cos(fle) + psi * k * sam1)

    cis2 = 2 * np.sqrt(1 - p*p - q*q)

    tp = 1 - 2 * p*p
    tq = 1 - 2 * q*q
    dg = 2 * p * q
    
    pos = np.empty(3)
    vit = np.empty(3)

    pos[0]=x1*tp+y1*dg
    pos[1]=x1*dg+y1*tq
    pos[2]=(-x1*p+y1*q)*cis2
    vit[0]=vx1*tp+vy1*dg
    vit[1]=vx1*dg+vy1*tq
    vit[2]=(-vx1*p+vy1*q)*cis2

    return pos, vit

@njit
def cartesian_to_orbital(pos, vit, GM):
    
    """
    Converts cartesian coordinates to orbital elements
    :param pos: cartesian coordinates [x, y, z]
    :param vit: cartesian speeds [vx, vy, vz]
    :param GM: gravitational parameter
    :return: orbital elements [a, e, i, Ω, ω+Ω, λ]
    """

    out = np.zeros((len(pos),6))

    for j in range(len(pos)):
        rayon = np.sqrt(np.sum(pos[j]**2))
        v2 = np.sum(vit[j]**2)
        a = GM * rayon / (2 * GM - rayon * v2)
        gx = pos[j,1] * vit[j,2] - pos[j,2] * vit[j,1]
        gy = pos[j,2] * vit[j,0] - pos[j,0] * vit[j,2]
        gz = pos[j,0] * vit[j,1] - pos[j,1] * vit[j,0]
        gg = np.sqrt(gx*gx + gy*gy + gz*gz)

        cis2 = np.sqrt(0.5 * (1 + gz / gg))

        q = -gy / (2 * gg * cis2)
        p = gx / (2 * gg * cis2)

        tp = 1 - 2 * p*p
        tq = 1 - 2 * q*q
        dg = 2 * p * q

        x1 = tp * pos[j,0] + dg * pos[j,1] - 2 * p * cis2 * pos[j,2]
        y1 = dg * pos[j,0] + tq * pos[j,1] + 2 * q * cis2 * pos[j,2]
        vx1 = tp * vit[j,0] + dg * vit[j,1] - 2 * p * cis2 * vit[j,2]
        vy1 = dg * vit[j,0] + tq * vit[j,1] + 2 * q * cis2 * vit[j,2]

        k = gg * vy1 / GM - x1 / rayon
        h = -gg * vx1 / GM - y1 / rayon

        psi = 1/(1+np.sqrt(1-k*k-h*h))

        ach = 1 - psi * h*h
        ack = 1 - psi * k*k

        adg = psi * h * k
        det = ach * ack - adg*adg
        sm1 = x1/a + k
        sm2 = y1/a + h

        cf = (sm1*ack - sm2*adg) / det
        sf = (ach*sm2 - adg*sm1) / det

        fle = np.arctan2(sf, cf)

        l = fle - k * sf + h * cf

        # Parsing the original function outputs to the orbital elements
        e = np.sqrt(k*k + h*h)
        i = 2 * np.arcsin(np.sqrt(p**2 + q**2))
        omega = np.arctan2(p, q)
        varpi = np.arctan2(h, k)

        out[j] = [a, e, i, omega, varpi, l]

    return out

if __name__ == "__main__":

    GM = 1
    ell = np.array([100, 0.5, 1.5, 1.5, 1.5, 1.5])

    print('Orbital elements :')
    print(ell)

    pos, vit = orbital_to_cartesian(ell, GM)

    print('Conversion to cartesian coordinates :')
    print(pos,vit)
    
    pos2 = np.zeros((2,3))
    vit2 = np.zeros((2,3))
    pos2[0],pos2[1] = pos,pos
    vit2[0],vit2[1] = vit,vit

    ell = cartesian_to_orbital(pos2,vit2, GM)

    print('Conversion to orbital elements :')
    for i in ell[0]:
        print(i, end="\t")
    print()