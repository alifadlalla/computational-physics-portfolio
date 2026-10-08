import numpy as np
from numba import njit
from convert import cartesian_to_orbital

def none(*args):
    """Placeholder function for optional parameters."""
    return 0

class Body:
    def __init__(self, name, gm, radius, position, velocity, color):
        self.name = name
        self.gm = gm
        self.radius = radius
        self.position = position
        self.velocity = velocity
        self.color = color

class Acceleration:
    """Handles acceleration calculations including perturbations and oblateness."""
    
    #@staticmethod
    #@njit
    def compute(X, gm, perturb=none, oblat=none):
        r = X[:3].T
        v = X[3:].T
        a = np.zeros_like(v)
        n = len(gm)

        for i in range(1, n):
            a[i] = -(gm[0] + gm[i]) * r[i] / np.linalg.norm(r[i])**3
            a[i] += perturb(X, gm, i)
            a[i] += oblat(X, gm, i)

        X_dot = np.concatenate((v.T, a.T))
        return X_dot

    #@staticmethod
    #@njit
    def perturb(X, gm, i):
        r = X[:3].T
        a = np.zeros(3)

        for j in range(1, len(gm)):
            if i != j:
                r_ij = r[j] - r[i]
                a += gm[j] * (r_ij / np.linalg.norm(r_ij)**3 - r[j] / np.linalg.norm(r[j])**3)

        return a

    #@staticmethod
    #@njit
    def oblat(X, gm, i):
        j2 = 1.4696e-2
        R = 60268 # equatorial radius
        r = X[:3].T
        xi, yi, zi = r[i]
        ri = np.linalg.norm(r[i])

        ax = -3 / 2 * (gm[0] + gm[i]) * j2 * R**2 * (xi / ri**7) * (xi**2 + yi**2 - 4 * zi**2)
        ay = -3 / 2 * (gm[0] + gm[i]) * j2 * R**2 * (yi / ri**7) * (xi**2 + yi**2 - 4 * zi**2)
        az = -3 / 2 * (gm[0] + gm[i]) * j2 * R**2 * (zi / ri**7) * (3 * xi**2 + 3 * yi**2 - 2 * zi**2)

        return np.array([ax, ay, az])


class NumericalMethod:
    """Implements various numerical integration methods."""
    
    #@staticmethod
    #@njit
    def runge_kutta(X, gm, dt, perturb=none, oblat=none):
        k1 = Acceleration.compute(X, gm, perturb, oblat)
        k2 = Acceleration.compute(X + 0.5 * dt * k1, gm, perturb, oblat)
        k3 = Acceleration.compute(X + 0.5 * dt * k2, gm, perturb, oblat)
        k4 = Acceleration.compute(X + dt * k3, gm, perturb, oblat)

        return X + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    #@staticmethod
    #@njit
    def leapfrog(X, gm, dt, perturb=none, oblat=none):
        Xdot = Acceleration.compute(X, gm, perturb, oblat)
        v, a = Xdot[:3], Xdot[3:]

        r = X[:3] + v * dt + 0.5 * a * dt**2

        Xdot = Acceleration.compute(np.concatenate((r, v)), gm, perturb, oblat)
        ap = Xdot[3:]
        v = v + 0.5 * (a + ap) * dt

        return np.concatenate((r, v))


class Simulation:
    """Handles the overall simulation of celestial systems."""
    
    def __init__(self, bodies, dt, nt, method, perturb=none, oblat=none):
        self.bodies = bodies
        self.dt = dt
        self.nt = nt
        self.method = method
        self.perturb = perturb
        self.oblat = oblat

        self.gm = np.array([body.gm for body in bodies])
        self.pos = np.zeros((nt, len(bodies), 3))
        self.vel = np.zeros((nt, len(bodies), 3))
        self.X = np.zeros((6, len(bodies)))

        self.run()

    def run(self):
        """Initializes and runs the simulation."""
        # Initialize positions and velocities
        for i, body in enumerate(self.bodies):
            self.X[:3, i] = body.position
            self.X[3:, i] = body.velocity

        # Run the simulation over the specified number of timesteps
        for t in range(self.nt):
            self.pos[t] = self.X[:3].T
            self.vel[t] = self.X[3:].T
            self.X = self.method(self.X, self.gm, self.dt, self.perturb, self.oblat)

    def get_orbital_elements(self):
        """Computes orbital elements for all bodies."""
        p = np.transpose(self.pos, (1, 0, 2))
        v = np.transpose(self.vel, (1, 0, 2))

        orb = np.array([cartesian_to_orbital(p[i], v[i], self.gm[0] + self.gm[i]).T for i in range(1, len(self.bodies))])
        return orb
