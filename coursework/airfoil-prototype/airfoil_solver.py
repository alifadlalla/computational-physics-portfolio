import numpy as np
from airfoils import Airfoil
import matplotlib.pyplot as plt
from matplotlib import cm



class AirfoilSolver:
    def __init__(self, nx, ny, length, width, time_step, num_steps, nit, u_freestream, v_freestream, density, viscosity, airfoil_code):
        self.nx = nx
        self.ny = ny
        self.length = length
        self.width = width
        self.dt = time_step
        self.num_steps = num_steps
        self.nit = nit
        self.u_freestream = u_freestream
        self.v_freestream = v_freestream
        self.rho = density  # density
        self.nu = viscosity  # viscosity
        self.nit = nit
        


        # Ensure airfoil_code is a string
        if not isinstance(airfoil_code, str):
            raise ValueError("airfoil_code must be a string.")
        self.airfoil_code = airfoil_code

        # pressure
        self.P0 = 0.5*self.rho* self.u_freestream**2

        # Initialize mesh        
        self.X, self.Y = np.meshgrid(np.linspace(0, self.length, self.nx), np.linspace(0, self.width, self.nx))
        # Calculate grid spacing
        self.dx = self.length / self.nx
        self.dy = self.width / self.ny
        
        self.coeff_upper, self.coeff_lower = self.fit_curves()


        self.airfoil_x, self.airfoil_y = self.load_airfoil_coordinates()
        self.initialize_flow_field()
        self.inside_indices = self.indices_airfoil()
        
       
    def initialize_flow_field(self):
        # Initialize flow field variables
        self.u = np.zeros((self.nx, self.ny))
        self.v = np.zeros((self.nx, self.ny))
        self.p = np.zeros((self.nx, self.ny))
        self.b = np.zeros((self.nx, self.ny))

  

    def load_airfoil_coordinates(self, combined_output=True):
        """
        Load airfoil coordinates using the NACA4 series and scale them to the domain.

        Returns:
        - airfoil_x: X-coordinates of airfoil points
        - airfoil_y: Y-coordinates of airfoil points
        """

        def scale_x(x):
            return (x-0.5) * (0.5*self.length) +(0.5*self.length)
    
        def scale_y(y):
            return (y-0.5) * (0.5*self.width) +(0.75*self.width)

        
        foil = Airfoil.NACA4(self.airfoil_code, n_points=int(self.nx/2))
        x_upper, y_upper, x_lower, y_lower = foil._x_upper, foil._y_upper, foil._x_lower, foil._y_lower
        x_upper, y_upper, x_lower, y_lower = scale_x(x_upper), scale_y(y_upper), scale_x(x_lower), scale_y(y_lower)
        all_points = foil.all_points

        airfoil_x = all_points[0, :]
        airfoil_y = all_points[1, :]

        # Scale airfoil coordinates to match simulation domain
        airfoil_x_scaled = scale_x(airfoil_x)
        airfoil_y_scaled = scale_y(airfoil_y)
        

        if combined_output is True:
            return airfoil_x_scaled, airfoil_y_scaled
        elif combined_output is False:
            return x_upper, y_upper, x_lower, y_lower # useful for the plot
    
    def indices_airfoil(self):

        # Airfoil surface indices

        inside_indices = [] 
        for i in range(1,self.ny-1):
            for j in range(1,self.nx-1):
                x = self.dx*j
                y = self.dy*i
                if min(self.airfoil_x)<=x<=max(self.airfoil_x):
                    if min(self.airfoil_y)<=y<=max(self.airfoil_y):    
                        if self.is_inside_airfoil(self.dx*j, self.dy*i):
                            inside_indices.append((i, j))
        return inside_indices

    def fit_curves(self):
        # Fit polynomials to the data
        x_upper, y_upper, x_lower, y_lower = self.load_airfoil_coordinates(combined_output=False)
        degree = 3  
        self.coeff_upper = np.polyfit(x_upper, y_upper, degree)
        self.coeff_lower = np.polyfit(x_lower, y_lower, degree)

        return self.coeff_upper, self.coeff_lower

    # Function to check if a point (x, y) lies between the upper and lower curves
    def is_inside_airfoil(self, x, y):
        upper_bound = np.polyval(self.coeff_upper, x)
        lower_bound = np.polyval(self.coeff_lower, x)
        return lower_bound <= y <= upper_bound

    
    
    def apply_boundary_velocity(self):

         # Left Boundary
        self.u[:, 0] = self.u_freestream
        self.v[:, 0] = self.v_freestream 

         # Right Boundary
        self.u[:, -1] = self.u[:,-2]  
        self.v[:, -1] = self.v[:,-2]
        
        # Top boundary
        self.u[0, :] = self.u[1,:]
        self.v[0, :] = self.v[1,:] 

        # Bottom Boundary
        self.u[-1, :] = self.u[-2,:]
        self.v[-1, :] = self.v[-2,:] 
        
        # Airfoil surface boundary
        for i, j in self.inside_indices:
            # Set velocity to 0 inside the airfoil
            self.u[i, j] = 0.0
            self.v[i, j] = 0.0
    def apply_boundary_pressure(self):
        #left 
        self.p[:,0] = self.P0
        #self.p[:,0] = self.p[:,1]
        #right
        #self.p[:, -1] = self.p[:,-2]
        self.p[:, -1] = self.P0
        # up
        #self.p[0, :] = self.p[1,:] 
        self.p[0, :] = self.P0    
        # bottom
        #self.p[-1, :] = self.p[-2,:]
        self.p[-1, :] = self.P0

        for i, j in self.inside_indices:
            self.p[i,j] = self.P0

    def airfoil_influence(self, i, j, strength=1.0):
        # Calculate distance between point and airfoil
        x = self.dx*j
        y = self.dy*i
        distance = np.sqrt((x - self.airfoil_x[i])**2 + (y - self.airfoil_y[j])**2)

        # Influence of airfoil on velocity
        influence = strength / (1 + 0.2 * np.exp(-10 * distance))

        return influence

    def apply_airfoil_influence(self, strength=1.0):
        for i in range(self.ny):
            for j in range(self.nx):
                if (i, j) not in self.inside_indices:
                    # Calculate influence of airfoil on velocity at each grid point
                    influence = self.airfoil_influence(i, j, strength)
                    self.u[i, j] += influence*self.u[i, j]
                    self.v[i, j] += influence*self.v[i, j]
    

    def build_up_b(self):
        # velocity term used in the pressure computation
        self.b[1:-1, 1:-1] = (self.rho * (1 / self.dt * 
                        ((self.u[1:-1, 2:] - self.u[1:-1, 0:-2]) / 
                        (2 * self.dx) + (self.v[2:, 1:-1] - self.v[0:-2, 1:-1]) / (2 * self.dy)) -
                        ((self.u[1:-1, 2:] - self.u[1:-1, 0:-2]) / (2 * self.dx))**2 -
                        2 * ((self.u[2:, 1:-1] - self.u[0:-2, 1:-1]) / (2 * self.dy) *
                            (self.v[1:-1, 2:] - self.v[1:-1, 0:-2]) / (2 * self.dx))-
                            ((self.v[2:, 1:-1] - self.v[0:-2, 1:-1]) / (2 * self.dy))**2))


    
    def calculate_pressure(self):
        # compute pressure using pressure Poisson equation
        pn = np.empty_like(self.p)
        pn = self.p.copy()
        
        self.build_up_b()
        for _ in range(self.nit):
            pn = self.p.copy()
            self.p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * self.dy**2 + 
                            (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * self.dx**2) /
                            (2 * (self.dx**2 + self.dy**2)) -
                            self.dx**2 * self.dy**2 / (2 * (self.dx**2 + self.dy**2)) * 
                            self.b[1:-1,1:-1])
            self.apply_boundary_pressure()          
        return self.p

    def calculate_u_velocity(self, un, vn):
        # Update x-velocity using explicit finite difference method
        self.u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * self.dt / self.dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * self.dt / self.dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         self.dt / (2 * self.rho * self.dx) * (self.p[1:-1, 2:] - self.p[1:-1, 0:-2]) +
                         self.nu * (self.dt / self.dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         self.dt / self.dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

    def calculate_v_velocity(self, un, vn):
        # Update y-velocity using explicit finite difference method
        self.v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * self.dt / self.dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * self.dt / self.dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        self.dt / (2 * self.rho * self.dy) * (self.p[2:, 1:-1] - self.p[0:-2, 1:-1]) +
                        self.nu * (self.dt / self.dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        self.dt / self.dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
    

    def airfoil_flow(self):
        # iterate over the number od steps and compute u,v,p
        un = np.empty_like(self.u)
        vn = np.empty_like(self.v)

        for t in range(self.num_steps):
            un = self.u.copy()
            vn = self.v.copy()

            self.calculate_pressure()

            self.calculate_u_velocity(un, vn)
            self.calculate_v_velocity(un, vn)

            self.apply_boundary_velocity()

            #self.apply_airfoil_influence() #This is not working for the moment it needs some adjustment, since it gives incorrect pressure
            
        return self.u, self.v, self.p
    
    def plot_airfoil(self):
        # Plot airfoil
        foil = Airfoil.NACA4(self.airfoil_code, n_points=int(self.nx/2))
        foil.plot(show=True, save=False, settings={'camber': True})

        
    def visualize(self, show=True, save=False):
        self.airfoil_flow()
        
        # Plot airfoil
        x_upper_scaled, y_upper_scaled, x_lower_scaled, y_lower_scaled = self.load_airfoil_coordinates(combined_output=False)
        plt.plot(x_upper_scaled, y_upper_scaled, '-', color='blue')
        plt.plot(x_lower_scaled, y_lower_scaled, '-', color='green')

        plt.contourf(self.X, self.Y, self.p, alpha=0.5, cmap=cm.viridis)
        cbar = plt.colorbar()
        cbar.set_label('Pressure')
        plt.streamplot(self.X, self.Y, self.u, self.v)
        plt.xlabel('X')
        plt.ylabel('Y')
        
        if save:
            plt.savefig('airfoil_simulation.pdf', dpi=300)

        if show:
            plt.show()


