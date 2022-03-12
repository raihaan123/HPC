import numpy as np
from matplotlib import pyplot as plt


# Solve the Barkley PDE using numerical methods including central finite difference

class ReactionDiffusion:
    '''
    Class to solve the Barkley PDE
    '''

    def _set_params(self, dt=1e-3, T=1,
                          Nx=11, Ny=11, 
                          a=0.75, b=0.06, 
                          mu1=5.0, mu2=0.0, 
                          epsilon=50):

        self.a = a
        self.b = b
        self.mu1 = mu1
        self.mu2 = mu2
        self.epsilon = epsilon
        self.Nx = Nx
        self.Ny = Ny
        self.dt = dt
        self.T = T
        

        self.dx = self.dy = 1.0


        # Create the mesh - x is from 0 to Nx in intervals of dx
        self.X = np.arange(0, self.Nx, self.dx)
        self.Y = np.arange(0, self.Ny, self.dy)

        print(f"X = {self.X}, Y = {self.Y}")


    def _set_initial_conditions(self):

        '''
        Initial conditions - u = 1 for y > Ny*dy/2 and u = 0 for y <= Ny*dy/2
                             v = 1 for x < Nx*dx/2 and v = 0 for x >= Nx*dx/2
        '''

        # Print the boundary conditions self.Ny*self.dy/2 and self.Nx*self.dx/2
        print(f"Ny*dy/2 = {self.Ny*self.dy/2}")
        print(f"Nx*dx/2 = {self.Nx*self.dx/2}")


        # Using np.where
        self.u = np.where(self.Y > self.Ny*self.dy/2, 1, 0)
        self.u = np.tile(self.u, (self.Nx, 1))
        self.u = self.u.T

        # Again for v
        self.v = np.where(self.X < self.Nx*self.dx/2, self.a/2, 0)
        self.v = np.tile(self.v, (self.Ny, 1))
        
        # Checking the inital conditions
        print(f"u = {self.u}, v = {self.v}")


    def _f1(self, u, v):
        return self.epsilon * u * (1-u) * (u - (v+self.b)/self.a)


    def _f2(self, u, v):
        return u**3 - v

    
    def _enforce_neumann_bc(self):
        '''
        Enforce Neumann boundary conditions on u and v
        '''

        u = self.u
        v = self.v
        Nx = self.Nx
        Ny = self.Ny

        # Enforce Neumann boundary conditions on u
        u[0, :] = u[1, :]
        u[Nx-1, :] = u[Nx-2, :]
        u[:, 0] = u[:, 1]
        u[:, Ny-1] = u[:, Ny-2]

        # Enforce Neumann boundary conditions on v
        v[0, :] = v[1, :]
        v[Nx-1, :] = v[Nx-2, :]
        v[:, 0] = v[:, 1]
        v[:, Ny-1] = v[:, Ny-2]

        # Check the new u and v
        print(f"u = {self.u}, v = {self.v}")

    
    def _time_integrate(self):

        u = self.u
        v = self.v
        Nx = self.Nx
        Ny = self.Ny
        mu1 = self.mu1
        mu2 = self.mu2
        dt = self.dt

        # Add the boundary conditions for u and v matrices
        
        for j in range(1, Ny-1):
            for i in range(1, Nx-1):
                u[i, j] += dt*((mu1/1) * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 * u[i, j]) + self._f1(u[i, j], v[i, j]))
                v[i, j] += dt*((mu2/1) * (v[i+1, j] + v[i-1, j] + v[i, j+1] + v[i, j-1] - 4 * v[i, j]) + self._f2(u[i, j], v[i, j]))


    def plot(self):
        fig, ax = plt.subplots()
        


    def solve(self):

        u = self.u
        v = self.v

        # Create time distribution from T and dt
        t = np.arange(0, self.T, self.dt)

        for _ in t:
            self.time_integrate()
        

        print(f"u = {self.u}, v = {self.v}")
        # self.plot()


x = ReactionDiffusion()
x._set_params()
x._set_initial_conditions()
x.solve()