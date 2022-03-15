import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import ArtistAnimation


# Solve the Barkley PDE using numerical methods including central finite difference

class ReactionDiffusion:
    '''
    Class to solve the Barkley PDE
    '''

    def _set_params(self, dt=1e-3, T=10,
                          Nx=101, Ny=101, 
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
        Initial conditions
        ------------------
        u = 1 for y > Ny*dy/2 and u = 0 for y <= Ny*dy/2
        v = 1 for x < Nx*dx/2 and v = 0 for x >= Nx*dx/2
        '''

        # Print the initial conditions self.Ny*self.dy/2 and self.Nx*self.dx/2
        print("Initial conditions:")
        print(f"Ny*dy/2 = {self.Ny*self.dy/2}")
        print(f"Nx*dx/2 = {self.Nx*self.dx/2}")


        # Using np.where() to find the nodes that satisfy the conditions
        self.u = np.where(self.Y > self.Ny*self.dy/2, 1, 0)
        self.u = np.tile(self.u, (self.Nx, 1))
        self.u = self.u.T

        # Again for v
        self.v = np.where(self.X < self.Nx*self.dx/2, self.a/2, 0)
        self.v = np.tile(self.v, (self.Ny, 1))

        # Add to u and v histories
        self.u_hist = []
        self.v_hist = []
        self.u_hist.append(self.u)
        self.v_hist.append(self.v)
        
        # Checking the inital conditions
        print(f"u = {self.u}, v = {self.v}")

        # Animated figure store
        self.fig, self.ax = plt.subplots()
        ax = self.ax
        ax.imshow(self.u, cmap='jet', interpolation='nearest', origin='lower', extent=[0, self.Nx, 0, self.Ny])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('u')

        self.figs = []


    def _f1(self, u, v):
        '''Reaction term 1'''
        return self.epsilon * u * (1-u) * (u - (v+self.b)/self.a)


    def _f2(self, u, v):
        '''Reaction term 2'''
        return u**3 - v



    def laplacian_stencil_interior(self):
        '''
        Laplacian contribution for interior nodes
        --> Uses central finite differencing
        '''

        # Import the necessary variables
        u = self.u
        v = self.v
        dx = self.dx        # Knowing dx=dy always
        mu1 = self.mu1
        mu2 = self.mu2

        h = dx


        # Create the spatial derivative matrices
        u_xx = np.zeros(self.u.shape)
        u_yy = np.zeros(self.u.shape)
        
        v_xx = np.zeros(self.v.shape)
        v_yy = np.zeros(self.v.shape)

        
        # Find the spatial contributions to u and v
        for i in range(1, self.Nx-1):
            u_yy[1:-2, i] = mu1/h**2 * (u[2:-1, i] + u[:-3, i] -2*u[1:-2, i])       # Vectorized operation
            v_yy[1:-2, i] = mu2/h**2 * (v[2:-1, i] + v[:-3, i] -2*v[1:-2, i])       # Vectorized operation

        for j in range(1, self.Ny-1):
            u_xx[j, 1:-2] = mu1/h**2 * (u[j, 2:-1] + u[j, :-3] -2*u[j, 1:-2])       # Vectorized operation
            v_xx[j, 1:-2] = mu2/h**2 * (v[j, 2:-1] + v[j, :-3] -2*v[j, 1:-2])       # Vectorized operation

        # Sum the derivatives
        self.u_laplacian = u_xx + u_yy
        self.v_laplacian = v_xx + v_yy

    

    def _time_integrate(self):

        u = self.u
        v = self.v
        dt = self.dt

        # Calculate the temporal derivative matrices
        u_t = self.u_t = np.zeros(u.shape)
        v_t = self.v_t = np.zeros(v.shape)

        # Add the laplacian contribution
        self.laplacian_stencil_interior()
        u_t += self.u_laplacian
        v_t += self.v_laplacian

        # Add the reaction terms
        u_t += self._f1(u, v)
        v_t += self._f2(u, v)

        # Update the solution with forward Euler with numpy
        self.u = u + dt*u_t
        self.v = v + dt*v_t
        
        # for j in range(1, Ny-1):
        #     for i in range(1, Nx-1):
        #         u[i, j] += dt*((mu1/1) * (u[i+1, j] + u[i-1, j] + u[i, j+1] + u[i, j-1] - 4 * u[i, j]) + self._f1(u[i, j], v[i, j]))
        #         v[i, j] += dt*((mu2/1) * (v[i+1, j] + v[i-1, j] + v[i, j+1] + v[i, j-1] - 4 * v[i, j]) + self._f2(u[i, j], v[i, j]))



    def plot(self):
        # fig, ax = plt.subplots()
        ax = self.ax
        # plot u over the entire domain x and y based on coordinates
        im = ax.imshow(self.u, cmap='jet', interpolation='nearest', origin='lower', extent=[0, self.Nx, 0, self.Ny])
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('u')

        # show the plot
        plt.show()

        # Save the image as png file - name is the current time step
        self.fig.savefig(f"Python/figs/{self.ts}.png")

        # add_arts = [im, ]

        # # Store the figure to animate it at the end of the simulation in self.fig
        # self.figs.append(add_arts)

        


    def solve(self):
        u = self.u
        v = self.v

        # Create time distribution from T and dt
        t = np.arange(0, self.T, self.dt)

        for i, self.ts in enumerate(t):
            self._time_integrate()

            # Every 100 iterations, print the solution
            if i % 100 == 0:
            #     print(f"\nAt time {ts}")
            #     print(f"u = {self.u}")
            #     print(f"v = {self.v}")

                # Save the solution to the history
                self.u_hist.append(self.u)
                self.v_hist.append(self.v)

        # Plot the solution
        self.plot()
            
        # Convert the list of figures to a movie
        # movie = plt.figure(figsize=(10, 10))
        # self.anim = ArtistAnimation(self.fig, self.figs, interval=10, blit=True)

        # Save the movie to a file
        # self.anim.save("movie.mp4", writer='ffmpeg')



x = ReactionDiffusion()
x._set_params()
x._set_initial_conditions()
x.solve()