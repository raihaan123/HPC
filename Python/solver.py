
class ReactionDiffusion:
    '''
    Class for solving the reaction diffusion equation.
    '''

    def spatial_bounded(self):
        '''
        Spatial derivative contribution for bounded nodes - uses central finite differencing
        '''

        u_xx = np.zeros(self.u.shape)
        u_yy = np.zeros(self.u.shape)
        v_xx = np.zeros(self.v.shape)
        v_yy = np.zeros(self.v.shape)

        
        # Find the spatial contributions to u and v
        for i in range(1, self.Nx-1):
            u_xx[1:-2, i] = dt/h**2 * mu1 * (u[2:-1, i] + u[:-3, i] -2*u[1:-2, i])
            v_xx[1:-2, i] = dt/h**2 * mu2 * (v[2:-1, i] + v[:-3, i] -2*v[1:-2, i])

        for j in range(1, self.Ny-1):
            u_yy[j, 1:-2] = dt/h**2 * mu1 * (u[j, 2:-1] + u[j, :-3] -2*u[j, 1:-2])
            v_yy[j, 1:-2] = dt/h**2 * mu2 * (v[j, 2:-1] + v[j, :-3] -2*v[j, 1:-2])

        # Sum the derivatives
        u_spatial = u_xx + u_yy
        v_spatial = v_xx + v_yy