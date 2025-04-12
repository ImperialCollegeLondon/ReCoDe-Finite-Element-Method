import pygmsh
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.tri import Triangulation


class Run_Cantilever:
    """ Class containing all the steps to run a Finite Element Method for
    2D Cantilever example
    """

    def __init__(self, domain_dim):
        """ To set up the class the key input is the dimension of the domain,
        whether the domain is in 2D or 3D"""
        self.solution = np.empty((0, domain_dim))
        self.domain_dim = domain_dim

        self.Youngs_modulus = 1e6
        self.Poisson_ratio = 0.25

        self.degree = 2

        self.element_list_2D = []
        self.global_node_coords = []
        self.x_start = 0
        self.x_end = 10
        self.c=3

        self.num_nodes = 0.
        self.A_matrix = np.empty((2,2))  ## note that for each node there are 2 d.o.f
        self.b = np.zeros((2,))
        self.x = np.zeros((2,))

        self.nodes_on_left_ = []
        self.nodes_on_right_ = []
        self.elements_left_ =  []
        self.elements_right_ = []

        self.P = 80
        
        


    def set_material_properties(self,Young_modulus, Poisson_ratio):
        self.Youngs_modulus = Young_modulus
        self.Poisson_ratio = Poisson_ratio
        

    def domain_mesh(self,  element_degree, refinement, x_limits, y_limits,
                    z_limits=[0, 0]):
        """This function creates the domain boundary and meshes the domain
        based on the input arguments.The domain is assumed to be rectangular
        shape.

        Parameters
        ----------
        element_degree (int) : 1 - linear elements, 2- quadratic and so on.
        refinement (float) : the size of the elements in units.
                                Larger number leads to coarser mesh.
        x_limits (list) : [x_start, x_end] are the lower and upper boundary
                                limits of the length along x-axis
        y_limits (list) : [y_start, y_end] are the lower and upper boundary
                                limits of the height along y-axis
        z_limits (list) : [z_start, z_end] are the lower and upper boundary
                                limits of the width along z-axis.
                              For 2D problem the z_limits are to be left
                              at default [0,0]

        Returns
        -------
        self.mesh : mygmsh.geo.Geometry that contains the list of elements
        with the corresponding nodes and node coordinates

        Example:
            >>> domain_mesh(2, 1., [0,10], [1,5])


        """
        if element_degree not in (1, 2):
            raise ValueError("For this exercise degree must be linear (1) or quadratic (2)")
        self.degree = element_degree

        # Create the pygmsh mesh object with the specified input arguments
        with pygmsh.geo.Geometry() as geom:
            geom.add_rectangle(xmin=x_limits[0], xmax=x_limits[1],
                               ymin=y_limits[0], ymax=y_limits[1],
                               z=z_limits[0], mesh_size=refinement)

            # Generate the mesh
            self.mesh = geom.generate_mesh(dim=self.domain_dim, algorithm=6,
                                           order=element_degree)

        # Define the name of elements in the mesh to be used in subsequent
        # functions
        if self.degree == 1:
            self.element_name = 'triangle'
        elif self.degree == 2:
            self.element_name = 'triangle6'

        # List of 2D elements to be used in accumulation
        self.element_list_2D = self.mesh.cells_dict[self.element_name]
        self.global_node_coords = self.mesh.points
        self.x_start = x_limits[0]
        self.x_end = x_limits[1]
        self.c=(y_limits[1]-y_limits[0])/2.
        
        
        return self.mesh

    def visualise_mesh(self):
        """Plot the mesh created in function domain_mesh() to visualise the
        domain, the refinement level and the nodes. The function requires the
        domain is created using domain_mesh() first.

        TO DO : Adjust the number of corners if extended to 3D

        Parameters
        ----------
        None

        Returns
        -------
        None

        Example:
            >>> visualise_mesh()

        """
        # Physical coordinates of the mesh nodes and the list of elements
        points = self.mesh.points
        cells = self.mesh.cells_dict

        # The Triangulation function needs only corners of elements, which
        # are the first three corners for triangulared mesh
        # For 2D elements this will be the first 3 nodes
        if self.element_name not in cells:
            raise ValueError("The mesh does not contain triangular elements.")
        triangles = np.array(cells[self.element_name])[:, :3]

        # Prepare data for Triangulation
        x, y = points[:, 0], points[:, 1]
        triangulation = Triangulation(x, y, triangles)

        # Plot the mesh
        plt.figure()
        plt.triplot(triangulation, color='blue', lw=0.8)
        plt.scatter(x, y, color='red', s=10, zorder=5)
        plt.title("Mesh Visualization")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis("equal")
        plt.grid(True)
        plt.show()

    def integration_points(self):
        """
        The Gaussian integration points for isoparametric 2D triangles.
        Currently supports only quadratic and cubic elements.

        TO DO: change to 3D points if extended to 3D.

        Parameters
        ----------
        None

        Returns
        -------
        Points : numpy.array of shape (n,2), where n depends on element degree.
        weights : numpy.array of shape (1,n)

        """
        Points = np.array([[1 / 6., 2 / 3.], [1 / 6., 1 / 6],
                               [2 / 3., 1 / 6.]])

        weights = np.array([1., 1., 1.]) * 1 / 6.


        return Points, weights

    def basis_functions_der(self, point):
        """
        Calculates the derivatives of the basis functions are the specified
        point. This point must lie inside isoparametric triangle using local
        coordinate system. Currently only 2D basis functions are considered.

        Parameters
        ----------
        point (list of floats) : a list of two coordinates  in local coordinate
                                system in range [0,1]

        Returns
        -------
        basis_functions : numpy.array of shape (2,n) where n is determined by
                        the elmeent degree

        """

        # linear elements
        if self.degree == 1:
            dNr = self.basis_functions_dNr(point)
            dNs = self.basis_functions_dNs(point)

        # quadratic elements
        if self.degree == 2:
            dNr = self.basis_functions_dNr(point)
            dNs = self.basis_functions_dNs(point)

        # stack the derivatives in two directions into a matrix
        basis_functions = np.vstack((dNr, dNs))

        return basis_functions

    def basis_functions_dNr(self, point):
        """
        Calculate the derivatives of basis functions in r direction (local
        coordinates).

        Parameters
        ----------
        point (list of floats) : a list of two coordinates  in local coordinate
                                system in range [0,1]

        Returns
        -------
        dNr : numpy.array of shape (1,n), where n is determined by the element
            degree

        """

        r = point[0]
        s = point[1]

        # linear element
        if self.degree == 1:
            dNr = np.zeros((3,))
            dNr[0] = -1
            dNr[1] = 1
            dNr[2] = 0.

        # quadratic elements
        elif self.degree == 2:
            dNr = np.zeros((6,))
            # Corner nodes first
            dNr[0] = -3. + 4. * r + 4. * s
            dNr[1] = -1. + 4. * r
            dNr[2] = 0.

            # Midside nodes
            dNr[3] = 4. - 8. * r - 4. * s
            dNr[4] = 4. * s
            dNr[5] = -4. * s

        return dNr

    def basis_functions_dNs(self, point):
        """
        Calculate the derivatives of basis functions in s direction (local
        coordinates).

        Parameters
        ----------
        point (list of floats) : a list of two coordinates  in local coordinate
                                system in range [0,1]

        Returns
        -------
        dNr : numpy.array of shape (1,n), where n is determined by the element
            degree

        """
        r = point[0]
        s = point[1]

        # linear elements
        if self.degree == 1:
            dNs = np.zeros((3,))
            dNs[0] = -1
            dNs[1] = 0.
            dNs[2] = 1.

        # quadratic elements
        if self.degree == 2:
            dNs = np.zeros((6,))
            # Corner nodes first
            dNs[0] = -3. + 4. * r + 4. * s
            dNs[1] = 0.
            dNs[2] = -1. + 4. * s

            # Midside nodes
            dNs[3] = -4. * r
            dNs[4] = 4. * r
            dNs[5] = 4. - 4. * r - 8. * s

        return dNs

    def basis_fun_line_quad(self,point):
        r = point
        N=np.zeros(3)
        N[0] = 0.5*(r*r -r)
        N[1] = 0.5*(r*r +r)
        N[2] = 1-r*r
        return N

    def der_basis_fun_line_quad(self,point):
        r = point
        dN=np.zeros(3)
        dN[0] = r-0.5
        dN[1] = r+0.5
        dN[2] = -2*r
        return dN


    def integration_points_1D(self):
        IP_line = [-1/np.sqrt(3.),  1/np.sqrt(3.)]
        IP_w = [1.,1.]
        return IP_line, IP_w


    def jacobian_1D(self,pt, element_xy_coord):
        dN = self.der_basis_fun_line_quad(pt)
        x_vals=element_xy_coord[:,1]
        return np.dot(dN,x_vals )
       
       
       
       
    def jacobian(self, point, corner_nodes):
        """
        Calculate the Jacobian matrix for mapping between the global and local
        coordinate systems. The Jacobian matrix is calculated at point in local
        coordinate system.

        TO DO: Extend the function to 3D
        Parameters
        ----------
        point (list of floats) : a list of two coordinates  in local coordinate
                                system in range [0,1]

        corner_nodes (list of list of floats): a list of coordinates of corner
                    nodes in the element in global coordinate system

        Returns
        -------
        det (float) : scale of volumetric differences between isoparametric
                      element and element in global coordinate system

        Jacobian_mat (numpy.array) : square matrix of size (dim,dim), where dim
                    is the spatial dimension of the domain

        """

        Jacobian_mat = np.zeros((self.domain_dim, self.domain_dim))

        # Get the derivatives of the basis functions at point
        dnr = self.basis_functions_dNr(point)
        dns = self.basis_functions_dNs(point)

        # take dot product between the derivatives and the corner nodes for
        # each axis
        for dim in range(self.domain_dim):
            Jacobian_mat[0, dim] = np.dot(dnr, corner_nodes[:, dim])
            Jacobian_mat[1, dim] = np.dot(dns, corner_nodes[:, dim])

        # get the determinant of the Jacobian matrix
        det = np.linalg.det(Jacobian_mat)
        return det, Jacobian_mat

    def dN_To2DOF(self, B):
        """
        This functions splits the derivative matrix B into the strain
        components. Depending on the dimensions of thedomain. The first row is
        the strain component in x-direction and the second row is the strain in
        y-direction.
        If original B matrix for element with n nodes is

            dN1/dx dN2/dx ... dNn/dx
            dN1/dy dN2/dy ... dNn/dy

        The converted strain matrix is
            dN1/dx    0     dN2/dx   0   ...  dNn/dx    0
              0     dN1/dy     0   dN2/dy...   0     dNn/dy
            dN1/dy  dN1/dx  dN2/dy dN2/dx... dNn/dy  dNn/dx

        Parameters
        ----------
        B (numpy matrix) : matrix containing the shape function derivatives at
        a point. The expected shape is (dim,n), where dim is domain dimension
        and n is the number of nodes in the element

        Returns
        -------
        new_B (numpy matrix) : converted matrix to strain displacement matrix
        in shape (dim C (dim-1), n * dim), where
        C is the choose function, so for dim 2 the shape will be (3,n)


        """
        shape = np.shape(B)
        new_B = np.zeros((3, shape[1] * 2))

        for i in range(shape[1]):
            new_B[0][i * 2] = B[0][i]
            new_B[1][i * 2 + 1] = B[1][i]

            new_B[2][i * 2] = B[1][i]  # /dy
            new_B[2][i * 2 + 1] = B[0][i]  # /dx

        return new_B

    def stiffness_matrix_2D(self, E, nu):
        """
        This function calculates the stiffness matrix for each element based on
        the material properties.

        Parameters
        ----------
        E (float) : the Young's modulus of the domain, this measure the
        stiffness of the material in response to uniaxial stress
        nu (float) : the Poisson ratio of the domain, this measures the
        lateral strain relative to axial strain,when domain is subjected to
        uniaxial stress.

        Returns
        -------
        D (numpy matrix) : the stiffness matrix for 2D domain of shape (3,3)


        """
        # 2D case
        D = np.zeros((3, 3))

        # calculate the coefficients of the matrix
        frac = nu / (1. - nu)
        coeff = (E * (1 - nu)) / ((1. + nu) * (1. - 2. * nu))

        # fill the diagonals of the matrix
        np.fill_diagonal(D, 1)
        D[2][2] = (1. - 2. * nu) / (2. * (1. - nu))

        # fill the off diagonals of the matrix
        D[0, 1] = frac
        D[1, 0] = frac

        D = D * coeff
        return D
        
    def accumulate_global_matrix(self, Youngs_modulus, Poisson_ratio):
        
        self.set_material_properties(Youngs_modulus, Poisson_ratio)
        
        
        self.num_nodes = len(self.mesh.points)
        self.A_matrix = np.zeros((self.num_nodes*self.domain_dim, self.num_nodes*self.domain_dim)) ## note that for each node there are 2 d.o.f
        self.b = np.zeros((self.num_nodes*self.domain_dim,))

        ## Now we iterate through each element in the mesh in the same way before
        for element in self.element_list_2D:
            
            # ------------ User Input required ------------------------------------------#
            # -- Step 1. Set the integration points IP and their weights 
            IPs, IP_weights = self.integration_points()

            # -- Step 2. Create empty LHS matrix for the element stiffness matrix ------------#
            # -- recall that variable element contains the list of the nodes and that each node has 2 d.o.f
            LHS = np.zeros((len(element)*self.domain_dim, len(element)*self.domain_dim)) 
                            
            # -- Step 3. Set up for loop to iterate over each integration point ------------#
            for IP, weight in zip(IPs,IP_weights):
                            
                # -- Step 4. Get the list of corner nodes in the element and calculate
                # --------   the Jacobian matrix using function Jacobian in FEM_module
                e_nodes_xy = self.global_node_coords[element][:,:3]
                detJ, JMat = self.jacobian(IP, e_nodes_xy)
                
                # -- Step 5. Calculate the derivative basis functions at the integration point 
                # --------   use the function basis_functions_der() and multiply it by the inverse of the Jacobian matrix
                der_matrix = self.basis_functions_der(IP)
                der_matrix_JInv = np.matmul(np.linalg.inv(JMat), der_matrix)
                
                # -- Step 6. Transform the derivative matrix into the B described above
                # --------  use the function dN_To2DOF() in FEM_module
                B = self.dN_To2DOF(der_matrix_JInv) 
                # -- Step 7. Take the transpose of the matrix B
                BT =  B.transpose()
                # -- Step 8. Calculate matrix D. You can use the function stiffness_matrix_2D in FEM module
                # --------   use the material properties Youngs_modulus and Poisson_ratio defined above
                D = self.stiffness_matrix_2D(self.Youngs_modulus,self.Poisson_ratio)        
                
                # -- Step 9. Multiply BT by D
                BT = np.matmul(BT, D)
                # -- Step 10. Multiply the matrix from Step 9 by B
                BT = np.matmul(BT, B)
                # -- Step 11. Multiply the final matrix by Jacobian determinant and weight of the integration point
                BT = BT * weight*(detJ) ## currently assume all the weight of IP are the same
                # -- Step 12. Add the resulting matrix to the element stiffness matrix LHS
                LHS += BT     
            # -- Step 13. Accumulate the element stiffness matrix LHS to global matrix A
    # ----------  The accumulation of the first d.o.f for the node is given as a hint
            for i in range(len(element)):
                for j in range(len(element)):
                    col = int(element[i])
                    row = int(element[j])
                    self.A_matrix[row*self.domain_dim,col*self.domain_dim]+= LHS[int(j*self.domain_dim)][int(i*self.domain_dim)] # xx
                    self.A_matrix[row*self.domain_dim,col*self.domain_dim+1]+= LHS[int(j*self.domain_dim)][int(i*self.domain_dim+1)] # xy
                    self.A_matrix[row*self.domain_dim+1,col*self.domain_dim]+= LHS[int(j*self.domain_dim+1)][int(i*self.domain_dim)] # y
                    self.A_matrix[row*self.domain_dim+1,col*self.domain_dim+1]+= LHS[int(j*self.domain_dim+1)][int(i*self.domain_dim+1)] # y

    def apply_BC_NeumannDirich(self, P):
        # ------------ User Input required ------------------------------------------#
        self.nodes_on_left_ = np.where(np.isclose(self.mesh.points[:,0], self.x_start))[0]
        self.nodes_on_right_ =np.where(np.isclose(self.mesh.points[:,0], self.x_end))[0]
        

        # ------------ End of User Input ------------------------------------------#
            
        # Find line elements that have all nodes the list of nodes_on_left_ boundary
        result = np.all(np.isin(self.mesh.cells_dict['line3'], self.nodes_on_left_), axis=1)
        self.elements_left_=self.mesh.cells_dict['line3'][np.where(result)[0]]

        # Find line elements that have all nodes the list of nodes_on_left_ boundary
        result = np.all(np.isin(self.mesh.cells_dict['line3'], self.nodes_on_right_), axis=1)
        self.elements_right_=self.mesh.cells_dict['line3'][np.where(result)[0]]

        # set force value
        self.P = P
        # Iterate through all the elements on the left boundary
        for element in self.elements_left_: 

            # ------------ User Input required ------------------------------------------#
            # -- Step 1. Set the integration points IP and their weights for 1D element
            IPs, IP_weights = self.integration_points_1D()

            # -- Step 2. Create empty RHS vector for element level force accumulation ------------#
            # -- recall that variable element contains the list of the nodes and that each node has 2 d.o.f
            RHS = np.zeros((len(element)*self.domain_dim)) 

            # -- Step 3. Set up for loop to iterate over each integration point ------------#
            for IP, weight in zip(IPs,IP_weights):

                # take the y-coordinate of the middle node as the average of the element
                y_value_ = self.mesh.points[element[2]][1] 
                # initialise the empty traction vector in x- and y-directions respectively
                traction = [0,0]
                # set the traction in the x-direction
                traction[0]= 0.
                 # calculate the traction in y-direction
                traction[1] = (3. / 4.) * (self.P / self.c) * (1. - pow((y_value_ / self.c), 2)) *(-1)
            
                # -- Step 4. Get the list of corner nodes in the element and calculate
                # --------   the Jacobian matrix in 1D using function Jacobian_1D
                e_nodes_xy = self.global_node_coords[element][:,:3]
                detJ = self.jacobian_1D(IP, e_nodes_xy)

                # -- Step 5. Calculate the basis functions of 1D element
                N = self.basis_fun_line_quad(IP)
                
                # -- Step 6. Create an empty vector for adjusted shape functions
                Nnc = np.zeros(len(N)*self.domain_dim) 
                
                # -- Step 7. Adjust the basis function vector for 2 d.o.f 
                # ---------   and multiply by the traction force in that direction
                # iterate through all the nodes in element
                for node_i in range(len(N)):
                    # for each node consider two degrees of freedom
                    for dof in range(self.domain_dim):  # for every dimension
                            Nnc[node_i * 2 + dof] += N[node_i] * traction[dof]
                # -- Step 8. Multiply the vector by determinant and the integration point weight
                Nnc = Nnc* weight * (detJ) 
                
                # add to the RHS vector for accumulation
                RHS+=Nnc
            
            # add the element level RHS vector to the global RHS vector b
            for i in range(len(element)):
                node = int(element[i])
                self.b[node*self.domain_dim] += RHS[i*self.domain_dim]
                self.b[node*self.domain_dim+1] += RHS[(i*self.domain_dim)+1]

        # check all the nodes along the right boundary
        for node in self.nodes_on_right_:
            
            # get the global coordinate of each node
            x_val = self.mesh.points[node][0]
            y_val = self.mesh.points[node][1]
            # ------------ User Input required ------------------------------------------#
            
            # -- Step 1. Check if y coordinate is equal to 0
            if(np.isclose(y_val,0)):
                
                # -- Step 2. Set the b vector to 0 at both (x- and y-) degrees of freedom of that node
                # ---------  Hint. Do not forget to set the A_matrix rows to identity row
                self.A_matrix[node*2,:] = 0
                self.A_matrix[node*2, node*2] = 1
                self.b[node*2 ] = 0

                self.A_matrix[node*2+1,:] = 0
                self.A_matrix[node*2+1, node*2+1] = 1
                self.b[node*2 ] = 0


            # -- Step 3. Check if y coordinate is equal to either c or -c
            if(np.isclose(abs(y_val),self.c)):
                
                # -- Step 4. Set the b vector to 0 at degrees of freedom  in x-direction for that node
                # ---------  Hint. Do not forget to set the A_matrix rows to identity row
                self.A_matrix[node*self.domain_dim,:] = 0
                self.A_matrix[node*self.domain_dim, node*self.domain_dim] = 1
                self.b[node*self.domain_dim] = 0


    def apply_BC_Dirich(self, P):
        # ------------ User Input required ------------------------------------------#
        self.nodes_on_left_ = np.where(np.isclose(self.mesh.points[:, 0], self.x_start))[0]
        self.nodes_on_right_ = np.where(np.isclose(self.mesh.points[:, 0], self.x_end))[0]

        # ------------ End of User Input ------------------------------------------#

        # Find line elements that have all nodes the list of nodes_on_left_ boundary
        result = np.all(np.isin(self.mesh.cells_dict['line3'], self.nodes_on_left_), axis=1)
        self.elements_left_ = self.mesh.cells_dict['line3'][np.where(result)[0]]

        # Find line elements that have all nodes the list of nodes_on_left_ boundary
        result = np.all(np.isin(self.mesh.cells_dict['line3'], self.nodes_on_right_), axis=1)
        self.elements_right_ = self.mesh.cells_dict['line3'][np.where(result)[0]]

        # set force value
        self.P = P
        L = self.x_end - self.x_start

        # Iterate through all the elements on the left boundary
        # check all the nodes along the right boundary
        for node in self.nodes_on_right_:

            # get the global coordinate of each node
            x_val = self.mesh.points[node][0]
            y_val = self.mesh.points[node][1]
            # ------------ User Input required ------------------------------------------#

            # -- Step 1. Check if y coordinate is equal to 0
            if (np.isclose(y_val, 0)):
                # -- Step 2. Set the b vector to 0 at both (x- and y-) degrees of freedom of that node
                # ---------  Hint. Do not forget to set the A_matrix rows to identity row
                self.A_matrix[node * 2, :] = 0
                self.A_matrix[node * 2, node * 2] = 1
                self.b[node * 2] = 0

                self.A_matrix[node * 2 + 1, :] = 0
                self.A_matrix[node * 2 + 1, node * 2 + 1] = 1
                self.b[node * 2] = 0

            # -- Step 3. Check if y coordinate is equal to either c or -c
            if (np.isclose(abs(y_val), self.c)):
                # -- Step 4. Set the b vector to 0 at degrees of freedom  in x-direction for that node
                # ---------  Hint. Do not forget to set the A_matrix rows to identity row
                self.A_matrix[node * self.domain_dim, :] = 0
                self.A_matrix[node * self.domain_dim, node * self.domain_dim] = 1
                self.b[node * self.domain_dim] = 0
        for node in self.nodes_on_left_:
            # get the global coordinate of each node
            x_val = self.mesh.points[node][0]
            y_val = self.mesh.points[node][1]

            value_at_x = self.analytical_solution_x(x_val, y_val, P, L, self.Youngs_modulus, self.Poisson_ratio, self.c,
                                                    1.)
            value_at_y = self.analytical_solution_y(x_val, y_val, P, L, self.Youngs_modulus, self.Poisson_ratio, self.c,
                                                    1.)

            self.A_matrix[node * self.domain_dim, :] = 0
            self.A_matrix[node * self.domain_dim, node * self.domain_dim] = 1.
            self.b[node * self.domain_dim] = value_at_x

            self.A_matrix[node * self.domain_dim + 1, :] = 0
            self.A_matrix[node * self.domain_dim + 1, node * self.domain_dim + 1] = 1.
            self.b[node * self.domain_dim + 1] = value_at_y



    def solve(self, A_matrix, b):
        """
        Solve the linear system of equations that is accumulated element by
        element during the assembly step

        Parameters
        ----------
        A_matrix (numpy matrix) : a square matrix contains all the element
        contributions.

        b (numpy row vector) : the right-hand-side vector of the linear system
        of equations. This contains the boundary conditions information of the
        system.

        Returns
        -------
        x (numpy vector) :  the solution vector that contains the nodal
        solutions of the differential system of equations.

        """
        self.x = np.linalg.solve(A_matrix, b)
        return self.x

    def analytical_solution_x(self, x, y, P, L, E, Poisson, c, t):
        """ Analytical solution in x -direction at point (x,y)

        Parameters
        -----------
        x,y (floats) : global coordinates of the point
        P (float) : load applied at the end of the beam
        L (float) : length of the beam
        E (float) : Youngs modulus
        Poisson (float) : Poisson ratio of material
        c (float) : Half of the height of the beam
        t (float) :
        """
        I = (2 * t * c**3) / 3.
        G = E / (2 * (1 + Poisson))
        part0 = (P * (x**2 - L**2) * y) / (2 * E * I)
        part1 = (Poisson * P * y * (y**2 - c**2)) / (6 * E * I)
        part2 = (P * y * (y**2 - c**2)) / (6 * G * I)

        return -part0-part1+part2

    def analytical_solution_y(self, x, y, P, L, E, Poisson, c, t):
        """ Analytical solution in y -direction at point (x,y)

        Parameters
        -----------
        x,y (floats) : global coordinates of the point
        P (float) : load applied at the end of the beam
        L (float) : length of the beam
        E (float) : Youngs modulus
        Poisson (float) : Poisson ratio of material
        c (float) : Half of the height of the beam
        t (float) :
        """
        I = (2 * t * c * c**2) / 3
        G = E/(2*(1+Poisson))
        part0 = (Poisson * P * x * y**2) / (2 * E * I)
        part1 = (P * (x**3 - L**3)) / (6 * E * I)
        part2 = (P * L**2) / (2 * E * I)
        part3 = (Poisson * P * c**2) / (6 * E * I)
        part4 = (P * c**2) / (3 * G * I)

        return part0+part1-(part2+part3+part4)*(x-L)

    def calculate_errors(self, P, L,c):

        # Set up empty containers for errors
        error_abs = np.empty((0, 2))
        error_mse = np.empty((0, 2))

        for node_ID in range(len(self.mesh.points)):
            x_val, y_val = self.mesh.points[node_ID][0], self.mesh.points[node_ID][1]

            # Get analytical solution for the node
            analytical_disp_x = self.analytical_solution_x(x_val, y_val, P, L, self.Youngs_modulus, self.Poisson_ratio, c, 1)
            analytical_disp_y = self.analytical_solution_y(x_val, y_val, P, L, self.Youngs_modulus, self.Poisson_ratio, c, 1)

            # Get the numerical displacement at the node from solution vector x
            numerical_disp_x = self.x[node_ID * self.domain_dim]
            numerical_disp_y = self.x[node_ID * self.domain_dim + 1]

            # Absolute error:
            error_x = abs(numerical_disp_x - analytical_disp_x)
            error_y = abs(numerical_disp_y - analytical_disp_y)

            error_abs = np.append(error_abs, [[error_x, error_y]], axis=0)

            # MSE error:
            error_x = pow((numerical_disp_x - analytical_disp_x), 2)
            error_y = pow((numerical_disp_y - analytical_disp_y), 2)

            error_mse = np.append(error_mse, [[error_x, error_y]], axis=0)


        errors_summary_x =np.array([len(self.mesh.cells_dict['triangle6']), np.max(error_abs[:, 0]), \
                                         np.mean(error_abs[:, 0]), np.mean(error_mse[:, 0])])

        errors_summary_y = np.array([len(self.mesh.cells_dict['triangle6']), np.max(error_abs[:, 1]), \
                                         np.mean(error_abs[:, 1]), np.mean(error_mse[:, 1])])

        return errors_summary_x,errors_summary_y