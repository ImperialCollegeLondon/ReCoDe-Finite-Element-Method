import pygmsh
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.tri import Triangulation


class FEM_model:
    """ Class containing the key functions required to run Finite Element
    Method model
    """

    def __init__(self, domain_dim):
        """ To set up the class the key input is the dimension of the domain,
        whether the domain is in 2D or 3D"""
        self.solution = np.empty((0, domain_dim))
        self.domain_dim = domain_dim
        

    def domain_mesh(self,  element_degree, refinement, x_limits, y_limits,
                    z_limits=[0, 0]):
        """This function creates the domain boundary and meshes the domain
        based on the input argumnts.The domain is assumed to be rectangular
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

        self.degree = element_degree

        # Create the pygmsh mesh objec with the specified input arguments
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
        return self.mesh

    def visualise_mesh(self):
        """Plot the mesh created in function domain_mesh() to visualise the
        domain, the refinement level and the nodes. The function requires the
        domain is created using domain_mesh() first.

        TO DO : Adjust the number of corners if exteded to 3D

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
        if self.element_name in cells:
            triangles = np.array(cells[self.element_name])[:, :3]
        else:
            raise ValueError("The mesh does not contain triangular elements.")

        # Prepare data for Triangulation
        x, y = points[:, 0], points[:, 1]
        triangulation = Triangulation(x, y, triangles)

        # Plot the mesh
        plt.figure(figsize=(10, 5))
        plt.triplot(triangulation, color='blue', lw=0.8)
        plt.scatter(x, y, color='red', s=10, zorder=5)
        plt.title("Mesh Visualization")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.axis("equal")
        plt.grid(True)
        plt.show()

    def integrationPoints(self):
        """
        The Gaussian integeration points for isoparametric 2D triangles.
        Currently supports only quadratic and cubic elements.

        TO DO: change to 3D points if extended to 3D.

        Parameters
        ----------
        None

        Returns
        -------
        Points : numpy.array of shape (n,2), where n depends on element degree.
        weights : numpy.array of shape (1,n)

        Example:
            >>> integrationPoints()
            np.array([[1 / 6., 2 / 3.], [1 / 6., 1 / 6],[2 / 3., 1 / 6.]]),
            np.array([1., 1., 1., 1.]) * 1 / 6.

        """
        if self.degree == 2:
            Points = np.array([[1 / 6., 2 / 3.], [1 / 6., 1 / 6],
                               [2 / 3., 1 / 6.]])

            weights = np.array([1., 1., 1., 1.]) * 1 / 6.

        # elif self.degree == 3:
        Points = np.array([
                [0.4459484909, 0.4459484909],
                [0.4459484909, 0.1081030182],
                [0.1081030182, 0.4459484909],
                [0.0915762135, 0.0915762135],
                [0.0915762135, 0.8168475730],
                [0.8168475730, 0.0915762135]
        ])
        weights = np.array([
                0.2233815897,
                0.2233815897,
                0.2233815897,
                0.1099517437,
                0.1099517437,
                0.1099517437
        ])

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

        # linear eleemnt
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

    def Jacobian(self, point, corner_nodes):
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
                    nodes in the elment in global coordinate system

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

    def solve(self, A_matrix, b):
        """
        Solve the linear system of equations that is accumulated element by
        element during the assembly step

        Parameters
        ----------
        A_matrix (numpy matrix) : a square matrix contains all the element
        contributions.

        b (numpy row vector) : the right-hand-side vector of the linear system
        of equations. This contains the boundary condtions information of the
        system.

        Returns
        -------
        x (numpy vector) :  the solution vector that contains the nodal
        solutions of the differential system of equations.

        """
        x = np.linalg.solve(A_matrix, b)
        return x

