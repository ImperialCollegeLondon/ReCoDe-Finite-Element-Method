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

        elif self.degree == 3:
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
