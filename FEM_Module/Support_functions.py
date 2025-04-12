from plotly.subplots import make_subplots
import numpy as np
import plotly.graph_objects as go
import os
import pygmsh


def write_vtk(degree, mesh, solution_x, domain_dim, variable_name, variable_type, filename):
    """ Function to create a vtk file for the solution field to be viewed in Paraview

        Parameters
        -----------
        degree : degree of elements in the mesh
        mesh : pygmsh.geo.Geometry() object that contain the information about the elements and coordinates of the domain
        solution_x : the solution array x
        domain_dim : dimension of the domain space. In this case 2
        variable_name : name of the solution field in string format  that will be displayed in Paraview
        variable_type : "scalar" or "vector" for the type of solution field
        filename : name of the new vtk file to be saved in string format with no extension. Automatic extension of "vtk"
                    will be added
    """

    solution = solution_x

    filename = filename + ".vtk"
    num_nodes = len(mesh.points)

    if os.path.exists(filename):
        raise FileExistsError(f"File '{filename}' already exists. Please provide a different name")

    file = open(filename, "w")  # creae an empty file
    file.write("# vtk DataFile Version 3.0 ")
    file.write("\nFinite-element dataset: variable: " + variable_name + ", timestep: 100")
    file.write("\nASCII \n\nDATASET UNSTRUCTURED_GRID \nPOINTS " + str(num_nodes) + " float\n")
    ## add points to the file
    for point in mesh.points:
        file.write(str(round(point[0], 6)) + " " + str(round(point[1], 6)) + " " + str(round(point[2], 6)) + "\n")

    if degree == 1:
        element_types = [("line", 3), ("triangle", 5)]
        num_nodes_per_el = [2, 3]
    elif degree == 2:
        element_types = [("line3", 21), ("triangle6", 22)]
        num_nodes_per_el = [3, 6]

    num_line_elements = len(mesh.cells_dict[element_types[0][0]])
    num_triangular_elements = len(mesh.cells_dict[element_types[1][0]])

    file.write("\nCELLS " + str(
        num_line_elements + num_triangular_elements))  ## add total number of elements (line & triangular)
    len_cell_list = num_line_elements * (num_nodes_per_el[0] + 1) + num_triangular_elements * (num_nodes_per_el[1] + 1)
    file.write(" " + str(len_cell_list) + "\n")  ## add number of entries in this section
    ## (number of nodes in line element + 1 for element type ) * number of line elements

    ## add all the line elements
    for element in mesh.cells_dict[element_types[0][0]]:
        file.write(str(len(element)) + " ")
        for node in element:
            file.write(str(node) + " ")
        file.write("\n")

    ## add the triangular elements
    for element in mesh.cells_dict[element_types[1][0]]:
        file.write(str(len(element)) + " ")
        for node in element:
            file.write(str(node) + " ")
        file.write("\n")

    ## add element type (numerical value of the type)
    file.write("\nCELL_TYPES " + str(num_line_elements + num_triangular_elements) + "\n")
    for element in mesh.cells_dict[element_types[0][0]]:
        file.write(str(element_types[0][1]) + "\n")
    for element in mesh.cells_dict[element_types[1][0]]:
        file.write(str(element_types[1][1]) + "\n")

    ## add the value of the solution fiel for each node
    ## note the difference in set up based on whether the solution is scalar or vector
    file.write("\nPOINT_DATA " + str(len(mesh.points)) + "\n")

    if variable_type == "scalar":
        file.write("SCALARS " + variable_name + " float\n")
        file.write("LOOKUP_TABLE default\n")
        for value in solution:
            file.write(str(round(value, 10)) + "\n")

    if variable_type == "vector":
        file.write("VECTORS " + variable_name + " float\n")
        for value_ind in range(int(len(solution) / 2)):
            file.write(
                str(round(solution[value_ind * domain_dim], 10)) + " "
                + str(round(solution[value_ind * domain_dim + 1], 10)) + " "
                + str(round(0, 10)) + "\n")
    file.close()

def visualise_basis_fn(show_basis_fn, basis_functions):
    """ Function to visualise how the basis functions vary inside the quadratic isoparametric element
        -----------
        show_basis_fn : list of integers indicating the local node ID (up to 6). This will plot the corresponding basis functions
        basis_functions : function that calculates the basis function for a given point
        
    """
    degree = 2
    # Define the element nodes for unit triangle of degree 2
    nodes = np.array([
        [0, 0, 0],
        [0.5, 0, 0],
        [1, 0, 0],
        [0.5, 0.5, 0],
        [0, 1, 0],
        [0, 0.5, 0],
        [0, 0, 0],
    ])

    # Create a mesh grid of the coordinate system
    r = np.linspace(0, 1, 30)
    s = np.linspace(0, 1, 30)
    r, s = np.meshgrid(r, s)

    # Filter points to stay inside the triangle (r+s <= 1)
    x = r
    y = s
    mask = (x + y <= 1)
    x = x[mask]
    y = y[mask]

    # Create subplot figure with 1 row and 3 columns
    fig = make_subplots(
        rows=1, cols=3,
        specs=[[{'type': 'scatter3d'}, {'type': 'scatter3d'}, {'type': 'scatter3d'}]],
        subplot_titles=["Basis function " + str(show_basis_fn[0]), \
                        "Basis function" + str(show_basis_fn[1]), \
                        "Basis function " + str(show_basis_fn[2])]
    )

    # Iterate through every basis functions ID to calculate the basis function value and add it to the plot

    for i in range(len(show_basis_fn)):
        z = []
        basis_fn_id = show_basis_fn[i]
        # Get the values of basis function at points x_,y_
        for x_, y_ in zip(x, y):
            phi_all = basis_functions(degree, [x_, y_])
            z.append(phi_all[basis_fn_id])
        # Add the values to the plot
        fig.add_trace(go.Scatter3d(
            x=nodes[:, 0], y=nodes[:, 1], z=nodes[:, 2],
            mode='markers+lines',
            name=f'Unit Triangle {i}',
            marker=dict(color='red', size=7),
            line=dict(color='blue', width=4, dash='solid')
        ), row=1, col=i + 1)

        fig.add_trace(go.Mesh3d(
            x=x, y=y, z=z,
            color='lightblue', opacity=0.5,
            name=f'Basis Function {basis_fn_id}',
            flatshading=True
        ), row=1, col=i + 1)

    fig.update_layout(
        title="Basis Function Surfaces inside Unit Triangle",
        height=600, width=1200,
        scene=dict(
            xaxis_title="r",
            yaxis_title="s",
            zaxis_title="$\phi$",
            aspectmode="cube"
        ),
        # zoom out in each plot
        margin=dict(l=0, r=0, b=0, t=40),
        scene1=dict(camera=dict(
            eye=dict(x=2, y=2, z=2))),
        scene2=dict(camera=dict(
            eye=dict(x=2, y=2, z=2))),
        scene3=dict(camera=dict(
            eye=dict(x=2, y=2, z=2))),

        showlegend=False
    )

    # Show figure
    fig.show()
