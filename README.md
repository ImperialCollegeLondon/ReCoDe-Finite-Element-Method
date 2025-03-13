<!-- Your Project title, make it sound catchy! -->

# Numerical Simulation using Finite Element Method 

<!-- Provide a short description to your project -->

## Description

Finite Element Method (FEM) is a powerful computational method used across a wide variety of disciplines, including engineering, mathematics and earth science. It is used to solve complex differential equations and simulate real-world scenarios. This project focuses on application of FEM in solid mechanics, but the teachings can be extended to other areas of study.
It provides a comprehensive, step-by-step tutorial designed to introduce users to the fundamentals of FEM and its practical implementation in Python. The tutorial takes users through all the key stages of numerical simulations, starting from domain meshing and problem discretisation, moving to solving the differential equations and finally analysing and understanding the results.  
To make the learning process intuitive, the project begins with a straightforward example that demonstrates the methodology and code implementation. As users progress, they are gradually introduced to more advanced aspects of FEM, ensuring a solid foundation and deeper understanding of this method.

<!-- What should the students going through your exemplar learn -->

## Learning Outcomes

- Develop understanding of the fundaments of FEM
- Learn how to create custom FEM code and adapt it to solve partial differential equations
- Gain familiarity with existing libraries and tools available for FEM implementation
- Develop skills to interpret solutions and effectively analyse the simulation results

<!-- How long should they spend reading and practising using your Code.
Provide your best estimate -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

## Requirements

#### Academic
To sucessfully complete this examplar the user should have basic understanding of key mathematical concepts, including linear algebra (such as matrices, vectors, and determinants) and differential equations (both ordinary and partial differential equations, linear system of equations)

Familiarity with Python programming, using existing libraries, writing classes, functions and graph plotting.

<!-- List the system requirements and how to obtain them, that can be as simple
as adding a hyperlink to as detailed as writting step-by-step instructions.
How detailed the instructions should be will vary on a case-by-case basis.

Here are some examples:

- 50 GB of disk space to hold Dataset X
- Anaconda
- Python 3.11 or newer
- Access to the HPC
- PETSc v3.16
- gfortran compiler
- Paraview
-->

### System

System requirements:
- Python 3.11 or newer
- Paraview 4.0 or newer
- Jupyter Notebooks

Python libraries used:
- pygmsh
- numpy
- matplotlib
  
<!-- Instructions on how the student should start going through the exemplar.

Structure this section as you see fit but try to be clear, concise and accurate
when writing your instructions.

For example:
Start by watching the introduction video,
then study Jupyter notebooks 1-3 in the `intro` folder
and attempt to complete exercise 1a and 1b.

Once done, start going through through the PDF in the `main` folder.
By the end of it you should be able to solve exercises 2 to 4.

A final exercise can be found in the `final` folder.

Solutions to the above can be found in `solutions`.
-->

## Getting Started

<!-- An overview of the files and folder in the exemplar.
Not all files and directories need to be listed, just the important
sections of your project, like the learning material, the code, the tests, etc.

A good starting point is using the command `tree` in a terminal(Unix),
copying its output and then removing the unimportant parts.

You can use ellipsis (...) to suggest that there are more files or folders
in a tree node.

-->
A significant part of our understanding of the real world comes from differential equations, which describe the behaviour of various quantities. These equations are fundamental across many scientific fields, including physics, geoscience, and medicine, as they explain how systems behave under specific conditions, how they respond to external forces and change with time.

Solving differential equations allows us to predict and analyse real world phenomena. In many cases, computer simulations are used to solve these equations for specific scenarios. Engineers use computer simulations to optimize product designs, ensuring they can withstand real world conditions, while geoscientists use them to study complex natural processes or our interactions with them, such as energy storage. These simulations enable us to test hypotheses and evaluate outcomes without the need for expensive real life experiments. However, their accuracy and reliability depend on our ability to effectively solve the underlying systems of differential equations.

In certain simplified cases, differential equations can be solved analytically, providing exact solutions. However, these cases are often over simplified and may not be flexible enough for practical applications. When analytical solutions are not feasible, we use numerical methods to approximate solutions. A wide range of numerical techniques exists, each designed for solving specific  problems and work to improve accuracy and efficiency of the solutions..

One of the most widely used numerical methods is the Finite Element Method (FEM). It is particularly popular in mechanics, including structural and rock mechanics, fluid dynamics, and heat transfer. Since its introduction in the 1950s, FEM has become an essential tool in numerical modeling, with extensive research dedicated to refining and advancing its capabilities.  

The core principle of Finite Element Method is discretising the problem domain into smaller interconnected elements. This process, known as discretisation, allows the problem to be described at the element level. Each element has an associated number of nodes, which are shared by adjacent elements. These nodes form relationships between elements that can be aggregated to a global system of equations.  We By solving this system, we obtain an approximate solution for the entire domain at the nodes. Using interpolation we can then calculate the solution at any point in the domain.

Therefore the key stages of Finite Element Methods are:

- Domain discretisation  -  dividing the domain into finite elements
- Element level problem formulation - defining the equations on each element using basis functions
- Assembly to global system - combining the element equations into a global system
- Application of boundary conditions - applying constrains and conditions to the problem to solve for unique solutions.
- Solving linear system of equations - using a numerical solver to find the solution
- Post processing of solution - interpreting, visualising and analysing the solution

This exercise will guide you through the key stages of constructing a Finite Element code. It will focus on a simple problem of heat transfer to build the basics of the method. At the end of this exercise you should have a good understanding of how Finite Element Method works. The next exercise will build on these concepts and look at solving a slightly more complex problem. 

## Project Structure

```log
.
C:.
|   .gitignore
|   LICENSE.md
|   mkdocs.yml
|   README.md
|   requirements.txt
|
+---docs
|   |   index.md
|   |
|   +---.icons
|   |   \---logos
|   |           iclogo.svg
|   |
|   \---assets
|           iclogo.png
|
+---FEM_Module
|       FEM_Module_class.py
|
\---notebooks
    |   Part_1.ipynb
    |   Part_2.ipynb
    |
    \---img
            Cantilever_diagram.png
            element_in_mesh.png
            element_types.png
            hat_functions_as_basis.png
            isoparametric_triangle.png
            mapping_to_isoparametric.png
            rectangularmesh.png
            rectangular_mesh_simple.png
            temp_result.png
            triangular_mesh.png

```

<!-- Change this to your License. Make sure you have added the file on GitHub -->

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
