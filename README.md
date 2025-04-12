# Numerical Simulation using Finite Element Method 

## Description

Finite Element Method (FEM) is a powerful computational method used across a wide variety of disciplines, including engineering, mathematics and earth science. It is used to solve complex differential equations and simulate real-world scenarios. This project focuses on application of FEM in solid mechanics, but the teachings can be extended to other areas of study.
It provides a comprehensive, step-by-step tutorial designed to introduce users to the fundamentals of FEM and its practical implementation in Python. The tutorial takes users through all the key stages of numerical simulations, starting from domain meshing and problem discretisation, moving to solving the differential equations and finally analysing and understanding the results.  
To make the learning process intuitive, the project works through two examples of increasing difficulty. The first example is straightforward and demonstrates the methodology and code implementation. As users progress, they are gradually introduced to more advanced aspects of FEM, ensuring a solid foundation and deeper understanding of this method.

## Learning Outcomes

- Develop understanding of the fundaments of FEM
- Learn how to create custom FEM code and adapt it to solve partial differential equations
- Gain familiarity with existing libraries and tools available for FEM implementation
- Develop skills to interpret solutions and effectively analyse the simulation results

In this project, the exercises are integrated with the theoretical components.

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

## Requirements

#### Academic
To successfully complete this examplar the user should have basic understanding of key mathematical concepts, including linear algebra (such as matrices, vectors, and determinants) and differential equations (both ordinary and partial differential equations, linear system of equations)

Familiarity with Python programming, using existing libraries, writing classes, functions and graph plotting.


### System

System requirements:
- Python 3.11 or newer
- Paraview 4.0 or newer
- Jupyter Notebooks

Python libraries used:
- pygmsh
- numpy
- matplotlib
  
  
  
## Getting Started

### Finite Element Method Overview
A significant part of our understanding of the real world comes from differential equations, which describe the behaviour of various quantities. These equations are fundamental across many scientific fields, including physics, geoscience, and medicine, as they explain how systems behave under specific conditions, how they respond to external forces and change with time.

Solving differential equations allows us to predict and analyse real world phenomena. In many cases, computer simulations are used to solve these equations for specific scenarios. Engineers use computer simulations to optimize product designs, ensuring they can withstand real world conditions, while geoscientists use them to study complex natural processes or our interactions with them, such as energy storage. These simulations enable us to test hypotheses and evaluate outcomes without the need for expensive real life experiments. However, their accuracy and reliability depend on our ability to effectively solve the underlying systems of differential equations.

In certain simplified cases, differential equations can be solved analytically, providing exact solutions. However, these cases are often over simplified and may not be flexible enough for practical applications. When analytical solutions are not feasible, we use numerical methods to approximate solutions. A wide range of numerical techniques exists, each designed to solve specific problems and improve the accuracy and efficiency of the solutions.

One of the most widely used numerical methods is the Finite Element Method (FEM). It is particularly popular in mechanics, including structural and rock mechanics, fluid dynamics, and heat transfer. Since its introduction in the 1950s, FEM has become an essential tool in numerical modeling, with extensive research dedicated to refining and advancing its capabilities.  

The core principle of FEM is discretising the problem domain into smaller interconnected elements. This process, known as discretisation, allows the problem to be described at the element level. Each element has an associated number of nodes, which are shared by adjacent elements. These nodes form relationships between elements that can be aggregated to a global system of equations. Solving this system yields an approximate solution at the nodes across the domain. Interpolation then allows us to calculate the solution at any point.

Therefore the key stages of Finite Element Methods are:

- Domain discretisation  -  dividing the domain into finite elements
- Element level problem formulation - defining the equations on each element using basis functions
- Assembly to global system - combining the element equations into a global system
- Application of boundary conditions - applying constrains and conditions to the problem to solve for unique solutions.
- Solving linear system of equations - using a numerical solver to find the solution
- Post processing of solution - interpreting, visualising and analysing the solution

### Project Overview
The project is split into two main examples. Inside each examples there are smaller exercises that allow the user to apply and practice the theoretical components presented. For some exercises the user is expected to complete the code. The missing parts are indicated with `...` and the solutions are provided in the hidden cells marked "Solution" in red. By double-clicking on the hidden cells, the user can reference solutiona and copy and past into the code cells.

[This exercise](notebooks/Part_1.ipynb) will guide you through the key stages of constructing a Finite Element code. It will focus on a simple problem of heat transfer to build the basics of the method. At the end of this exercise you should have a good understanding of how FEM works.

[The next exercise](notebooks/Part_2.ipynb) will build on these concepts and look at solving a slightly more complex problem: an engineering problem where we solve for displacement vector. In this exercise we focus on applying different type of boundary conditions and analysing solution convergence to analytical solution.


## Project Structure

```
│   .gitignore
│   LICENSE.md
│   mkdocs.yml
│   README.md
│   requirements.txt
│
├───.devcontainer
│       devcontainer.json
│       Dockerfile
│
├───.github
│   └───workflows
│           docs.yml
│           link_checker.yml
│
├───.ipynb_checkpoints
├───docs
│   │   index.md
│   │
│   ├───.icons
│   │   └───logos
│   │           iclogo.svg
│   │
│   └───assets
│           iclogo.png
│
├───FEM_Module
│   │   FEM_Module_class.py
│   │   Run_Cantilever_Example.py
│   │   Support_functions.py
│   │
│   └───__pycache__
│           FEM_Module_class.cpython-311.pyc
│           Run_Cantilever_Example.cpython-311.pyc
│           Support_functions.cpython-311.pyc
│
└───notebooks
    │   .placeholder
    │   Part_1.ipynb
    │   Part_2.ipynb
    │
    ├───.ipynb_checkpoints
    │       Part_1-checkpoint.ipynb
    │       Part_2-checkpoint.ipynb
    │       Part_2_withlinear-checkpoint.ipynb
    │
    └───img
            Cantilever_diagram.png
            cantilever_result.png
            element_in_mesh.png
            element_types.png
            hat_functions_as_basis.png
            isoparametric_triangle.png
            line_element_diag.png
            mapping_to_isoparametric.png
            rectangularmesh.png
            rectangular_mesh_simple.png
            temp_result.png
            triangular_mesh.png

```

## License

This project is licensed under the [BSD-3-Clause license](LICENSE.md)
