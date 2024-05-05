# Vibrations in Random Spring Mass Systems
Recently I have been introduced to the theory of random matrices. A friend wrote his honours thesis on using random matrix theory to analyze the energy spectrum of a black hole in an Anti-de Sitter space. I would like to learn more about the theory of random matrices, but in a simpler context. I settled on the dynamics of systems of masses and springs, since this provides a simple context in which random matrices can be investigated. Additionally over this summer I need to brush up on my C++ skills for grad school, so I will aim to write simulations in C++ of the vibrations in random networks of masses and springs. I will likely be using BLAS and LAPACK in the code.

## Theory: Spring Mass System
Consider a system with $n$ masses and springs that can connect any pair of masses together, so there is at most $\frac{n(n - 1)}{2}$ springs. Label the masses with an integer $i \in [1, \cdots n]$. For the $i$th mass, denote its position as $\vec{x}_i$ and its mass as $m_i$. The velocity is $\vec{v}_i = \dot{\vec{x}}_i$.

### Springs
Springs connect pairs of masses and so can be labelled by a pair of indices, $(i, j)$, indicating the masses connected by the spring. Note the indices $(i, j)$ and $(j, i)$ refer to the same spring.

The spring constant for any spring, $(i, j)$ in the system is denoted as $K_{ij}$. If two masses $i$ and $j$ are not connected by a spring this this is equivalent to being connected by a spring with a spring constant of zero, $K_{ij} = 0$. Also, note $K_{ij} = K_{ji}$ since both refer to the same spring. It is nonsensical for a spring to connect a mass to itself so let $K_{ii} = 0$. The spring constant $K_{ij}$ is defined for all values of $i$ and $j$, and so the set of all spring constants define a square matrix. Note that the matrix $K$ with elements $K_{ij}$ is a non-negative symmetric matrix with zeros on the diagonal.

Each spring has a natural length, $L_{ij}$, corresponding to the distance between the masses if the system was in a steady state and all other spring constants where zero. Using similar arguments as with the spring constant $L_{ij} = L_{ji}$, and $L_{ii} = 0$. Like with the spring constants, these lengths determine a matrix $L$ with elements $L_{ij}$ and is a non-negative symmetric matrix with zeros on the diagonal.

Ignoring the physics for a moment. The matrix $K$ defines the network of springs between the masses which forms a graph, in the mathematical sense. The Adjacency matrix for this graph is the matrix $A$,

$$A_{ij} = 
\begin{cases}
1 & K_{ij} \neq 0 \\
0 & K_{ij} = 0
\end{cases}$$

The matrix $K$ can be thought of as the weights of the graph. Since $K$ and $A$ are symmetric the graph is undirected. Note, isolated subgraphs in the network can be identified from the matrix $(I + A)^n$. Any two masses $i$ and $j$ are isolated from each other if the element ${(I + A)^n}_{ij} = 0$.

## Lagrangian
Each mass has kinetic energy $T_i = \frac{1}{2}m_i \vec{v}_i \cdot \vec{v}_i$ and each spring has potential energy $V_{ij} = \frac{1}{2} K_{ij} \left( \left| \vec{x}_i - \vec{x}_j \right| - L_{ij} \right)^2$. The Lagrangian for the system is then $\mathcal{L} = \sum_i \left(T_i - \sum_j^{j < i} V_{ij}\right)$. Adding a factor to account for double counting of the potential energy of the springs, the Lagrangian can be given as,

$$\mathcal{L} = \sum_i \frac{1}{2}m_i \vec{v}_i \cdot \vec{v}_i - \sum_{i, j}\frac{1}{4} K_{ij} \left( \left| \vec{x}_i - \vec{x}_j \right| - L_{ij} \right)^2$$

Computing terms for the Euler-Lagrange equations. Starting with the kinetic term,

$$\frac{d}{dt} \frac{\partial \mathcal{L}}{\partial \vec{v}_k} = \frac{d}{dt}\sum_i \frac{1}{2}m_i \frac{\partial}{\partial \vec{v}_k} \left( \vec{v}_i \cdot \vec{v}_i \right) + 0 = \frac{d}{dt} m_k\vec{v}_k = m_k \dot{\vec{v}}_k$$

Moving on to the potential term,

$$\frac{\partial \mathcal{L}}{\partial \vec{x}_k} = 0 - \sum_{i, j}\frac{1}{4} K_{ij} \frac{\partial}{\partial \vec{x}_k} \left( \left| \vec{x}_i - \vec{x}_j \right| - L_{ij} \right)^2 = -\sum_{i, j}\frac{1}{2} K_{ij} \left( \left| \vec{x}_i - \vec{x}_j \right| - L_{ij} \right) \frac{\vec{x}_i - \vec{x}_j}{\left| \vec{x}_i - \vec{x}_j \right|} \frac{\partial}{\partial \vec{x}_k} \left( \vec{x}_i - \vec{x}_j \right)$$

After evaluating this last partial derivative the sum simplifies to,

$$\frac{\partial \mathcal{L}}{\partial \vec{x}_k} = \sum_{i} K_{ik} \left( \left| \vec{x}_i - \vec{x}_k \right| - L_{ik} \right) \frac{\vec{x}_i - \vec{x}_k}{\left| \vec{x}_i - \vec{x}_k \right|}$$

Using the Euler-Lagrange equations, the equation of motion for the $k$th mass in the system is,

$$m_k \dot{\vec{v}}_k - \sum_{i} K_{ik} \left( \left| \vec{x}_i - \vec{x}_k \right| - L_{ik} \right) \frac{\vec{x}_i - \vec{x}_k}{\left| \vec{x}_i - \vec{x}_k \right|} = 0$$

If system exists in one dimension, where $\vec{x}_i$ is simply a number, then the equations of motion for this system is a linear equation, provided that no two particles ever have the same position. This will be true for perturbations of the system when near equilibrium, so vibrations in $1$-D are described by this matrix equation. However, if the space the particles are in is not one dimensional, then this equation is non-linear.

## Vibrations
Lets consider vibrations in the mass spring system. Given $\vec{x_0}_i$ is an equilibrium solution to the equations of motion, then denote the perturbed system as $\vec{x}_i = \vec{x_0}_i + \epsilon \vec{y}_i$. We will expand the equations of motion in terms of $\epsilon$ while neglecting terms of order two or higher. Expanding $\left| \vec{x}_i - \vec{x}_k \right|$ gives, 

$$\left| \vec{x}_i - \vec{x}_k \right| = \left| \vec{x_0}_i - \vec{x_0}_k \right| + \epsilon \frac{\vec{x_0}_i - \vec{x_0}_k}{\left| \vec{x_0}_i - \vec{x_0}_k \right|} \cdot \left(\vec{y}_i - \vec{y}_k\right) + \mathcal{O}(\epsilon^2)$$

Then the expansion of the term $\frac{\vec{x}_i - \vec{x}_k}{\left| \vec{x}_i - \vec{x}_k \right|}$, is,

$$\frac{\vec{x}_i - \vec{x}_k}{\left| \vec{x}_i - \vec{x}_k \right|} = \frac{\left( \vec{x_0}_i - \vec{x_0}_k \right) + \epsilon \left(\vec{y}_i - \vec{y}_k\right)}{\left| \vec{x_0}_i - \vec{x_0}_k \right| + \epsilon \frac{\vec{x_0}_i - \vec{x_0}_k}{\left| \vec{x_0}_i - \vec{x_0}_k \right|} \cdot \left(\vec{y}_i - \vec{y}_k\right) + \mathcal{O}(\epsilon^2)} = \frac{\vec{x_0}_i - \vec{x_0}_k}{\left| \vec{x_0}_i - \vec{x_0}_k \right|} + \mathcal{O}(\epsilon^2)$$

So, expanding the full equations of motion give,

$$m_k \left( \ddot{\vec{x_0}}_k + \epsilon \ddot{\vec{y}}_k \right) - \sum_{i} K_{ik} \left( \left(\vec{x_0}_i - \vec{x_0}_k \right) - L_{ik} \frac{\vec{x_0}_i - \vec{x_0}_k}{\left| \vec{x_0}_i - \vec{x_0}_k \right|} + \epsilon \left(\vec{y}_i - \vec{y}_k\right) \right) + \mathcal{O}(\epsilon^2)= 0$$

The dynamical equation for the perturbation of the system is,

$$m_k \ddot{\vec{y}}_k - \sum_i K_{ik} \left( \vec{y}_i - \vec{y}_k \right) = 0$$

So the perturbations from equilibrium, i.e. vibrations, are governed by a simple linear matrix equation. Also, note that full equations of motion in terms of $\vec{x}_i$ has the same form as this equation in the case where $L_{ij} = 0$. The equations of motion are the same for as spring system where each spring has a natural length of zero and the vibrations of an arbitrary spring system.

## Vibration Solution
To solve the equations of motion for vibrations in an arbitrary spring system, lets translate the equation from being expressed in terms of indicies to being expressed as a matrix equation. The sum $\sum_i K_{ik} \left( \vec{y}_i - \vec{y}_k \right)$ can not be directly interpreted as a matrix equation in terms of $K$, but it can be rewritten as $\sum_{i} K_{ik} \vec{y}_i - \left( \sum_i K_{ik} \right) \vec{y}_k$. These two terms can be identified as the diagonal and off diagonal components of a matrix. Denote this matrix as $K'$ with the elements,

$$K'_{ij} =
\begin{cases}
\sum_k K_{ik} & i = j \\
-K_{ij} & i \neq j
\end{cases}$$

Also, define the vector $\vec{Y}$ such its elements are the vectors $\vec{y}_i$ and the diagonal matrix $M$, where the elements on the diagonal are $m_i$. Then the vibrational equation of motion can be expressed as the matrix equation,

$$M \ddot{\vec{Y}} + K' \vec{Y} = 0$$

Finaly, rearanging terms puts the equation in a form analogus to the equation $a = -\frac{k}{m}x$, for a simple mass on a spring,

$$\ddot{\vec{Y}} = -M^{-1} K' \vec{Y}$$

The matrix $M^{-1} K'$ is symmetric, so it is diagonalizable. There are two matrices $\pm \Omega$ such that ${\Omega}^2 = M^{-1} K'$, which are trivial compute when in diagonal form. The solutions to the vibrational equation of motion are given by $\vec{Y}(t) = e^{\pm i \Omega t}\vec{A}$. The general solution is then,

$$\vec{Y}(t) = e^{i \Omega t} \vec{A} + e^{-i \Omega t} \vec{B}$$

Restricting to real solutions and using initial conditions,

$$\vec{Y}(t) = \frac{e^{i \Omega t} + e^{-i \Omega t}}{2} \vec{Y}_0 + \frac{e^{i \Omega t} - e^{-i \Omega t}}{2i} \Omega^{-1} \dot{\vec{Y}}_0 = \cos(\Omega t) \vec{Y}_0 + \Omega^{-1} \sin(\Omega t) \dot{\vec{Y}}_0$$

where $\vec{Y}_0 = \vec{Y}(0)$ and $\dot{\vec{Y}}_0 = \dot{\vec{Y}}(0)$.

