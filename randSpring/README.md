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
To solve the equations of motion for vibrations in an arbitrary spring system, lets translate the equation from being expressed in terms of indicies to being expressed as a matrix equation. The sum $\sum_i K_{ik} \left( \vec{y}_i - \vec{y}_k \right)$ can not be directly interpreted as a matrix equation in terms of $K$, but it can be rewritten as $\sum_{i} K_{ik} \vec{y}_i - \left( \sum_i K_{ik} \right) \vec{y}_k$. These two terms can be identified as the diagonal and off diagonal components of a matrix. Denote this matrix as $K'$ with the elements (Note this is like a generalization or weighted version of the graph Laplacian),

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

## Error and stability
So far I have made a 1D linear simulation of the general problem, a 2D nonlinear simulation and I have not yet written the 2D linear simulation of vibrational solutions. In these two programs I have used forward Euler integration and semi-implicit Euler integration, where the nonlinear code is restricted to only using the forward Euler integrator. I started by looking into stability analysis of the linear system, but since this is a linear system the dynamics of the error is equivelent to that of the system as a whole. Since the matrix defining the linear system has purely imaginary eigenvalues, with some eigenvalues of zero, then $||1 + \Delta t \lambda|| \geq 1$ while the region of stability is $||1 + \Delta t \lambda|| \leq 1$. Since this system does not include any dampening to bring it into the region of stability, it can only be stable in the mathematical limit when $\Delta t \to 0$, the integrator is unstable for any non-zero $\Delta t$.

### Target System and Analytic Solution
Before we continue, let's define a general system to consider. Neglecting some of the spesific symmetry constraints for the spring systems that we have been considering leaves, a vector in $\mathbb{R}^n$ of posistions $\vec{x}$ and velocities $\vec{v}$ subject to the system of equations,

$$\begin{align}
\frac{d}{dt}\vec{x} & = \vec{v} \\
\frac{d}{dt}\vec{v} & = A\vec{x}
\end{align}$$

This can be represented as a linear system in $\mathbb{R^{2n}}$ space. Define a vector $\vec{w}$ where the first $n$ elements are $\vec{x}$ and the last $n$ elements are $\vec{v}$ and a matrix $B$ to reproduce the linear system. The resulting system is,

$$\vec{w} =
\begin{bmatrix}
\vec{x} \\
\vec{v}
\end{bmatrix}$$
$$B = 
\begin{bmatrix}
0 & I \\
A & 0
\end{bmatrix}
$$
$$\frac{d}{dt}\vec{w} = B\vec{w}$$

The solution to this system, for some inital values $\vec{w}_0$, is the vector $\vec{w}_f$ given by the matrix exponential,

$$\vec{w}_f = e^{tB}\vec{w}_0 = \left[ \sum_{k=0}^\infty \frac{(tB)^k}{k!} \right] \vec{w}_0$$

Provided that $B$ is diagonalizable, the system reduces to $2n$ decoupled equations. Lableing the eigenvalues and eigenvectors as $\lambda_k$ and $\vec{w}_{\lambda_k}$, with $\vec{w}_0 = \sum_{k=1}^{2n} a_k \vec{w}_{\lambda_k}$ and $|\vec{w}_{\lambda_k}| = 1$, then the solution to the system is

$$\vec{w}_f = \sum_{k=1}^{2n} a_k e^{t\lambda_k} \vec{w}_{\lambda_k}$$

### Forward Euler
Discritize the system such that there are $N$ intemederiary steps from $\vec{w}_0$ at $t_0=0$ to $\vec{w}_f$ at $t_f=t$. The intermediary vectors are labeld as $\vec{w}_i$ at $t_i = i \Delta t$, where $\Delta t = \frac{t}{N}$. Using a forward Euler integrator the system becomes,

$$\vec{w}_{i + 1} = \vec{w}_i + \Delta t B \vec{w}_i = \left( I + \Delta t B \right) \vec{w}_i$$

Colapsing this recursive equation gives the approximate solution for $\vec{w}_f$,

$$\vec{w}_N = \left( I + \Delta t B \right)^{N} \vec{w}_0 = \left( I + \frac{tB}{N} \right)^{N} \vec{w}_0$$

Note, in the limit where $N \to \infty$ this reproduces the analytic solution,

$$\begin{align}
\lim_{N \to \infty} \left( I + \frac{tB}{N} \right)^{N} \vec{w}_0 & = \lim_{N \to \infty} \sum_{k = 0}^N \binom{N}{k}I^{N - k}\left( \frac{tB}{N} \right)^k \vec{w}_0 \\
& = \lim_{N \to \infty} \sum_{k = 0}^N \frac{N!}{(N - k)!N^k}\frac{(tB)^k}{k!} \vec{w}_0 \\
& = \left[ \lim_{N \to \infty} I + tB + (1 - \frac{1}{N}) \frac{(tB)^2}{2!}+ (1 - \frac{1}{N})(1 - \frac{2}{N}) \frac{(tB)^3}{3!} + \cdots \right] \vec{w}_0 \\
& = \sum_{k=0}^\infty \frac{(tB)^k}{k!} \vec{w}_0 \\
& = e^{tB} \vec{w}_0
\end{align}$$

Define the error $\varepsilon_f$ as the difference between the discritized and the analytic solutions, so

$$\varepsilon_f = \vec{w}_N - \vec{w}_f = \left(I + \Delta t B\right)^N \vec{w}_0 - e^{tB} \vec{w}_0 = \left[ \left( I + \Delta t B \right)^N - e^{tB} \right] \vec{w}_0$$

Lets solve for this error as an asymptotic series of the number of intemediary steps $N$. Make the substitution $u = \frac{1}{N}$ and expand $\varepsilon_f$ about $u = 0$

$$\begin{align}
\varepsilon_f & = \left[ \left( I + \Delta t B \right)^N - e^{tB} \right] \vec{w}_0 \\
& = \left[ \left( I + \frac{tB}{N} \right)^N - e^{tB} \right] \vec{w}_0 \\
& = \left[ \left( I + utB \right)^{1/u} - e^{tB} \right] \vec{w}_0 \\
& = \left[ e^{\frac{1}{u}\ln(I + utB}) - e^{tB} \right] \vec{w}_0 \\
& = \left[ e^{\frac{1}{u}\left( utB - \frac{1}{2}(utB)^2 + \frac{1}{3}(utB)^3 + \mathcal{O}(u^4) \right)} - e^{tB} \right] \vec{w}_0 \\
& = \left[ e^{tB - \frac{u}{2}(tB)^2 + \frac{u^2}{3}(tB)^3 + \mathcal{O}(u^3)} - e^{tB} \right] \vec{w}_0 \\
& = \left[ e^{-\frac{u}{2}(tB)^2 + \frac{u^2}{3}(tB)^3 + \mathcal{O}(u^3)} - I \right] e^{tB} \vec{w}_0 \\
& = \left[ I - \frac{u}{2}(tB)^2 + \frac{u^2}{3}(tB)^3 + \frac{u^2}{8}(tB)^4 + \mathcal{O}(u^3) - I \right] \vec{w}_f \\
& = \left[ -\frac{(tB)^2}{2N} + \frac{(tB)^3}{24N^2}(8 + 3tB) + \mathcal{O}(1/N^3) \right] \vec{w}_f
\end{align}$$

Given the inital value is a unit eigenvector, $\vec{w}_0 = \vec{w}_{\lambda_k}$, then to first order the magnitude of the error will be,

$$||\varepsilon_f|| = \frac{t^2 \left( \Re(\lambda_k)^2 + \Im(\lambda_k)^2 \right)}{2N} e^{t\Re(\lambda_k)}$$

Applying this to the mass spring systems where $\lambda_k = \pm i 2 \pi f$ the number of steps needed keep the final relative error to $||\varepsilon_f||$ is,

$$N = \frac{2\left( \pi tf \right)^2}{||\varepsilon_f||}$$

