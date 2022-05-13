# Polytrope #
Code for finding solutions to the Lane-Emden equation which are known as a polytrope. The four equations of stellar ($P$, $M$, $L$, $T$) structure can be reduced down to the two equations for $P$ and $M$ by enforcing a relationship between the density and known varaibles $r$, $P$ and $M$. Using the equation $P = K\rho^{\frac{n + 1}{n}}$, where $K$ is a constant and $n$ is the polytropic index, leads to the polytrope solution described by the Lane-Emden equation.


## Lane-Emden Equation ##
Solve the stellar structur equations,
$$\frac{dP(r)}{dr} = -\frac{GM(r)\rho(r)}{r^2}$$
$$\frac{dM(r)}{dr} = 4\pi r^2 \rho(r)$$
given the equation $P = K\rho^{\frac{n + 1}{n}}$ where $K$ is a constant and $n$ is the polytropic index. In this case the equations do not have an explicit tempurature dependence and can be conbined to the form,

$$\frac{1}{r^2}\frac{d}{dr}\left( \frac{r^2}{\rho(r)}\frac{dP(r)}{dr} \right) = -4\pi G\rho(r)$$

Written in dimensionless form with $\rho(r) = \rho_c \theta^n(r)$, $P(r) = P_c \theta^{n + 1}(r)$, $\alpha^2 = \frac{(n + 1)P_c}{4\pi G\rho_c^2}$, $\xi = \frac{r}{\alpha}$, where $P_c = K\rho_c^{\frac{n+1}{n}}$, then the result is the dimensionless Lane-Emden equation

$$\frac{1}{\xi^2}\frac{d}{d\xi}\left( \xi^2\frac{d\theta(\xi)}{d\xi} \right) = -\theta^n(\xi)$$

Solutions should satisfy the boundary conditions $\left. \frac{d\theta(\xi)}{\xi} \right|_{\xi = 0} = 0$ and $\left. \theta(\xi) \right|_{\xi = 0} = 1$

In this model a star has radius $R_* = \alpha \xi_0$ and mass $M_* = 4\pi \alpha^3 \rho_c \int_0^{\xi_0} \xi^2 \theta^n(\xi)d\xi$, where $\xi_0$ is the first zero of the function $\theta(\xi)$. Notes from LeBlanc 5.4


# Mass Coordinate #
Changeing the free varaible from $r$ to $M$ in the structure equations leads to the pair of equations,

$$\frac{dP(M)}{dM} = -\frac{GM}{4\pi r^4(M)}$$

$$\frac{dr(M)}{dM} = \frac{1}{4\pi r^2(M) \rho(M)}$$

given the equation $P = K\rho^{\frac{n + 1}{n}}$ where $K$ is a constant and $n$ is the polytropic index. Since the goal is to a system of first order equation instead of trying to combine these equations let us directly insert the equation for $P$ and change to dimensionless varaibles.

using the variables $\rho = \rho_c \theta^n$, $r = R_\odot \chi$ and $M = M_\odot \phi$ and let $P_c = k\rho_c^{\frac{n + 1}{n}}$, note $P = P_c\theta^{n + 1}$, then the equations become

$$P_c\frac{d\theta^{n+1}}{d\phi} = (n + 1)P_c\theta^n \frac{d\theta}{d\phi} = -\frac{G M_\odot^2}{4\pi R_\odot^4} \frac{\phi}{\chi^4}$$

$$\frac{d\chi}{d\phi} = \frac{M_\odot}{4\pi R_\odot^3 \rho_c}\frac{1}{\chi^2\theta^n}$$

Finaly simplifying the equations yields the set of first order equations

$$\frac{d\theta}{d\phi} = -\frac{\alpha}{n + 1}\frac{\phi}{\chi^4\theta^n}$$

$$\frac{d\chi}{d\phi} = \beta\frac{1}{\chi^2\theta^n}$$

where $\alpha = \frac{GM_\odot^2}{4\pi R_\odot^4P_c}$ and $\beta = \frac{M_\odot}{4\pi R_\odot^3 \rho_c}$
