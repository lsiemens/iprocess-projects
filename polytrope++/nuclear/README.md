# Nuclear #
Code and utilities for managing nuclear reactions and calculateing nuclear reaction rates. Nuclear reaction rates can be downloaded from [JINA Reaclib]("https://reaclib.jinaweb.org"). Each reaclib 1/2 file records metadata, the energy released in the reaction and parameters for a rate parameterization function [see](#reaclib_database). The rate parameterization function is

$$\lambda = exp\left[ a_0 + \sum_{i=1}^5 a_iT_9^{\frac{2i-5}{3}} + a_6\ln T_9 \right]$$

Where $T_9$ is the tempurature in Gigakelvin. Given a single nuclear reaction of the form $A(B, C)D$ the couple differential equations describing the reaction is

$$-\frac{1}{\mathcal{N_A}}\partial_tY_A = -\frac{1}{\mathcal{N_B}}\partial_tY_B = \frac{1}{\mathcal{N_C}}\partial_tY_C = \frac{1}{\mathcal{N_D}}\partial_tY_D = \frac{Y_A^{\mathcal{N_A}}Y_B^{\mathcal{N_B}}}{\mathcal{N_A}!\mathcal{N_B}!}\rho^\nu\lambda(T_9)$$

Where $rho$ is the density, $\mathcal{N_i}$ is the number of particle of isotope $i$, $Y_i$ is the molar abundance per gram of isotope $i$ and $\nu = \mathcal{N_A} + \mathcal{N_B} - 1$. Note for an isotope $i$ with molar mass $m_i$ and mass fraction $X_i$ its molar abundance is $Y_i = \frac{X_i}{m_i}$. For a $n$-ary reaction the units of $\lambda$ are $[cm^{3 \cdot n} s^{-1} mol^{-n}]$. The energy realeased / Q value of the reactions are given in $[MeV]$. For systems involving many nuclear reactions then $\frac{dY_i}{dt} = \left. \sum \partial_t Y_i \right|_{A(B, C)D}$

# To Do #
Calculate reverse reaction, Ex SkyNet: A Modular Nuclear Reaction Network Library

<a name="reaclib_database"></a>Cyburt, Richard H., et al. "The JINA REACLIB database: its recent updates and impact on type-I X-ray bursts." The Astrophysical Journal Supplement Series 189.1 (2010): 240.
