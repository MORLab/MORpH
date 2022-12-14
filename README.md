# MORpH: Model reduction of linear port-Hamiltonian systems

## Description

**MORpH** is a MATLAB toolbox designed for
- efficient storage, 
- system analysis, 
- interconnection and 
- structure-preserving model reduction (MOR) 

of large-scale port-Hamiltonian (pH) models. The model class of pH systems enables energy-based modeling and a flexible coupling of models across different physical domains. This makes them well-suited for the simulation and control of complex technical systems. 
    
## Port-Hamiltonian systems

In **MORpH**, we work with pH models of the following form:

$$
\begin{align*}
Ex(t) &= (J-R)Qx(t) + (G-P)u(t), \\
y(t) &= (G+P)^TQx(t) + (S+N)u(t),
\end{align*}
$$

where $E,Q,J,R \in \mathbb{R}^{n \times n}$, $G,P \in \mathbb{R}^{n \times m}$, $S,N \in \mathbb{R}^{m \times m}$ fulfill the following constraints:

(i) The structure matrix 

$$ \Gamma := \left\lbrack \matrix{Q^T J Q & Q^T G \cr -G^T Q & N} \right\rbrack $$ 

is skew-symmetric, i.e., $\Gamma=-\Gamma^T$.

(ii) The dissipation matrix

$$ W := \left\lbrack \matrix{Q^T R Q & Q^T P \cr P^T Q & S} \right\rbrack $$ 

is positive semi-definite, i.e., $W=W^T\geq 0$.

(iii) The quadratic Hamiltonian is stated as:

$$ \mathcal{H}(x) = \frac{1}{2}x^T Q^T E x $$

with $Q^T E = E^T Q \geq 0$.

PH models are represented as objects of the *phs* class. They can be created in the following way:
```
sys = phs(J,R,Q,G,E,P,S,N);
```
Upon object creation, the *phs* class checks the pH structural constraints. The system matrices are saved in MATLAB's *sparse* format to reduce storage demands for large-scale systems. 

With the *phs* class comes a set of wrappers for system analysis. For example, you can create a Bode plot for your system via
```
bode(sys)
```

### Structure-preserving MOR

The goal of structure-preserving MOR is to approximate a large-scale pH model with a much smaller reduced model 

$$
\begin{align*}
\hat{E}\hat{x}(t) &= (\hat{J}-\hat{R})\hat{Q}\hat{x}(t) + (\hat{G}-\hat{P})u(t), \\
\hat{y}(t) &= (\hat{G}+\hat{P})^T\hat{Q}\hat{x}(t) + (\hat{S}+\hat{N})u(t),
\end{align*}
$$

with reduced state vector $\hat{x} \in \mathbb{R}^{r}$ such that (i) $r \ll n$, (ii) $\hat{y} \approx y$ for certain $u$ and (iii) the reduced model fulfills the pH structural constraints.

**MORpH** offers different algorithms to reduce large-scale pH models with an intuitive user interface.
For example, we can reduce our created pH model ```sys``` with the IRKA-PH algorithm (see [[1]](https://doi.org/10.1016/j.automatica.2012.05.052)) to a system with dimension ```redOrder``` via
```
redSys = irkaPH(sys, redOrder, Opts);
```
The ```Opts``` struct can be used to configure the algorithm's parameters.
An overview of the implemented algorithms is given [here](src/MOR/README.md). For a more detailed description of the different algorithms, pleases have a look at the [demos](/demos). 


## Installation

After downloading **MORpH**, change to the installation directory in MATLAB and run the script `setup_morph.m`.
It also assists with the installation of third-party software that may be required for some functionalities within **MORpH**.

## Getting started
You successfully installed **MORpH**? Great! The best way to get familiar with the toolbox is the [demo for creating pH models](demos/demo_phs_class.mlx). After you created your first pH model, you can analyze your model with the functions provided [here](src/@phs). Next, you can have a look [here](src/MOR/README.md) for an overview of the different MOR methods.

The toolbox is structured as follows:

**MORpH** (main folder)
- **src**
    - **@phs**: Class definition and functions to store and analyze large-scale, sparse pH models 
    - **MOR**: Algorithms for structure-preserving MOR of pH models
        - **@phsRed**: Class to store reduced pH models
        - **pa-ENF**: Passivity enforcement algorithms
        - **pa-MOR**: Passivity-preserving MOR algorithms
        - **pH-MOR**: PH-preserving MOR algorithms
    - **thirdParty**: Third-party code
    - **Utility**: Helper functions
- **demos**: Test systems and demos for different algorithms
- **test**: Unit tests for core functionalities

## Dependencies
Additional to the Control System Toolbox and Optimization Toolbox of MATLAB, **MORpH** partially relies on  (or provides interfaces to) the following third-party open-source software packages which are greatfully acknowledged:
- [Manopt](https://github.com/NicolasBoumal/manopt) for optimization on manifolds
- [M-M.E.S.S.](https://github.com/mpimd-csc/mmess) for solving matrix equations
- [GRANSO](https://gitlab.com/timmitchell/GRANSO) for optimization
- [CVX](https://github.com/cvxr/CVX) for optimization
- [YALMIP](https://github.com/yalmip/YALMIP) for optimization
- [SADPA/SAMDP](https://sites.google.com/site/rommes/software) for computing dominant spectral zeros
- [LINORM_SUBSP](http://www.math.tu-berlin.de/index.php?id=186267&L=1) for H-infinity norm computations
- [HINORM](http://www.mpi-magdeburg.mpg.de/mpcsc/software/infnorm) for H-infinity norm computations

> Note: If any of these software packages is used, please comply with the licensing and citation rules of each package! 

## Contributing
Pull requests are welcome! For major changes, please open an issue first to discuss what you would like to change.
For more information on how to contribute, please look [HERE](/CONTRIBUTING.md).

## Copyright
This toolbox is developed by [MORLab](https://www.epc.ed.tum.de/en/rt/research/model-order-reduction/), the model reduction lab at the [Chair of Automatic Control](https://www.epc.ed.tum.de/en/rt/home/) at [TUM](https://www.tum.de/en/). 
For more information, see [HERE](/LICENSE.md).
