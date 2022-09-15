# MORpH: Model reduction of linear port-Hamiltonian systems

## Description

**MORpH** is a MATLAB toolbox designed for
- efficient storage, 
- system analysis, 
- interconnection and 
- structure-preserving model reduction (MOR) 

of large-scale port-Hamiltonian (pH) models. The model class of pH systems enables energy-based modeling and a flexible coupling of models across different physical domains. This makes them well-suited for the simulation and control of complex technical systems. 
    
### The *phs* class

In **MORpH**, pH models are represented as objects of the *phs* class. They can be created in the following way:
```
sys = phs(J,R,Q,G,E,P,S,N);
```
Upon object creation, the *phs* class checks the pH structural constraints. The system matrices are saved in MATLAB's *sparse* format to reduce storage demands for large-scale systems. 

With the *phs* class comes a set of wrappers for system analysis. For example, you can create a Bode plot for your system via
```
bode(sys)
```

### Structure-preserving MOR

**MORpH** offers different ways to reduce large-scale pH models in a structure-preserving way, meaning that the reduced model also has pH form.
For example, we can reduce our created pH model ```sys``` with the IRKA-PH algorithm (see [[1]](https://doi.org/10.1016/j.automatica.2012.05.052)) to a system with dimension ```redOrder``` via
```
redSys = irkaPH(sys, redOrder, Opts);
```
The ```Opts``` struct can be used to configure the algorithm's parameters.
An overview of the implemented algorithms is given [HERE](src/MOR/README.md). For a more detailed description of the different algorithms, pleases have a look at the [DEMO](/demos) files. 


## Installation

After downloading **MORpH**, change to the installation directory in MATLAB and run the script `setup_morph.m`.
It also assists with the installation of third-party software that may be required for some functionalities within **MORpH**.

### Toolbox structure
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