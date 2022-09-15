# Model Order Reduction (MOR)

The **MORpH** toolbox enables three different strategies to reduce large-scale pH models in a structure-preserving way.

## (1) PH-preserving MOR (pH-MOR)
The algorithms in the folder **pH-MOR** directly enforce the pH structure on the reduced-order model (ROM). 
The algorithms in this category follow this syntax:
```
ph_rom = <pHMOR_function>(ph_fom, redOrder, Opts);
```
This creates a reduced pH model ```ph_rom``` (represented as a *phsRed* object) with dimension ```redOrder``` for the full-order pH model ```ph_fom```. The ```Opts``` struct can be used to configure the algorithm's parameters.

## (2) Passivity-preserving MOR (pa-MOR) 
Since every pH model is inherently passive, one may also apply passivity-preserving MOR methods to pH models. Algorithms for this purpose are in the folder **pa-MOR**.
To use them, we first convert the *phs* object to a sparse *sss* object (the class *sss* is defined in the **sss** toolbox):

```
pa_fom = phs2sss(ph_fom);
```

Then we can use any function in the **pa-MOR** folder to reduce the passive ```pa_fom``` as follows:

```
pa_rom = <paMOR_function>(sss_fom, redOrder, Opts);
```
This creates a reduced, passive model ```pa_rom``` (represented as a *ss* object) with dimension ```redOrder```. Note that the pH form is generally lost in this process. However, since the reduced passive model ```pa_rom``` can be assumed to be minimal, it may be transformed back to pH form via
```
ph_rom = ss2phs(pa_rom);
```
The function ```ss2phs``` solves the Kalman-Yakubovich-Popov (KYP) inequality for this transformation.

## (3) Traditional MOR + Passivity enforcement (pa-ENF)

One can also use traditional, unstructured MOR methods to reduce pH models. The **MORpH** toolbox includes the **sssMOR** toolbox for this purpose. 
With **sssMOR**, you can reduce pH systems in the following way:

```
fom = phs2sss(ph_fom);
rom = <sssMOR_function>(fom, redOrder, Opts)
```

Note that the created ROM  ```rom``` is generally neither passive nor in pH form - its structure is destroyed! However, we can restore passivity of ```rom``` by applying passivity enforcement techniques which are located in the folder **pa-ENF**. They are called as follows:

```
pa_rom = <paENF_function>(rom, Opts);
```

Since the ROM ```pa_rom``` is passive, we can again transform it back to pH form via 

```
ph_rom = ss2phs(pa_rom);
```