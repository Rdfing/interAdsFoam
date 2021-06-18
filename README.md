# interAdsFoam
Multi-phase CFD solver for chemical transport and adsorption with dynamic mass transfer zone induced by varying fluid interface

## Description
This solver is developed to simulate turbulence mixing and chemical adsorption in urban water treatment systems. The solver integrates VOF with chemical adsorption.
- Fluid interface capture inherited from interFoam
- Darcy and Forchheimer momentum sink for reactor
- Turbulence modeling is suppressed in reactor (granular Reynold number too small, laminar flow)
- Non-equilibrium (Langmuir-like) adsorption model 
- Hydrodynamic dispersion coefficient tensor can be defined by longitudinal and transverse dispersivity
- First-oder time splitting for transport and reaction (for now)
- Stiff ODE solver for reaction
- Example of inverse modeling of adsorption and headloss are provided in the reference

Here is an example 
![Alt Text](https://github.com/Rdfing/interAdsFoam/blob/main/example.gif)

## Prerequisites
OpenFOAM-v1912 

## Authors
Haochen Li, PhD <br />
Department of Environmental Engineering Sciences <br />
University of Florida

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD®  trade marks.

## Reference
Li, H., & Sansalone, J. (under review). CFD modeling of adsorption reactor with dynamic reaction zone under uncontrolled environmental loading. 
