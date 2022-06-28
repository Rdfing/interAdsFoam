# interAdsFoam [![Github All Releases](https://img.shields.io/github/downloads/Rdfing/interAdsFoam/total.svg)]()

Multi-phase CFD solver for chemical transport and adsorption with dynamic mass transfer zone induced by varying fluid interface. The solver is developed to simulate system where bulk fluid flow and porous media flow both present.

## Description
This solver is developed to simulate turbulence mixing and chemical adsorption in urban water treatment systems. The solver integrates VOF with chemical adsorption.
- Fluid interface capturing by VOF (inherited from interFoam)
- Darcy and Forchheimer momentum sink for reactor
- Turbulence modeling is suppressed in reactor (granular Reynold number too small, laminar flow)
- Non-equilibrium (Langmuir-like) adsorption model 
- Hydrodynamic dispersion coefficient tensor can be defined by longitudinal and transverse dispersivity
- First-order time splitting for transport and reaction (for now)
- Stiff ODE solver for reaction
- Examples for inverse modeling of adsorption and headloss are provided in the reference

Here is an example 
![Alt Text](https://github.com/Rdfing/interAdsFoam/blob/main/example.gif)

## Prerequisites
OpenFOAM-v1912 

## Authors
Haochen Li, PhD <br />
Department of Environmental Engineering Sciences <br />
University of Florida

## Reference
Li, H., Sansalone, J., 2022. InterAdsFoam: An Open-Source CFD Model for Granular Media–Adsorption Systems with Dynamic Reaction Zones Subject to Uncontrolled Urban Water Fluxes. J. Environ. Eng. 148, 04022049. https://doi.org/10.1061/(ASCE)EE.1943-7870.0002027

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD®  trade marks.
