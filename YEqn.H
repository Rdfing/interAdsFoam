{
    /*************Step 1 solve the ODE over deltaT*********************/
    Info << "Step 1 ODE: Solving for reaction" << endl;
    ODEdeltaT = deltaT.value();

    Info << "ODEdeltaT = "<< ODEdeltaT << endl;

    #include "solveAdsorption.H"

    forAll(dYdt,celli)
    {
        dYdt[celli] = (Y_ODE[celli]-Y[celli])/ODEdeltaT;
    }
    //dYdt.write();

    /*******Step 2 solve the transport equation over deltaT************/
    Info << "Step 2 PDE: Solving for transport" << endl;

    // Calculate the diffusivity
    Dt = turbulence->nut()/Sct;
    Dt.correctBoundaryConditions();
    volScalarField Deff("Deff", turbulence->nu()/Sc + Dt); //Diffusion at the bulk fluid region

    // construct the flux
    Info << "Executing the flux correction" << endl;

    volScalarField alphaPorosity ("alphaPorosity", alpha1*porosityField);

    const volScalarField& alpha = alphaPorosity; // consider the influence of the porosity

    // Make a crude guess of the phase flux using default interpolation
    surfaceScalarField alphaPhi
    (
        "alphaPhi",
        phi*fvc::interpolate(alpha)
    );

    /**********************Get the potential field*********************/
    autoPtr<volScalarField> PhiPtr_;

    const volScalarField& p_ = p_rgh;

    wordList PhiPatchFieldTypes(mesh.boundaryMesh().size());
    
    forAll(p_.boundaryField(), patchi)
    {
        PhiPatchFieldTypes[patchi] =
            p_.boundaryField()[patchi].fixesValue()
        ? fixedValueFvPatchField<scalar>::typeName
        : zeroGradientFvPatchField<scalar>::typeName;
    }

    PhiPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "PhiY",
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(phi.dimensions()/dimLength, Zero),
            PhiPatchFieldTypes
        )
    );
    mesh.setFluxRequired(PhiPtr_->name());
    const volScalarField& Phi = PhiPtr_;

    /******************************************************************/
    // Solve for the potential and correct alphaPhi with the resulting flux
    // Lookup the non-orthogonal solution control!
    // note that the nonOrthogonalSolutionControl.H is not implement in v1912
    // So I use the pimple NonOrthogonal control since the model has to correct the NonOrthogonal anyhow in other equation.
    // It must be notice that this NonOrthogonal control is for solving the laplacian equation, which has nothing to do the mesh orthanonality.
    // Even for the perfect orthogonal mesh, this computation is still needed!!!!
    if (phi.dimensions() == dimVolume/dimTime)
    {
        while (pimple.correctNonOrthogonal()) 
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(Phi)
              + fvc::ddt(alpha)
              + fvc::div(alphaPhi)
            );

            PhiEqn.solve(); //this equation solve for Phi!!!

            if (pimple.finalNonOrthogonalIter())
            {
                alphaPhi += PhiEqn.flux();
            }
        }
    }
    else
    {
        SeriousErrorInFunction();
    }
    /******************************************************************/

    // Updates the hydrodynamic dispersion tensor
    dimensionedScalar RTSMALL("",dimLength/dimTime,ROOTSMALL); //prevent for small velocity 
    volTensorField Dhh("Dhh", (D*tau_p+alpha_T*mag(U)/porosity)*tensor::I+(alpha_L-alpha_T)*(U*U)/(porosity*mag(U)+RTSMALL));

    fvScalarMatrix alphaYEqn
    (
          fvm::ddt(alpha, Y)
        + fvm::div(alphaPhi, Y)
        - fvm::laplacian(fvc::interpolate(alpha)*fvc::interpolate(Dhh*iD*porosity), Y)
        - fvm::laplacian(fvc::interpolate(alpha)*fvc::interpolate(Deff*(1.0-iD)), Y)
        ==
          dYdt*iD*alpha
        + fvOptions(Y)
        - fvm::ddt(residualAlpha_, Y)
        + fvc::ddt(residualAlpha_, Y)
    );

    alphaYEqn.relax();
    fvOptions.constrain(alphaYEqn);
    alphaYEqn.solve();
    fvOptions.correct(Y);
  
  
    Y.max(0.0); //clip off the negtive values

    // update the results
    alphaY = alpha*Y/porosityField; // should use this guys as Y.

}



