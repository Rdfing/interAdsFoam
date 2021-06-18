//solve the ODE over ODEdeltaT

// update the varible
q_ODE = q; //current q
Y_ODE = Y; //current Y

// Start time and end time
scalar xStart = 0;
scalar xEnd = xStart + ODEdeltaT;

// loop over all the cells
//Info << "ODE: Solving for reaction " << endl;
forAll(Y_ODE,celli)
{
    if (iD[celli] == 1 && alpha1[celli] > 0.99) // only solve for the porous media and water phase
    {

        // ****************invoke the ODE solver********************
        yStart[0] = Y[celli]; //current Y
        yStart[1] = q[celli]; //current q
        scalarField y(yStart);

        // Initial ODE step size (does not affect the final solution)
        // but will affect the speed to some extend
        scalar dxEst =  1e-7;

        //- Solve the ODE system from xStart to xEnd, update the state
        //  and return an estimate for the next step in dxTry
        ode.derivatives(xStart, yStart, dyStart);
        odeSolver->solve(xStart, xEnd, y, dxEst);

        //Info<< "Numerical:  y = " << y << ", dxEst = " << dxEst << endl;

        // store the final results of ODE
        Y_ODE[celli] = y[0]; 
        q_ODE[celli] = y[1];
        // *****************end of the ODE solver*******************
        
    }

}

//updated the q to keep the tab
q_ODE.max(0.0); //clip


forAll(q,celli)
{
    q[celli] = q_ODE[celli];
}
