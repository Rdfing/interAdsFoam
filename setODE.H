    // constant parameter to simplify the expression
    const dimensionedScalar beta = ((1-porosity)/porosity)*rhop;

    // Create the selected ODE system solver
    // ODE system parameter
    const scalar beta_1 = beta.value();

    const scalar kf_1 = kf.value(); 
    const scalar kr_1 = kr.value();

    // langmuir isotherm parameters, pay attention to the unit!
    const scalar qMax_1 = Langmuir_qMax.value();

    // Initial condition
    //const scalar y0 = 2e-20;

    // Create the ODE system
    adsODE ode(beta_1, kf_1, qMax_1, kr_1);

    // construct the solver
    word ODESolverName ("rodas23");
    dictionary dict;
    dict.add("solver", ODESolverName);
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    // Required to store initial condition 
    scalarField yStart(ode.nEqns());
    // Required to store dydx
    scalarField dyStart(ode.nEqns());

    // some solving controls
    
    //odeSolver->maxSteps() = 20000; //10000 is the default
    //odeSolver->absTol() = 1e-12; //
    odeSolver->relTol() = 1e-4; //1e-4; is the default