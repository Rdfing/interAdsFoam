// singlePhaseTransportModel laminarTransport(U, phi);

// Laminar Schmidt number
dimensionedScalar Sc("Sc", dimless, mixture);

// Turbulent Schmidt number
dimensionedScalar Sct("Sct", dimless, mixture);


// molecular diffusivity
const dimensionedScalar D("D", dimViscosity, mixture);

// adsorbent density
const dimensionedScalar rhop("rhop", dimDensity, mixture);

//porosity
const dimensionedScalar porosity("porosity", dimless, mixture);


// Langmuir isotherm parameters
// the isotherm plot must be provided in the SI unit, C (kg/m3) and Qe (g/g)
const dimensionedScalar Langmuir_qMax("Langmuir_qMax", dimless, mixture);

// effective diffusivity based on Millington-Quirk correlation
const dimensionedScalar Dm = D*porosity/(Foam::pow(porosity,-1.0/3.0));

// This is the Hydrodyanmic dispersion tensor 
const dimensionedScalar alpha_L("alpha_L", dimLength, mixture); //longitudinal dispersivity (1/m)
const dimensionedScalar alpha_T("alpha_T", dimLength, mixture); //transverse dispersivity (1/m)
const dimensionedScalar tau_p("tau_p", dimless, mixture); // tortuosity factor < 1.0 (squared ratio of sample length to flow-path length)


// Forward and backward reaction rate, both has unit of 1/s
const dimensionedScalar kf("kf", dimensionSet(0,0,-1,0,0,0,0), mixture);
const dimensionedScalar kr("kr", dimensionSet(0,0,-1,0,0,0,0), mixture);


//- Residual volume-fraction to stablize the solution, choose a small value
const dimensionedScalar residualAlpha_ (dimensionedScalar::lookupOrAddToDict("residualAlpha",mixture,ROOTSMALL));
Info << "residualAlpha: " << residualAlpha_ << endl;

// swithc for freezing momentum calculation
const Switch limitY_ (Switch::lookupOrAddToDict("limitY",mixture,false));
