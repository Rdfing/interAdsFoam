// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "argList.H"
#include "IOmanip.H"
#include "ODESystem.H"
#include "ODESolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class adsODE
:
public ODESystem
{
    // some parameters
    const scalar beta_;
    const scalar kf_;
    const scalar qMax_;
    const scalar kr_;

    public:

    adsODE(const scalar& beta_1, const scalar& kf_1, const scalar& qMax_1, const scalar& kr_1)
    :
    ODESystem(),
    beta_(beta_1),
    kf_(kf_1),
    qMax_(qMax_1),
    kr_(kr_1)
    {}

    label nEqns() const
    {
        return 2;
    }

    void derivatives
    (
        const scalar x,
        const scalarField& y,
        scalarField& dydx
    ) const
    {
        dydx[0] = -kf_*y[0]*(1-(y[1]/qMax_))+beta_*kr_*y[1];           //dydt
        dydx[1] = (1/beta_)*kf_*y[0]*(1-(y[1]/qMax_))-kr_*y[1];        //dq/dt
    }

    void jacobian
    (
        const scalar x,
        const scalarField& y,
        scalarField& dfdx,
        scalarSquareMatrix& dfdy
    ) const
    {
        dfdx[0] = 0.0;
        dfdx[1] = 0.0;

        dfdy(0, 0) = -kf_*(1-(y[1]/qMax_));
        dfdy(0, 1) = kf_*y[0]/qMax_+beta_*kr_;
        dfdy(1, 0) = (1/beta_)*kf_*(1-(y[1]/qMax_));
        dfdy(1, 1) =  -(1/beta_)*kf_*y[0]/qMax_-kr_;
    }
    
};

// ************************************************************************* //
