// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

START(MatrixAssemblyU);
fvVectorMatrix UEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
STOP(MatrixAssemblyU);

UEqn.relax();

fvOptions.constrain(UEqn);

if (piso.momentumPredictor())
{
    START(MomentumPredictor);
    solve(UEqn == -fvc::grad(p));
    STOP(MomentumPredictor);

    fvOptions.correct(U);
}
