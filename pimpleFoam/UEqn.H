// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

START(MatrixAssemblyU);
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();
STOP(MatrixAssemblyU);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    START(MomentumPredictor);
    solve(UEqn == -fvc::grad(p));
    STOP(MomentumPredictor);

    fvOptions.correct(U);
}
