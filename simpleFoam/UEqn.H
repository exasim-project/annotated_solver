    // Momentum predictor

    MRF.correctBoundaryVelocity(U);

    START(MatrixAssemblyU);
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    STOP(MatrixAssemblyU);
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    START(MomentumPredictor);
    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
    STOP(MomentumPredictor);
