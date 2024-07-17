    // Solve the momentum equation

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

    if (pimple.momentumPredictor())
    {
	START(MomentumPredictor);
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
	STOP(MomentumPredictor);

        fvOptions.correct(U);
    }
