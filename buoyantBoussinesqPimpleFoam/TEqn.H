{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    START(MatrixAssemblyT);
    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );
    STOP(MatrixAssemblyT);

    TEqn.relax();

    fvOptions.constrain(TEqn);

    START(SolveT);
    TEqn.solve();
    STOP(SolveT);

    radiation->correct();

    fvOptions.correct(T);

    rhok = 1.0 - beta*(T - TRef);
}
