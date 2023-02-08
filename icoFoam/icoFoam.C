/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    icoFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids.

    \heading Solver details
    The solver uses the PISO algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}}
          + \div \left( \vec{U} \vec{U} \right)
          - \div \left(\nu \grad \vec{U} \right)
          = - \grad p
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define START(NAME) auto start_##NAME = std::chrono::steady_clock::now()

#define STOP(NAME)                                                       \
	auto end_##NAME = std::chrono::steady_clock::now();              \
	Info << "[INFO] " #NAME ": " <<                                  \
	std::chrono::duration_cast<std::chrono::microseconds>            \
	(end_##NAME-start_##NAME).count()/1000.0 << " [ms]" << endl 

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, laminar flow"
        " of Newtonian fluids."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
	START(TimeStep);
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // Momentum predictor

	START(MatrixAssemblyU);
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );
	STOP(MatrixAssemblyU);


        if (piso.momentumPredictor())
        {
            START(MomentumPredictor);
            solve(UEqn == -fvc::grad(p));
	    STOP(MomentumPredictor);
        }

        // --- PISO loop
        while (piso.correct())
        {
            START(PISOStep);
	    START(MatrixAssemblyPI);
            volScalarField rAU(1.0/UEqn.A());
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
            );
	    STOP(MatrixAssemblyPI);

	    START(AdjustPhi);
            adjustPhi(phiHbyA, U, p);
	    STOP(AdjustPhi);

            // Update the pressure BCs to ensure flux consistency
	    START(ConstrainPressure);
            constrainPressure(p, U, phiHbyA, rAU);
	    STOP(ConstrainPressure);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                // Pressure corrector

		START(MatrixAssemblyPII);
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                );
		STOP(MatrixAssemblyPII);

                pEqn.setReference(pRefCell, pRefValue);

		START(SolveP);
                pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
		STOP(SolveP);

                if (piso.finalNonOrthogonalIter())
                {
                    phi = phiHbyA - pEqn.flux();
                }

            }

	    START(ContinuityError);
            #include "continuityErrs.H"
	    STOP(ContinuityError);

            U = HbyA - rAU*fvc::grad(p);
            START(CorrectBoundary);
            U.correctBoundaryConditions();
	    STOP(CorrectBoundary);
            STOP(PISOStep);
        }

	STOP(TimeStep);
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
