
#include "mixed.hxx"

#include <bout/fv_ops.hxx>
#include <boutexception.hxx>
#include <options.hxx>

#include "div_ops.hxx"

using bout::globals::mesh;

NeutralMixed::NeutralMixed(Solver *solver, Mesh *UNUSED(mesh), Options &options)
    : NeutralModel(options) {
  solver->add(Nn, "Nn");
  solver->add(Pn, "Pn");
  solver->add(NVn, "NVn");

  sheath_ydown = options["sheath_ydown"]
                     .doc("Enable wall boundary conditions at ydown")
                     .withDefault<bool>(true);

  sheath_yup = options["sheath_yup"]
                   .doc("Enable wall boundary conditions at yup")
                   .withDefault<bool>(true);

  neutral_gamma = options["neutral_gamma"]
          .doc("Heat flux to the wall q = neutral_gamma * n * T * cs")
          .withDefault(5. / 4);

  nn_floor = options["nn_floor"]
                 .doc("A minimum density used when dividing NVn by Nn. "
                      "Normalised units.")
                 .withDefault(1e-4);

  precondition = options["precondition"]
                     .doc("Enable preconditioning in neutral model?")
                     .withDefault<bool>(true);

  if (precondition) {
    inv = std::unique_ptr<Laplacian>(Laplacian::create(&options["precon_laplace"]));

    inv->setInnerBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);
    inv->setOuterBoundaryFlags(INVERT_DC_GRAD | INVERT_AC_GRAD);

    inv->setCoefA(1.0);
  }

  // Optionally output time derivatives
  if (options["output_ddt"]
          .doc("Save derivatives to output?")
          .withDefault<bool>(false)) {
    SAVE_REPEAT(ddt(Nn), ddt(Pn), ddt(NVn));
  }

  if (options["diagnose"]
          .doc("Save additional diagnostic fields (Dnn) ?")
          .withDefault<bool>(false)) {
    SAVE_REPEAT(Dnn);
  }

  Dnn = 0.0; // Neutral gas diffusion

  S = 0;
  F = 0;
  Qi = 0;
  Rp = 0;
}

void NeutralMixed::update(const Field3D &Ne, const Field3D &Te,
                          const Field3D &Ti, const Field3D &Vi) {
  TRACE("NeutralMixed::update");

  mesh->communicate(Nn, Pn, NVn);

  Nn.clearParallelSlices();
  Pn.clearParallelSlices();
  NVn.clearParallelSlices();

  Coordinates *coord = mesh->getCoordinates();

  //////////////////////////////////////////////////////
  // 3D model, diffusive in X-Z, fluid in Y
  //
  // Evolves neutral density Nn, pressure Pn, and
  // parallel momentum NVn
  //

  Nn = floor(Nn, 1e-8);
  Pn = floor(Pn, 1e-10);
  // Nnlim Used where division by neutral density is needed
  Field3D Nnlim = floor(Nn, nn_floor);
  Field3D Tn = Pn / Nn;
  // Tn = floor(Tn, 0.01 / Tnorm);
  Vn = NVn / Nn;
  Field3D Vnlim = NVn / Nnlim; // Neutral parallel velocity

  Pnlim = Nn * Tn;
  Pnlim.applyBoundary("neumann");

  /////////////////////////////////////////////////////
  // Boundary conditions
  TRACE("Neutral boundary conditions");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->ystart, jz) -
                                 Nn(r.ind, mesh->ystart + 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->ystart, jz);

        Nn(r.ind, mesh->ystart - 1, jz) =
            2 * nnwall - Nn(r.ind, mesh->ystart, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->ystart - 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->ystart - 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->ystart, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->ystart - 1, jz) = Pn(r.ind, mesh->ystart, jz);

        // No flow into wall
        Vn(r.ind, mesh->ystart - 1, jz) = -Vn(r.ind, mesh->ystart, jz);
        Vnlim(r.ind, mesh->ystart - 1, jz) = -Vnlim(r.ind, mesh->ystart, jz);
        NVn(r.ind, mesh->ystart - 1, jz) = -NVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        // Free boundary (constant gradient) density
        BoutReal nnwall = 0.5 * (3. * Nn(r.ind, mesh->yend, jz) -
                                 Nn(r.ind, mesh->yend - 1, jz));
        if (nnwall < 0.0)
          nnwall = 0.0;

        BoutReal tnwall = Tn(r.ind, mesh->yend, jz);

        Nn(r.ind, mesh->yend + 1, jz) = 2 * nnwall - Nn(r.ind, mesh->yend, jz);

        // Zero gradient temperature, heat flux added later
        Tn(r.ind, mesh->yend + 1, jz) = tnwall;

        // Set pressure consistent at the boundary
        // Pn(r.ind, mesh->yend + 1, jz) =
        //     2. * nnwall * tnwall - Pn(r.ind, mesh->yend, jz);

        // Zero-gradient pressure
        Pn(r.ind, mesh->yend + 1, jz) = Pn(r.ind, mesh->yend, jz);

        // No flow into wall
        Vn(r.ind, mesh->yend + 1, jz) = -Vn(r.ind, mesh->yend, jz);
        Vnlim(r.ind, mesh->yend + 1, jz) = -Vnlim(r.ind, mesh->yend, jz);
        NVn(r.ind, mesh->yend + 1, jz) = -NVn(r.ind, mesh->yend, jz);
      }
    }
  }

  /////////////////////////////////////////////////////
  // Atomic processes
  TRACE("Atomic processes");

  Field3D Riz, Rrc, Rcx;
  neutral_rates(Ne, Te, Ti, Vi, Nn, Tn, Vn, S, F, Qi, Rp, Riz, Rrc, Rcx);

  // Neutral cross-field diffusion coefficient
  BoutReal neutral_lmax = 0.1 / Lnorm;
  Field3D Rnn = Nn * sqrt(Tn) / neutral_lmax; // Neutral-neutral collisions
  Dnn = Pnlim / (Riz + Rcx + Rnn);

  mesh->communicate(Dnn);
  Dnn.clearParallelSlices();
  Dnn.applyBoundary("dirichlet_o2");

  // Apply a Dirichlet boundary condition to all the coefficients
  // used in diffusion operators. This is to ensure that the flux through
  // the boundary is zero.
  Field3D DnnPn = Dnn * Pn;
  DnnPn.applyBoundary("dirichlet_o2");
  Field3D DnnNn = Dnn * Nn;
  DnnNn.applyBoundary("dirichlet_o2");
  Field3D DnnNVn = Dnn * NVn;
  DnnNVn.applyBoundary("dirichlet_o2");

  if (sheath_ydown) {
    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->ystart - 1, jz) = -Dnn(r.ind, mesh->ystart, jz);
        DnnPn(r.ind, mesh->ystart - 1, jz) = -DnnPn(r.ind, mesh->ystart, jz);
        DnnNn(r.ind, mesh->ystart - 1, jz) = -DnnNn(r.ind, mesh->ystart, jz);
        DnnNVn(r.ind, mesh->ystart - 1, jz) = -DnnNVn(r.ind, mesh->ystart, jz);
      }
    }
  }

  if (sheath_yup) {
    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {
        Dnn(r.ind, mesh->yend + 1, jz) = -Dnn(r.ind, mesh->yend, jz);
        DnnPn(r.ind, mesh->yend + 1, jz) = -DnnPn(r.ind, mesh->yend, jz);
        DnnNn(r.ind, mesh->yend + 1, jz) = -DnnNn(r.ind, mesh->yend, jz);
        DnnNVn(r.ind, mesh->yend + 1, jz) = -DnnNVn(r.ind, mesh->yend, jz);
      }
    }
  }

  // Logarithms used to calculate perpendicular velocity
  // V_perp = -Dnn * ( Grad_perp(Nn)/Nn + Grad_perp(Tn)/Tn )
  //
  // Grad(Pn) / Pn = Grad(Tn)/Tn + Grad(Nn)/Nn
  //               = Grad(logTn + logNn)
  // Field3D logNn = log(Nn);
  // Field3D logTn = log(Tn);

  Field3D logPnlim = log(Pnlim);
  logPnlim.applyBoundary("neumann");

  // Sound speed appearing in Lax flux for advection terms
  Field3D sound_speed = sqrt(Tn * (5. / 3));

  /////////////////////////////////////////////////////
  // Neutral density
  TRACE("Neutral density");

  ddt(Nn) = -FV::Div_par(Nn, Vnlim, sound_speed) // Advection
            + S // Source from recombining plasma
            + FV::Div_a_Laplace_perp(DnnNn, logPnlim) // Perpendicular diffusion
      ;

  /////////////////////////////////////////////////////
  // Neutral momentum
  TRACE("Neutral momentum");

  // Relationship between heat conduction and viscosity for neutral
  // gas Chapman, Cowling "The Mathematical Theory of Non-Uniform
  // Gases", CUP 1952 Ferziger, Kaper "Mathematical Theory of
  // Transport Processes in Gases", 1972
  // eta_n = (2. / 5) * kappa_n;

  ddt(NVn) =
      -FV::Div_par(NVn, Vnlim, sound_speed)          // Momentum flow
      + F                                            // Friction with plasma
      - Grad_par(Pnlim)                              // Pressure gradient
      + FV::Div_a_Laplace_perp(DnnNVn, logPnlim)     // Perpendicular diffusion
      + FV::Div_a_Laplace_perp((2. / 5) * DnnNn, Vn) // Perpendicular viscosity
      + FV::Div_par_K_Grad_par((2. / 5) * DnnNn, Vn) // Parallel viscosity
      ;

  Fperp = Rcx + Rrc;

  /////////////////////////////////////////////////////
  // Neutral pressure
  TRACE("Neutral pressure");

  ddt(Pn) = -FV::Div_par(Pn, Vnlim, sound_speed) // Advection
            - (2. / 3) * Pnlim * Div_par(Vn)     // Compression
            + (2. / 3) * Qi                      // Energy exchange with ions
            + FV::Div_a_Laplace_perp(DnnPn, logPnlim) // Perpendicular diffusion
            + FV::Div_a_Laplace_perp(DnnNn, Tn)       // Conduction
            + FV::Div_par_K_Grad_par(DnnNn, Tn)       // Parallel conduction
      ;

  //////////////////////////////////////////////////////
  // Boundary condition on fluxes
  TRACE("Neutral boundary fluxes");

  if (neutral_gamma > 0.0) {
    if (sheath_ydown) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          // Loss of thermal energy to the target.
          // This depends on the reflection coefficient
          // and is controlled by the option neutral_gamma
          //         q = neutral_gamma * n * T * cs

          // Density at the target
          BoutReal Nnout = 0.5 * (Nn(r.ind, mesh->ystart, jz) +
                                  Nn(r.ind, mesh->ystart - 1, jz));
          if (Nnout < 0.0)
            Nnout = 0.0;
          // Temperature at the target
          BoutReal Tnout = 0.5 * (Tn(r.ind, mesh->ystart, jz) +
                                  Tn(r.ind, mesh->ystart - 1, jz));
          if (Tnout < 0.0)
            Tnout = 0.0;

          // gamma * n * T * cs
          BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
          // Multiply by cell area to get power
          BoutReal heatflux = q *
                              (coord->J(r.ind, mesh->ystart) +
                               coord->J(r.ind, mesh->ystart - 1)) /
                              (sqrt(coord->g_22(r.ind, mesh->ystart)) +
                               sqrt(coord->g_22(r.ind, mesh->ystart - 1)));

          // Divide by volume of cell, and multiply by 2/3 to get pressure
          ddt(Pn)(r.ind, mesh->ystart, jz) -=
              (2. / 3) * heatflux /
              (coord->dy(r.ind, mesh->ystart) * coord->J(r.ind, mesh->ystart));
        }
      }
    }

    if (sheath_yup) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          // Loss of thermal energy to the target.
          // This depends on the reflection coefficient
          // and is controlled by the option neutral_gamma
          //         q = neutral_gamma * n * T * cs

          // Density at the target
          BoutReal Nnout =
              0.5 * (Nn(r.ind, mesh->yend, jz) + Nn(r.ind, mesh->yend + 1, jz));
          if (Nnout < 0.0)
            Nnout = 0.0;
          // Temperature at the target
          BoutReal Tnout =
              0.5 * (Tn(r.ind, mesh->yend, jz) + Tn(r.ind, mesh->yend + 1, jz));
          if (Tnout < 0.0)
            Tnout = 0.0;

          // gamma * n * T * cs
          BoutReal q = neutral_gamma * Nnout * Tnout * sqrt(Tnout);
          // Multiply by cell area to get power
          BoutReal heatflux =
              q *
              (coord->J(r.ind, mesh->yend) + coord->J(r.ind, mesh->yend + 1)) /
              (sqrt(coord->g_22(r.ind, mesh->yend)) +
               sqrt(coord->g_22(r.ind, mesh->yend + 1)));

          // Divide by volume of cell, and multiply by 2/3 to get pressure
          ddt(Pn)(r.ind, mesh->yend, jz) -=
              (2. / 3) * heatflux /
              (coord->dy(r.ind, mesh->yend) * coord->J(r.ind, mesh->yend));
        }
      }
    }
  }

  BOUT_FOR(i, Pn.getRegion("RGN_ALL")) {
    if ((Pn[i] < 1e-9) && (ddt(Pn)[i] < 0.0)) {
      ddt(Pn)[i] = 0.0;
    }
    if ((Nn[i] < 1e-7) && (ddt(Nn)[i] < 0.0)) {
      ddt(Nn)[i] = 0.0;
    }
  }
}

void NeutralMixed::precon(BoutReal UNUSED(t), BoutReal gamma,
                          BoutReal UNUSED(delta)) {
  if (!precondition) {
    return;
  }

  // Neutral gas diffusion
  // Solve (1 - gamma*Dnn*Delp2)^{-1}

  inv->setCoefD(-gamma * Dnn);

  ddt(Nn) = inv->solve(ddt(Nn));
  ddt(NVn) = inv->solve(ddt(NVn));
  ddt(Pn) = inv->solve(ddt(Pn));
}
