
#include "neutral-background.hxx"

NeutralBackground::NeutralBackground(Solver* UNUSED(solver), Mesh* mesh,
                                     Options& options) : NeutralModel(options) {
  auto& optroot = Options::root();
  // Get the fixed density, pressure and velocity
  Nn = optroot["Nn"]["function"].as(Field3D(0.0, mesh)); // Must set the density
  Pn = optroot["Pn"]["function"].as(Field3D(0.0, mesh)); // Must set pressure
  Vn = optroot["Vn"]["function"].withDefault(Field3D(0.0, mesh)); // Optional velocity

  Tn = Pn / Nn;

  SAVE_ONCE(Nn, Pn, Vn, Tn);
}

void NeutralBackground::update(const Field3D& Ne, const Field3D& Te, const Field3D& Ti,
                               const Field3D& Vi) {
  // Atomic processes
  TRACE("Atomic processes");

  Field3D Riz, Rrc, Rcx;
  neutral_rates(Ne, Te, Ti, Vi, Nn, Tn, Vn, S, F, Qi, Rp, Riz, Rrc, Rcx);

  // Collision rate which appears in the vorticity equation
  Fperp = Rcx + Rrc;

  // Assume that the neutrals radiate any energy from the ions
  Rn = Qi;
}
