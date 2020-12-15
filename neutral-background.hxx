/*
 * Fixed neutral gas background, which doesn't change with the plasma
 * This is appropriate for cases where the neutral mean free path is long
 * in comparison to the plasma
 */

#ifndef NEUTRAL_BACKGROUND_H__
#define NEUTRAL_BACKGROUND_H__

#include "neutral-model.hxx"

class NeutralBackground : public NeutralModel {
public:
  NeutralBackground(Solver* solver, Mesh* mesh, Options& options);
  ~NeutralBackground() {}

  /// Update plasma quantities
  void update(const Field3D& Ne, const Field3D& Te, const Field3D& Ti, const Field3D& Vi);

  Field3D getDensity() override { return Nn; }
private:
  Field3D Nn, Pn, Vn; ///< Density, pressure and parallel velocity
  Field3D Tn; ///< neutral temperature
};

#endif // NEUTRAL_BACKGROUND_H__
