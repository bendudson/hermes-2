
#pragma once
#ifndef HERMES_SPECIES_H
#define HERMES_SPECIES_H

#include "state.hxx"

/// Class representing a particle species
/// This is the interface that the main Hermes code uses to access state
struct Species {
  SpeciesState getState(BoutReal time);
  void setTimeDerivatives(const SpeciesState &state);
};

#endif // HERMES_SPECIES_H
