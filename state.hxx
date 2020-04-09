
#pragma once
#ifndef HERMES_STATE_H
#define HERMES_STATE_H

#include <map>
#include <string>

/// Normalisation factors. A normalised quantity multiplied
/// by these values has the corresponding units
///
/// e.g. Ni : Dimensionless density, evolving variable
///      Ni * cubic_meters : Density in cubic meters
struct Normalisation {
  BoutReal cubic_meters;  // density normalisation
  BoutReal eV;  // temperature normalisation
  BoutReal meters_per_second; // velocity normalisation
  BoutReal seconds;  // time normalisation
};

struct SpeciesState {
  BoutReal AA; // Atomic mass
  BoutReal ZZ; // Charge
  
  Field3D N;  // Density
  Field3D P;  // Pressure
  Field3D NV; // Momentum, includes factor of AA
  Field3D T;  // Temperature
  Field3D V;  // Velocity
}

/// Represents the state of the simulation at a given time
struct SimulationState {
  BoutReal time;
  Normalisation units; // Normalisation factors

  std::map<std::string, Field3D> fields;       ///< Electromagnetic fields
  std::map<std::string, SpeciesState> species; ///< Particle species
};

#endif // HERMES_STATE_H
