#ifndef ODIN_GRAVITY_H
#define ODIN_GRAVITY_H

#include "gravitational/gravitational_base.hpp"
#include "gravitational/gravitational_concept.hpp"

#ifdef USE_GSL_ELLIPSOIDAL_GRAVITY
#include "gravitational/ellipsoidal.hpp"
#endif//USE_GSL_ELLIPSOIDAL_GRAVITY

#ifdef USE_ESA_POLYHEDRAL_GRAVITY
#include "gravitational/polyhedral.hpp"
#endif//USE_ESA_POLYHEDRAL_GRAVITY

#include "gravitational/spherical.hpp"

#endif// ODIN_GRAVITY_H