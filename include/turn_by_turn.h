#ifndef TURN_BY_TURN_H_INCLUDED
#define TURN_BY_TURN_H_INCLUDED

#include <vector>
#include "beam.h"
#include "cooler.h"
#include "ring.h"

void initialize_turn_by_turn_model(Beam &ion, Ring &ring);
void turn_by_turn_move_particles(Beam &ion, Ring &ring, Cooler &cooler);

#endif // TURN_BY_TURN_H_INCLUDED
