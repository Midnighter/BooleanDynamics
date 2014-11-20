
#ifndef DYNAMICS_H
#define DYNAMICS_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


unsigned char majority_activation(unsigned char current, unsigned char *states,
        int *inc_adj, int *regulation, const int num_pred);
void attractor_finder(unsigned char *states, int *inc_adj, int *regulation, int *inc_ptr,
        const int num_nodes, size_t *attractor_id, const size_t num_states);
void time_series(unsigned char *states, int *inc_adj, int *regulation, int *inc_ptr,
        const int num_nodes, unsigned char *series, const size_t time);

#endif // DYNAMICS_H

