
#include "dynamics.h"


static void
set_states(size_t state_id, unsigned char *states, const int num_states)
{
    int i = 0;
    for (i = 0; i < num_states; ++i) {
        states[i] = (unsigned char)(state_id & 1);
        state_id >>= 1;
    }
}

static size_t
get_state_id(unsigned char *states, const int num_states)
{
    int i = 0;
    size_t id = states[num_states - 1];
    for (i = num_states - 2; i >= 0 ; --i) {
        id <<= 1;
        id |= states[i];
    }
    return id;
}

unsigned char
majority_activation(unsigned char current, unsigned char *states,
        int *inc_adj, int *regulation, const int num_pred)
{
    int i = 0;
    int pred = 0;
    int state_sum = 0;
    for (i = 0; i < num_pred; ++i) {
        state_sum += states[inc_adj[i]] * regulation[i];
    }
    if (state_sum > 0) {
        return (unsigned char)1;
    }
    else if (state_sum < 0) {
        return (unsigned char)0;
    }
    else {
        return current;
    }
}

void
attractor_finder(unsigned char *states, int *inc_adj, int *regulation, int *inc_ptr,
        const int num_nodes,
        size_t *attractor_id, const size_t num_states)
{
    size_t i = 0;
    int n = 0;
    int begin = 0;
    int end = 0;
    size_t id = 0;
    size_t old_id = 0;
    unsigned char *old_states = NULL;
    old_states = malloc(num_nodes * sizeof(unsigned char));
    if (old_states == NULL) {
        printf("memory allocation failed");
        return;
    }
    for (i = 0; i < num_states; ++i) {
        memcpy(old_states, states, num_nodes * sizeof(unsigned char));
        old_id = id;
        for (n = 0; n < num_nodes; ++n) {
            begin = inc_ptr[n];
            end = inc_ptr[n + 1];
            states[n] = majority_activation(old_states[n], old_states,
                    &inc_adj[begin], &regulation[begin], end - begin);
            id = get_state_id(states, num_nodes);
        }
        attractor_id[i] = old_id;
    }
    free(states);
    free(old_states);
}

void
time_series(unsigned char *states, int *inc_adj, int *regulation, int *inc_ptr, const int num_nodes,
        unsigned char *series, const size_t time)
{
    size_t t = 0;
    int n = 0;
    int begin = 0;
    int end = 0;
    unsigned char *old_states = NULL;
    old_states = malloc(num_nodes * sizeof(unsigned char));
    if (old_states == NULL) {
        printf("memory allocation failed\n");
        return;
    }
    memcpy(&series[0], states, num_nodes * sizeof(unsigned char));
    for (t = 1; t < time; ++t) {
        memcpy(old_states, states, num_nodes * sizeof(unsigned char));
        for (n = 0; n < num_nodes; ++n) {
            begin = inc_ptr[n];
            end = inc_ptr[n + 1];
            states[n] = majority_activation(old_states[n], old_states,
                    &inc_adj[begin], &regulation[begin], end - begin);
            memcpy(&series[t * num_nodes], states, num_nodes * sizeof(unsigned char));
        }
    }
    free(old_states);
}

