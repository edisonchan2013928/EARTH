#pragma once
#ifndef COUNT_H
#define COUNT_H

#include "Network.h"

double point_edge_count_same_edge(model& our_model, int p_index, int e_index);
double point_edge_count(model& our_model, int e_index, int node_u, int node_v);

#endif