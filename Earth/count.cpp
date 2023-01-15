#include "count.h"

double point_edge_count_same_edge(model& our_model, int p_index, int e_index)
{
	double edge_K_value = 0;
	double delta_spatial;
	double delta_temporal;
	vector<Point>& PS = our_model.edge_set[e_index].PS;

	if (our_model.method == 1 || our_model.method == 2) //RQS and SPS
	{
		for (int p = 0; p < (int)PS.size(); p++)
		{
			if (p_index == p)
				continue;

			delta_spatial = fabs(PS[p_index].dist_n1 - PS[p].dist_n1);
			delta_temporal = fabs(PS[p_index].time - PS[p].time);
			if (delta_spatial <= our_model.s_threshold && delta_temporal <= our_model.t_threshold)
				edge_K_value++;
		}
	}

	if (our_model.method == 3) //EAR
	{
		int count_value;
		Range_Tree& t = our_model.edge_set[e_index].t;
		if ((int)PS.size() == 0)
			return 0;

		double dist_L = max(PS[p_index].dist_n1 - our_model.s_threshold, 0.0);
		double dist_U = min(our_model.edge_set[e_index].length - (PS[p_index].dist_n2 - our_model.s_threshold), 
			our_model.edge_set[e_index].length);

		count_value = t.two_D_range_counting(t, dist_L, dist_U, PS[p_index].time - our_model.t_threshold, PS[p_index].time + our_model.t_threshold);

		if (count_value > 0) //Remove the count value of itself.
			count_value--;

		return count_value;
	}

	return edge_K_value;
}

double point_edge_count(model& our_model, int e_index, int node_u, int node_v)
{
	double delta_spatial;
	double delta_temporal;
	double K_value = 0;
	vector<Point>& PS = our_model.edge_set[e_index].PS;
	vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS;

	if (our_model.method == 1) //RQS
	{
		for (int p = 0; p < (int)PS.size(); p++)
		{
			delta_spatial = min(our_model.sp_node_vec_first[node_u].cur_sp_value + PS[p].dist_n1,
				our_model.sp_node_vec_first[node_v].cur_sp_value + PS[p].dist_n2);
			delta_temporal = fabs(PS_wide_e[our_model.cur_point_index].time - PS[p].time);
			if (delta_spatial <= our_model.s_threshold && delta_temporal <= our_model.t_threshold)
				K_value++;
		}
	}

	if (our_model.method == 2) //SPS
	{
		double delta_a;
		double delta_b;
		vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS;
		int p_wide_e = our_model.cur_point_index;

		for (int p = 0; p < (int)PS.size(); p++)
		{
			delta_a = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_u].cur_sp_value + PS[p].dist_n1,
				PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_v].cur_sp_value + PS[p].dist_n2);
			delta_b = min(PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_u].cur_sp_value + PS[p].dist_n1,
				PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_v].cur_sp_value + PS[p].dist_n2);

			delta_spatial = min(delta_a, delta_b);
			delta_temporal = fabs(PS_wide_e[our_model.cur_point_index].time - PS[p].time);

			if (delta_spatial <= our_model.s_threshold && delta_temporal <= our_model.t_threshold)
				K_value++;
		}
	}

	if (our_model.method == 3) //EAR
	{
		vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS;
		int p_wide_e = our_model.cur_point_index;
		double dist_p_i_u; //Equation 10
		double dist_p_i_v; //Equation 11

		if ((int)PS.size() == 0)
			return 0;

		dist_p_i_u = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_u].cur_sp_value,
			PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_u].cur_sp_value);
		dist_p_i_v = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_v].cur_sp_value,
			PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_v].cur_sp_value);

		if (dist_p_i_u > our_model.s_threshold && dist_p_i_v > our_model.s_threshold)
			return 0;

		if (dist_p_i_u <= our_model.s_threshold && dist_p_i_v > our_model.s_threshold)
		{
			Range_Tree& t = our_model.edge_set[e_index].t;

			return t.two_D_range_counting(t, 0, our_model.s_threshold - dist_p_i_u,
				PS_wide_e[our_model.cur_point_index].time - our_model.t_threshold, PS_wide_e[our_model.cur_point_index].time + our_model.t_threshold);
		}
		
		if (dist_p_i_u > our_model.s_threshold && dist_p_i_v <= our_model.s_threshold)
		{
			Edge& edge = our_model.edge_set[e_index];
			Range_Tree& t = edge.t;

			return t.two_D_range_counting(t, edge.length - (our_model.s_threshold - dist_p_i_v), edge.length,
				PS_wide_e[our_model.cur_point_index].time - our_model.t_threshold, PS_wide_e[our_model.cur_point_index].time + our_model.t_threshold);
		}

		if (dist_p_i_u <= our_model.s_threshold && dist_p_i_v <= our_model.s_threshold)
		{
			Edge& edge = our_model.edge_set[e_index];
			Range_Tree& t = edge.t;
			double sigma_u = our_model.s_threshold - dist_p_i_u;
			double sigma_v = our_model.s_threshold - dist_p_i_v;

			if (sigma_u + sigma_v > edge.length)
				return t.two_D_range_counting(t, 0, edge.length, 
					PS_wide_e[our_model.cur_point_index].time - our_model.t_threshold, PS_wide_e[our_model.cur_point_index].time + our_model.t_threshold);
			else
				return (t.two_D_range_counting(t, 0, sigma_u, PS_wide_e[our_model.cur_point_index].time - our_model.t_threshold, PS_wide_e[our_model.cur_point_index].time + our_model.t_threshold) +
					t.two_D_range_counting(t, edge.length - sigma_v, edge.length, PS_wide_e[our_model.cur_point_index].time - our_model.t_threshold, PS_wide_e[our_model.cur_point_index].time + our_model.t_threshold));
		}
	}

	return K_value;
}