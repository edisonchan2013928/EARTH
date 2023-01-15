#include "count.h"

double point_edge_count_same_edge(model& our_model, int p_index, int e_index)
{
	double edge_K_value = 0;
	double delta_spatial;
	double delta_temporal;
	vector<Point>& PS = our_model.edge_set[e_index].PS_vec[our_model.cur_data_index];

	if (our_model.method == 1 || our_model.method == 2) //RQS and SPS
	{
		for (int p = 0; p < (int)PS.size(); p++)
		{
			if (p_index == p)
				continue;

			delta_spatial = fabs(PS[p_index].dist_n1 - PS[p].dist_n1);
			delta_temporal = fabs(PS[p_index].time - PS[p].time);
			if (delta_spatial <= our_model.s_tau_vec[our_model.cur_s_tau_index] && delta_temporal <= our_model.t_tau_vec[our_model.cur_t_tau_index])
				edge_K_value++;
		}
	}

	if (our_model.method == 3 || our_model.method == 5) //EAR or EATH methods
	{
		int count_value;
		Range_Tree& t = our_model.edge_set[e_index].t_vec[our_model.cur_data_index];
		if ((int)PS.size() == 0)
			return 0;

		double dist_L = max(PS[p_index].dist_n1 - our_model.s_tau_vec[our_model.cur_s_tau_index], 0.0);
		double dist_U = min(our_model.edge_set[e_index].length - (PS[p_index].dist_n2 - our_model.s_tau_vec[our_model.cur_s_tau_index]),
			our_model.edge_set[e_index].length);

		count_value = t.two_D_range_counting(t, dist_L, dist_U, PS[p_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS[p_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]);

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
	vector<Point>& PS = our_model.edge_set[e_index].PS_vec[our_model.cur_data_index];
	vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS_vec[our_model.cur_data_index];

	if (our_model.method == 1) //RQS
	{
		for (int p = 0; p < (int)PS.size(); p++)
		{
			delta_spatial = min(our_model.sp_node_vec_first[node_u].cur_sp_value + PS[p].dist_n1,
				our_model.sp_node_vec_first[node_v].cur_sp_value + PS[p].dist_n2);
			delta_temporal = fabs(PS_wide_e[our_model.cur_point_index].time - PS[p].time);
			if (delta_spatial <= our_model.s_tau_vec[our_model.cur_s_tau_index] && delta_temporal <= our_model.t_tau_vec[our_model.cur_t_tau_index])
				K_value++;
		}
	}

	if (our_model.method == 2) //SPS
	{
		double delta_a;
		double delta_b;
		vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS_vec[our_model.cur_data_index];
		int p_wide_e = our_model.cur_point_index;

		for (int p = 0; p < (int)PS.size(); p++)
		{
			delta_a = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_u].cur_sp_value + PS[p].dist_n1,
				PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_v].cur_sp_value + PS[p].dist_n2);
			delta_b = min(PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_u].cur_sp_value + PS[p].dist_n1,
				PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_v].cur_sp_value + PS[p].dist_n2);

			delta_spatial = min(delta_a, delta_b);
			delta_temporal = fabs(PS_wide_e[our_model.cur_point_index].time - PS[p].time);

			if (delta_spatial <= our_model.s_tau_vec[our_model.cur_s_tau_index] && delta_temporal <= our_model.t_tau_vec[our_model.cur_t_tau_index])
				K_value++;
		}
	}

	if (our_model.method == 3 || our_model.method == 5) //EAR
	{
		vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS_vec[our_model.cur_data_index];
		int p_wide_e = our_model.cur_point_index;
		double dist_p_i_u; //Equation 10
		double dist_p_i_v; //Equation 11

		if ((int)PS.size() == 0)
			return 0;

		dist_p_i_u = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_u].cur_sp_value,
			PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_u].cur_sp_value);
		dist_p_i_v = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_v].cur_sp_value,
			PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_v].cur_sp_value);

		if (dist_p_i_u > our_model.s_tau_vec[our_model.cur_s_tau_index] && dist_p_i_v > our_model.s_tau_vec[our_model.cur_s_tau_index])
			return 0;

		if (dist_p_i_u <= our_model.s_tau_vec[our_model.cur_s_tau_index] && dist_p_i_v > our_model.s_tau_vec[our_model.cur_s_tau_index])
		{
			Range_Tree& t = our_model.edge_set[e_index].t_vec[our_model.cur_data_index];

			return t.two_D_range_counting(t, 0, our_model.s_tau_vec[our_model.cur_s_tau_index] - dist_p_i_u,
				PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS_wide_e[our_model.cur_point_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]);
		}

		if (dist_p_i_u > our_model.s_tau_vec[our_model.cur_s_tau_index] && dist_p_i_v <= our_model.s_tau_vec[our_model.cur_s_tau_index])
		{
			Edge& edge = our_model.edge_set[e_index];
			Range_Tree& t = edge.t_vec[our_model.cur_data_index];

			return t.two_D_range_counting(t, edge.length - (our_model.s_tau_vec[our_model.cur_s_tau_index] - dist_p_i_v), edge.length,
				PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS_wide_e[our_model.cur_point_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]);
		}

		if (dist_p_i_u <= our_model.s_tau_vec[our_model.cur_s_tau_index] && dist_p_i_v <= our_model.s_tau_vec[our_model.cur_s_tau_index])
		{
			Edge& edge = our_model.edge_set[e_index];
			Range_Tree& t = edge.t_vec[our_model.cur_data_index];
			double sigma_u = our_model.s_tau_vec[our_model.cur_s_tau_index] - dist_p_i_u;
			double sigma_v = our_model.s_tau_vec[our_model.cur_s_tau_index] - dist_p_i_v;

			if (sigma_u + sigma_v > edge.length)
				return t.two_D_range_counting(t, 0, edge.length,
					PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS_wide_e[our_model.cur_point_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]);
			else
				return (t.two_D_range_counting(t, 0, sigma_u, PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS_wide_e[our_model.cur_point_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]) +
					t.two_D_range_counting(t, edge.length - sigma_v, edge.length, PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.cur_t_tau_index], PS_wide_e[our_model.cur_point_index].time + our_model.t_tau_vec[our_model.cur_t_tau_index]));
		}
	}

	return K_value;
}

void point_edge_count_same_edge_MTS(model& our_model, int e_index, int node_u, int node_v)
{
	vector<Point>& PS = our_model.edge_set[e_index].PS_vec[our_model.cur_data_index];
	int time_index;
	int dist_index;
	int t_min, t_max;
	double time_fix;
	//double time_min = PS[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];

	for (int p = 0; p < (int)PS.size(); p++)
	{
		if (our_model.cur_point_index == p)
			continue;

		if (PS[p].time > PS[our_model.cur_point_index].time)
		{
			time_fix = PS[our_model.cur_point_index].time + our_model.start_t_tau;
			time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold - eps) + 1;
			time_index = max(time_index, 0) + our_model.T;
		}
		else
		{
			time_fix = PS[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];
			time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold + eps);
			time_index = min(time_index, our_model.T - 1);
		}

		if (time_index < 0 || time_index > 2 * our_model.T - 1)
			continue;

		dist_index = floor((fabs(PS[p].dist_n1 - PS[our_model.cur_point_index].dist_n1) - our_model.s_tau_vec[0]) / our_model.delta_s_threshold - eps);
		dist_index = max(dist_index, -1) + 1;

		if (dist_index < 0 || dist_index > our_model.M - 1)
			continue;

		our_model.prefix_count[time_index][dist_index]++;
	}
	
	construct_prefix_count_array(our_model, our_model.prefix_count, our_model.M);
	for (int s = 0; s < our_model.M; s++)
	{
		for (int t = 0; t < our_model.T; t++)
		{
			t_min = our_model.T - t - 2;
			t_max = our_model.T + t;

			if (t_min < 0)
				our_model.multiple_K_values[our_model.cur_data_index][s][t] += our_model.prefix_count[t_max][s];
			else
				our_model.multiple_K_values[our_model.cur_data_index][s][t] += (our_model.prefix_count[t_max][s] - our_model.prefix_count[t_min][s]);
		}
	}
	clear_prefix_count_array(our_model, our_model.prefix_count, our_model.M);
}

void CASE(vector<Point>& PS_wide_e, vector<Point>& PS, model& our_model, double dist_p_i_node, int d_max_index, bool is_Case2)
{
	int time_index;
	int dist_index;
	int t_min, t_max;
	double time_fix;

	double sigma_node_1 = our_model.s_tau_vec[0] - dist_p_i_node;
	//double time_min = PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];

	for (int p = 0; p < (int)PS.size(); p++)
	{
		if (PS[p].time > PS_wide_e[our_model.cur_point_index].time)
		{
			time_fix = PS_wide_e[our_model.cur_point_index].time + our_model.start_t_tau;
			time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold - eps) + 1;
			time_index = max(time_index, 0) + our_model.T;
		}
		else
		{
			time_fix = PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];
			time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold + eps);
			time_index = min(time_index, our_model.T - 1);
		}

		if (time_index < 0 || time_index > 2 * our_model.T - 1)
			continue;

		/*time_index = (int)floor((PS[p].time - time_min) / our_model.delta_t_threshold);
		if (time_index < 0 || time_index > 2 * our_model.T - 1)
			continue;*/

		if (is_Case2 == true) //Pass through the node u
			dist_index = (int)max(floor((PS[p].dist_n1 - sigma_node_1) / our_model.delta_s_threshold - eps), -1.0) + 1;
		else //Pass through the node v
			dist_index = (int)max(floor((PS[p].dist_n2 - sigma_node_1) / our_model.delta_s_threshold - eps), -1.0) + 1;

		//if (dist_index < 0 || dist_index > d_max_index)
		if (dist_index < 0 || dist_index >= d_max_index)
			continue;

		our_model.prefix_count[time_index][dist_index]++;
	}

	construct_prefix_count_array(our_model, our_model.prefix_count, d_max_index);
	for (int s = 0; s < d_max_index; s++)
	{
		for (int t = 0; t < our_model.T; t++)
		{
			t_min = our_model.T - t - 2;
			t_max = our_model.T + t;

			if (t_min < 0)
				our_model.multiple_K_values[our_model.cur_data_index][s][t] += our_model.prefix_count[t_max][s];
			else
				our_model.multiple_K_values[our_model.cur_data_index][s][t] += (our_model.prefix_count[t_max][s] - our_model.prefix_count[t_min][s]);
		}
	}
	clear_prefix_count_array(our_model, our_model.prefix_count, d_max_index);
}

void point_edge_count_MTS(model& our_model, int e_index, int node_u, int node_v)
{
	vector<Point>& PS_wide_e = our_model.edge_set[our_model.cur_edge_index].PS_vec[our_model.cur_data_index];
	vector<Point>& PS = our_model.edge_set[e_index].PS_vec[our_model.cur_data_index];
	int p_wide_e = our_model.cur_point_index;
	double dist_p_i_u; //Equation 10
	double dist_p_i_v; //Equation 11

	if ((int)PS.size() == 0)
		return;

	dist_p_i_u = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_u].cur_sp_value,
		PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_u].cur_sp_value);
	dist_p_i_v = min(PS_wide_e[p_wide_e].dist_n1 + our_model.sp_node_vec_first[node_v].cur_sp_value,
		PS_wide_e[p_wide_e].dist_n2 + our_model.sp_node_vec_second[node_v].cur_sp_value);

	if (dist_p_i_u > our_model.s_tau_vec[our_model.M - 1] && dist_p_i_v > our_model.s_tau_vec[our_model.M - 1])
		return;
	if (dist_p_i_u <= our_model.s_tau_vec[our_model.M - 1] && dist_p_i_v > our_model.s_tau_vec[our_model.M - 1])
		CASE(PS_wide_e, PS, our_model, dist_p_i_u, our_model.M, true);
	if (dist_p_i_u > our_model.s_tau_vec[our_model.M - 1] && dist_p_i_v <= our_model.s_tau_vec[our_model.M - 1])
		CASE(PS_wide_e, PS, our_model, dist_p_i_v, our_model.M, false);

	if (dist_p_i_u <= our_model.s_tau_vec[our_model.M - 1] && dist_p_i_v <= our_model.s_tau_vec[our_model.M - 1])
	{
		double sigma_node_u; 
		double sigma_node_v;
		int max_s = -1;
		int time_index;
		//double time_min;
		double time_fix;
		int t_min; 
		int t_max;

		for (int s = our_model.M - 1; s >= 0; s--)
		{
			sigma_node_u = our_model.s_tau_vec[s] - dist_p_i_u;
			sigma_node_v = our_model.s_tau_vec[s] - dist_p_i_v;

			if (sigma_node_u + sigma_node_v <= our_model.edge_set[e_index].length)
			{
				max_s = s;
				break;
			}
		}

		if (max_s != -1)
		{
			CASE(PS_wide_e, PS, our_model, dist_p_i_u, max_s + 1, true);
			CASE(PS_wide_e, PS, our_model, dist_p_i_v, max_s + 1, false);
		}

		for (int p = 0; p < (int)PS.size(); p++)
		{
			if (PS[p].time > PS_wide_e[our_model.cur_point_index].time)
			{
				time_fix = PS_wide_e[our_model.cur_point_index].time + our_model.start_t_tau;
				time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold - eps) + 1;
				time_index = max(time_index, 0) + our_model.T;
			}
			else
			{
				time_fix = PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];
				time_index = (int)floor((PS[p].time - time_fix) / our_model.delta_t_threshold + eps);
				time_index = min(time_index, our_model.T - 1);
			}

			if (time_index < 0 || time_index > 2 * our_model.T - 1)
				continue;

			our_model.prefix_count_one_D[time_index]++;
		}

		/*time_min = PS_wide_e[our_model.cur_point_index].time - our_model.t_tau_vec[our_model.T - 1];
		for (int p = 0; p < (int)PS.size(); p++)
		{
			time_index = (int)floor((PS[p].time - time_min) / our_model.delta_t_threshold);
			if (time_index < 0 || time_index > 2 * our_model.T - 1)
				continue;

			our_model.prefix_count_one_D[time_index]++;
		}*/

		construct_prefix_count_one_D_array(our_model);

		for (int t = 0; t < our_model.T; t++)
		{
			t_min = our_model.T - t - 2;
			t_max = our_model.T + t;

			for (int s = max_s + 1; s < our_model.M; s++)
			{
				if (t_min < 0)
					our_model.multiple_K_values[our_model.cur_data_index][s][t] += our_model.prefix_count_one_D[t_max];
				else
					our_model.multiple_K_values[our_model.cur_data_index][s][t] += (our_model.prefix_count_one_D[t_max] - our_model.prefix_count_one_D[t_min]);
			}
		}

		clear_prefix_count_one_D_array(our_model);
	}
}