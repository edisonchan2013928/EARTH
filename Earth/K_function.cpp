#include "K_function.h"

void call_RQS(model& our_model)
{
	//Preprocessing
	int e_index;
	int node_u, node_v;
	our_model.K_value = 0;
	Edge*edge_set = our_model.edge_set;
	Dijkstra_preprocess(our_model);

	for (int wide_e_index = 0; wide_e_index < our_model.num_edges; wide_e_index++)
	{
		our_model.cur_edge_index = wide_e_index;
		vector<Point>& PS = edge_set[wide_e_index].PS;
		for (int p = 0; p < (int)PS.size(); p++)
		{
			our_model.cur_point_index = p;
			Dijkstra(our_model, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
			for (int e = 0; e < (int)our_model.access_edge_list_first.size(); e++)
			{
				e_index = our_model.access_edge_list_first[e];

				if (wide_e_index != e_index)
				{
					node_u = edge_set[e_index].n1;
					node_v = edge_set[e_index].n2;
					if (our_model.sp_node_vec_first[node_u].cur_sp_value > our_model.s_threshold &&
						our_model.sp_node_vec_first[node_v].cur_sp_value > our_model.s_threshold)
						continue;

					//call point-edge_count function
					our_model.K_value += point_edge_count(our_model, e_index, node_u, node_v);
				}
				else
					our_model.K_value += point_edge_count_same_edge(our_model, p, e_index);
			}

			clear_Dijkstra(edge_set, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		}
	}
}

void call_SPS(model& our_model)
{
	//Preprocessing
	int e_index;
	int node_u, node_v;
	our_model.K_value = 0;
	Edge*edge_set = our_model.edge_set;
	//int edges_totalNum;
	int num_edges_two_sets;
	Dijkstra_preprocess(our_model);

	if (our_model.method == 3) //EAR method (preprocessing step)
	{
		for (int e = 0; e < our_model.num_edges; e++)
		{
			init_tree(our_model.edge_set[e]);
			our_model.edge_set[e].t.build_range_tree();
		}
		
		//cout << our_model.edge_set[0].t.two_D_range_counting(our_model.edge_set[0].t, 3, 5, 0.5, 2.5) << endl;
		//cout << our_model.edge_set[0].t.two_D_range_counting(our_model.edge_set[0].t, 0, 10, 0, 10) << endl;

		//return; //check the range-tree
	}

	for (int wide_e_index = 0; wide_e_index < our_model.num_edges; wide_e_index++)
	{
		our_model.cur_edge_index = wide_e_index;
		vector<Point>& PS = edge_set[wide_e_index].PS;
		if (PS.size() == 0)
			continue;

		our_model.is_first = true;
		Dijkstra(our_model, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		our_model.is_first = false;
		Dijkstra(our_model, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
		num_edges_two_sets = (int)our_model.access_edge_list_first.size() + (int)our_model.access_edge_list_second.size();

		for (int e = 0; e < num_edges_two_sets; e++)
		{
			if (e < (int)our_model.access_edge_list_first.size())
			{
				e_index = our_model.access_edge_list_first[e];
				our_model.detect_DoubleScan_edge_list[e_index] = true;
			}
			else
			{
				e_index = our_model.access_edge_list_second[e - our_model.access_edge_list_first.size()];
				if (our_model.detect_DoubleScan_edge_list[e_index] == true)
					continue;
				else
					our_model.detect_DoubleScan_edge_list[e_index] = true;
			}

			if (wide_e_index != e_index)
			{
				node_u = edge_set[e_index].n1;
				node_v = edge_set[e_index].n2;
				if (our_model.sp_node_vec_first[node_u].cur_sp_value > our_model.s_threshold &&
					our_model.sp_node_vec_first[node_v].cur_sp_value > our_model.s_threshold &&
					our_model.sp_node_vec_second[node_u].cur_sp_value > our_model.s_threshold &&
					our_model.sp_node_vec_second[node_v].cur_sp_value > our_model.s_threshold)
					continue;

				for (int p = 0; p < (int)PS.size(); p++)
				{
					our_model.cur_point_index = p;
					our_model.K_value += point_edge_count(our_model, e_index, node_u, node_v);
				}
			}
			else
			{
				for (int p = 0; p < (int)PS.size(); p++)
					our_model.K_value += point_edge_count_same_edge(our_model, p, e_index);
			}
		}

		clear_double_scan_edgeList(our_model.access_edge_list_first, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_first, our_model.access_edge_list_first, our_model.is_scan_edge_list_first);
		clear_double_scan_edgeList(our_model.access_edge_list_second, our_model.detect_DoubleScan_edge_list);
		clear_Dijkstra(edge_set, our_model.sp_node_vec_second, our_model.access_edge_list_second, our_model.is_scan_edge_list_second);
	}
}

void compute_K_function(model& our_model)
{
	double run_time;

	sort_data_points(our_model);

	auto start_s = chrono::high_resolution_clock::now();

	//Step 1: Compute shortest path distance
	//Step 2: Compute C_P^{(wide_e,e)}
	if (our_model.method == 1) //RQS
		call_RQS(our_model);
	if (our_model.method == 2) //SPS
		call_SPS(our_model);
	if (our_model.method == 3) //EAR
		call_SPS(our_model);

	auto end_s = chrono::high_resolution_clock::now();

	run_time = (chrono::duration_cast<chrono::nanoseconds>(end_s - start_s).count()) / 1000000000.0;
	cout << "method " << our_model.method << " (run time):" << run_time << endl;
}