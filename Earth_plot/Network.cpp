#include "Network.h"

void init_paremeters(int argc, char**argv, model& our_model)
{
	//debug
	/*our_model.network_fileName = (char*)"../../../Datasets/Testing/Testing_network";
	our_model.out_fileName = (char*)"./Results/Testing_M5";
	our_model.method = 5;
	our_model.L = 3;
	our_model.M = 3;
	our_model.T = 3;
	our_model.start_s_tau = 4;
	our_model.incr_s_tau = 2;
	our_model.start_t_tau = 5;
	our_model.incr_t_tau = 2;*/

	/*our_model.network_fileName = (char*)"../../../Datasets/London/London_network";
	our_model.out_fileName = (char*)"./Results/London_M3";
	our_model.method = 3;
	our_model.L = 2;
	our_model.M = 4;
	our_model.T = 4;
	our_model.start_s_tau = 200;
	our_model.incr_s_tau = 200;
	our_model.start_t_tau = 3;
	our_model.incr_t_tau = 2;*/

	our_model.network_fileName = argv[1];
	our_model.out_fileName = argv[2];
	our_model.method = atoi(argv[3]);
	our_model.L = atoi(argv[4]);
	our_model.M = atoi(argv[5]);
	our_model.T = atoi(argv[6]);
	our_model.start_s_tau = atof(argv[7]);
	our_model.incr_s_tau = atof(argv[8]);
	our_model.start_t_tau = atof(argv[9]);
	our_model.incr_t_tau = atof(argv[10]);

	our_model.s_tau_vec = new double[our_model.M];
	our_model.t_tau_vec = new double[our_model.T];
	our_model.multiple_K_values = new double**[our_model.L];
	for (int l = 0; l < our_model.L; l++)
	{
		our_model.multiple_K_values[l] = new double*[our_model.M];
		for (int s = 0; s < our_model.M; s++)
			our_model.multiple_K_values[l][s] = new double[our_model.T];
	}

	for (int l = 0; l < our_model.L; l++)
		for (int s = 0; s < our_model.M; s++)
			for (int t = 0; t < our_model.T; t++)
				our_model.multiple_K_values[l][s][t] = 0;

	for (int s = 0; s < our_model.M; s++)
		our_model.s_tau_vec[s] = our_model.start_s_tau + s * our_model.incr_s_tau;
	for (int t = 0; t < our_model.T; t++)
		our_model.t_tau_vec[t] = our_model.start_t_tau + t * our_model.incr_t_tau;
}

void load_Network(model& our_model)
{
	fstream network_file;
	int num_nodes;
	int num_edges;
	int num_points;
	Point pt;

	network_file.open(our_model.network_fileName);
	if (network_file.is_open() == false)
	{
		cout << "Cannot open this network file." << endl;
		exit(0);
	}

	network_file >> num_edges;
	our_model.num_edges = num_edges;
	our_model.num_points = 0;

	//Process edge_set
	our_model.edge_set = new Edge[num_edges];
	init_PS_memory_edges(our_model, num_edges);
	our_model.min_time = inf; our_model.max_time = -inf;
	for (int e = 0; e < num_edges; e++)
	{
		network_file >> our_model.edge_set[e].n1;
		network_file >> our_model.edge_set[e].n2;
		network_file >> our_model.edge_set[e].length;
		network_file >> num_points;

		our_model.num_points += num_points;
		for (int p = 0; p < num_points; p++)
		{
			pt.edge_index = e;
			network_file >> pt.dist_n1;
			network_file >> pt.dist_n2;
			network_file >> pt.time;

			if (pt.time < our_model.min_time)
				our_model.min_time = pt.time;
			if (pt.time > our_model.max_time)
				our_model.max_time = pt.time;

			our_model.edge_set[e].PS_vec[0].push_back(pt);
		}
	}

	generate_points(our_model);

	//Process the road network
	vector<int> edge_index_vec;
	int degree;
	int edge_index;
	network_file >> num_nodes;
	our_model.num_nodes = num_nodes;
	for (int n = 0; n < num_nodes; n++)
	{
		edge_index_vec.clear();

		network_file >> degree;
		for (int d = 0; d < degree; d++)
		{
			network_file >> edge_index;
			edge_index_vec.push_back(edge_index);
		}

		our_model.Network.push_back(edge_index_vec);
	}

	network_file.close();
}

void output_multiple_results(model& our_model)
{
	fstream output_file;

	output_file.open(our_model.out_fileName, ios::in | ios::out | ios::trunc);
	if (output_file.is_open() == false)
	{
		cout << "Cannot open the output file!" << endl;
		exit(0);
	}

	//The K-function values for the original dataset
	output_file << "The K-function values for the original dataset: " << endl;
	for (int s = 0; s < our_model.M; s++)
	{
		for (int t = 0; t < our_model.T; t++)
			output_file << our_model.multiple_K_values[0][s][t] << " ";
		output_file << endl;
	}

	output_file << endl << "The K-function values for the " << our_model.L - 1 << " random datasets: " << endl;
	for (int l = 1; l < our_model.L; l++)
	{
		for (int s = 0; s < our_model.M; s++)
		{
			for (int t = 0; t < our_model.T; t++)
				output_file << our_model.multiple_K_values[l][s][t] << " ";
			output_file << endl;
		}
		output_file << endl;
	}
}

void init_PS_memory_edges(model& our_model, int num_edges)
{
	vector<Point> temp_vec;

	for (int e = 0; e < num_edges; e++)
		for (int s = 0; s < our_model.L; s++)
			our_model.edge_set[e].PS_vec.push_back(temp_vec);
}

void generate_points(model& our_model)
{
	default_random_engine generator_int;
	default_random_engine generator_double;
	default_random_engine generator_double_time;
	int edge_index;
	Point pt;

	uniform_int_distribution<int> distribution(0, our_model.num_edges - 1);

	for (int l = 1; l < our_model.L; l++)
	{
		for (int i = 0; i < our_model.num_points; i++)
		{
			edge_index = distribution(generator_int);
			vector<Point>& PS_cur = our_model.edge_set[edge_index].PS_vec[l];
			uniform_real_distribution<double> distribution_double(0.0, our_model.edge_set[edge_index].length);
			uniform_real_distribution<double> distribution_double_time(our_model.min_time, our_model.max_time);
			pt.dist_n1 = distribution_double(generator_double);
			pt.dist_n2 = our_model.edge_set[edge_index].length - pt.dist_n1;
			pt.time = distribution_double_time(generator_double_time);
			pt.edge_index = edge_index;

			PS_cur.push_back(pt);
		}
	}
}

void sort_data_points(model& our_model)
{
	//Sorting the data points
	for (int e = 0; e < our_model.num_edges; e++)
		for (int l = 0; l < our_model.L; l++)
			sort(our_model.edge_set[e].PS_vec[l].begin(), our_model.edge_set[e].PS_vec[l].end());
}

void init_tree(Edge& e)
{
	int num_points;
	int L = (int)e.PS_vec.size();
	Range_Tree t;

	for (int l = 0; l < L; l++)
	{
		num_points = (int)e.PS_vec[l].size();
		e.t_vec.push_back(t);
		e.t_vec[l].num_points = num_points;
		e.t_vec[l].Point_Matrix_distance = new double*[num_points];
		for (int i = 0; i < num_points; i++)
		{
			e.t_vec[l].Point_Matrix_distance[i] = new double[2];
			e.t_vec[l].Point_Matrix_distance[i][0] = e.PS_vec[l][i].time;
			e.t_vec[l].Point_Matrix_distance[i][1] = e.PS_vec[l][i].dist_n1;
		}
	}
}

void init_prefix_count_array(model& our_model)
{
	our_model.prefix_count_one_D = new double[2 * our_model.T];
	our_model.prefix_count = new double*[2 * our_model.T];

	for (int t = 0; t < 2 * our_model.T; t++)
		our_model.prefix_count[t] = new double[our_model.M];

	clear_prefix_count_array(our_model, our_model.prefix_count, our_model.M);
	clear_prefix_count_one_D_array(our_model);
}

void construct_prefix_count_array(model& our_model, double**prefix_count, int d_max_index)
{
	double current_value;
	for (int t = 0; t < 2 * our_model.T; t++)
	{
		for (int d = 0; d < d_max_index; d++)
		{
			current_value = prefix_count[t][d];

			if (t != 0 && d != 0)
				prefix_count[t][d] = prefix_count[t - 1][d] + prefix_count[t][d - 1] - prefix_count[t - 1][d - 1] + current_value;
			if (t == 0 && d != 0)
				prefix_count[t][d] = prefix_count[t][d - 1] + current_value;
			if (t != 0 && d == 0)
				prefix_count[t][d] = prefix_count[t - 1][d] + current_value;
		}
	}
}

void clear_prefix_count_array(model& our_model, double**prefix_count, int d_max_index)
{
	for (int t = 0; t < 2 * our_model.T; t++)
		for (int d = 0; d < d_max_index; d++)
			prefix_count[t][d] = 0;
}

void construct_prefix_count_one_D_array(model& our_model)
{
	double current_value;
	for (int t = 0; t < 2 * our_model.T; t++)
	{
		current_value = our_model.prefix_count_one_D[t];
		if (t != 0)
			our_model.prefix_count_one_D[t] = our_model.prefix_count_one_D[t - 1] + current_value;
	}
}

void clear_prefix_count_one_D_array(model& our_model)
{
	for (int t = 0; t < 2 * our_model.T; t++)
		our_model.prefix_count_one_D[t] = 0;
}