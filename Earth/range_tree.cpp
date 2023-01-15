#include "range_tree.h"

void Range_Tree::build_range_tree_Recur(Range_Node*node, bool is_base)
{
	if (node->index_u < node->index_l) //Base case (no data point)
	{
		cout << "Empty node!" << endl;
		exit(0);
	}

	if (node->index_u == node->index_l) //Base case (leaf node)
	{
		if (is_base == true)
			build_asso_structure(node);
		return;
	}

	Range_Node*left_child;
	Range_Node*right_child;

	node->left_child = new Range_Node();
	node->right_child = new Range_Node();

	left_child = node->left_child;
	right_child = node->right_child;

	left_child->index_l = node->index_l;
	left_child->index_u = node->index_m;
	left_child->index_m = (int)floor((left_child->index_l + left_child->index_u) / 2.0);
	right_child->index_l = node->index_m + 1;
	right_child->index_u = node->index_u;
	right_child->index_m = (int)floor((right_child->index_l + right_child->index_u) / 2.0);

	if (is_base == true)
	{
		left_child->value_l = Point_Matrix_distance[left_child->index_l][1];
		left_child->value_u = Point_Matrix_distance[left_child->index_u][1];
		left_child->value_m = Point_Matrix_distance[left_child->index_m][1];
		right_child->value_l = Point_Matrix_distance[right_child->index_l][1];
		right_child->value_u = Point_Matrix_distance[right_child->index_u][1];
		right_child->value_m = Point_Matrix_distance[right_child->index_m][1];

		build_range_tree_Recur(left_child, true);
		build_range_tree_Recur(right_child, true);
		build_asso_structure(node);
	}
	else
	{
		left_child->value_l = Point_Matrix_time[left_child->index_l][0];
		left_child->value_u = Point_Matrix_time[left_child->index_u][0];
		left_child->value_m = Point_Matrix_time[left_child->index_m][0];
		right_child->value_l = Point_Matrix_time[right_child->index_l][0];
		right_child->value_u = Point_Matrix_time[right_child->index_u][0];
		right_child->value_m = Point_Matrix_time[right_child->index_m][0];

		build_range_tree_Recur(left_child, false);
		build_range_tree_Recur(right_child, false);
	}
}

void Range_Tree::build_range_tree()
{
	rootNode = new Range_Node();

	if (num_points == 0)
		return;

	rootNode->index_l = 0;
	rootNode->index_u = num_points - 1;
	rootNode->index_m = (int)floor((rootNode->index_l + rootNode->index_u) / 2.0);
	rootNode->value_l = Point_Matrix_distance[0][1];
	rootNode->value_u = Point_Matrix_distance[num_points - 1][1];
	rootNode->value_m = Point_Matrix_distance[rootNode->index_m][1];

	build_range_tree_Recur(rootNode, true);
}

void Range_Tree::merge_time_Matrix(double**lc_Point_Matrix_time, double**rc_Point_Matrix_time, double**& Point_Matrix_time, int lc_num_points, int rc_num_points, int num_points)
{
	if (lc_num_points + rc_num_points != num_points)
	{
		cout << "ERROR!" << endl;
		exit(0);
	}

	int cur_l = 0;
	int cur_r = 0;
	Point_Matrix_time = new double*[num_points];
	for (int i = 0; i < num_points; i++)
		Point_Matrix_time[i] = new double[2];

	//code here
	for (int i = 0; i < num_points; i++)
	{
		if (cur_l == lc_num_points)
		{
			Point_Matrix_time[i][0] = rc_Point_Matrix_time[cur_r][0];
			Point_Matrix_time[i][1] = rc_Point_Matrix_time[cur_r][1];
			cur_r++;
			continue;
		}

		if (cur_r == rc_num_points)
		{
			Point_Matrix_time[i][0] = lc_Point_Matrix_time[cur_l][0];
			Point_Matrix_time[i][1] = lc_Point_Matrix_time[cur_l][1];
			cur_l++;
			continue;
		}

		if (lc_Point_Matrix_time[cur_l][0] <= rc_Point_Matrix_time[cur_r][0])
		{
			Point_Matrix_time[i][0] = lc_Point_Matrix_time[cur_l][0];
			Point_Matrix_time[i][1] = lc_Point_Matrix_time[cur_l][1];
			cur_l++;
		}
		else
		{
			Point_Matrix_time[i][0] = rc_Point_Matrix_time[cur_r][0];
			Point_Matrix_time[i][1] = rc_Point_Matrix_time[cur_r][1];
			cur_r++;
		}
	}
}

void Range_Tree::build_asso_structure(Range_Node*node)
{
	Range_Tree& tt = node->time_tree;
	tt.rootNode = new Range_Node();
	tt.rootNode->index_l = 0;
	tt.rootNode->index_u = node->index_u - node->index_l;
	//tt.rootNode->index_u = num_points - 1;
	tt.rootNode->index_m = (int)floor((tt.rootNode->index_l + tt.rootNode->index_u) / 2.0);
	tt.num_points = node->index_u - node->index_l + 1;

	if (node->index_u == node->index_l) //Base case (leaf node)
	{
		//tt.num_points = 1;
		tt.Point_Matrix_time = new double*[1];
		tt.Point_Matrix_time[0] = new double[2];
		tt.Point_Matrix_time[0][0] = Point_Matrix_distance[node->index_l][0];
		tt.Point_Matrix_time[0][1] = Point_Matrix_distance[node->index_l][1];
		tt.rootNode->value_l = tt.Point_Matrix_time[0][0];
		tt.rootNode->value_u = tt.Point_Matrix_time[0][0];
		tt.rootNode->value_m = tt.Point_Matrix_time[0][0];
		return;
	}

	//tt.num_points = num_points;
	merge_time_Matrix(node->left_child->time_tree.Point_Matrix_time, node->right_child->time_tree.Point_Matrix_time, tt.Point_Matrix_time,
		node->left_child->time_tree.num_points, node->right_child->time_tree.num_points, tt.num_points);
	tt.rootNode->value_l = tt.Point_Matrix_time[0][0];
	tt.rootNode->value_u = tt.Point_Matrix_time[tt.num_points - 1][0];
	tt.rootNode->value_m = tt.Point_Matrix_time[tt.rootNode->index_m][0];

	tt.build_range_tree_Recur(tt.rootNode, false);
}

void Range_Tree::find_split_node(Range_Tree& t, double l, double u, Range_Node*& node)
{
	node = t.rootNode;
	while (node->index_u > node->index_l && (u < node->value_m || l > node->value_m))
	{
		if (u < node->value_m)
			node = node->left_child;
		else
			node = node->right_child;
	}
}

double Range_Tree::one_D_range_counting(Range_Tree& t, double l, double u)
{
	double count = 0;
	Range_Node*splitNode;
	Range_Node*node;

	find_split_node(t, l, u, splitNode);
	if (splitNode->index_u == splitNode->index_l) //leaf node
	{
		if (splitNode->value_m >= l && splitNode->value_m <= u)
			return 1;
		else
			return 0;
	}

	//Left child
	node = splitNode->left_child;
	while (node->index_u > node->index_l)
	{
		if (l <= node->value_m)
		{
			count += (node->right_child->index_u - node->right_child->index_l + 1);
			node = node->left_child;
		}
		else
			node = node->right_child;
	}

	if (node->value_m >= l && node->value_m <= u) //leaf node
		count++;

	//Right child
	node = splitNode->right_child;
	while (node->index_u > node->index_l)
	{
		if (u > node->value_m)
		{
			count += (node->left_child->index_u - node->left_child->index_l + 1);
			node = node->right_child;
		}
		else
			node = node->left_child;
	}

	if (node->value_m >= l && node->value_m <= u) //leaf node
		count++;

	return count;
}

double Range_Tree::two_D_range_counting(Range_Tree& t, double dist_L, double dist_U, double time_L, double time_U)
{
	double count = 0;
	Range_Node*splitNode;
	Range_Node*node;

	find_split_node(t, dist_L, dist_U, splitNode);
	if (splitNode->index_u == splitNode->index_l) //leaf node
	{
		if ((t.Point_Matrix_distance[splitNode->index_l][0] >= time_L && t.Point_Matrix_distance[splitNode->index_l][0] <= time_U) &&
			(splitNode->value_m >= dist_L && splitNode->value_m <= dist_U))
			return 1;
		else
			return 0;
	}

	node = splitNode->left_child;
	while (node->index_u > node->index_l)
	{
		if (dist_L <= node->value_m)
		{
			count += one_D_range_counting(node->right_child->time_tree, time_L, time_U);
			node = node->left_child;
		}
		else
			node = node->right_child;
	}

	if ((t.Point_Matrix_distance[node->index_l][0] >= time_L && t.Point_Matrix_distance[node->index_l][0] <= time_U)
		&& (t.Point_Matrix_distance[node->index_l][1] >= dist_L && t.Point_Matrix_distance[node->index_l][1] <= dist_U))
		count++;

	//Right child
	node = splitNode->right_child;
	while (node->index_u > node->index_l)
	{
		if (dist_U > node->value_m)
		{
			count += one_D_range_counting(node->left_child->time_tree, time_L, time_U);
			node = node->right_child;
		}
		else
			node = node->left_child;
	}

	if ((t.Point_Matrix_distance[node->index_l][0] >= time_L && t.Point_Matrix_distance[node->index_l][0] <= time_U)
		&& (t.Point_Matrix_distance[node->index_l][1] >= dist_L && t.Point_Matrix_distance[node->index_l][1] <= dist_U))
		count++;

	return count;
}