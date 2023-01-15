#pragma once
#ifndef RANGE_TREE_H
#define RANGE_TREE_H

#include "Library.h"

class Range_Node;
class Range_Tree;

class Range_Tree
{
public:
	//x-coordinate: time 
	//y-coordinate: distance
	double**Point_Matrix_time; //sorted according to time
	double**Point_Matrix_distance;  //sorted according to distance
	Range_Node*rootNode;
	int num_points;

	//Build the range tree
	void build_range_tree_Recur(Range_Node*node, bool is_base);
	void build_range_tree();
	void merge_time_Matrix(double**lc_Point_Matrix_time, double**rc_Point_Matrix_time, double**& Point_Matrix_time, int lc_num_points, int rc_num_points, int num_points);
	//void build_asso_structure_Recur(Range_Node*asso_node);
	void build_asso_structure(Range_Node*node);

	//Search the range tree
	void find_split_node(Range_Tree& t, double l, double u, Range_Node*& node);
	double one_D_range_counting(Range_Tree& t, double l, double u);
	double two_D_range_counting(Range_Tree& t, double dist_L, double dist_U, double time_L, double time_U);
};

class Range_Node
{
public:
	int index_l; //smallest index
	int index_u; //largest index
	int index_m; //middle index
	double value_l; //value of the smallest index
	double value_u; //value of the largest index
	double value_m; //value of the middle index
	Range_Node*left_child;
	Range_Node*right_child;
	Range_Tree time_tree;
};

#endif