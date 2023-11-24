#pragma once

#include"base.h"

enum Border {
	TOP = 0,
	BOTTOM = 1,
	LEFT = 2,
	RIGHT = 3
};

enum SeamDirection {
	SEAM_VERTICAL = 0,
	SEAM_HORIZONTAL = 1
};

vector<vector<Coordinate>> Local_warp(const CVMat src, CVMat& wrap_img, CVMat mask);
CVMat Insert_local_seam(CVMat src, CVMat& seam_img, CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend);
int* Get_local_seam_improved(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end);
vector<vector<Coordinate>> Get_Local_warp_displacement(CVMat& warp_img, CVMat mask);
pair<int, int> find_longest_boundary_segment(CVMat mask, Border& direction);
void mesh_warp_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config);
vector<vector<CoordinateDouble>> get_rectangle_mesh( Config config);

