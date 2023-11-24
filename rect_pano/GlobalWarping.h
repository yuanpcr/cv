#pragma once

#include"base.h"

struct Line_rotate {
	Vector2d pstart = Vector2d::Zero();
	Vector2d pend = Vector2d::Zero();
	double angle = 0;
	Line_rotate(Vector2d pstart, Vector2d pend, double angle) 
	{
		this->pstart = pstart;
		this->pend = pend;
		this->angle = angle;
	}
};

struct InvBilinearWeights 
{
	double u;
	double v;
	InvBilinearWeights(double u_, double v_) { u = u_; v = v_; }
};

SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config);
SpareseMatrixD_Row meshvector_to_Vq_factor(vector<vector<CoordinateDouble>> mesh, Config config);
pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(vector<vector<CoordinateDouble>> mesh, Config config);
VectorXd get_vertice(int row, int col, vector<vector<CoordinateDouble>> mesh);
vector<vector<vector<Line_seg>>> init_line_seg(const CVMat src, CVMat mask, Config config, vector < Line_seg > &lineSeg_flatten,vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>&id_theta, vector<double> &rotate_theta);
SpareseMatrixD_Row get_line_mat(CVMat mask, vector<double>rotate_theta, vector<vector<vector<Line_seg>>> lineSeg,
	vector<pair<MatrixXd, MatrixXd>> BilinearVec, Config config, int& linenum);
void get_bilinearVec(vector<vector<vector<Line_seg>>>lineSeg, Config config,
	const vector<vector<CoordinateDouble>>& mesh, vector<pair<MatrixXd, MatrixXd>>& BilinearVec);
double cross(CoordinateDouble a, CoordinateDouble b);
InvBilinearWeights get_ibilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, const vector<vector<CoordinateDouble>>& mesh);
MatrixXd IBilinearWeightsToMatrix(InvBilinearWeights w);
bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight,
	CoordinateDouble bottomLeft, CoordinateDouble bottomRight);
CVMat fill_missing_pixel(CVMat &img, const CVMat mask);