#pragma once

#define INF 1e8
#define PI 3.14159265358979323846

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include<iostream>
#include<vector>
#include<algorithm>
#include"lsd.h"
#include <Eigen/Sparse>
#include<Eigen/Dense>
#include<cmath>
#include "GL/glut.h"

typedef cv::Mat CVMat;
typedef cv::Vec3b colorPixel;

typedef Eigen::SparseMatrix<double> SparseMatrixD;//������
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpareseMatrixD_Row;//������
typedef Eigen::VectorXd VectorXd;
typedef Eigen::MatrixXd MatrixXd;
typedef Eigen::Vector2d Vector2d;
typedef Eigen::Matrix2d Matrix2d;

using namespace std;

//ͼ�β�����mesh����
struct Config {//����ǰ����
	int rows;
	int cols;
	int mesh_Rownum;
	int mesh_Colnum;
	int Quad_Rownum;
	int Quad_Colnum;
	double row_spacing;
	double col_spacing;
	Config(int rows, int cols, int mesh_Rownum, int mesh_Colnum) {
		this->rows = rows;
		this->cols = cols;
		this->mesh_Rownum = mesh_Rownum;
		this->mesh_Colnum = mesh_Colnum;
		this->Quad_Colnum = mesh_Colnum - 1;
		this->Quad_Rownum = mesh_Rownum - 1;
		this->row_spacing = double(rows - 1) / (mesh_Rownum - 1);//����Ǿ�׼���㰡...
		this->col_spacing = double(cols - 1) / (mesh_Colnum - 1);
	}
};

//��������
struct Coordinate {
	int row;
	int col;
	bool operator==(const Coordinate& rhs) const {return (row == rhs.row && col == rhs.col);}
	bool operator<(const Coordinate& rhs) const 
	{
		// this operator is used to determine equality, so it must use both x and y
		if (row < rhs.row) {return true;}
		if (row > rhs.row) {return false;}
		return col < rhs.col;
	}
	Coordinate() :row(0), col(0) { };
	Coordinate(double setRow, double setCol) :row(setRow), col(setCol) { };
};

//��������
struct CoordinateDouble {
	double row;
	double col;

	bool operator==(const CoordinateDouble& rhs) const {
		return (row == rhs.row && col == rhs.col);
	}
	bool operator<(const CoordinateDouble& rhs) const {
		// this operator is used to determine equality, so it must use both x and y
		if (row < rhs.row) {return true;}
		if (row > rhs.row) {return false;}
		return col < rhs.col;
	}

	CoordinateDouble operator+(const CoordinateDouble& b)
	{
		CoordinateDouble temp;
		temp.row = row + b.row;
		temp.col = col + b.col;
		return temp;
	}
	CoordinateDouble operator-(const CoordinateDouble& b)
	{
		CoordinateDouble temp;
		temp.row = row - b.row;
		temp.col = col - b.col;
		return temp;
	}

	friend ostream &operator<<(ostream &stream, const CoordinateDouble &p){
		stream << "("<<p.col<<","<<p.row<<")";
		return stream;
	}
	CoordinateDouble():row(0),col(0) { };
	CoordinateDouble(double setRow, double setCol) :row(setRow), col(setCol) { };
};

//�߶�
struct Line_seg {
	double row1, col1;
	double row2, col2;
	Line_seg(double _row1, double _col1, double _row2, double _col2):row1(_row1),row2(_row2),col1(_col1),col2(_col2) {}
	Line_seg() :row1(0), row2(0), col1(0), col2(0) {}
	Line_seg(CoordinateDouble p1, CoordinateDouble p2) :row1(p1.row), row2(p2.row), col1(p1.col), col2(p2.col) {} 
};

CVMat Mask_contour(const CVMat src);
vector<vector<CoordinateDouble>> vector_to_mesh(VectorXd x, Config config);
SpareseMatrixD_Row row_stack(SparseMatrixD origin, SpareseMatrixD_Row diag);
SpareseMatrixD_Row row_stack(SpareseMatrixD_Row origin, SpareseMatrixD_Row diag);
MatrixXd row_stack(MatrixXd mat1, MatrixXd mat2);
MatrixXd col_stack(MatrixXd mat1, MatrixXd mat2);
void DrawLine(CVMat& img, CoordinateDouble coordstart, CoordinateDouble coordend);
void DrawLine(CVMat& img, Line_seg line);
void draw_savemesh(const CVMat src, string filename, vector<vector<CoordinateDouble>> mesh, Config config);
void enlarge_mesh(vector<vector<CoordinateDouble>>& mesh, double enlarge_x, double enlarge_y, Config config);

void compute_scaling(double &sx_avg, double& sy_avg, const vector<vector<CoordinateDouble>> mesh,
	const vector<vector<CoordinateDouble>> outputmesh, const Config config);

GLuint matToTexture(cv::Mat mat, GLenum minFilter = GL_LINEAR, GLenum magFilter = GL_LINEAR, GLenum wrapFilter = GL_CLAMP);