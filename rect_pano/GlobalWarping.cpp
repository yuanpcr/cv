#include"GlobalWarping.h"

//计算逆双线性插值计算出来的权重的对应矩阵T（2，8）维度, P = T * Vq，该矩阵用于乘以Vq=[x0,y0,...,y3,y4] V的顺序都是（左上 右上 左下 右下）
MatrixXd IBilinearWeightsToMatrix(InvBilinearWeights w)
{
	MatrixXd mat(2, 8);
	//a-左上 b-右上 d-左下 c-右下
	double a_w = 1 - w.u - w.v + w.u*w.v;
	double b_w = w.u - w.u*w.v;
	double d_w = w.v - w.u*w.v;
	double c_w = w.u*w.v;
	mat << a_w, 0, b_w, 0, d_w, 0, c_w, 0,
		   0, a_w, 0, b_w, 0, d_w, 0, c_w;
	return mat;
}

//计算a.x*b.y-a.y-b.x
double cross(CoordinateDouble a, CoordinateDouble b){ return a.col*b.row - a.row*b.col; }

//计算逆双线性插值的u,v
InvBilinearWeights get_ibilinear_weights(CoordinateDouble point, Coordinate upperLeftIndices, const vector<vector<CoordinateDouble>>& mesh)
{
	//row为y轴，col为x轴
	//通过mesh和左上顶点的索引计算四个顶点坐标(顺序！)
	CoordinateDouble a = mesh[upperLeftIndices.row][upperLeftIndices.col]; // topLeft
	CoordinateDouble b = mesh[upperLeftIndices.row][upperLeftIndices.col + 1]; // topRight
	CoordinateDouble d = mesh[upperLeftIndices.row + 1][upperLeftIndices.col]; // bottomLeft
	CoordinateDouble c = mesh[upperLeftIndices.row + 1][upperLeftIndices.col + 1]; // bottomRight

	//E、F、G、H
	CoordinateDouble e = b - a;
	CoordinateDouble f = d - a;
	CoordinateDouble g = a - b + c - d;
	CoordinateDouble h = point - a;

	//k2,k1,k0
	double k2 = cross(g, f);
	double k1 = cross(e, f) + cross(h, g);
	double k0 = cross(h, e);

	double u, v;

	if ((int)k2 == 0)
	{
		//平行四边形
		v = -k0 / k1;
		u = (h.col - f.col*v) / (e.col + g.col*v);
	}
	else
	{
		double w = k1 * k1 - 4.0*k0*k2;
		assert(w >= 0.0);//如果point在quad内部，必然不会小于0
		w = sqrt(w);

		double v1 = (-k1 - w) / (2.0*k2);
		double u1 = (h.col - f.col*v1) / (e.col + g.col*v1);

		double v2 = (-k1 + w) / (2.0*k2);
		double u2 = (h.col - f.col*v2) / (e.col + g.col*v2);

		u = u1;
		v = v1;

		if (v<0.0 || v>1.0 || u<0.0 || u>1.0) { u = u2;   v = v2; }
		if (v-0.0 || v>1.0 || u<0.0 || u>1.0) 
		{
			//这里基本都是因为误差导致的，先忽略吧~
			//u = -1.0; 
			//v = -1.0; 
		}
	}
	return InvBilinearWeights(u, v);
}


pair<SpareseMatrixD_Row, VectorXd> get_boundary_mat(vector<vector<CoordinateDouble>> mesh, Config config) {
	
	int rows = config.rows;
	int cols = config.cols;
	int Mesh_Rownum = config.mesh_Rownum;
	int Mesh_Colnum = config.mesh_Colnum;
	int vertexnum = Mesh_Rownum * Mesh_Colnum;
	VectorXd dvec = VectorXd::Zero(vertexnum * 2);//dvec标记需要约束的坐标
	VectorXd B = VectorXd::Zero(vertexnum * 2);//B表示约束的位置
	for (int i = 0; i < vertexnum * 2; i += Mesh_Colnum * 2) {dvec(i) = 1;B(i) = 0;}//左边界，对左侧的顶点x坐标进行约束
	for (int i = Mesh_Colnum * 2 - 2; i < vertexnum * 2; i += Mesh_Colnum * 2) {dvec(i) = 1;B(i) = cols - 1;}//右边界，对右侧顶点的x坐标进行约束
	for (int i = 1; i < 2 * Mesh_Colnum; i += 2) {dvec(i) = 1;B(i) = 0;}//上边界，对上侧顶点的y坐标进行约束
	for (int i = 2 * vertexnum - 2 * Mesh_Colnum + 1; i < vertexnum * 2; i += 2) {dvec(i) = 1;B(i) = rows - 1;}//下边界，对下侧顶点的y坐标进行约束
	SpareseMatrixD_Row diag(dvec.size(), dvec.size());
	for (int i = 0; i < dvec.size(); i++) 
		diag.insert(i, i) = dvec(i);
	diag.makeCompressed();
	return make_pair(diag, B);
};

//获取顶点坐标向量
VectorXd get_vertice(int row, int col, vector<vector<CoordinateDouble>> mesh) 
{   
	//x0,y0,x1,y1...（左上 右上 左下 右下）,row为y轴，col为x轴
	VectorXd Vq = VectorXd::Zero(8);
	CoordinateDouble p0 = mesh[row][col];//左上
	CoordinateDouble p1 = mesh[row][col + 1];//右上
	CoordinateDouble p2 = mesh[row + 1][col];//左下
	CoordinateDouble p3 = mesh[row + 1][col + 1];//右下
	Vq << p0.col, p0.row, p1.col, p1.row, p2.col, p2.row, p3.col, p3.row;
	return Vq;
}


SpareseMatrixD_Row get_shape_mat(vector<vector<CoordinateDouble>> mesh, Config config) 
{
	int Mesh_Rownum = config.mesh_Rownum;
	int Mesh_Colnum = config.mesh_Colnum;
	int Quad_Rownum = config.Quad_Rownum;
	int Quad_Colnum = config.Quad_Rownum;
	SpareseMatrixD_Row Shape_energy(8 * Quad_Rownum*Quad_Colnum, 8 * Quad_Rownum*Quad_Colnum);
	for (int row = 0; row < Quad_Rownum; row++) 
	{
		for (int col = 0; col < Quad_Colnum; col++) 
		{
			CoordinateDouble p0 = mesh[row][col];
			CoordinateDouble p1 = mesh[row][col + 1];
			CoordinateDouble p2 = mesh[row + 1][col];
			CoordinateDouble p3 = mesh[row + 1][col + 1];
			MatrixXd Aq(8, 4);
			Aq << p0.col, -p0.row, 1, 0,
				  p0.row,  p0.col, 0, 1,
				  p1.col, -p1.row, 1, 0,
				  p1.row,  p1.col, 0, 1,
				  p2.col, -p2.row, 1, 0,
				  p2.row,  p2.col, 0, 1,
				  p3.col, -p3.row, 1, 0,
				  p3.row,  p3.col, 0, 1;
			MatrixXd Aq_trans = Aq.transpose(); //Aq^T
			MatrixXd Aq_trans_mul_Aq_reverse = (Aq_trans * Aq).inverse();//(Aq^TAq)^-1
			MatrixXd I = MatrixXd::Identity(8, 8);
			MatrixXd coeff = (Aq*(Aq_trans_mul_Aq_reverse)*Aq_trans - I);

			int left_top_x = (row*Quad_Colnum + col) * 8;
			for (int i = 0; i < 8; i++) {
				for (int j = 0; j < 8; j++) {
					Shape_energy.insert(left_top_x + i, left_top_x + j) = coeff(i, j);
				}
			}
		}
	}
	Shape_energy.makeCompressed();
	return Shape_energy;
}


SpareseMatrixD_Row meshvector_to_Vq_factor(vector<vector<CoordinateDouble>> mesh, Config config) {
	int Mesh_Rownum = config.mesh_Rownum;
	int Mesh_Colnum = config.mesh_Colnum;
	int Quad_Rownum = config.Quad_Rownum;
	int Quad_Colnum = config.Quad_Colnum;
	SpareseMatrixD_Row Q(8 * Quad_Rownum*Quad_Colnum, 2 * Mesh_Rownum*Mesh_Colnum);
	for (int row = 0; row < Quad_Rownum; row++) {
		for (int col = 0; col < Quad_Colnum; col++) {
			int quadid = 8 * (row*Quad_Colnum + col);
			int topleftvertexId = 2 * (row*Mesh_Colnum + col);
			Q.insert(quadid, topleftvertexId) = 1;
			Q.insert(quadid+1, topleftvertexId+1) = 1;
			Q.insert(quadid+2, topleftvertexId+2) = 1;
			Q.insert(quadid+3, topleftvertexId+3) = 1;
			Q.insert(quadid + 4, topleftvertexId + 2 * Mesh_Colnum) = 1;
			Q.insert(quadid + 5, topleftvertexId + 2 * Mesh_Colnum+1) = 1;
			Q.insert(quadid + 6, topleftvertexId + 2 * Mesh_Colnum+2) = 1;
			Q.insert(quadid + 7, topleftvertexId + 2 * Mesh_Colnum+3) = 1;
		}
	}
	Q.makeCompressed();
	return Q;
}

//判断点是否在quad中
bool is_in_quad(CoordinateDouble point, CoordinateDouble topLeft, CoordinateDouble topRight,
	CoordinateDouble bottomLeft, CoordinateDouble bottomRight)
{
	
	//左边的线的右边
	if (topLeft.col == bottomLeft.col)
	{
		if (point.col < topLeft.col)
			return false;
	}
	else 
	{		
		double leftSlope = (topLeft.col - bottomLeft.col) / (topLeft.row - bottomLeft.row);
		double yOnLineX = leftSlope * (point.row - bottomLeft.row) + bottomLeft.col;
		if (point.col < yOnLineX) 
			return false;
	}
	// 右边的线的左边
	if (topRight.col == bottomRight.col)
	{
		if (point.col > topRight.col)
			return false;
	}
	else 
	{
		double rightSlope = (topRight.col - bottomRight.col) / (topRight.row - bottomRight.row);//k=(y2-y1)/(x2-x1)
		double yOnLineX = rightSlope * (point.row - bottomRight.row) + bottomRight.col;//y=k*(x-x1)+y1
		if (point.col > yOnLineX)
			return false;	
	}
	// 上边的线的下边
	if (topLeft.row == topRight.row) 
	{
		if (point.row < topRight.row)
			return false;
	}
	else 
	{
		double topSlope = (topRight.col - topLeft.col) / (topRight.row - topLeft.row);//k=(y2-y1)/(x2-x1)
		double xOnLineY = 1/topSlope * (point.col - topLeft.col) + topLeft.row;//x=1/k*(y-y1)+x1
		if (point.row < xOnLineY)
			return false;
	}
	// 下边线的上边
	if (bottomLeft.row == bottomRight.row) 
	{
		if (point.row > bottomRight.row)
			return false;
	}
	else 
	{

		double bottomSlope = (bottomRight.col - bottomLeft.col) / (bottomRight.row - bottomLeft.row);//k=(y2-y1)/(x2-x1)
		double xOnLineY = 1 / bottomSlope * (point.col - bottomLeft.col) + bottomLeft.row;//x=1/k*(y-y1)+x1
		if (point.row > xOnLineY) 
			return false;
	}
	return true;
}

//判断这个线是否有效（沿着原图边缘是个必须处理掉的情况）
bool line_in_mask(CVMat mask, Line_seg line)
{
	int row1 = round(line.row1), row2 = round(line.row2), 
		col1 = round(line.col1), col2 = round(line.col2);
	//不在内部不算
	if (mask.at<uchar>(row1, col1) != 0 && mask.at<uchar>(row2, col2) != 0)
		return false;
	//边界上不算
	if ((col1 == mask.cols - 1 && col2 == mask.cols - 1) ||(col1 == 0 && col2 == 0))
		return false;
	if ((row1 == mask.rows - 1 && row2 == mask.rows - 1) ||(row1 == 0 && row2 == 0))
		return false;

	//单点边界的时候
	if (row1 == 0 || row1 == mask.rows - 1 || col1 == 0 || col1 == mask.cols - 1)
	{
		try
		{
			if (mask.at<uchar>(row2 + 1, col2) == 255 || mask.at<uchar>(row2 - 1, col2) == 255|| mask.at<uchar>(line.row2, line.col2 + 1) == 255 || mask.at<uchar>(line.row2, line.col2 - 1) == 255)
				return false;
		}
		catch (std::exception) {}
		return true;
	}
	if (row2 == 0 || row2 == mask.rows - 1 || col2 == 0 || col2 == mask.cols - 1)
	{
		try
		{
			if (mask.at<uchar>(row1 + 1, col1) == 255 || mask.at<uchar>(row1 - 1, col1) == 255|| mask.at<uchar>(row1, col1 + 1) == 255 || mask.at<uchar>(row1, col1 - 1) == 255)
				return false;
		}
		catch (std::exception) {}
		return true;
	}

	//一般情况
	try
	{
		if (mask.at<uchar>(row1 + 1, col1) == 255 || mask.at<uchar>(row1 - 1, col1) == 255
			|| mask.at<uchar>(row1, col1 + 1) == 255 || mask.at<uchar>(row1, col1 - 1) == 255)
			return false;
		else
		{
			if (mask.at<uchar>(row2 + 1, col2) == 255 || mask.at<uchar>(row2 - 1, col2) == 255
				|| mask.at<uchar>(line.row2, line.col2 + 1) == 255 || mask.at<uchar>(line.row2, line.col2 - 1) == 255)
				return false;
			else return true;
		}
	}
	catch(std::exception){throw "line 的判断异常";}
}

//判断quad的某一边界与线段的交点
bool does_segment_intersect_line(Line_seg lineSegment, double slope, double intersect, CoordinateDouble& intersectPoint)
{
	/*
    row为y轴，col为x轴，斜率k=(y2-y1)/(x2-x1) 截距b=y1-k*x1=y2-k*x2 从而y=k*x+b
	*/
	double lineSegmentSlope = INF;
	if (lineSegment.col1 != lineSegment.col2)
		lineSegmentSlope = (lineSegment.row2 - lineSegment.row1) / (lineSegment.col2 - lineSegment.col1);//k
	double lineSegmentIntersect = lineSegment.row1 - lineSegmentSlope * lineSegment.col1;//b
	// calculate intersection
	if (lineSegmentSlope == slope) 
	{
		//如果重叠的话！会被另外的情况下捕捉到与另外相交两边的两交点，私以为，可不管先
		if (lineSegmentIntersect == intersect) {return false;}
		else return false;//平行无交点
	}
	//y=k1*x+b1 y=k2*x+b2 → 交点：x0=(b2-b1)/(k1-k2),y0=k1*x+b1
	double intersectX = (intersect - lineSegmentIntersect) / (lineSegmentSlope - slope);
	double intersectY = lineSegmentSlope * intersectX + lineSegmentIntersect;
	// 检查交点是否在线段内，用行坐标或者列坐标一个就够
	if ((intersectY <= lineSegment.row1 && intersectY >= lineSegment.row2) ||(intersectY <= lineSegment.row2 && intersectY >= lineSegment.row1))
	{
		intersectPoint.col = intersectX;
		intersectPoint.row = intersectY;
		return true;
	}
	else return false;

}

//当一条线的两个端点不在同一个quad中的时候，计算该线与这个quad的交点！（该quad用四个坐标表示）
vector<CoordinateDouble> intersections_with_quad(Line_seg lineSegment, CoordinateDouble topLeft,
	CoordinateDouble topRight, CoordinateDouble bottomLeft, CoordinateDouble bottomRight)
{
	/*
	这里可以这样理解：row为y轴，col为x轴，斜率k=(y2-y1)/(x2-x1) 截距b=y1-k*x1=y2-k*x2 从而y=k*x+b 绝了
	下述代码 **Slope=k,**Intersect=b
	*/
	vector<CoordinateDouble> intersections;

	// 左边
	double leftSlope = INF;
	if (topLeft.col != bottomLeft.col)
		leftSlope = (topLeft.row - bottomLeft.row) / (topLeft.col - bottomLeft.col);

	double leftIntersect = topLeft.row - leftSlope * topLeft.col;
	// 检查是否有交点
	CoordinateDouble leftIntersectPoint;
	if (does_segment_intersect_line(lineSegment, leftSlope, leftIntersect, leftIntersectPoint))
		if (leftIntersectPoint.row >= topLeft.row && leftIntersectPoint.row <= bottomLeft.row)
			intersections.push_back(leftIntersectPoint);

	// right
	double rightSlope = INF;
	if (topRight.col != bottomRight.col)
		rightSlope = (topRight.row - bottomRight.row) / (topRight.col - bottomRight.col);
	double rightIntersect = topRight.row - rightSlope * topRight.col;
	// check
	CoordinateDouble rightIntersectPoint;
	if (does_segment_intersect_line(lineSegment, rightSlope, rightIntersect, rightIntersectPoint))
		if (rightIntersectPoint.row >= topRight.row && rightIntersectPoint.row <= bottomRight.row)
			intersections.push_back(rightIntersectPoint);

	// top
	double topSlope = INF;
	if (topLeft.col != topRight.col)
		topSlope = (topRight.row - topLeft.row) / (topRight.col - topLeft.col);
	double topIntersect = topLeft.row - topSlope * topLeft.col;
	// check
	CoordinateDouble topIntersectPoint;
	if (does_segment_intersect_line(lineSegment, topSlope, topIntersect, topIntersectPoint))
		if (topIntersectPoint.col >= topLeft.col && topIntersectPoint.col <= topRight.col)
			intersections.push_back(topIntersectPoint);

	// bottom
	double bottomSlope = INF;
	if (bottomLeft.col != bottomRight.col)
		bottomSlope = (bottomRight.row - bottomLeft.row) / (bottomRight.col - bottomLeft.col);
	double bottomIntersect = bottomLeft.row - bottomSlope * bottomLeft.col;
	// check
	CoordinateDouble bottomIntersectPoint;
	if (does_segment_intersect_line(lineSegment, bottomSlope, bottomIntersect, bottomIntersectPoint))
		if (bottomIntersectPoint.col >= bottomLeft.col && bottomIntersectPoint.col <= bottomRight.col)
			intersections.push_back(bottomIntersectPoint);

	return intersections;
}

//用Quad切割line
vector<vector<vector<Line_seg>>> segment_line_in_quad(vector<Line_seg> lines, vector<vector<CoordinateDouble>> mesh,Config config) 
{
	int QuadnumRow = config.Quad_Rownum;
	int QuadnumCol = config.Quad_Colnum;
	vector<vector<vector<Line_seg>>> quad_line_seg;
	
	for (int row = 0; row < QuadnumRow; row++) {
		vector<vector<Line_seg>> vec_row;
		for (int col = 0; col < QuadnumCol; col++) {
			CoordinateDouble lefttop = mesh[row][col];
			CoordinateDouble righttop = mesh[row][col+1];
			CoordinateDouble leftbottom = mesh[row+1][col];
			CoordinateDouble rightbottom = mesh[row+1][col+1];
			vector<Line_seg> lineInQuad;
			for (int i = 0; i < lines.size(); i++) 
			{
				Line_seg line = lines[i];
				CoordinateDouble point1(line.row1, line.col1);
				CoordinateDouble point2(line.row2, line.col2);
				bool p1InQuad = is_in_quad(point1, lefttop, righttop, leftbottom, rightbottom);
				bool p2InQuad = is_in_quad(point2, lefttop, righttop, leftbottom, rightbottom);
				if (p1InQuad && p2InQuad) 
					lineInQuad.push_back(line);
				else if (p1InQuad) 
				{
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);				
					if (intersections.size() != 0) 
					{
						Line_seg cutLine(point1,intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else if (p2InQuad) 
				{
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() != 0) 
					{
						Line_seg cutLine(point2, intersections[0]);
						lineInQuad.push_back(cutLine);
					}
				}
				else 
				{
					vector<CoordinateDouble> intersections = intersections_with_quad(line, lefttop, righttop, leftbottom, rightbottom);
					if (intersections.size() ==2) {
						Line_seg cutLine(intersections[0], intersections[1]);
						lineInQuad.push_back(cutLine);
					}
				}
			}
			vec_row.push_back(lineInQuad);
		}
		quad_line_seg.push_back(vec_row);
	}

	return quad_line_seg;
}



//以对角的形式存储C*T 之后直接乘V就可以得到结果
SpareseMatrixD_Row block_diag(SpareseMatrixD_Row origin, MatrixXd addin, int QuadID, Config config) 
{
	int cols_total = 8 * config.Quad_Rownum*config.Quad_Colnum;
	SpareseMatrixD_Row res(origin.rows() + addin.rows(), cols_total);
	res.topRows(origin.rows()) = origin;
	
	int lefttop_row = origin.rows();
	int lefttop_col = 8 * QuadID;
	for (int row = 0; row < addin.rows(); row++) {
		for (int col = 0; col < addin.cols(); col++) {
			res.insert(lefttop_row+row, lefttop_col+col) = addin(row,col);
		}
	}
	res.makeCompressed();
	return res;
	
}

//利用lsd的代码进行line的检测
vector<Line_seg> lsd_detect(const CVMat src, CVMat mask)
{
	CVMat temp;
	src.copyTo(temp);
	int rows = temp.rows;
	int cols = temp.cols;
	//转化为灰度图
	CVMat gray_img;
	cv::cvtColor(temp, gray_img, cv::COLOR_BGR2GRAY);
	//构造数组
	double *image = new double[gray_img.rows*gray_img.cols];
	for (int row = 0; row < gray_img.rows; row++)
		for (int col = 0; col < gray_img.cols; col++)
			image[row*gray_img.cols + col] = gray_img.at<uchar>(row, col);
	vector<Line_seg> lines;
	double * out;
	int num_lines;
	//检测线段
	out = lsd(&num_lines, image, gray_img.cols, gray_img.rows);
	for (int i = 0; i < num_lines; i++)
	{
		Line_seg line(out[i * 7 + 1], out[i * 7 + 0], out[i * 7 + 3], out[i * 7 + 2]);
		if (line_in_mask(mask, line))
		{
			lines.push_back(line);
			DrawLine(temp, line);
		}
	}
	cv::imwrite("line_detect.png", temp);
	delete[] image;
	return lines;
}

//把3维的line的vector数组拉平成一维的
void flatten(vector<vector<vector<Line_seg>>> lineSeg, vector<Line_seg>& line_vec, Config config) {
	int numQuadRow = config.Quad_Rownum;
	int numQuadCol = config.Quad_Colnum;
	int linetmpnum = -1;
	for (int row = 0; row < numQuadRow; row++) {
		for (int col = 0; col < numQuadCol; col++) {
			vector<Line_seg>linesegInquad=lineSeg[row][col];
			Coordinate topleft(row,col);
			for (int k = 0; k < linesegInquad.size(); k++) {
				line_vec.push_back(lineSeg[row][col][k]);
			}
		}
	}
}
void get_bilinearVec(vector<vector<vector<Line_seg>>>lineSeg, Config config, 
	const vector<vector<CoordinateDouble>>& mesh, vector<pair<MatrixXd, MatrixXd>>& BilinearVec)
{
	int numQuadRow = config.Quad_Rownum;
	int numQuadCol = config.Quad_Colnum;
	int linetmpnum = -1;
	for (int row = 0; row < numQuadRow; row++) 
	{
		for (int col = 0; col < numQuadCol; col++) 
		{
			vector<Line_seg>linesegInquad = lineSeg[row][col];
			Coordinate topleft(row, col);
			for (int k = 0; k < linesegInquad.size(); k++) 
			{
				linetmpnum++;
				Line_seg line = linesegInquad[k];
				CoordinateDouble linestart(line.row1, line.col1);
				CoordinateDouble lineend(line.row2, line.col2);
				InvBilinearWeights startWeight = get_ibilinear_weights(linestart, topleft, mesh);
				MatrixXd start_W_mat = IBilinearWeightsToMatrix(startWeight);
				InvBilinearWeights endWeight = get_ibilinear_weights(lineend, topleft, mesh);
				MatrixXd end_W_mat = IBilinearWeightsToMatrix(endWeight);
				BilinearVec.push_back(make_pair(start_W_mat, end_W_mat));
			}
		}
	}
}
//初始化计算三维line的vector数组及线段角度
vector<vector<vector<Line_seg>>> init_line_seg(const CVMat src, const CVMat mask,Config config, vector < Line_seg > &lineSeg_flatten,
	vector<vector<CoordinateDouble>> mesh, vector<pair<int, double>>&bin_id_theta,vector<double> &rotate_theta ) 
{
	double thetaPerbin = PI / 49;//将PI分成了50份
	//找到图中所有的线段
	vector<Line_seg> lines = lsd_detect(src, mask);
	//用Quad分割线段
	vector<vector<vector<Line_seg>>> lineSeg = segment_line_in_quad(lines, mesh, config);
	//用一维vector储存线段
	flatten(lineSeg,lineSeg_flatten, config);
	
	for (int i = 0; i < lineSeg_flatten.size(); i++) {
		Line_seg line = lineSeg_flatten[i];
		double theta = atan((line.row1 - line.row2) / (line.col1 - line.col2));//-PI/2~PI/2
		int lineSegmentBucket =(int) round((theta + PI / 2) / thetaPerbin);
		assert(lineSegmentBucket < 50);
		bin_id_theta.push_back(make_pair(lineSegmentBucket, theta));
		rotate_theta.push_back(0);
	}
	return lineSeg;
}

//获取计算线的能量的矩阵 对角上存储（2，8）—— C*T，稀疏存储
SpareseMatrixD_Row get_line_mat(CVMat mask,vector<double>rotate_theta, vector<vector<vector<Line_seg>>> lineSeg,
	vector<pair<MatrixXd,MatrixXd>> BilinearVec,Config config,int &linenum)
{
	int linetmpnum = -1;
	int QuadnumRow = config.Quad_Rownum;
	int QuadnumCol = config.Quad_Colnum;

	SpareseMatrixD_Row energy_line;
	for (int row = 0; row < QuadnumRow; row++) 
	{
		for (int col = 0; col < QuadnumCol; col++) 
		{
			vector<Line_seg> linesegInquad = lineSeg[row][col];
			int QuadID = row * QuadnumCol + col;
			if (linesegInquad.size() == 0) continue;
			else {
				Coordinate topleft(row, col);
				MatrixXd C_row_stack(0,8);
				for (int k = 0; k < linesegInquad.size(); k++) 
				{
					linetmpnum++;
					Line_seg line = linesegInquad[k];
					double theta = rotate_theta[linetmpnum];	
					Matrix2d R;//旋转矩阵R
					R << cos(theta), -sin(theta), 
						 sin(theta), cos(theta);
					MatrixXd ehat(2,1);
					ehat << line.col1 - line.col2, line.row1-line.row2;
					MatrixXd tmp = (ehat.transpose()*ehat).inverse();
					Matrix2d I = Matrix2d::Identity();
					MatrixXd C = R * ehat*tmp*(ehat.transpose())*(R.transpose()) - I;
					MatrixXd start_W_mat = BilinearVec[linetmpnum].first; MatrixXd end_W_mat = BilinearVec[linetmpnum].second;
					MatrixXd CT = C * (start_W_mat - end_W_mat);
					C_row_stack = row_stack(C_row_stack, CT);
				}
				energy_line = block_diag(energy_line, C_row_stack, QuadID,config);
			}
		}
	}
	linenum = linetmpnum;
	return energy_line;
}


//用最蠢得办法来修补。。
CVMat fill_missing_pixel(CVMat &img, const CVMat mask)
{
	//row为y轴，col为x轴
	assert(img.rows == mask.rows);
	assert(img.cols == mask.cols);
	CVMat mask_2;
	int size_erode = 9;
	CVMat element = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(size_erode, size_erode));
	cv::erode(mask, mask_2, element);//255是非图的部分，，腐蚀掉一些这个
	/*cv::imshow("temp", mask_2);
	cv::imshow("temp2", mask);
	cv::waitKey(0);	*/
	for (int row = 0; row < mask.rows; row++)
	{
		for (int col = 0; col < mask.cols; col++)
		{
			if (mask.at<uchar>(row, col) == 255 && mask_2.at<uchar>(row,col)==0)
			{
				for (int i = 0; i < size_erode; i++)
				{
					int temp_y = row - 2 + i / size_erode;
					int temp_x = col - 2 + i % size_erode;
					if (temp_y >= 0 && temp_y <= mask.rows&&temp_x >= 0 && temp_x <= mask.cols)
					{
						if (mask.at<uchar>(temp_y,temp_x) == 0)
						{
							img.at<cv::Vec3b>(row, col) = img.at<cv::Vec3b>(temp_y,temp_x);
							break;
						}
					}
				}
			}
		}
	}
	return mask_2;
}