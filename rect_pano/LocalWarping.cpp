#include"LocalWarping.h"
//给sort用的比较函数
bool cmp(const pair<int, float> a, const pair<int, float> b) {
	return a.second<b.second;
}

//初始化
void init_displacement(vector<vector<Coordinate>>& displacement, int rows, int cols) {
	for (int row = 0; row < rows; row++) 
	{
		vector<Coordinate> displacement_row;
		for (int col = 0; col < cols; col++) 
		{
			Coordinate c;
			displacement_row.push_back(c);
		}
		displacement.push_back(displacement_row);
	}
}

//warp_back操作
void mesh_warp_back(vector<vector<CoordinateDouble>>& mesh, vector<vector<Coordinate>> displacementMap, Config config) {
	int meshrow_num = config.mesh_Rownum;
	int meshcol_num = config.mesh_Colnum;
	for (int meshrow = 0; meshrow< meshrow_num; meshrow++) {
		for (int meshcol = 0; meshcol < meshcol_num; meshcol++) {		
			CoordinateDouble& vetex_coor = mesh[meshrow][meshcol];
			Coordinate vertexDisplacement;
			vertexDisplacement = displacementMap[(int)floor(vetex_coor.row)][(int)floor(vetex_coor.col)];//取出结点的位移
			vetex_coor.row += vertexDisplacement.row;//计算出结点横坐标的原位置
			vetex_coor.col += vertexDisplacement.col;//计算出结点纵坐标的原位置
		}
	}
}

//找到具有最长的边界
pair<int, int> find_longest_boundary_segment(CVMat mask, Border& direction) 
{

	int rows = mask.rows;
	int cols = mask.cols;

	int maxLength = 0;
	int max_startIndex = 0;
	int max_endIndex = 0;

	int tmp_maxLength, tmp_startIndex, tmp_endIndex;
	
	//left	
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	for (int row = 0; row < rows; row++) 
	{
		if (mask.at<uchar>(row, 0) == 0) //遇到有像素的点
		{
			if (tmp_maxLength > maxLength) 
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = LEFT;
			}
			tmp_startIndex = tmp_endIndex = row + 1;
			tmp_maxLength = 0;
		}
		else if (row == rows - 1)//最后一个像素
		{
			if (mask.at<uchar>(row, 0) != 0) //判断最后一个像素
			{
				tmp_endIndex++;
				tmp_maxLength++;
			}
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = LEFT;
			}
		}
		else 
		{//无像素，继续计数
			tmp_endIndex++;
			tmp_maxLength++;
		}
	}
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	for (int row = 0; row < rows; row++)
	{
		if (mask.at<uchar>(row, cols-1) == 0)
		{
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = RIGHT;
			}
			tmp_startIndex = tmp_endIndex = row + 1;
			tmp_maxLength = 0;
		}
		else if (row == rows - 1)
		{
			if (mask.at<uchar>(row, cols-1) != 0) 
			{
				tmp_endIndex++;
				tmp_maxLength++;
			}
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = RIGHT;
			}
		}
		else
		{
			tmp_endIndex++;
			tmp_maxLength++;
		}
	}
	//top
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	for (int col = 0; col < cols; col++)
	{
		if (mask.at<uchar>(0, col) == 0)
		{
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = TOP;
			}
			tmp_startIndex = tmp_endIndex = col + 1;//从下一个开始的
			tmp_maxLength = 0;
		}
		else if (col == cols - 1)
		{
			if (mask.at<uchar>(0, cols - 1) != 0) //判断最后一个元素的
			{
				tmp_endIndex++;
				tmp_maxLength++;
			}
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = TOP;
			}
		}
		else
		{//计数
			tmp_endIndex++;
			tmp_maxLength++;
		}
	}
	//bottom
	tmp_maxLength = tmp_startIndex = tmp_endIndex = 0;
	for (int col = 0; col < cols; col++)
	{
		if (mask.at<uchar>(rows-1, col) == 0)
		{
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = BOTTOM;
			}
			tmp_startIndex = tmp_endIndex = col + 1;//从下一个开始的
			tmp_maxLength = 0;
		}
		else if (col == cols - 1)
		{
			if (mask.at<uchar>(rows-1, cols - 1) != 0) //判断最后一个元素的
			{
				tmp_endIndex++;
				tmp_maxLength++;
			}
			if (tmp_maxLength > maxLength)
			{
				maxLength = tmp_maxLength;
				max_startIndex = tmp_startIndex;
				max_endIndex = tmp_endIndex;
				direction = BOTTOM;
			}
		}
		else
		{//继续计数
			tmp_endIndex++;
			tmp_maxLength++;
		}
	}
	if(maxLength==0)return make_pair(0, 0);
	return make_pair(max_startIndex, max_endIndex - 1);
}

//就是把边界画个红线呢
CVMat Show_longest_border(CVMat src, pair<int, int>begin_end, Border direction) {
	CVMat tmpsrc;
	src.copyTo(tmpsrc);
	int rows = src.rows;
	int cols = src.cols;
	switch (direction) {
	case LEFT:
		for (int row = begin_end.first; row < begin_end.second; row++)
			tmpsrc.at<colorPixel>(row, 0) = colorPixel(0, 0, 255);
		break;
	case RIGHT:
		for (int row = begin_end.first; row < begin_end.second; row++)
			tmpsrc.at<colorPixel>(row, cols - 1) = colorPixel(0, 0, 255);
		break;
	case TOP:
		for (int col = begin_end.first; col < begin_end.second; col++)
			tmpsrc.at<colorPixel>(0, col) = colorPixel(0, 0, 255);
		break;
	case BOTTOM:
		for (int col = begin_end.first; col < begin_end.second; col++)
			tmpsrc.at<colorPixel>(rows - 1, col) = colorPixel(0, 0, 255);
		break;
	default:
		break;
	}

	cv::namedWindow("Border", cv::WINDOW_AUTOSIZE);
	cv::imshow("Border", tmpsrc);
	cv::waitKey(0);

	return tmpsrc;
}

//执行local warp
vector<vector<Coordinate>> Local_warp(const CVMat src, CVMat& warp_img, CVMat mask) 
{
	src.copyTo(warp_img);
	vector<vector<Coordinate>> displacementMap = Get_Local_warp_displacement(warp_img, mask);
	return displacementMap;
}

//通过xxx诸多过程最终获得local warp的位移矩阵,并且通过引用将warp_img插成矩形
vector<vector<Coordinate>> Get_Local_warp_displacement(CVMat& warp_img, CVMat mask) {

	CVMat seam_img;
	warp_img.copyTo(seam_img);

	int rows = warp_img.rows;
	int cols = warp_img.cols;

	vector<vector<Coordinate>> displacementMap;//这个用来存储最终位移矩阵
	init_displacement(displacementMap, rows, cols);

	while (true) {
		Border direction;
		pair<int, int> begin_end = find_longest_boundary_segment(mask, direction);
		if (begin_end.first == begin_end.second)
		{
			cv::imwrite("local_warping.png", warp_img);
			cv::imwrite("seam_img.png", seam_img);
			return displacementMap;
		}
		else
		{
			bool shift_to_end = false;
			SeamDirection seamdirection;
			switch (direction)
			{
			case LEFT:seamdirection = SEAM_VERTICAL;shift_to_end = false;break;
			case RIGHT:seamdirection = SEAM_VERTICAL;shift_to_end = true;break;
			case TOP:seamdirection = SEAM_HORIZONTAL;shift_to_end = false;break;
			case BOTTOM:seamdirection = SEAM_HORIZONTAL;shift_to_end = true;break;
			default:throw "direction值异常";break;}
			int* seam = Get_local_seam_improved(warp_img, mask, seamdirection, begin_end);
			warp_img = Insert_local_seam(warp_img, seam_img, mask, seam, seamdirection, begin_end, shift_to_end);
			//更新位移矩阵
			if (seamdirection == SEAM_VERTICAL)
			{
				for (int row = begin_end.first; row <= begin_end.second; row++)
				{
					int local_row = row - begin_end.first;
					if (shift_to_end)
					{
						for (int col = cols-1; col>seam[local_row]; col--)
						{
							displacementMap[row][col].col = displacementMap[row][col-1].col-1;
							displacementMap[row][col].row = displacementMap[row][col - 1].row;
						}
					}
					else
					{
						for(int col=0;col<seam[local_row];col++)
						{
							displacementMap[row][col].col = displacementMap[row][col+1].col+1;
							displacementMap[row][col].row = displacementMap[row][col + 1].row;
						}
					}
				}
			}
			else if (seamdirection == SEAM_HORIZONTAL)
			{
				for (int col = begin_end.first; col <= begin_end.second; col++)
				{
					int local_col = col - begin_end.first;
					if (shift_to_end)
					{
						for (int row = rows - 1; row > seam[local_col]; row--)
						{
							displacementMap[row][col].row = displacementMap[row - 1][col].row -1;
							displacementMap[row][col].col = displacementMap[row-1][col].col;
						}
					}
					else
					{
						for(int row=0;row<seam[local_col];row++)
						{
							displacementMap[row][col].row =displacementMap[row+1][col].row+1;
							displacementMap[row][col].col = displacementMap[row+1][col ].col;
						}
					}
				}
			}
		}
		
	}
	throw "seam carving 异常终止";
}
//对seam执行插入操作
CVMat Insert_local_seam(CVMat src, CVMat& seam_img,CVMat& mask, int* seam, SeamDirection seamdirection, pair<int, int> begin_end, bool shiftToend) 
{
	if (seamdirection == SEAM_HORIZONTAL) 
	{
		cv::transpose(src, src);
		cv::transpose(seam_img, seam_img);
		cv::transpose(mask, mask);
	}

	CVMat resimg;
	src.copyTo(resimg);

	int begin = begin_end.first;
	int end = begin_end.second;

	int rows = src.rows;
	int cols = src.cols;
	//处理seam边界线段一侧的像素
	for (int row = begin; row <= end; row++) {
		int local_row = row - begin;		
		if (!shiftToend) 
			for (int col = 0; col < seam[local_row]; col++)
			{
				resimg.at<colorPixel>(row, col) = src.at<colorPixel>(row, col + 1);
				seam_img.at<colorPixel>(row, col) = seam_img.at<colorPixel>(row, col + 1);
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col + 1);
			}
		else
			for (int col = cols - 1; col > seam[local_row]; col--)
			{
				resimg.at<colorPixel>(row, col) = src.at<colorPixel>(row, col - 1);
				seam_img.at<colorPixel>(row, col) = seam_img.at<colorPixel>(row, col - 1);
				mask.at<uchar>(row, col) = mask.at<uchar>(row, col - 1);
			}
		
		//seam上的值为两边取平均，边界情况单独处理
		mask.at<uchar>(row, seam[local_row]) = 0;
		if (seam[local_row] == 0)
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] + 1);
		else if (seam[local_row] == cols - 1)
			resimg.at<colorPixel>(row, seam[local_row]) = src.at<colorPixel>(row, seam[local_row] - 1);
		else 
			resimg.at<colorPixel>(row, seam[local_row]) = 0.5* src.at<colorPixel>(row, seam[local_row] + 1) + 0.5* src.at<colorPixel>(row, seam[local_row] - 1);
		seam_img.at<colorPixel>(row, seam[local_row]) = colorPixel(0, 255, 0);
	}

	if (seamdirection == SEAM_HORIZONTAL)
	{
		cv::transpose(resimg, resimg);
		cv::transpose(seam_img, seam_img);
		cv::transpose(mask, mask);
	}
	
	/*cv::namedWindow("insert_seam", cv::WINDOW_AUTOSIZE);
	cv::imshow("insert_seam", resimg);
	cv::waitKey(0);*/

	return resimg;
}
int* Get_local_seam_improved(CVMat src, CVMat mask, SeamDirection seamdirection, pair<int, int> begin_end)
{
	//都转换成处理垂直的情况
	if (seamdirection == SEAM_HORIZONTAL)
	{
		cv::transpose(src, src);
		cv::transpose(mask, mask);
	}
	
	int rows = src.rows;
	int cols = src.cols;

	int row_start = begin_end.first;
	int row_end = begin_end.second;
	int range = row_end - row_start + 1;

	int col_start = 0;
	int col_end = cols - 1;

	CVMat local_img = src(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));
	CVMat local_mask = mask(cv::Range::Range(row_start, row_end + 1), cv::Range::Range(col_start, col_end + 1));

	CVMat gray,U, U2, L, L2, R, R2;
	cv::cvtColor(local_img, gray, cv::COLOR_BGR2GRAY);
	U = gray(cv::Range(1, range), cv::Range(0, cols));
	U2 = gray(cv::Range(0, 1), cv::Range(0, cols));
	cv::vconcat(U, U2, U);
	L = gray(cv::Range(0, range), cv::Range(1, cols));
	L2 = gray(cv::Range(0, range), cv::Range(0, 1));
	cv::hconcat(L, L2, L);
	R = gray(cv::Range(0, range), cv::Range(0, cols - 1));
	R2 = gray(cv::Range(0, range), cv::Range(cols - 1, cols));
	cv::hconcat(R2,R, R);

	CVMat cU,cL,cR;
	cU = abs(R - L);
	cL = abs(U - L) + cU;
	cR = abs(U - R) + cU;
	CVMat M(range, cols, CV_64FC1, cv::Scalar(0));
	
	//非图中部分能量置为无穷大
	for (int row = 0; row < range; row++)
		for (int col = col_start; col <= col_end; col++)
			if ((int)local_mask.at<uchar>(row, col) == 255) M.at<double>(row, col) = INF;
	for (int i = 1; i < range; i++)
		for (int j = col_start; j <= col_end; j++)
		{
			double MU = M.at<double>(i - 1, j) + cU.at<uchar>(i, j);
			double ML = M.at<double>(i - 1, j - 1)+ cL.at<uchar>(i, j);
			double MR = M.at<double>(i - 1, j + 1)+ cR.at<uchar>(i, j);
			if ((int)local_mask.at<uchar>(i, j) == 255)continue;
			M.at<double>(i, j) = min(MU, min(ML, MR));
		}
	


	//DP计算最小seam
	CVMat tmpenergy;//这个指的是seam的累计能量
	M.copyTo(tmpenergy);

	//找最小的前向能量
	vector<pair<int, double>> last_row;
	for (int col = col_start; col <= col_end; col++)
		last_row.push_back(make_pair(col, tmpenergy.at<double>(range - 1, col)));
	sort(last_row.begin(), last_row.end(), cmp);
	int* seam = new int[range];
	seam[range - 1] = last_row[0].first;

	//回溯
	for (int row = range - 2; row >= 0; row--)
	{
		if (seam[row + 1] == col_start)
		{
			if (tmpenergy.at<double>(row, seam[row + 1] + 1) < tmpenergy.at<double>(row, seam[row + 1]))
				seam[row] = seam[row + 1] + 1;
			else 
				seam[row] = seam[row + 1];
		}
		else if (seam[row + 1] == col_end)
		{
			if (tmpenergy.at<double>(row, seam[row + 1] - 1) < tmpenergy.at<double>(row, seam[row + 1]))
				seam[row] = seam[row + 1] - 1;
			else
				seam[row] = seam[row + 1];
		}
		else
		{
			double min_energy = min(tmpenergy.at<double>(row, seam[row + 1] - 1), min(tmpenergy.at<double>(row, seam[row + 1]), tmpenergy.at<double>(row, seam[row + 1] + 1)));
			if (min_energy == tmpenergy.at<double>(row, seam[row + 1] - 1))
				seam[row] = seam[row + 1] - 1;
			else if (min_energy == tmpenergy.at<double>(row, seam[row + 1] + 1))
				seam[row] = seam[row + 1] + 1;
			else 
				seam[row] = seam[row + 1];
		}
		
	}

	for (int row = 0; row < range; row++)
		local_img.at<colorPixel>(row, seam[row]) = colorPixel(255, 0, 0);

	/*cv::namedWindow("seam", cv::WINDOW_AUTOSIZE);
	cv::imshow("seam", local_img);
	cv::waitKey(0);*/

	return seam;//对于水平方向的情况不需要转置回来，插入过程会转置以对应
}

//计算mesh顶点坐标（0~row-1和0~col-1的额）
vector<vector<CoordinateDouble>> get_rectangle_mesh(Config config) {
	int rows = config.rows;
	int cols = config.cols;
	int meshnum_row = config.mesh_Rownum;
	int meshnum_col = config.mesh_Colnum;
	double row_per_mesh = config.row_spacing;
	double col_per_mesh = config.col_spacing;
	vector<vector<CoordinateDouble>> mesh;
	for (int row_mesh = 0; row_mesh < meshnum_row; row_mesh++) {
		vector<CoordinateDouble> meshrow;
		for (int col_mesh = 0; col_mesh < meshnum_col; col_mesh++) {
			CoordinateDouble coord;
			coord.row = row_mesh * row_per_mesh;
			coord.col = col_mesh * col_per_mesh;
			meshrow.push_back(coord);
		}
		mesh.push_back(meshrow);
	}
	return mesh;
}