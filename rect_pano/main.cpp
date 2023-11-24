#include"base.h"
#include"LocalWarping.h"
#include"GlobalWarping.h"

CVMat img;
vector<vector<CoordinateDouble>> outputmesh;
vector<vector<CoordinateDouble>> mesh;
double sx_avg = 1, sy_avg = 1;
bool flag_display = true;

GLuint matToTexture(cv::Mat mat, GLenum minFilter, GLenum magFilter, GLenum wrapFilter)
{	
	// Generate a number for our textureID's unique handle
	GLuint textureID;
	glGenTextures(1, &textureID);
	// Bind to our texture handle
	glBindTexture(GL_TEXTURE_2D, textureID);
	// Catch silly-mistake texture interpolation method for magnification
	if (magFilter == GL_LINEAR_MIPMAP_LINEAR ||
		magFilter == GL_LINEAR_MIPMAP_NEAREST ||
		magFilter == GL_NEAREST_MIPMAP_LINEAR ||
		magFilter == GL_NEAREST_MIPMAP_NEAREST)
	{
		magFilter = GL_LINEAR;
	}
	// Set texture interpolation methods for minification and magnification
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, minFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, magFilter);
	// Set texture clamping method
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapFilter);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapFilter);
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	GLenum inputColourFormat = GL_BGR_EXT;
	if (mat.channels() == 1)
	{
		inputColourFormat = GL_LUMINANCE;
	}
	// Create the texture
	glTexImage2D(
		GL_TEXTURE_2D,     // Type of texture
		0,                 // Pyramid level (for mip-mapping) - 0 is the top level
		GL_RGB,            // Internal colour format to convert to
		mat.cols,          // Image width  i.e. 640 for Kinect in standard mode
		mat.rows,          // Image height i.e. 480 for Kinect in standard mode
		0,                 // Border width in pixels (can either be 1 or 0)
		inputColourFormat, // Input image format (i.e. GL_RGB, GL_RGBA, GL_BGR etc.)
		GL_UNSIGNED_BYTE,  // Image data type
		mat.ptr());        // The actual image data itself

	return textureID;
}
void display() 
{	
	// 清除屏幕
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GLuint texGround = matToTexture(img);
    glViewport(0, 0, (GLsizei)img.cols, (GLsizei)img.rows);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, img.cols, img.rows, 0);

	glEnable(GL_TEXTURE_2D);    // 启用纹理	
	glBindTexture(GL_TEXTURE_2D, texGround);	
	if (flag_display)
	{
		for (int row = 0; row < 20; row++)
		{
			for (int col = 0; col < 20; col++)
			{
				CoordinateDouble &coord = outputmesh[row][col];
				CoordinateDouble &localcoord = mesh[row][col];
				localcoord.row /= img.rows;
				localcoord.col /= img.cols;
			}
		}
		flag_display = false;
	}
	for (int row = 0; row < 19; row++) {
		for (int col = 0; col < 19; col++) {
			CoordinateDouble local_left_top = mesh[row][col];
			CoordinateDouble local_right_top = mesh[row][col + 1];
			CoordinateDouble local_left_bottom = mesh[row + 1][col];
			CoordinateDouble local_right_bottom = mesh[row + 1][col + 1];

			CoordinateDouble global_left_top = outputmesh[row][col];
			CoordinateDouble global_right_top = outputmesh[row][col + 1];
			CoordinateDouble global_left_bottom = outputmesh[row + 1][col];
			CoordinateDouble global_right_bottom = outputmesh[row + 1][col + 1];
			
			glBegin(GL_QUADS);
			glTexCoord2d(local_left_top.col, local_left_top.row); glVertex2d(global_left_top.col, global_left_top.row);
			glTexCoord2d(local_right_top.col, local_right_top.row); glVertex2d(global_right_top.col,  global_right_top.row);
			glTexCoord2d(local_right_bottom.col, local_right_bottom.row); glVertex2d(global_right_bottom.col,  global_right_bottom.row);
			glTexCoord2d(local_left_bottom.col, local_left_bottom.row);	glVertex2d(global_left_bottom.col,  global_left_bottom.row);		
			glEnd();			
		}
	}
	glDisable(GL_TEXTURE_2D);
	glutSwapBuffers();

}

int main(int argc, char* argv[]) 
{
    img = cv::imread("D:\\nk2\\exp\\rect_pano\\data\\p2.jpg");
	double Time = (double)cv::getTickCount();
	double fator_scale = 1.0;
	CVMat scaled_img;
	cv::resize(img, scaled_img, cv::Size(0, 0), 1.0 / fator_scale, 1.0 / fator_scale);
	Config config(img.rows, img.cols, 20, 20);
	CVMat mask = Mask_contour(scaled_img);
	CVMat tmpmask;
	mask.copyTo(tmpmask);

	CVMat warpped_img = CVMat::zeros(scaled_img.size(), CV_8UC3);
	vector<vector<Coordinate>> displacementMap = Local_warp(scaled_img, warpped_img, tmpmask);
	mesh = get_rectangle_mesh(config);
	mesh_warp_back(mesh, displacementMap, config);
	cout << "mesh warp back 完成" << endl;
	double temp_time;



	SpareseMatrixD_Row Q = meshvector_to_Vq_factor(mesh, config);
	SpareseMatrixD_Row shape_energy = get_shape_mat(mesh, config);
	cout << "获得保形能量参数" << endl;
	pair<SpareseMatrixD_Row, VectorXd> pair_dvec_B = get_boundary_mat(mesh, config);
	cout << "获得边界约束参数" << endl;
	vector<pair<int, double>>binId_originTheta;
	vector <Line_seg> line_flatten;
	vector<double> rotate_theta;
	vector<pair<MatrixXd, MatrixXd>> BilinearVec;
	vector<vector<vector<Line_seg>>> LineSeg = init_line_seg(scaled_img, mask, config, line_flatten, mesh, binId_originTheta, rotate_theta);
	get_bilinearVec(LineSeg,config,mesh, BilinearVec);




	//10 iteration
	for (int iter = 1; iter <= 10; iter++) 
	{
		cout << "迭代次数：" << iter << endl;
		int Nl = 0;

		SpareseMatrixD_Row line_energy = get_line_mat(mask, rotate_theta, LineSeg, BilinearVec, config, Nl);
		cout << "获得直线约束能量参数 "<< endl;

		double Nq = config.Quad_Rownum*config.Quad_Colnum;
		double lambdaB = 1e8;
		double lambdaL = 100;
		SpareseMatrixD_Row shape = (1 / sqrt(Nq))*(shape_energy*Q);
		SpareseMatrixD_Row boundary = sqrt(lambdaB) * pair_dvec_B.first;
		SpareseMatrixD_Row line = sqrt((lambdaL / Nl))*(line_energy*Q);
		SpareseMatrixD_Row K = row_stack(shape, line);
		SpareseMatrixD_Row K2 = row_stack(K, boundary);
		VectorXd B = pair_dvec_B.second;
		VectorXd b = VectorXd::Zero(K2.rows());
		b.tail(B.size()) = sqrt(lambdaB) * B;
		SparseMatrixD K2_trans = K2.transpose();
		MatrixXd temp = K2_trans * K2;
		temp = temp.inverse();
		VectorXd x;	
		x = temp * K2_trans*b;

	    //update theta
		outputmesh = vector_to_mesh(x, config);
		int tmplinenum = -1;
		VectorXd thetagroup = VectorXd::Zero(50);
		VectorXd thetagroupcnt = VectorXd::Zero(50);
		for (int row = 0; row < config.Quad_Rownum; row++) 
		{
			for (int col = 0; col < config.Quad_Colnum; col++)
			{
				vector<Line_seg> linesegInquad = LineSeg[row][col];
				if (linesegInquad.size() == 0) continue;
				else {
					VectorXd Vq = get_vertice(row, col, outputmesh);
					for (int k = 0; k < linesegInquad.size(); k++) {
						tmplinenum++;
						pair<MatrixXd, MatrixXd> Bstartend = BilinearVec[tmplinenum];					
						MatrixXd start_W_mat = Bstartend.first;
						MatrixXd end_W_mat = Bstartend.second;
						Vector2d newstart = start_W_mat * Vq;
						Vector2d newend = end_W_mat * Vq;

						double theta = atan((newstart(1) - newend(1)) / (newstart(0) - newend(0)));
						double deltatheta = theta - binId_originTheta[tmplinenum].second;
						if (isnan(binId_originTheta[tmplinenum].second) || isnan(deltatheta))continue;
						if (deltatheta > (PI / 2))deltatheta -= PI;					
						if (deltatheta < (-PI / 2))deltatheta += PI;
						thetagroup(binId_originTheta[tmplinenum].first) += deltatheta;
						thetagroupcnt(binId_originTheta[tmplinenum].first) += 1;
					}
				}
			}
		}
		for (int i = 0; i < thetagroup.size(); i++)
			thetagroup(i) /= thetagroupcnt(i);
		for (int i = 0; i < rotate_theta.size(); i++)
			rotate_theta[i] = thetagroup[binId_originTheta[i].first];
	}


	enlarge_mesh(mesh, fator_scale, fator_scale, config);//因为开始缩小了
	enlarge_mesh(outputmesh, fator_scale, fator_scale, config);

	CVMat erode_mask;
	draw_savemesh(img,"quad1.png", mesh, config);//保存，注释掉了画图的代码
	draw_savemesh(img,"quad2.png", outputmesh, config);
	
	//erode_mask=fill_missing_pixel(img, mask);//修补像素
	//erode_mask=fill_missing_pixel(img, erode_mask);//修补像素
	
	compute_scaling(sx_avg, sy_avg, mesh, outputmesh, config);//计算缩放因子

    enlarge_mesh(mesh, 1/sx_avg, 1/sy_avg, config);
	enlarge_mesh(outputmesh, 1 / sx_avg, 1 / sy_avg, config);

	cv::resize(img, img, cv::Size(0, 0), 1 / sx_avg, 1 / sy_avg);

	draw_savemesh(img, "quad1_scaling.png", mesh, config);//保存，注释掉了画图的代码
	draw_savemesh(img, "quad2_scaling.png", outputmesh, config);

	temp_time = (double)cv::getTickCount() - Time;
	std::printf("run time = %f s\n", temp_time / (cv::getTickFrequency()));
	//glut
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);	
	glutInitWindowSize(img.cols, img.rows);
	glutInitWindowPosition(100, 100);	
	glutCreateWindow("Panoramic_image");
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(&display);   //注册函数 

	temp_time = (double)cv::getTickCount() - Time;
	std::printf("run time = %f s\n", temp_time / (cv::getTickFrequency()));
	glutMainLoop();
	
	return 0;
}