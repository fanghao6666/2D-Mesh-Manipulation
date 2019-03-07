#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <glm.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;

// *********点基类类************
class Point
{
public:
	Point() {};
	Point(float x, float y, bool out_in)
	{
		position = glm::vec2(x, y);
		out_or_in = out_in;
	}
	~Point() {};

	// 判等运算符
	bool operator==(const Point &p)
	{
		return (this->position == p.position);
	}
public:
	// 点的位置
	glm::vec2 position;
	// 点在外部还是内部 0/内部 1/外部
	bool out_or_in;		

};
//********外部点类**********
class OPoint : public Point
{
public:
	OPoint() {};
	OPoint(float x, float y,bool _is_control,float _t)
	{
		position = glm::vec2(x, y);
		is_control = _is_control;
		t = _t;
	}
	~OPoint() {};

public:
	// 0表示不是关键点
	bool is_control;

	// 曲线或者直线的参数信息
	float t;			
};

// ********内部点类***********
class IPoint : public Point
{
public:
	IPoint() {};
	~IPoint() {};

public:
	// 邻接点的索引和求得的权重
	vector<pair<int, float>> adj_points;

};

// **********线段类***********
class Line
{
public:
	Line() {};
	~Line() {};

	// 根据控制点计算内部点
	void calculatePoints();

	// 计算一个点是否在线段内部离某些点很近
	bool isCloseToLine(glm::vec2 point);
public:
	// 线段中间的点
	vector<OPoint> point_list;

	// 线段中包含点的数量
	//static float point_distance;
	int point_size;

	// 线段的端点
	pair<OPoint, OPoint> end_points;

	// 线段的中间控制点
	OPoint curve_control_point;

	// 0表示直线，1表示曲线
	bool is_curve;		

};

// ***********轮廓线类*********
class Contour
{
public:
	Contour() {};
	// 从文件中读取
	Contour(string filepath);
	~Contour() {};

	struct Triangle {
		int v0;
		int v1;
		int v2;
	};

public:
	// Contour的初始化
	void Init();

	// 获取轮廓线控制点和轮廓线信息
	void getOPints();

	// 获取所有点的系数矩阵
	void getCoefficientMatrix();

	// MVC变形，内部点根据外部点变化
	void MVCDeform();

	// handle_index为变化的点，new_positon为变化之后的点
	void Deform(int handle_index, glm::vec2 new_position);	

	// 传入一个点，返回该点所在的线段编号，以及点在线段中的位置 -1/起点 0/中间控制点 1/终点
	vector<pair<int, int>> controlToLines(int control_index);

	// 传入一个控制点，返回该点在一个点集中的索引
	int getControlPointIndex(OPoint p);

	// 根据轮廓线进行德劳累三角化得到mesh点和三角形
	void delaunayTrianglation();

	// 面片重新三角化
	void reMesh();

	// 从OBJ中读取点和三角形
	void readOBJ(string path);

	// 清除操作
	void clear();



public:
	// 轮廓的所有线段
	vector<Line>  lines;	

	// 轮廓的所有控制点
	vector<OPoint> control_points;

	// 轮廓的所有轮廓点
	vector<OPoint> contour_points;

	// 所有的点包括内部点和外部点
	vector<Point> vertices;

	// 所有的三角形
	vector<Triangle> triangles;

	// 所有点的系数矩阵
	MatrixXd coe_matrix;
	LLT<MatrixXd> llt;

	// 轮廓点和总点的对应表
	vector<int> out_all_map;
	
};

