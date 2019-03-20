#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <glm.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>


#define POINT_DISTANCE 15.0f


using namespace std;
using namespace Eigen;

// *********�������************
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

	// �е������
	bool operator==(const Point &p)
	{
		return (this->position == p.position);
	}
public:
	// ���λ��
	glm::vec2 position;
	// �����ⲿ�����ڲ� 0/�ڲ� 1/�ⲿ
	bool out_or_in;		

};
//********�ⲿ����**********
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
	// 0��ʾ���ǹؼ���
	bool is_control;

	// ���߻���ֱ�ߵĲ�����Ϣ
	float t;			
};

// ********�ڲ�����***********
class IPoint : public Point
{
public:
	IPoint() {};
	~IPoint() {};

public:
	// �ڽӵ����������õ�Ȩ��
	vector<pair<int, float>> adj_points;

};

// **********�߶���***********
class Line
{
public:
	Line() { is_dart_line = false; };
	Line(bool _is_curve,bool _is_dart_line) :is_curve(_is_curve),is_dart_line(_is_dart_line) {};
	~Line() {};

	// ���ݿ��Ƶ�����ڲ���
	void calculatePoints();

	// ����һ�����Ƿ����߶��ڲ���ĳЩ��ܽ�
	bool isCloseToLine(glm::vec2 point);
public:
	// �߶��м�ĵ�
	vector<OPoint> point_list;

	// �߶��а����������
	//static float point_distance;
	int point_size;

	// �߶εĶ˵�
	pair<OPoint, OPoint> end_points;

	// �߶ε��м���Ƶ�
	OPoint curve_control_point;

	// 0��ʾֱ�ߣ�1��ʾ����
	bool is_curve;		
	
	// �Ƿ���dart���߶�
	bool is_dart_line;

};

// ***********��������*********
class Contour
{
public:
	Contour() {};
	// ���ļ��ж�ȡ
	Contour(string filepath);
	~Contour() {};

	struct Triangle {
		int v0;
		int v1;
		int v2;
	};

	struct Dart
	{
		int l0;
		int l1;
		int l2;
		int l3;
	};

public:
	// Contour�ĳ�ʼ��
	void Init();

	// ��ȡ�����߿��Ƶ����������Ϣ
	void getOPints();

	// ��ȡ���е��ϵ������
	void getCoefficientMatrix();

	// MVC���Σ��ڲ�������ⲿ��仯
	void MVCDeform();

	// handle_indexΪ�仯�ĵ㣬new_positonΪ�仯֮��ĵ�
	void Deform(int handle_index, glm::vec2 new_position);	

	// ����һ���㣬���ظõ����ڵ��߶α�ţ��Լ������߶��е�λ�� -1/��� 0/�м���Ƶ� 1/�յ�
	vector<pair<int, int>> controlToLines(int control_index);

	// ����һ�����Ƶ㣬���ظõ���һ���㼯�е�����
	int getControlPointIndex(OPoint p);

	// ���������߽��е��������ǻ��õ�mesh���������
	void delaunayTrianglation();

	// ��Ƭ�������ǻ�
	void reMesh();

	// ��OBJ�ж�ȡ���������
	void readOBJ(string path);

	// �������
	void clear();

	// ����dart
	void createDart(glm::vec2 dart_begin,glm::vec2 dart_end);

	//������������ true��ʾ��Ҫ�������ǻ� false��ʾ����Ҫ�������ǻ�
	bool evaluateMesh();



public:
	// �����������߶�
	vector<Line>  lines;	

	// �������е�dart,��¼line�ı��
	vector<Dart> darts;

	// ���������п��Ƶ�
	vector<OPoint> control_points;

	// ����������������
	vector<OPoint> contour_points;

	// ���������ǻ������е�Լ��
	vector<vector<OPoint>> segments;

	// ���еĵ�����ڲ�����ⲿ��
	vector<Point> vertices;

	// ���е�������
	vector<Triangle> triangles;

	// ���е��ϵ������
	MatrixXd coe_matrix;
	LLT<MatrixXd> llt;

	// ��������ܵ�Ķ�Ӧ��
	vector<int> out_all_map;
	
};

