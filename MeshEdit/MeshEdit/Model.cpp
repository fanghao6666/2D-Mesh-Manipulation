#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <glm.hpp>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Fade_2D.h>

#include "Model.h"

using namespace std;
using namespace Eigen;
using namespace GEOM_FADE25D;


//******************LINE*************************//

// ͨ���߶εĶ˵�����߶ε��ڲ������
void Line::calculatePoints()
{
	point_list.clear();

	if (is_curve == false)
	{
		for (int i = 1; i <= point_size; ++i)
		{
			// ÿ�������Ĳ���t
			float _t = float(point_size + 1 - i) / float(point_size + 1);
			// ÿ��������λ��
			glm::vec2 insert_position = end_points.first.position * _t + end_points.second.position*(1 - _t);
			point_list.push_back(OPoint(insert_position.x, insert_position.y, false, _t));
		}
	}
	else
	{
		for (int i = 1; i <= point_size; ++i)
		{
			// ÿ�������Ĳ���t
			float _t = float(point_size + 1 - i) / float(point_size + 1);
			// ÿ��������λ��
			glm::vec2 insert_position = _t * _t * end_points.first.position + 2 * _t * (1 - _t) * curve_control_point.position + (1 - _t) * (1 - _t) * end_points.second.position;
			point_list.push_back(OPoint(insert_position.x, insert_position.y, false, _t));
		}
	}
}

// �ж�һ�����Ƿ����߶κܽ�
bool Line::isCloseToLine(glm::vec2 point)
{
	for (int i = 1; i < point_list.size()-1; ++i)
	{
		if (glm::distance(point_list[i].position, point) < (glm::distance(end_points.first.position, end_points.second.position) / (point_list.size() + 1)))
		{
			return true;
		}
	}
	return false;
}

//********************CONTOUR*******************//

// Contour��ʼ���������ȶ����߶�Ȼ���ٴ���
void Contour::Init()
{
	// ��ȡ���Ƶ���Ϣ����������Ϣ
	getOPints();

	// ��ȡOBJ�ļ��õ����������
	//readOBJ("C:/Users/ai/Desktop/mesh.obj");

	// ���������ǻ�
	delaunayTrianglation();

	// ��ȡϵ������
	getCoefficientMatrix();
}

// ��ȡ���Ƶ����������Ϣ
void Contour::getOPints()
{
	// �����һ�εĵ�
	control_points.clear();
	contour_points.clear();

	// ��ȡ���º�ĵ�
	for (int i = 0; i < lines.size(); ++i)
	{
		if (lines[i].is_curve == true)
		{
			// ��ȡ���Ƶ�
			control_points.push_back(lines[i].end_points.first);
			control_points.push_back(lines[i].curve_control_point);
			contour_points.push_back(lines[i].end_points.first);

			// ��ȡ�ڲ���
			lines[i].calculatePoints();
			contour_points.insert(contour_points.end(), lines[i].point_list.begin(), lines[i].point_list.end());
		}
		else if (lines[i].is_curve == false)
		{
			// ��ȡ���Ƶ�
			control_points.push_back(lines[i].end_points.first);
			contour_points.push_back(lines[i].end_points.first);

			// ��ȡ�ڲ���
			lines[i].calculatePoints();
			contour_points.insert(contour_points.end(), lines[i].point_list.begin(), lines[i].point_list.end());
		}
	}
}

// ��ȡ���е��ϵ������
void Contour::getCoefficientMatrix()
{
	// ��ʼ��һ��n*n�ľ���nΪ���е�Ĵ�С
	coe_matrix = MatrixXd::Identity(vertices.size(), vertices.size());

	// ����ϵ������
	for (int i = 0; i < triangles.size(); ++i)
	{
		// ��ȡ�����ε��������������
		int v0 = triangles[i].v0;
		int v1 = triangles[i].v1;
		int v2 = triangles[i].v2;
		// �����һ����
		if (vertices[v0].out_or_in == false)
		{
			// ������v0��ļн�
			glm::vec3 v01 = glm::vec3(vertices[v1].position - vertices[v0].position, 0.0f);
			glm::vec3 v02 = glm::vec3(vertices[v2].position - vertices[v0].position, 0.0f);

			// tan(��/2) = 1-cos(��)/ sin(��) = sin(��)/1 + cos(��)
			float cos0 = glm::dot(glm::normalize(v01), glm::normalize(v02));
			float sin0 = sqrtf(1.0f - powf(cos0, 2));
			float tan0 = (1 - cos0) / sin0;

			float w01 = tan0 / glm::length(v01);
			float w02 = tan0 / glm::length(v02);

			// ����ϵ��ֵ
			coe_matrix(v0, v1) += (-w01);
			coe_matrix(v0, v2) += (-w02);
		}
		// ����ڶ�����
		if (vertices[v1].out_or_in == false)
		{
			// ������v1��ļн�
			glm::vec3 v10 = glm::vec3(vertices[v0].position - vertices[v1].position, 0.0f);
			glm::vec3 v12 = glm::vec3(vertices[v2].position - vertices[v1].position, 0.0f);

			// tan(��/2) = 1-cos(��)/ sin(��) = sin(��)/1 + cos(��)
			float cos1 = glm::dot(glm::normalize(v10), glm::normalize(v12));
			float sin1 = sqrtf(1.0f - powf(cos1, 2));
			float tan1 = (1 - cos1) / sin1;

			float w10 = tan1 / glm::length(v10);
			float w12 = tan1 / glm::length(v12);

			// ����ϵ��ֵ
			coe_matrix(v1, v0) += (-w10);
			coe_matrix(v1, v2) += (-w12);
		}
		// �����������
		if (vertices[v2].out_or_in == false)
		{
			// ������v2��ļн�
			glm::vec3 v20 = glm::vec3(vertices[v0].position - vertices[v2].position, 0.0f);
			glm::vec3 v21 = glm::vec3(vertices[v1].position - vertices[v2].position, 0.0f);

			// tan(��/2) = 1-cos(��)/ sin(��) = sin(��)/1 + cos(��)
			float cos2 = glm::dot(glm::normalize(v20), glm::normalize(v21));
			float sin2 = sqrtf(1.0f - powf(cos2, 2));
			float tan2 = (1 - cos2) / sin2;

			float w20 = tan2 / glm::length(v20);
			float w21 = tan2 / glm::length(v21);

			// ����ϵ��ֵ
			coe_matrix(v2, v0) += (-w20);
			coe_matrix(v2, v1) += (-w21);
		}
	}

	//��һ������
	MatrixXd row_matrix = coe_matrix.rowwise().sum();
	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = 0; j < vertices.size(); ++j)
		{
			if ((coe_matrix(i, j) != 1) && (coe_matrix(i, j) != 0))
			{
				coe_matrix(i, j) /= -(row_matrix(i, 0) - 1.0f);
			}
		}
	}

	llt = (coe_matrix.transpose()*coe_matrix).llt();

}


// ����һ���㣬���ظõ����ڵ��߶α�ţ��Լ������߶��е�λ�� -1/��� 0/�м���Ƶ� 1/�յ�
vector<pair<int, int>> Contour::controlToLines(int control_index)
{
	vector<pair<int, int>> result;

	for (int i = 0; i < lines.size(); ++i)
	{
		if (lines[i].end_points.first == control_points[control_index])
		{
			result.push_back(pair<int, int>(i, -1));
		}
		if (lines[i].curve_control_point == control_points[control_index])
		{
			result.push_back(pair<int, int>(i, 0));
		}
		if (lines[i].end_points.second == control_points[control_index])
		{
			result.push_back(pair<int, int>(i, 1));
		}
	}
	return result;
}

// ����һ�����Ƶ㣬���ظõ���һ���㼯�е�����
int Contour::getControlPointIndex(OPoint p)
{
	for (int i = 0; i < control_points.size(); ++i)
	{
		if (p == control_points[i])
		{
			return i;
		}
	}
	cout << "Can not find P" << endl;
	return -1;
}

// Mean Value Coordinate����
void Contour::MVCDeform()
{
	// ����֮�����е���������n*2
	MatrixXd new_vertices_matrix;

	// �����������
	MatrixXd static_vertices_matrix = MatrixXd::Zero(vertices.size(), 2);

	// ����vertices�е��ⲿ������
	for (int i = 0; i < out_all_map.size(); ++i)
	{
		vertices[out_all_map[i]].position = contour_points[i].position;
		static_vertices_matrix(out_all_map[i], 0) = contour_points[i].position.x;
		static_vertices_matrix(out_all_map[i], 1) = contour_points[i].position.y;
	}
	
	// �ⷽ�� ATAx = ATb 
	
	// ����ATA
	//static LLT<MatrixXd> llt = (coe_matrix.transpose()*coe_matrix).llt();
	//auto llt = (coe_matrix.transpose()*coe_matrix).llt();

	clock_t start = clock();
	// ���ֲ�ͬ�����Է������
	// LDLT�ֽ�
	//new_vertices_matrix = (coe_matrix.transpose()*coe_matrix).ldlt().solve(coe_matrix.transpose()*static_vertices_matrix);
	// LLT�ֽ�
	new_vertices_matrix = llt.solve(coe_matrix.transpose()*static_vertices_matrix);
	// QR�ֽ�
	//new_vertices_matrix = coe_matrix.colPivHouseholderQr().solve(static_vertices_matrix);
	clock_t end = clock();

	cout << "Calculate Time : " << (end - start) << " ms" << endl;


	//���������ڲ����ֵ
	for (int i = 0; i < vertices.size(); ++i)
	{
		vertices[i].position = glm::vec2(new_vertices_matrix(i, 0), new_vertices_matrix(i, 1));
	}
}

// �����߱���
void Contour::Deform(int handle_index, glm::vec2 new_pos)
{
	vector<pair<int, int>> control_lines = controlToLines(handle_index);

	for (int i = 0; i < control_lines.size(); ++i)
	{
		if (control_lines[i].second == -1)
		{
			lines[control_lines[i].first].end_points.first.position = new_pos;
		}
		else if (control_lines[i].second == 0)
		{
			lines[control_lines[i].first].curve_control_point.position = new_pos;
		}
		else if (control_lines[i].second == 1)
		{
			lines[control_lines[i].first].end_points.second.position = new_pos;
		}
	}
	getOPints();
	MVCDeform();
}

// ���������߽��е��������ǻ�
void Contour::delaunayTrianglation()
{

	// ����������
	Fade_2D dt;
	vector<Point2> vPoints;
	for (int i = 0;i < contour_points.size(); ++i)
	{
		vPoints.push_back(Point2(contour_points[i].position.x, contour_points[i].position.y,0.0f));
	}
	dt.insert(vPoints);
	
	// ������Χ����
	vector<Segment2> vSegments;
	for (int i = 0; i < vPoints.size(); ++i)
	{
		Point2& p0(vPoints[i]);
		Point2& p1(vPoints[(i + 1) % vPoints.size()]);
		vSegments.push_back(Segment2(p0, p1));
	}
	ConstraintGraph2* pCG = dt.createConstraint(vSegments, CIS_CONSTRAINED_DELAUNAY);

	// ���Լ��
	dt.applyConstraintsAndZones();
	dt.show("../Mesh/contour.ps");

	// ������������
	// �Զ��������ӵ�
	float x(0.0), y(0.0);
	for (int i = 0; i < vPoints.size(); ++i)
	{
		x += vPoints[i].x();
		y += vPoints[i].y();
	}
	x /= vPoints.size();
	y /= vPoints.size();

	Point2 seedPoint(x, y, 0.0f);
	vector<ConstraintGraph2*> vCG;
	vCG.push_back(pCG);
	Zone2* pGrowZone = dt.createZone(vCG, ZL_GROW, seedPoint);

	Zone2* pBoundedZone(pGrowZone->convertToBoundedZone());

	// ����Ӧ�ļ���ѡ�������̱߳�
	float mean_length = 0;
	for (int i = 0; i < vPoints.size(); ++i)
	{
		mean_length += sqrt(sqDistance2D(vPoints[i], vPoints[(i + 1)%vPoints.size()]));
	}
	mean_length /= vPoints.size();
	// ���һ��������ʾԼ�����Ƿ���Ա��ָͨ���趨Ϊtrue
	// �����ڱ�������Ӧ���趨Ϊfalse
	// ��ΪԼ���������ǻ�֮ǰ���ϸ���ģ����ɽ����޸ģ�
	// ��Ϊһ���޸���ᵼ�±߽��ϳ��ֶ���ĵ�����ɷ�����ִ���
	//dt.refine(pBoundedZone, 27, 100, 50, false);
	dt.refine(pBoundedZone, 27, mean_length*3.0f, mean_length*1.5f, false);
	dt.writeObj("../Mesh/mesh.obj", pBoundedZone); //���õ���������Ƭ����ΪOBJ

	readOBJ("../Mesh/mesh.obj");
}

void Contour::reMesh()
{
	Init();
}

//��OBJ�ļ��ж�ȡ���������
void Contour::readOBJ(string path)
{
	// �������������ɵ��еĹ�ϵ
	vector<int>map(contour_points.size());
	// ���еĶ�ȡ�ĵ�����ڲ�����ⲿ��
	vector<Point> tmp_vertices;
	// ���ж�ȡ��������
	vector<Triangle> tmp_triangles;

	ifstream infile(path);
	float x, y, z;
	char c;
	string line;
	while(!infile.eof())
	{
		getline(infile, line);
		if (line == "")continue;
		stringstream stringin(line);
		stringin >> c >> x >> y>> z;
		// �ж϶�����ǵ㻹��������
		if (c == 'v')
		{
			Point insert_point(x, y, false);

			// ���ڼ�С�����²�����ȣ�����ֻ��ͨ����������
			int index = -1;
			for (int i = 0; i < contour_points.size(); ++i)
			{
				if (glm::distance(insert_point.position, contour_points[i].position) < 0.1)
				{
					index = i;
				}
			}

			if (index != -1)
			{
				insert_point.out_or_in = true;
				map[index] = tmp_vertices.size();
				tmp_vertices.push_back(insert_point);
			}
			else
			{
				tmp_vertices.push_back(insert_point);
			}
		}
		if (c == 'f')
		{
			// ����face�е�������1��ʼ
			Triangle insert_triangle = {x-1,y-1,z-1};
			tmp_triangles.push_back(insert_triangle);
		}
	}
	//�洢�㣬�����κͶ�Ӧ��ϵ
	vertices = tmp_vertices;
	triangles = tmp_triangles;
	out_all_map = map;

	cout << "Read " << tmp_vertices.size() << " vertices and " << tmp_triangles.size() << " triangles from OBJ" << endl;
}

// ���contour
void Contour::clear()
{
	lines.clear();
	control_points.clear();
	contour_points.clear();
	vertices.clear();
	triangles.clear();
	out_all_map.clear();
}
