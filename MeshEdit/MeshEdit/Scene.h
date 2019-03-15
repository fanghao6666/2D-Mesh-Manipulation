#pragma once
#include<iostream>
#include "Model.h"

using namespace std;

class Scene
{
public:
	Scene() {};
	~Scene() {};

public:
	// �򳡾������Contour
	static void addContour(Contour contour);

	// ������ʼ��
	static void Init(int argc, char* *argv);

	// ��Ⱦ����
	static void OnRender();

	// ���ڴ�С����
	static void OnReshape(int width, int height);

	// ���������
	static void OnMouseClick(int button,int state, int x, int y);

	// ����ƶ�����
	static void OnMouseMove(int x, int y);

	// ���̲���
	static void OnKey(unsigned char key, int, int);

	// ѭ����Ⱦ
	void Render();
public:
	static Contour contour;

	// ѡ��ĵ��������
	static int select_point_index;
	static int select_line_index;

	// ѡ��ĵ���ߵ�����
	static glm::vec2 click_point;
	static glm::vec2 click_line_s;
	static glm::vec2 click_line_m;
	static glm::vec2 click_line_e;

	// ������껭���������Ƶ�
	static vector<pair<glm::vec2,bool>> draw_contour_points;

	// ��ǰ�Ƿ������߿��Ƶ�
	static bool is_curve_control;

	// ��ǰ�Ƿ�dart
	static bool is_draw_dart;
	// dart����ʼ��
	static glm::vec2 dart_begin;
	static glm::vec2 dart_end;
	

};