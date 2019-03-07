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
	// 向场景中添加Contour
	static void addContour(Contour contour);

	// 场景初始化
	static void Init(int argc, char* *argv);

	// 渲染场景
	static void OnRender();

	// 窗口大小调整
	static void OnReshape(int width, int height);

	// 鼠标点击操作
	static void OnMouseClick(int button,int state, int x, int y);

	// 鼠标移动操作
	static void OnMouseMove(int x, int y);

	// 键盘操作
	static void OnKey(unsigned char key, int, int);

	// 循环渲染
	void Render();
public:
	static Contour contour;

	// 选择的点和线索引
	static int select_point_index;
	static int select_line_index;

	// 选择的点和线的坐标
	static glm::vec2 click_point;
	static glm::vec2 click_line_s;
	static glm::vec2 click_line_m;
	static glm::vec2 click_line_e;

	// 保存鼠标画的轮廓控制点
	static vector<pair<glm::vec2,bool>> draw_contour_points;

	// 当前是否是曲线控制点
	static bool is_curve_control;
	

};