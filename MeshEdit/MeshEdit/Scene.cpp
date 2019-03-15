#include<iostream>
#include <vector>
#include <math.h>
#include <assert.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
#include <glm.hpp>
#include <windows.h>

#include "Scene.h"

using namespace std;

Contour Scene::contour;
int Scene::select_point_index = -1;
int Scene::select_line_index = -1;
glm::vec2 Scene::click_point;
glm::vec2 Scene::click_line_s;
glm::vec2 Scene::click_line_m;
glm::vec2 Scene::click_line_e;
vector<pair<glm::vec2, bool>> Scene::draw_contour_points;
bool Scene::is_curve_control = false;
bool Scene::is_draw_dart = false;
glm::vec2 Scene::dart_begin = glm::vec2(0.0f, 0.0f);
glm::vec2 Scene::dart_end = glm::vec2(0.0f, 0.0f);


// �߶���������֮��ļ��
//float point_distance = 15.0f;


// ѡ����ʾ������
void selectFont(int size, int charset)
{
	HFONT hFont = CreateFont(size, 0, 0, 450, FW_MEDIUM, 0, 0, 0,
		charset, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
		DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS,NULL);
	HFONT hOldFont = (HFONT)SelectObject(wglGetCurrentDC(), hFont);
	DeleteObject(hOldFont);
}

// ��OPENGL����ʾ�ַ�
void drawString(const char*str)
{
	static int isFirstCall = 1;
	static GLuint lists;

	if (isFirstCall) { // ����ǵ�һ�ε��ã�ִ�г�ʼ��
					   // Ϊÿһ��ASCII�ַ�����һ����ʾ�б�
		isFirstCall = 0;

		// ����MAX_CHAR����������ʾ�б���
		lists = glGenLists(128);

		// ��ÿ���ַ��Ļ������װ����Ӧ����ʾ�б���
		wglUseFontBitmaps(wglGetCurrentDC(), 0, 128, lists);
	}
	// ����ÿ���ַ���Ӧ����ʾ�б�����ÿ���ַ�
	for (; *str != '\0'; ++str)
		glCallList(lists + *str);
}



// ������ʼ��
void Scene::Init(int argc,char* *argv)
{
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(1280, 960);
	glutInitWindowPosition(250, 30);
	glutCreateWindow("Mesh Edit");

	glClearColor(1.0, 1.0, 1.0, 0.0);
	glColor3f(0.5f, 0.4f, 0.9f);
	glPointSize(5.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluOrtho2D(0.0, 1280, 960, 0.0);
}

// �򳡾������contour
void Scene::addContour(Contour _contour)
{
	contour = _contour;
}

// ��ʾ
void Scene::OnRender()
{
	glClear(GL_COLOR_BUFFER_BIT);

	// ����ʵʱ����������
	glBegin(GL_POINTS);
	glColor3f(1.0f, 0.0f, 0.0f);
	for (int i = 0; i < draw_contour_points.size(); ++i)
	{
		glVertex2f(draw_contour_points[i].first.x, draw_contour_points[i].first.y);
	}
	glEnd();

	// �������Ƶ�
	glBegin(GL_POINTS);
	glColor3f(1.0f, 0.0f, 0.0f);
	for (int i = 0; i < contour.control_points.size(); ++i)
	{
		if (i != select_point_index)
		{
			glVertex2f(contour.control_points[i].position.x, contour.control_points[i].position.y);
		}
	}
	glEnd();

	// ������ѡ�еĿ��Ƶ�
	glBegin(GL_POINTS);
	glColor3f(0.0f, 0.0f, 1.0f);
	if (select_point_index != -1)
	{
		glVertex2f(contour.control_points[select_point_index].position.x, contour.control_points[select_point_index].position.y);
	}
	glEnd();

	// ���������Ƶ���ı߽��
	glBegin(GL_POINTS);
	glColor3f(0.0f, 1.0f, 0.0f);
	for (int i = 0; i < contour.lines.size(); ++i)
	{
		for (int j = 0; j < contour.lines[i].point_list.size(); ++j)
		{
			glVertex2f(contour.lines[i].point_list[j].position.x, contour.lines[i].point_list[j].position.y);
		}
	}
	glEnd();

	// �����ڲ���
	glBegin(GL_POINTS);
	glColor3f(0.0f, 0.0f, 1.0f);
	for (int i = 0; i < contour.vertices.size(); ++i)
	{
		if (contour.vertices[i].out_or_in == false)
		{
			glVertex2f(contour.vertices[i].position.x, contour.vertices[i].position.y);
		}
	}
	glEnd();

	// �������е��߶�,�����е��߶λ�������
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < contour.triangles.size(); ++i)
	{
		glVertex2f(contour.vertices[contour.triangles[i].v0].position.x, contour.vertices[contour.triangles[i].v0].position.y);
		glVertex2f(contour.vertices[contour.triangles[i].v1].position.x, contour.vertices[contour.triangles[i].v1].position.y);

		glVertex2f(contour.vertices[contour.triangles[i].v0].position.x, contour.vertices[contour.triangles[i].v0].position.y);
		glVertex2f(contour.vertices[contour.triangles[i].v2].position.x, contour.vertices[contour.triangles[i].v2].position.y);

		glVertex2f(contour.vertices[contour.triangles[i].v1].position.x, contour.vertices[contour.triangles[i].v1].position.y);
		glVertex2f(contour.vertices[contour.triangles[i].v2].position.x, contour.vertices[contour.triangles[i].v2].position.y);
	}
	glEnd();

	// ����dart��
	if (is_draw_dart)
	{
		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		if (dart_begin != glm::vec2(0.0f, 0.0f) && dart_end != glm::vec2(0.0f, 0.0f))
		{
			glVertex2f(dart_begin.x, dart_begin.y);
			glVertex2f(dart_end.x, dart_end.y);
		}
		glEnd();
	}

	glFlush();
}

// ���������
void Scene::OnMouseClick(int button,int state, int x, int y)
{
	// ����������
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		// ��ȡ�����
		 click_point = glm::vec2(x, y);
		 if (is_draw_dart == true)
		 {
			 dart_end = glm::vec2(0.0f, 0.0f);
			 dart_begin = click_point;
		 }
		 else
		 {
			 // �ж��Ƿ�ѡ�����Ƶ�
			 for (int i = 0; i < contour.control_points.size(); ++i)
			 {
				 if (glm::distance(click_point, contour.control_points[i].position) < 5)
				 {
					 select_point_index = i;
					 click_point = contour.control_points[i].position;
				 }
			 }
			 //�ж��Ƿ�ѡ��������
			 for (int i = 0; i < contour.lines.size(); ++i)
			 {
				 if (contour.lines[i].isCloseToLine(click_point))
				 {
					 select_line_index = i;
					 click_line_s = contour.lines[i].end_points.first.position;
					 click_line_m = contour.lines[i].curve_control_point.position;
					 click_line_e = contour.lines[i].end_points.second.position;
				 }
			 }
		 }
	 }
	// �������ɿ�
	else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
	{
		if (is_draw_dart)
		{
			contour.createDart(dart_begin, dart_end);
		}
		select_point_index = -1;
		select_line_index = -1;
	}

	// ����м�����
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
	{
		assert(draw_contour_points.size() > 2);

		Line new_line;
		if (draw_contour_points[draw_contour_points.size() - 1].second == false)
		{
			new_line.is_curve = false;
			new_line.end_points.first = OPoint(draw_contour_points[draw_contour_points.size() - 1].first.x, draw_contour_points[draw_contour_points.size() - 1].first.y, true, 1);
			new_line.end_points.second = OPoint(draw_contour_points[0].first.x, draw_contour_points[0].first.y, true, 0);
		}
		else
		{
			new_line.is_curve = true;
			new_line.end_points.first = OPoint(draw_contour_points[draw_contour_points.size() - 2].first.x, draw_contour_points[draw_contour_points.size() - 2].first.y, true, 1);
			new_line.end_points.second = OPoint(draw_contour_points[0].first.x, draw_contour_points[0].first.y, true, 0);
			new_line.curve_control_point = OPoint(draw_contour_points[draw_contour_points.size() - 1].first.x, draw_contour_points[draw_contour_points.size() - 1].first.y, true, -1);
		}
		new_line.point_size = floor(glm::distance(new_line.end_points.first.position, new_line.end_points.second.position) / POINT_DISTANCE);
		contour.lines.push_back(new_line);

		// ��contour���г�ʼ��
		contour.Init();
	}
	// ����м��ɿ�
	else if (button == GLUT_MIDDLE_BUTTON && state == GLUT_UP)
	{
		draw_contour_points.clear();
	}

	// ����Ҽ�����
	if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		draw_contour_points.push_back(pair<glm::vec2, bool>(glm::vec2(x,y), is_curve_control));

		if (is_curve_control == true)
		{
			is_curve_control = false;
		}
	}
	// ����Ҽ��ɿ�
	else if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
	{
		if (draw_contour_points.size() > 1 && draw_contour_points[draw_contour_points.size() - 1].second != true)
		{
			Line new_line;

			if (draw_contour_points[draw_contour_points.size() - 2].second == false)
			{
				new_line.is_curve = false;
				new_line.end_points.first = OPoint(draw_contour_points[draw_contour_points.size() - 2].first.x, draw_contour_points[draw_contour_points.size() - 2].first.y, true, 1);
				new_line.end_points.second = OPoint(draw_contour_points[draw_contour_points.size() - 1].first.x, draw_contour_points[draw_contour_points.size() - 1].first.y, true, 0);
			}
			else
			{
				new_line.is_curve = true;
				new_line.end_points.first = OPoint(draw_contour_points[draw_contour_points.size() - 3].first.x, draw_contour_points[draw_contour_points.size() - 3].first.y, true, 1);
				new_line.end_points.second = OPoint(draw_contour_points[draw_contour_points.size() - 1].first.x, draw_contour_points[draw_contour_points.size() - 1].first.y, true, 0);
				new_line.curve_control_point = OPoint(draw_contour_points[draw_contour_points.size() - 2].first.x, draw_contour_points[draw_contour_points.size() - 2].first.y, true, -1);
			}
			new_line.point_size = floor(glm::distance(new_line.end_points.first.position, new_line.end_points.second.position) / POINT_DISTANCE);
			contour.lines.push_back(new_line);
		}
		contour.getOPints();
	}

	OnRender();
}

// ����ƶ�����
void Scene::OnMouseMove(int x, int y)
{
	glm::vec2 new_position = glm::vec2(x, y);


	if (is_draw_dart == true)
	{
		dart_end = new_position;
	}
	else
	{
		// ѡ�е���е����
		if (select_point_index != -1)
		{
			contour.Deform(select_point_index, new_position);
		}
		// ѡ�б߽��б���
		if (select_line_index != -1)
		{
			OPoint first_point = contour.lines[select_line_index].end_points.first;
			OPoint middle_point = contour.lines[select_line_index].curve_control_point;
			OPoint second_point = contour.lines[select_line_index].end_points.second;
			// ��һ�������
			contour.Deform(contour.getControlPointIndex(first_point), click_line_s + (new_position - click_point));
			// ���ߵ��м���Ƶ����
			if (contour.lines[select_line_index].is_curve == true)
			{
				contour.Deform(contour.getControlPointIndex(middle_point), click_line_m + (new_position - click_point));
			}
			// �ڶ��������
			contour.Deform(contour.getControlPointIndex(second_point), click_line_e + (new_position - click_point));
		}
	}


	OnRender();
}

// ���̲���
void Scene::OnKey(unsigned char key, int, int)
{
	switch (key)
	{
	case 'q':
	case 'Q':
		is_curve_control = !(is_curve_control);
		break;
	case 'c':
	case 'C':
	{
		Contour *new_contour = new Contour();
		contour = *new_contour;
		draw_contour_points.clear();
	}
		break;
	case 'd':
	case 'D':
		is_draw_dart = !is_draw_dart;
		if (is_draw_dart)
		{
			cout << "Start draw dart" << endl;
		}
		else
		{
			cout << "Stop draw dart" << endl;
		}
		break;
	case 'r':
	case 'R':
		contour.reMesh();
		break;
	default:
		break;
	}
	OnRender();
}

// ��Ⱦѭ��
void Scene::Render()
{
	glutDisplayFunc(OnRender);
	glutMouseFunc(OnMouseClick);
	glutMotionFunc(OnMouseMove);
	glutKeyboardFunc(OnKey);

	glutMainLoop();
}