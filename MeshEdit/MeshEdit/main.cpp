#include <iostream>

#include "Model.h"
#include "Scene.h"


using namespace std;

int main(int argc, char* *argv)
{
	// 生成一个渲染场景
	Scene scene;
	scene.Init(argc, argv);
	scene.Render();

	return 0;
}