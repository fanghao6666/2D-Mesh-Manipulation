#include <iostream>

#include "Model.h"
#include "Scene.h"

using namespace std;



int main(int argc, char* *argv)
{
	//����һ����Ⱦ����
	Scene scene;
	scene.Init(argc, argv);
	scene.Render();

	return 0;
}