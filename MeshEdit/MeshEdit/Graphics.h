// Copyright (C) Common Graphics Algorithms e.U, Fang Hao
//
// This file is implementation of Common Graphics Algorithms.
//
// Please contact the author if any conditions of this file are
// not clear to you.
//
// Author: Fang Hao .Nanjing University ,VISG

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <glm.hpp>
#include <Fade_2D.h>

using namespace std;
using namespace GEOM_FADE25D;

/********************** Declare All Functions Here ***********************/

double Cos(glm::vec3 vec_1, glm::vec3 vec_2);
double Sin(glm::vec3 vec_1, glm::vec3 vec_2);
double Tan(glm::vec3 vec_1, glm::vec3 vec_2);

glm::vec2 getVerticalUnitVec(glm::vec2 vec, bool is_left);

double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);

glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2);
glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
glm::vec2 segToSegIntersection2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2);

bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2);
bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
bool isPointInPolygon2D(glm::vec2 p, vector<glm::vec2> polygon);
bool isSegmentInPolygon2D(glm::vec2 s1, glm::vec2 s2, vector<glm::vec2> polygon);
bool isSegmentIntersect2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2);

/********************** Declare All Functions Here ***********************/

