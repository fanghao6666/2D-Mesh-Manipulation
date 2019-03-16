
#include "Graphics.h"

using namespace std;
using namespace GEOM_FADE25D;


//*************************** Cout Test **********************************
/*
* Cout glm::vector in standard form: vec2,vec3,vec4
*/
void Out(glm::vec2 vec_2)
{
	cout << "Vec2 -- X : " << vec_2.x << "  Y : " << vec_2.y << endl;
}

void Out(glm::vec3 vec_3)
{
	cout << "Vec3 -- X : " << vec_3.x << "  Y : " << vec_3.y << "  Z : " << vec_3.z << endl;
}

void Out(glm::vec4 vec_4)
{
	cout << "Vec4 -- X : " << vec_4.x << "  Y : " << vec_4.y << "  Z : " << vec_4.z << " A : " << vec_4.a << endl;
}

//*************************** Delaunay ***********************************

/* \brief 2D Delaunay Triangulation
*
* Delaunay triangulation using exist contour points
*
* @param points are the exist contour points to delaunay
* @param segments are the constraint condition in delaunay
*  most of them are line segments
* @param is_loop specifies the line segments whether to loop
* @param max_length and min_length are the max and min limit
*  length of the segments that delaunay generated,default 0
* @param path is the path to save obj
* @param name is the name of obj
*
*	This method is to Delaunay Triangulation
*/

bool Delaunay(vector<glm::vec2> points, vector<glm::vec2> segments, bool is_loop,
	string name, double max_length = 0, double min_length = 0, string path = "defualt")
{
	// read the contour points
	Fade_2D dt;
	vector<Point2> vPoints;
	for (int i = 0; i < points.size(); ++i)
	{
		vPoints.push_back(Point2(points[i].x, points[i].y, 0.0f));
	}
	dt.insert(vPoints);
	// read the segment points
	vector<Point2> sPoints;
	for (int i = 0; i < segments.size(); ++i)
	{
		sPoints.push_back(Point2(segments[i].x, segments[i].y, 0.0f));
	}

	// read the constraint segments
	vector<Segment2> vSegments;
	if (is_loop == true)
	{
		for (int i = 0; i < sPoints.size(); ++i)
		{
			Point2& p0(sPoints[i]);
			Point2& p1(sPoints[(i + 1) % sPoints.size()]);
			vSegments.push_back(Segment2(p0, p1));
		}
	}
	else
	{
		for (int i = 1; i < sPoints.size(); ++i)
		{
			Point2& p0(sPoints[i - 1]);
			Point2& p1(sPoints[i]);
			vSegments.push_back(Segment2(p0, p1));
		}
	}
	ConstraintGraph2* pCG = dt.createConstraint(vSegments, CIS_CONSTRAINED_DELAUNAY);

	// apply the segment constraint
	dt.applyConstraintsAndZones();
	if (path == "default")
	{
		string _name = name + ".ps";
		dt.show(_name);
	}
	else
	{
		string _name = path + name + ".ps";
		dt.show(_name);
	}

	// generate the seed point
	// we use mean value point here
	double x(0.0), y(0.0);
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

	if (max_length == 0 && min_length == 0)
	{
		// Automatic calculation of the longest and shortest side lengths
		double mean_length = 0;
		for (int i = 0; i < vPoints.size(); ++i)
		{
			mean_length += sqrt(sqDistance2D(vPoints[i], vPoints[(i + 1) % vPoints.size()]));
		}
		mean_length /= vPoints.size();

		// 最后一个参数表示约束边是否可以被分割，通常设定为true
		// 但是在本程序中应该设定为false
		// 因为约束边在三角化之前是严格定义的，不可进行修改，
		// 因为一旦修改则会导致边界上出现多余的点则造成仿真出现错误
		dt.refine(pBoundedZone, 27, mean_length*3.0f, mean_length*1.5f, false);
	}
	else
	{
		dt.refine(pBoundedZone, 27, max_length, min_length, false);
	}

	// output mesh as .obj
	if (path == "default")
	{
		string _name = name + ".obj";
		dt.writeObj(_name, pBoundedZone);
	}
	else
	{
		string _name = path + name + ".obj";
		dt.writeObj(_name, pBoundedZone);
	}
	return true;
}

//*************************** Trigonometric ***********************************

/* \brief Calculate cosα
*
* Calculate cosine value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate cosine value of two vectors
*/
double Cos(glm::vec3 vec_1, glm::vec3 vec_2)
{
	return (glm::dot(glm::normalize(vec_1), glm::normalize(vec_2)));
}


/* \brief Calculate sinα
*
* Calculate sine value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate sine value of two vectors
*/
double Sin(glm::vec3 vec_1, glm::vec3 vec_2)
{
	return (sqrt(1.0f - pow(Cos(vec_1, vec_2), 2)));
}

/* \brief Calculate tanα
*
* Calculate tangent value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate tangent value of two vectors
*   Cosα could not be zero
*/
double Tan(glm::vec3 vec_1, glm::vec3 vec_2)
{
	double cos = Cos(vec_1, vec_2);
	double sin = Sin(vec_1, vec_2);

	assert(cos != 0);

	return (sin / cos);
}
//*************************** Vector ***********************************

/* \brief Vertical Unit Vector
*
* Get the Vertical Unit vector of a given vector
*
* @param vec is the given vector
* @param is_left mean whether the output vector lie on the left
*		of the given vector
*
*	This method is to Get the Vertical Unit vector of a given vector
*/
glm::vec2 getVerticalUnitVec(glm::vec2 vec, bool is_left)
{
	glm::vec2 unit_vec;

	unit_vec.x = vec.y / glm::length(vec);
	unit_vec.y = -1 * vec.x / glm::length(vec);

	if (is_left)
	{
		if (glm::cross(glm::vec3(vec, 0.0f), glm::vec3(unit_vec, 0.0f)).z > 0)
		{
			return unit_vec;
		}
		else
		{
			return (-1.0f * unit_vec);
		}
	}
	else
	{
		if (glm::cross(glm::vec3(vec, 0.0f), glm::vec3(unit_vec, 0.0f)).z > 0)
		{
			return (-1.0f * unit_vec);
		}
		else
		{
			return unit_vec;
		}
	}
}

//*************************** Geometry ***********************************

/* \brief Point to Line distance
*
* Calculate the distance of point and line in 3D space
*
* @param p is the point
* @param l1,l2 are points on the line
*
* This method is to Calculate the distance of point and line in 3D space
*
*/
double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	assert(p != l1 && p != l2);

	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;

	double d = glm::dot(line_vec, point_vec) / glm::length(line_vec);

	return (sqrt(pow(glm::length(point_vec), 2) - pow(d, 2)));
}

/* \brief Point to Plane distance
*
* Calculate the distance of point and plane in 3D space
*
* @param p is the point
* @param p1,p2,p3 are points on the plane
*
* This method is to Calculate the distance of point and plane in 3D space
*
*/
double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 project_point = pointToPlaneProjection(p, p1, p2, p3);

	return (glm::distance(p, project_point));
}

/* \brief Point to Plane Project
*
* Calculate the project point of point to plane in 3D space
*
* @param p is the point
* @param p1,p2,p3 are points on the plane
*
* This method is to Calculate the project point of point to plane in 3D space
*
*/
glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 plane_vec1 = p2 - p1;
	glm::vec3 plane_vec2 = p3 - p1;
	glm::vec3 point_vec = p - p1;
	glm::vec3 unit_normal_vec = glm::normalize(glm::cross(plane_vec1, plane_vec2));

	float dot = glm::dot(point_vec, unit_normal_vec);
	assert(dot != 0);

	return (p - dot * unit_normal_vec);
}

/* \brief Point to Line Project
*
* Calculate the project point of point to line in 3D space
*
* @param p is the point
* @param l1,l2 are points on the plane
*
* This method is to Calculate the project point of point to line in 3D space
*
*/
glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;
	glm::vec3 unit_line_vec = glm::normalize(line_vec);

	float dot = glm::dot(point_vec, unit_line_vec);

	return (l1 + dot * unit_line_vec);
}

/* \brief Line & Line interseciton point
*
* Calculate the intersection point of line and line in 3D space
*
* @param l1_1,l1_2 are points of line_1
* @param l2_1,l2_2 are points of line_2
*
* This method is to Calculate the intersection point of line and line in 3D space
*
*/
glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2)
{
	glm::vec3 line1 = l1_2 - l1_1;
	glm::vec3 line2 = l2_2 - l2_1;
	glm::vec3 norm1 = glm::normalize(line1);
	glm::vec3 norm2 = glm::normalize(line2);

	// whether two lines are collinear
	assert(glm::cross(line1.line2) != 0);

	// whether two lines are coplanar
	assert(glm::dot(glm::cross(l2_1 - l1_1, line1), line2) == 0);

	float t1 = glm::length(glm::cross(l1_1 - l2_1, norm2)) / glm::length(glm::cross(norm2, norm1));
	float t2 = glm::length(glm::cross(l2_1 - l1_1, norm1)) / glm::length(glm::cross(norm1, norm2));

	assert((l1_1 + t1 * norm1) == (l2_1 + t2 * norm2));

	return l1_1 + t1 * norm1;
}

/* \brief Segment & Segment interseciton point 2D
*
* Calculate the intersection point of Segment and Segment in 2D space
*
* @param s1_1,s1_2 are points of Segment1
* @param s2_1,s2_2 are points of Segment2
*
* This method is to Calculate the intersection point of line and line in 3D space
*
*/
glm::vec2 segToSegIntersection2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2)
{
	glm::vec3 l1_1 = glm::vec3(s1_1, 0.0f);
	glm::vec3 l1_2 = glm::vec3(s1_2, 0.0f);
	glm::vec3 l2_1 = glm::vec3(s2_1, 0.0f);
	glm::vec3 l2_2 = glm::vec3(s2_2, 0.0f);

	if (isSegmentIntersect2D(s1_1, s1_2, s2_1, s2_2))
	{
		glm::vec3 inter_point = lineToLineIntersection(l1_1, l1_2, l2_1, l2_2);
		return glm::vec2(inter_point);
	}
	else
	{
		cerr << "ERROR : two segments are not intersect!!!" << endl;
		return glm::vec2(0);
	}
}

/* \brief Plane & Line interseciton point
*
* Calculate the intersection point of line and plane in 3D space
*
* @param p1,p2,p3 are points in the plane
* @param l1,l2 are points on the line
*
* This method is to Calculate the intersection point of line and plane in 3D space
*
*/
glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 plane_norm_vec = glm::normalize(glm::cross(p2 - p1, p3 - p1));

	assert(glm::dot(line_vec, plane_norm_vec) != 0);

	glm::vec3 point_vec = p1 - l1;

	float t = (glm::dot(point_vec, plane_norm_vec)) / (glm::dot(plane_norm_vec, line_vec));

	return l1 + t * line_vec;
}

/* \brief Point on Line
*
* Judge whether point on line
*
* @param p is point
* @param l1,l2 are points on the line
*
* This method is to Judge whether point on line
*
*/
bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;

	return (glm::length(glm::cross(line_vec, point_vec)) == 0.0f);
}

/* \brief Point on Segment
*
* Judge whether point on segment
*
* @param p is point
* @param l1,l2 are points on the segment
*
* This method is to Judge whether point on segment
* Notice: end point is also belong to the segment
*
*/
bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2)
{
	if (p == s1 || p == s2) return true;

	glm::vec3 line_vec = s2 - s1;
	glm::vec3 point_vec1 = p - s1;
	glm::vec3 point_vec2 = p - s2;

	return (glm::length(glm::cross(line_vec, point_vec1)) == 0.0f && glm::dot(point_vec1, point_vec2) < 0.0f);
}

/* \brief Point in Plane
*
* Judge whether point in plane
*
* @param p is point
* @param p1,p2,p3 are points on the plane
*
* This method is to  Judge whether point in plane
*
*/
bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 p_vec1 = p2 - p1;
	glm::vec3 p_vec2 = p3 - p1;
	glm::vec3 point_vec = p - p1;

	return(glm::dot(glm::cross(p_vec1, p_vec2), point_vec) == 0.0f);
}

/* \brief Point in Triangle
*
* Judge whether point in Triangle
*
* @param p is point
* @param p1,p2,p3 are points on the Triangle
*
* This method is to  Judge whether point in Triangle
*
*/
bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 v0 = p3 - p1;
	glm::vec3 v1 = p2 - p1;
	glm::vec3 v2 = p - p1;

	float dot00 = glm::dot(v0, v0);
	float dot01 = glm::dot(v0, v1);
	float dot02 = glm::dot(v0, v2);
	float dot11 = glm::dot(v1, v1);
	float dot12 = glm::dot(v1, v2);

	float inverDeno = 1 / (dot00 * dot11 - dot01 * dot01);

	float u = (dot11 * dot02 - dot01 * dot12) * inverDeno;
	float v = (dot00 * dot12 - dot01 * dot02) * inverDeno;

	if (u < 0 || u > 1) return false;
	if (v < 0 || v > 1) return false;
	if (u + v < 1)
		return true;
	else
		return false;

}

/* \brief Point in Polygon
*
* Judge whether point in or out polygon
*
* @param p is point
* @param polygon are point vector of polygon contour
*
* This method is to  Judge whether point in or out polygon
* Draw a ray from the target point and see the number of
* intersections between the ray and all sides of the polygon.
* If there are odd intersections, it means inside, and if
* there are even intersections, it means outside.
*/
bool isPointInPolygon2D(glm::vec2 p, vector<glm::vec2> polygon)
{
	float min_x = 1000;
	for (int i = 0; i < polygon.size(); ++i)
	{
		if (polygon[i].x < min_x)
		{
			min_x = polygon[i].x;
		}
	}
	glm::vec2 out_p = glm::vec2(min_x - 10.0f, p.y);

	int cnt = 0;
	for (int i = 0; i < polygon.size(); ++i)
	{
		glm::vec2 first = polygon[i];
		glm::vec2 second = polygon[(i + 1) % polygon.size()];
		if (isSegmentIntersect2D(p, out_p, first, second))
		{
			cnt++;
		}
	}
	return (cnt % 2 == 1);
}

/* \brief Segment in Polygon
*
* Judge whether segment in or out polygon
*
* @param s1,s2 are endpoints of segment
* @param polygon are point vector of polygon contour
*
* This method is to  Judge whether segment in or out polygon
* First, determine whether the endpoint is inside the polygon
* Then determine whether the line segment intersects with any side.
*/
bool isSegmentInPolygon2D(glm::vec2 s1, glm::vec2 s2, vector<glm::vec2> polygon)
{
	bool is_s1_in = isPointInPolygon2D(s1, polygon);
	bool is_s2_in = isPointInPolygon2D(s2, polygon);

	if (is_s1_in && is_s2_in)
	{
		for (int i = 0; i < polygon.size(); ++i)
		{
			glm::vec2 first = polygon[i];
			glm::vec2 second = polygon[(i + 1) % polygon.size()];
			if (isSegmentIntersect2D(s1, s2, first, second))
			{
				return false;
			}
		}
		return true;
	}
	else
	{
		return false;
	}
}

/* \brief Two segments intersect
*
* Judge whether two segments intersect
*
* @param s1_1,s1_2 are endpoints of segment1
* @param s2_1,s2_2 are endpoints of segment2
*
* This method is to  Judge whether two segments intersect
*/
bool isSegmentIntersect2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2)
{
	glm::vec3 a = glm::vec3(s1_1, 0);
	glm::vec3 b = glm::vec3(s1_2, 0);
	glm::vec3 c = glm::vec3(s2_1, 0);
	glm::vec3 d = glm::vec3(s2_2, 0);
	glm::vec3 seg1 = b - a;
	glm::vec3 seg2 = d - c;

	float u = glm::dot(glm::cross(c - a, seg1), glm::cross(d - a, seg1));
	float v = glm::dot(glm::cross(a - c, seg2), glm::cross(b - c, seg2));

	return(u <= 0.00000001 && v <= 0.00000001);
}