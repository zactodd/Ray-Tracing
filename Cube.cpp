/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The Cube class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cube.h"
#include <math.h>
#include <algorithm> 

/**
* Cone's intersection method.  The input is a ray (pos, dir). 
*/
float Cube::intersect(Vector pos, Vector dir)
{


	float halfH = height / 2;
	float a = halfH * (dir.x * dir.x + dir.z * dir.z)
			+ height * dir.y * dir.y;
			
	float b = halfH * 2 * (pos.x * dir.x + pos.x * dir.x)
			+ height * 2 * pos.y * dir.y;
			
	float c = halfH * (pos.x * pos.x + pos.z * pos.z)
			+ height * pos.y * pos.y;
			
			
	float discriminant = sqrt(b*b - 4.0 * a * c) ;
	
	if(fabs(discriminant) < 0.001) return -1.0; 
    
    if(discriminant < 0.0) return -1.0;
    
	float t1 = (-b + discriminant) / (2 * a);
	float t2 = (-b - discriminant) / (2 * a);
	
	if(fabs(t1) < 0.001 )
    {
        if (t2 > 0) return t2;
        else t1 = -1.0;
    }
    if(fabs(t2) < 0.001 ) t2 = -1.0;
    float closet = std::min(t1, t2);
    float furest = std::max(t1, t2);
    float y1 = (pos.y + closet * dir.y) - center.y;
    float y2 = (pos.y + furest * dir.y) - center.y;
    if(! (y1 < 0 || y1 > height) && closet != -1.0) return closet;
    if(! (y2 < 0 || y2 > height) && furest != -1.0) return furest;			
	return -1.0;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the Cone.
*/
Vector Cube::normal(Vector p)
{

	Vector n = Vector(p.x, 0, p.z);
    return n;
}
