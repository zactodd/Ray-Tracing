/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The Cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include <math.h>
#include <algorithm> 

/**
* Sphere's intersection method.  The input is a ray (pos, dir). 
*/
float Cylinder::intersect(Vector pos, Vector dir)
{

	float a = pow(dir.x, 2) + pow(dir.z, 2);
	float b = 2 * (dir.x * (pos.x - center.x) + dir.z * (pos.z - center.z));
	float c = pow((pos.x - center.x), 2) + pow((pos.z - center.z), 2) - pow(radius, 2);
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
    
    float closest = std::min(t1, t2);
    float furthest = std::max(t1, t2);
    float y1 = (pos.y + closest * dir.y) - center.y;
    float y2 = (pos.y + furthest * dir.y) - center.y;
    if(! (y1 < 0 || y1 > height) && closest != -1.0) return closest;
    if(! (y2 < 0 || y2 > height) && furthest != -1.0) return furthest;			
	return -1.0;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the sphere.
*/
Vector Cylinder::normal(Vector p)
{
	Vector n = Vector((p.x - center.x), 0.0, (p.z - center.z));
    n.normalise();
    return n;
}
