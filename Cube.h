/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cone class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CUBE
#define H_CUBE

#include "Object.h"

/**
 * Defines a simple Cube located at 'center' 
 * with the specified radius and hight
 */
class Cube : public Object
{

private:
    Vector center;
    float height;

public:	
	Cube()
		: center(Vector()), height(1) //Default constructor creates a unit cube
	{
		color = Color::WHITE;
	};
	
    Cube(Vector c, float h, Color col)
		: center(c), height(h)
	{
		color = col;
	};


	float intersect(Vector pos, Vector dir);

	Vector normal(Vector p);

};

#endif //!H_CUBE
