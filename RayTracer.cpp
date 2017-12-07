// ========================================================================
// COSC 363  Computer Graphics  Lab07
// A simple ray tracer
// ========================================================================

#include <iostream>
#include <cmath>
#include <vector>
#include "Vector.h"
#include "Sphere.h"
#include "Cylinder.h"
#include "Plane.h"
#include "Cone.h"
#include "Cube.h"
#include "Color.h"
#include "Object.h"
#include <GL/glut.h>


using namespace std;

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int PPU = 60;     //Total 600x600 pixels
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;
const float GLOBAL_MEDIUM = 1.0;
vector<Object*> sceneObjects;
float zone = 0;

Vector light;
Color backgroundCol;

//A useful struct
struct PointBundle   
{
	Vector point;
	int index;
	float dist;
};

/*
* This function compares the given ray with all objects in the scene
* and computes the closest point  of intersection.
*/
PointBundle closestPt(Vector pos, Vector dir)
{
    Vector  point(0, 0, 0);
	float min = 10000.0;

	PointBundle out = {point, -1, 0.0};

    for(int i = 0;  i < sceneObjects.size();  i++)
	{
        float t = sceneObjects[i]->intersect(pos, dir);
		if(t > 0)        //Intersects the object
		{
			point = pos + dir*t;
			if(t < min)
			{
				out.point = point;
				out.index = i;
				out.dist = t;
				min = t;
			}
		}
	}

	return out;
}

void cubiod()
{
	Plane *front = new Plane(Vector(-30, 7, -67), Vector(-12, 7, -67), Vector(-12., 23, -67), Vector(-30, 23, -67), Color::YELLOW);
    Plane *down = new Plane(Vector(-30, 7, -67), Vector(-12, 7, -67), Vector(-12., 7, -87), Vector(-30, 7, -87), Color::YELLOW);
    Plane *right = new Plane(Vector(-12, 7, -67), Vector(-12., 7, -87), Vector(-12., 23, -87), Vector(-12., 23, -67), Color::YELLOW);
    Plane *up = new Plane(Vector(-12., 23, -67), Vector(-30, 23, -67), Vector(-30., 23, -87), Vector(-12, 23, -87), Color::YELLOW);
    Plane *left = new Plane(Vector(-30, 7, -67), Vector(-30, 7, -87), Vector(-30, 23, -87), Vector(-30, 23, -67), Color::YELLOW);
    Plane *back = new Plane(Vector(-30, 7, -87), Vector(-12., 7, -87), Vector(-12., 23, -87), Vector(-30, 23, -87), Color::YELLOW);
	sceneObjects.push_back(front);
	sceneObjects.push_back(down);
	sceneObjects.push_back(right);
	sceneObjects.push_back(up);
	sceneObjects.push_back(left);
	sceneObjects.push_back(back);
}


/*
* Computes the colour value obtained by tracing a ray.
* If reflections and refractions are to be included, then secondary rays will 
* have to be traced from the point, by converting this method to a recursive
* procedure.
*/

Color trace(Vector pos, Vector dir, int step)
{
	Color colorSum;
	zone++;
    PointBundle q = closestPt(pos, dir);
    if(q.index == -1) return backgroundCol;        //no intersection
    
	Vector n = sceneObjects[q.index]->normal(q.point);  //normal vector
	
	Color col = sceneObjects[q.index]->getColor(); //Object's colour
	if(q.index == 3){
		if(int(zone) % 2 && int(q.point.x * 2.5 - 100) % 3){
			col = col.phongLight(Color::YELLOW, 0.0, 0.0);
		}else{
			col = col.phongLight(Color::CYAN, 0.0, 0.0);
		}
	}
	
	if(q.index == 4){
		if(int(q.point.x - 1000) % 2 == int(q.point.z - 5040) % 2){
			col = col.phongLight(Color::BLACK, 0.0, 0.0);
		}else{
			col = col.phongLight(Color::WHITE, 0.0, 0.0);
		}
	}
	
	Vector l = light - q.point; //The light source vector
	
	float lightDist = l.length(); //Distance to light
	
	Vector v(-dir.x, -dir.y, -dir.z); //View vector;
	
	l.normalise(); //Normalise this vector, and compute the dot product l•n
	
	float lDotn = l.dot(n); 
	
	PointBundle s = closestPt(q.point, l);
	if(s.index >-1 && s.dist < lightDist){
		 colorSum = col.phongLight(backgroundCol, 0.0, 0.0);
	}
	else if(lDotn <= 0){
		colorSum = col.phongLight(backgroundCol, 0.0, 0.0);
	}
	else{
		// r = 2(L.n)n – L. ‘l’ = el
		Vector r = ((n * 2) * lDotn)- l; 
		r.normalise();
		float rDotv = r.dot(v);
		float spec;
		if(rDotv < 0) spec = 0.0;
		else spec = pow(rDotv, 10); //Phong exponent = 10
		colorSum = col.phongLight(backgroundCol, lDotn, spec); 
	}
	

	if(q.index == 2 && step < MAX_STEPS)
	{
		float n1 = 1.0;
		float n2 = 1.5;
		float nDotDir = dir.dot(n);
		float cosTheta = sqrt(1 - (pow(n1 / n2, 2.0)) * (1 - (pow(nDotDir, 2.0))));
		Vector gVector = (dir * (n1 / n2)) - n * (((n1 / n2) * nDotDir) + cosTheta);	
		Vector inverseN = sceneObjects[closestPt(q.point, gVector).index] -> normal(closestPt(q.point, gVector).point) * - 1;
		cosTheta = sqrt(1 - (pow(n2 / n1, 2.0)) * (1 - (pow(gVector.dot(inverseN), 2.0))));
		Vector refractionVector = (gVector * (n2 / n1)) - inverseN * (((n2 / n1) * gVector.dot(inverseN)) + cosTheta);
		refractionVector.normalise();
		Color refractionCol = trace(closestPt(q.point, gVector).point, refractionVector, step+1);
		colorSum.combineColor(refractionCol, 0.8);
	}

	//Transmisson
	if(q.index == 0 && step < MAX_STEPS){
		Color transmissionCol = trace(q.point, dir, step + 1);
		colorSum.combineColor(transmissionCol, 0.8);
	}

	
	//Generate reflection ray
	if((q.index == 1) && step < MAX_STEPS){
		Vector reflectionVector = (n * 2) * (n.dot(v)) - v;
		reflectionVector.normalise();
		
		Color reflectionCol = trace(q.point, reflectionVector, step + 1);
		colorSum.combineColor(reflectionCol, 1);
	}
	
	return colorSum;
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each pixel as quads.
//---------------------------------------------------------------------------------------
//~ void display()
//~ {
	//~ int widthInPixels = (int)(WIDTH * PPU);
	//~ int heightInPixels = (int)(HEIGHT * PPU);
	//~ float pixelSize = 1.0/PPU;
	//~ float halfPixelSize = pixelSize/2.0;
	//~ float quterPixelSize = pixelSize/4.0;
	//~ float x1, y1, xc, yc;
	//~ Vector eye(0., 0., 0.);
//~ 
	//~ glClear(GL_COLOR_BUFFER_BIT);
//~ 
	//~ glBegin(GL_QUADS);  //Each pixel is a quad.
//~ 
	//~ for(int i = 0; i < widthInPixels; i++)	//Scan every "pixel"
	//~ {
		//~ x1 = XMIN + i*pixelSize;
		//~ xc = x1 + halfPixelSize;
		//~ for(int j = 0; j < heightInPixels; j++)
		//~ {
			//~ y1 = YMIN + j*pixelSize;
			//~ yc = y1 + halfPixelSize;
//~ 
		    //~ Vector dir1(xc - quterPixelSize, yc - quterPixelSize, -EDIST);	//direction of the primary ray
		    //~ Vector dir2(xc + quterPixelSize, yc - quterPixelSize, -EDIST);	//direction of the primary ray
		    //~ Vector dir3(xc - quterPixelSize, yc + quterPixelSize, -EDIST);	//direction of the primary ray
		    //~ Vector dir4(xc + quterPixelSize, yc + quterPixelSize, -EDIST);	//direction of the primary ray
//~ 
		    //~ dir1.normalise();			//Normalise this direction
		    //~ dir2.normalise();			//Normalise this direction
		    //~ dir3.normalise();			//Normalise this direction
		    //~ dir4.normalise();			//Normalise this direction
//~ 
		    //~ Color col1 = trace (eye, dir1, 1); //Trace the primary ray and get the colour value
		    //~ Color col2 = trace (eye, dir2, 1); //Trace the primary ray and get the colour value
		    //~ Color col3 = trace (eye, dir3, 1); //Trace the primary ray and get the colour value
		    //~ Color col4 = trace (eye, dir4, 1); //Trace the primary ray and get the colour value
			//~ 
			//~ col1.scaleColor(0.25);
			//~ col1.combineColor(col2, 0.25);
			//~ col1.combineColor(col3, 0.25);
			//~ col1.combineColor(col4, 0.25);
			//~ 
//~ 
			//~ glColor3f(col1.r, col1.g, col1.b);
			//~ glVertex2f(x1, y1);				//Draw each pixel with its color value
			//~ glVertex2f(x1 + pixelSize, y1);
			//~ glVertex2f(x1 + pixelSize, y1 + pixelSize);
			//~ glVertex2f(x1, y1 + pixelSize);
        //~ }
    //~ }
//~ 
    //~ glEnd();
    //~ glFlush();
//~ }

//Without anti-alasing
void display(){
	int widthInPixels = (int)(WIDTH * PPU);
	int heightInPixels = (int)(HEIGHT * PPU);
	float pixelSize = 1.0/PPU;
	float halfPixelSize = pixelSize/2.0;
	float x1, y1, xc, yc;
	Vector eye(0., 0., 0.);

	glClear(GL_COLOR_BUFFER_BIT);

	glBegin(GL_QUADS);  //Each pixel is a quad.

	for(int i = 0; i < widthInPixels; i++)	//Scan every "pixel"
	{
		x1 = XMIN + i*pixelSize;
		xc = x1 + halfPixelSize;
		for(int j = 0; j < heightInPixels; j++)
		{
			y1 = YMIN + j*pixelSize;
			yc = y1 + halfPixelSize;

		    Vector dir(xc, yc, -EDIST);	//direction of the primary ray

		    dir.normalise();			//Normalise this direction

		    Color col = trace (eye, dir, 1); //Trace the primary ray and get the colour value
			glColor3f(col.r, col.g, col.b);
			glVertex2f(x1, y1);				//Draw each pixel with its color value
			glVertex2f(x1 + pixelSize, y1);
			glVertex2f(x1 + pixelSize, y1 + pixelSize);
			glVertex2f(x1, y1 + pixelSize);
        }
    }

    glEnd();
    glFlush();
}

void initialize()
{
	//Iniitialize background colour and light's position
	backgroundCol = Color::GRAY;
	light = Vector(10.0, 40.0, -5.0);
	//Add spheres to the list of scene objects here.	

	Sphere *sphere1 = new Sphere(Vector(6, 1, -60), 5.0, Color::CYAN);
	sceneObjects.push_back(sphere1);
	
	Sphere *sphere2 = new Sphere(Vector(-2, 12, -100), 16.0, Color::BLACK);
	sceneObjects.push_back(sphere2);
	
	Sphere *sphere3 = new Sphere(Vector(12, 12, -70), 5.0, Color::YELLOW);	
	sceneObjects.push_back(sphere3);
	
	Sphere *sphere4 = new Sphere(Vector(-8, -2, -50), 4.5, Color::WHITE);
	sceneObjects.push_back(sphere4);
	
	
	Plane *plane = new Plane(Vector(-1000, -10, 5000), Vector(1000, -10, 5000),
	Vector(1000., -10, -5000), Vector(-1000., -10, -5000), Color::WHITE);
	sceneObjects.push_back(plane);
	
	Cylinder *cylinder = new Cylinder(Vector(2, -10, -70), 3.0, 10.0, Color::YELLOW);
	sceneObjects.push_back(cylinder);
	
	Cone *cone = new Cone(Vector(20, -10, -65), 5.0, 16.0, Color::YELLOW);
	sceneObjects.push_back(cone);
	
	cubiod();
	
	sceneObjects.push_back(plane);
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glClearColor(0, 0, 0, 1);
}


int main(int argc, char *argv[]) 
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracing");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
