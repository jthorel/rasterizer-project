#include <iostream>
#include "glm/glm.hpp"
#include "glm/ext.hpp"
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <math.h>
#include "plane.h"
#include "ViewFrustum.h"


#define ANG2RAD 3.14159265358979323846/180.0

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::dot;
using glm::length;
using glm::fastNormalize;
using glm::fastLength;


// ----------------------------------------------------------------------------
// STRUCTS

struct Pixel {
    int x;
    int y;
    float zinv;
    //vec3 illumination;
    vec3 pos3d;
};

struct Vertex {
    vec3 position;
    /*vec3 normal;
    vec3 reflectance;*/
};



// ----------------------------------------------------------------------------
// GLOBAL VARIABLES

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;

// rotation matrix, identity matrix for start (no rotation, looking along z-axis).
mat3 R(1,0,0,
       0,1,0,
       0,0,1);

vec3 lightPosition = vec3(0,-0.5,-0.7);

// Fixed light coefficient
vec3 diffuseLightPower = 0.5f*vec3(1,1,1);
vec3 specularLightPower = 1.f*vec3(1,1,1);
vec3 indirectLightPowerPerArea = 0.1f*vec3( 1, 1, 1 );

bool blinn = true;
bool specular = false;
bool culling = false;
bool lambertianOff = false;

// Focal length same as width, will give FOV of 60 degrees.
vec3 cameraPos(0,0,-3.001);
float focalLength = SCREEN_WIDTH;
float camAngle = atan((SCREEN_WIDTH/2)/focalLength);
ViewFrustum viewFrustum;
glm::mat4 projectionMatrix(1.f, 0.f, float(SCREEN_WIDTH/2)/(focalLength), 0.f,
                           0.f, 1.f, float(SCREEN_HEIGHT/2)/(focalLength), 0.f,
                           0.f, 0.f, 1.f, 0.f,
                           0.f, 0.f, 1.f/(focalLength), 0.f);


float yaw = .0f;
float pitch = .0f;
float roll = .0f;


float moveSpeed = 1.f;
float x_sensitivity = 1.f; // Radians per second!
float y_sensitivity = x_sensitivity;

// mouse delta movements
int dx;
int dy;

vec3 currentColor = vec3(1,1,1);
vec3 currentNormal;
vec3 currentReflectance;
float currentSpecularHardness = 32.0f;

float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];

    // light
// VAR TVUNGEn FLYTTA LIGHTPOS FRÅN DENNA RADEN ANNARS ÄNDRADES DEN HELA TIDEN
// EVERY GLOBAL FLOAT DECLARED AFTER THIS LINE WILL CHANGE WHEN MOVING CAMERA
// FLOATING POINT MEMORY BUG?


// ----------------------------------------------------------------------------
// FUNCTION DECLARATIONS

void Update();
void Draw();
void DrawCulling();
void VertexShader(const Vertex& v, Pixel& p);
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result);
//void DrawPolygonEdges( const vector<vec3>& vertices );
//void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color );
void rotateCamera();
void computePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels);
void drawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels);
void drawPolygon( const vector<Vertex>& vertices );
void pixelShader( const Pixel p );
void lambertianLight(vec3 lightdir, float lightDistance, vec3& lambertLight);
void phongLight(vec3 lightdir, float lightDistance, vec3 viewDir, vec3& lambertLight, vec3& specularLight);
void blinnPhongLight(vec3 lightdir, float lightDistance, vec3 viewDir, vec3& lambertLight, vec3& specularLight);
void calcTriangleMidRad(Triangle triangle, vec3& middle, float& radius);


// ----------------------------------------------------------------------------
// MAIN PROGRAM

int main( int argc, char* argv[] )
{
	LoadTestModel( triangles );
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    viewFrustum.setDist(1.f,10.f);
    viewFrustum.setCamInternals(camAngle);
    
    
    // Main Loop
	while( NoQuitMessageSDL() )
	{
		Update();
        if(culling){
            DrawCulling();
        } else {
            Draw();
        }
	}

    // Save screenshot when quitting
	SDL_SaveBMP( screen, "screenshot.bmp" );
	return 0;
}


void calcTriangleMidRad(Triangle triangle, vec3& middle, float& radius){
    middle = (triangle.v0+triangle.v1+triangle.v2)/3.f;
    float l1 = fastLength(triangle.v0 - middle);
    float l2 = fastLength(triangle.v1 - middle);
    float l3 = fastLength(triangle.v2 - middle);
    
    // Find max of all 3?
    radius = max(l1,l2);
}


// Lambertian diffuse light only
void lambertianLight(vec3 lightdir, float lightDistance, vec3& lambertLight){
    float NdotL = dot(currentNormal,lightdir);
    // Don't really need 4*Pi as the r^2 is the important varying one
    //float piSquaredDistance = float(4*M_PI*lightDistance*lightDistance);
    float piSquaredDistance = float(lightDistance*lightDistance);
    
    if(NdotL > 0){
        lambertLight = (diffuseLightPower * NdotL) / piSquaredDistance;
    }
}


// Lambertian with phong specular
void phongLight(vec3 lightdir, float lightDistance, vec3 viewDir, vec3& lambertLight, vec3& specularLight){
    float NdotL = dot(currentNormal,lightdir);
    float piSquaredDistance = float(lightDistance*lightDistance);
    
    if(NdotL > 0){
        lambertLight = (diffuseLightPower * NdotL) / piSquaredDistance;
        vec3 reflectDir = lightdir - 2.f*dot(lightdir, currentNormal)*currentNormal;
        float RdotV = dot(viewDir, reflectDir);
        RdotV = pow(RdotV, currentSpecularHardness);
        specularLight = (specularLightPower * max(RdotV,0.f)) / piSquaredDistance;
    }

}

// Lambertian with blinn-phong specular
void blinnPhongLight(vec3 lightdir, float lightDistance, vec3 viewDir, vec3& lambertLight, vec3& specularLight){
    float NdotL = dot(currentNormal,lightdir);
    float piSquaredDistance = float(lightDistance*lightDistance);
    
    if(NdotL > 0){
        lambertLight = (diffuseLightPower * NdotL) / piSquaredDistance;
        vec3 H = fastNormalize(lightdir + viewDir);
        float HdotN = dot(currentNormal, H);
        HdotN = pow(HdotN, 2.f*currentSpecularHardness);
        specularLight = (specularLightPower * max(HdotN, 0.f)) / piSquaredDistance;
    }
}

// Draws a Pixel on the screen with it's x,y-coordinates if it's not occluded (z-buffer)
void pixelShader( const Pixel p ){
    
    // Is the pixel closer than what is currently stored in the z-buffer?
    if( p.zinv > depthBuffer[p.y][p.x] ){

        vec3 lambertLight;
        vec3 specularLight;
        
        //Set this pixel to the closest
        depthBuffer[p.y][p.x] = p.zinv;
        
        /* PART 3 */
        // Direction to light, 3D-position inverted back with z-inverse.
        vec3 lightdir = lightPosition - (p.pos3d/p.zinv);
        //vec3 viewDir = fastNormalize(vec3((p.x-SCREEN_WIDTH/2)/(focalLength*p.zinv), (p.y-SCREEN_HEIGHT/2)/(focalLength*p.zinv), (cameraPos[2] - 1/p.zinv)));
        vec3 viewDir = fastNormalize(cameraPos - p.pos3d/p.zinv);
        
        float lightDistance = fastLength(lightdir); // Inverse square root? Can be done faster
        lightdir = fastNormalize(lightdir); // Normalize
        
        if(specular){
            if(blinn){
                blinnPhongLight(lightdir, lightDistance, viewDir, lambertLight, specularLight);
            } else {
                phongLight(lightdir, lightDistance, viewDir, lambertLight, specularLight);
            }
        } else {
            lambertianLight(lightdir, lightDistance, lambertLight);
        }
        
        // just for test
        if(lambertianOff){
            lambertLight = vec3(0,0,0);
        }
        
        // Light at pixel using interpolated 3D-position
        vec3 illumination = currentReflectance * (specularLight + lambertLight + indirectLightPowerPerArea);
        /* */
        
        // Put pixel on screen/surface/canvas
        PutPixelSDL( screen, p.x, p.y, illumination );
    }
}


// Draws a polygon from vertices
// PRE: vector of Vertex-structs
void drawPolygon( const vector<Vertex>& vertices ){
    
    size_t V = vertices.size();
    vector<Pixel> vertexPixels( V );
    
    // Project vertices and save projected pixels
    for( size_t i=0; i<V; ++i )
        VertexShader( vertices[i], vertexPixels[i] );
    
    // LINE CLIPPING HERE
    // Make a frustrum from the camera with 6 planes
    // have a center of the frustrum planes
    // take a middle point of the triangle and make a sphere for border
    // check if that sphere is inside all planes
    // longest distance from center to a vertex is the radius?
    // dot(trianglecenterposition, normalofplane) + distance of plane + radius
    //      (normal of plane pointing out)
    // http://www.lighthouse3d.com/tutorials/view-frustum-culling/geometric-approach-extracting-the-planes/
    // Finding a sphere that contains all the vertices of the car is an easy task (the average of the vertices is the center of the sphere, and the radius is the distance to the farthest vertex), and testing a sphere is extremely fast as it will be shown next.
    
    // BOunding boxes can be tried for complete models or polygons. Our models are separate triangles
    // So bounding sphere for each polygon will be used.
    
    // _Signed_ distance
    // dist = A*rx + B*ry + C*rz + D = n . r  + D
    // D = – n . p0
    
    
    vector<Pixel> leftPixels;
    vector<Pixel> rightPixels;
    
    computePolygonRows( vertexPixels, leftPixels, rightPixels );
    drawPolygonRows( leftPixels, rightPixels );
}


// Draws the polygon-surface/face/inside its edges
void drawPolygonRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels){
    
    // Loop through rows, both vectors _should_ be the size of number of rows (add check?)
    for(int i = 0; i < leftPixels.size(); ++i){
        int numPixels = glm::abs(leftPixels[i].x - rightPixels[i].x) + 1;
        
        // Interpolate row on the polygon
        vector<Pixel> line(numPixels);
        Interpolate(leftPixels[i], rightPixels[i], line);
        
        // Draw the row
        for(Pixel p : line){
            // Need to be inside the screen or z-buffer will go out of range
            if(p.x > 0 && p.x < SCREEN_WIDTH && p.y > 0 && p.y < SCREEN_HEIGHT){
                pixelShader(p);
            }
        }
    }
}


// Computes the pixel-rows inside the polygon to be drawn
void computePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels){
    
    size_t V = vertexPixels.size();
    
    // 1. Find max and min y-value of the polygon
    // and compute the number of rows it occupies.
    
    // Set first value
    int maxY = vertexPixels[0].y;
    int minY = vertexPixels[0].y;
    
    // Check other vertices
    for(Pixel pixel : vertexPixels){
        minY = (pixel.y < minY) ? pixel.y : minY;
        maxY = (pixel.y > maxY) ? pixel.y : maxY;

    }

    // 2. Resize leftPixels and rightPixels
    // so that they have an element for each row.
    int ROWS = maxY - minY + 1;
    leftPixels.resize(ROWS);
    rightPixels.resize(ROWS);
    
    // 3. Initialize the x-coordinates in leftPixels
    // to some really large value and the x-coordinates
    // in rightPixels to some really small value.
    for( int i=0; i<ROWS; ++i ){
        leftPixels[i].x = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
    }
    
    // 4. Loop through all edges of the polygon and use
    // linear interpolation to find the x-coordinate for
    // each row it occupies. Update the corresponding
    // values in rightPixels and leftPixels.
    
    // Loop over all vertices and draw the edge from it to the next vertex:
    for( int i=0; i<V; ++i ) {
        
        // The next vertex
        Pixel a = vertexPixels[i];
        Pixel b = vertexPixels[(i+1)%V];
        
        // Need to find the delta 'separately', as Pixel is not a vector-type with overloaded operators
        ivec2 delta = glm::abs(ivec2(a.x - b.x, a.y - b.y));
        int numPixels = glm::max( delta.x, delta.y ) + 1;
        
        vector<Pixel> edge( numPixels );
        Interpolate( a, b, edge );
        
        
        for(Pixel pixel : edge){
            
            // get the intrinsic row of the polygon
            int row = pixel.y - minY;
            // check polygon boundaries
            if(pixel.x < leftPixels[row].x) {
                leftPixels[row] = pixel;
            }
            
            if(pixel.x > rightPixels[row].x) {
                rightPixels[row] = pixel;
            }
        }
        
    }

}



// EARLIER PART OF LAB
//void DrawPolygonEdges( const vector<vec3>& vertices )
//{
//	size_t V = vertices.size();
//    
//	// Transform each vertex from 3D world position to 2D image position: vector<ivec2> projectedVertices( V );
//	vector<ivec2> projectedVertices( V );
//	for( size_t i=0; i<V; ++i ){
//		VertexShader( vertices[i], projectedVertices[i] );
//	}
//	
//	// Loop over all vertices and draw the edge from it to the next vertex:
//	for( int i=0; i<V; ++i ) {
//		int j = (i+1)%V; // The next vertex
//		vec3 color( 1, 1, 1 );
//		DrawLineSDL( screen, projectedVertices[i], projectedVertices[j], currentColor);
//	}
//}


// Draws an interpolated line between PROJECTED POINTS WITH X,Y-coordinates
//void DrawLineSDL( SDL_Surface* surface, ivec2 a, ivec2 b, vec3 color ){
//	ivec2 delta = glm::abs( a - b );
//	int pixels = glm::max( delta.x, delta.y ) + 1;
//
//	// You can then get the pixel positions of the line by calling the Interpolation function:
//	vector<ivec2> line( pixels );
//	Interpolate( a, b, line );
//
//
//	for(int i=0; i<line.size(); ++i){
//		PutPixelSDL(surface, line[i].x, line[i].y, color);
//	}
//
//
//}




// Linear interpolation between 2 points
//
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ){

	size_t N = result.size();
    float n = float(max(int(N)-1,1));
    
    float xstep = float(b.x-a.x) / n;
    float ystep = float(b.y-a.y) / n;
    float zstep = float(b.zinv - a.zinv) / n;
    vec3 pos3dstep = (b.pos3d - a.pos3d) / n;
    
    
    /* PART2
    vec3 lightstep = (b.illumination - a.illumination)/float(n);
     */
    

    /* PART 1
	vec2 step = vec2(b-a) / float(max(int(N)-1,1));
    vec2 current(a);
     */
    
    Pixel current = a;

	for( size_t i=0; i<N; ++i ){
		result[i] = current;
		current.x = a.x + i*xstep;
        current.y = a.y + i*ystep;
        current.zinv += zstep;
        current.pos3d += pos3dstep;
        
        
        /* PART 2
         current.illumination += lightstep;
         */
        
	}
}

void VertexShader(const Vertex& v, Pixel& p){
	// Do transformations here on camera
    // Translate the coordinate system to origin at camera and rotation
	// P' = (P - C)*R (5)
	vec3 vTranslated = (v.position - cameraPos)*R;
    
    p.zinv = 1.f/vTranslated.z;
    p.pos3d = v.position/vTranslated.z;
    
    glm::vec4 d = glm::vec4(vTranslated.x, vTranslated.y, vTranslated.z, 1.f);
    glm::vec4 f = d*projectionMatrix;
    //cout<<vTranslated.z*1.f/focalLength<<endl;
    //cout<<f.w<<endl;
    p.x = f[0]/f[3];
    p.y = f[1]/f[3];
    //cout<<p.x<<endl;
    
    
	// Project on x,y screen with range [0,screensize]
	// x = fX/Z + W/2 (3)
	// y = fY/Z + H/2 (4)
//	p.x = (focalLength * (vTranslated.x/vTranslated.z)) + (SCREEN_WIDTH/2);
//	p.y = (focalLength * (vTranslated.y/vTranslated.z)) + (SCREEN_HEIGHT/2);
    
    // PART 3
    // SAVE INVERSE Z AND POSITION (z-inversed)

    
    
    
    /* PART 2
    vec3 lightdir = lightPosition - v.position;
    float lightLength = length(lightdir);
    p.illumination = v.reflectance * ((lightPower * max((dot(fastNormalize(lightdir),v.normal)), .0f))
    / float(4*M_PI*lightLength*lightLength)+indirectLightPowerPerArea);
     
     */

}

void rotateCamera(){
    
    // Yaw
    mat3 RotateY = mat3(cos(yaw), 0, sin(yaw),
                        0, 1, 0,
                        -sin(yaw), 0, cos(yaw));
    
    // Pitch
    mat3 RotateX = mat3(1, 0, 0,
                        0, cos(pitch), -sin(pitch),
                        0, sin(pitch), cos(pitch));
    
    // Roll
    mat3 RotateZ = mat3(cos(roll), -sin(roll), 0,
                        sin(roll), cos(roll), 0,
                        0, 0, 1);
    
    // Final rotation matrix, correct order
    R = RotateZ*RotateY*RotateX;
}




void Update() {
    bool updateCamera = false;
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t) / 1000;
	t = t2;
	//cout << "Render time: " << dt << " ms." << endl;
    cout<<"FPS: "<<int(1/dt)<<endl;
    
    vec3 right(R[0][0], R[0][1], R[0][2]);
    vec3 down(R[1][0], R[1][1], R[1][2]);
    vec3 up = -down;
    vec3 forward(R[2][0], R[2][1], R[2][2]);
	
    
    
    
    // MOUSE MOVEMENT
    //SDL_GetRelativeMouseState( &dx, &dy );
    //yaw -= dx*.01f;
    //pitch += dy*.01f;
    
    
    // GET KEYPRESSES
    Uint8* keystate = SDL_GetKeyState(0);
        if(keystate[SDLK_UP]){
            pitch -= y_sensitivity*dt;
        }
        
        if( keystate[SDLK_DOWN] ){
            pitch += y_sensitivity*dt;
        }
        
        if( keystate[SDLK_RIGHT] ){
            yaw -= x_sensitivity*dt;
            
        }
        
        if( keystate[SDLK_LEFT] ){
            yaw += x_sensitivity*dt;
        }
        
        if( keystate[SDLK_RSHIFT] )
            ;
        
        if( keystate[SDLK_RCTRL] )
            ;
        
        if( keystate[SDLK_w] )
            cameraPos += forward*moveSpeed*dt;
        
        if( keystate[SDLK_s] )
            cameraPos -= forward*moveSpeed*dt;
        
        if( keystate[SDLK_d] )
            cameraPos += right*moveSpeed*dt;
        
        if( keystate[SDLK_a] )
            cameraPos -= right*moveSpeed*dt;
        
        if( keystate[SDLK_e] )
            roll -= x_sensitivity*dt;
        
        if( keystate[SDLK_q] )
            roll += x_sensitivity*dt;
        
        if( keystate[SDLK_y] )
            lightPosition.z += 0.1f;
        
        if( keystate[SDLK_h] )
            lightPosition.z -= 0.1f;
        
        if( keystate[SDLK_j] )
            lightPosition.x += 0.1f;
        
        if( keystate[SDLK_g] )
            lightPosition.x -= 0.1f;
        
        
        // RESET POSITIONS AND ROTATIONS
        if( keystate[SDLK_SPACE] ){
            roll = 0;
            yaw = 0;
            pitch = 0;
            cameraPos = vec3(0,0,-3.0001);
            lightPosition = vec3(0, -0.5, -0.7);
        }
        
        if(keystate[SDLK_z]){
            if(specular){
                specular = false;
            } else {
                specular = true;
            }
        }
        
        if(keystate[SDLK_x]){
            if(blinn){
                blinn = false;
            } else {
                blinn = true;
            }
        }
        
        if(keystate[SDLK_c]){
            if(culling){
                culling = false;
                cout<<"culling off"<<endl;
            } else {
                culling = true;
                cout<<"culling on"<<endl;
            }
        }
        
        if(keystate[SDLK_v]){
            if(lambertianOff){
                lambertianOff = false;
                cout<<"lambertian on"<<endl;
            } else {
                lambertianOff = true;
                cout<<"lambertian off"<<endl;
            }
        }
    

    // Update camera rotation matrix

        rotateCamera();
        viewFrustum.computePlanes(cameraPos, forward, up, right);
    
    
}

void DrawCulling()
{
    // Clear depth buffer
    for( int y=0; y<SCREEN_HEIGHT; ++y )
        for( int x=0; x<SCREEN_WIDTH; ++x )
            depthBuffer[y][x] = 0;
    
    
	SDL_FillRect( screen, 0, 0 );

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);
	
	for( int i=0; i<triangles.size(); ++i )
	{
        vec3 middle;
        float radius;
        calcTriangleMidRad(triangles[i], middle, radius);
        
        if(viewFrustum.sphereInFrustum(middle, radius)){

            vector<Vertex> vertices(3);
            vec3 dirToTriangle = (cameraPos - middle);
            if(dot(dirToTriangle,triangles[i].normal) > 0){
            
                currentReflectance = triangles[i].color;
                currentNormal = triangles[i].normal;
                
                vertices[0].position = triangles[i].v0;
                vertices[1].position = triangles[i].v1;
                vertices[2].position = triangles[i].v2;
                
                
                /* PART 2
                 Save normal and color (reflectance) of the current triangle
                 to each vertex.
                 
                 for (int j = 0; j<3; ++j){
                 vertices[j].normal = triangles[i].normal;
                 vertices[j].reflectance = currentColor;
                 } */
                
                
                /* PART 1
                 Add drawing
                 
                 for(int v=0; v<3; ++v){
                 ivec2 projPos;
                 VertexShader(vertices[v], projPos);
                 vec3 color(1,1,1);
                 PutPixelSDL(screen, projPos.x, projPos.y, color);
                 }*/
                
                drawPolygon( vertices );
            }
        }
	}
	
	if ( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void Draw()
{
    // Clear depth buffer
    for( int y=0; y<SCREEN_HEIGHT; ++y )
        for( int x=0; x<SCREEN_WIDTH; ++x )
            depthBuffer[y][x] = 0;
    
    
    SDL_FillRect( screen, 0, 0 );
    
    if( SDL_MUSTLOCK(screen) )
        SDL_LockSurface(screen);
    
    for( int i=0; i<triangles.size(); ++i )
    {
        vec3 middle;
        float radius;
        calcTriangleMidRad(triangles[i], middle, radius);
        
        if(viewFrustum.sphereInFrustum(middle, radius)){
            vector<Vertex> vertices(3);
            
            currentReflectance = triangles[i].color;
            currentNormal = triangles[i].normal;
            
            vertices[0].position = triangles[i].v0;
            vertices[1].position = triangles[i].v1;
            vertices[2].position = triangles[i].v2;
            
            
            /* PART 2
             Save normal and color (reflectance) of the current triangle
             to each vertex.
             
             for (int j = 0; j<3; ++j){
             vertices[j].normal = triangles[i].normal;
             vertices[j].reflectance = currentColor;
             } */
            
            
            /* PART 1
             Add drawing
             
             for(int v=0; v<3; ++v){
             ivec2 projPos;
             VertexShader(vertices[v], projPos);
             vec3 color(1,1,1);
             PutPixelSDL(screen, projPos.x, projPos.y, color);
             }*/
            
            drawPolygon( vertices );
        }

    }
    
    if ( SDL_MUSTLOCK(screen) )
        SDL_UnlockSurface(screen);
    
    SDL_UpdateRect( screen, 0, 0, 0, 0 );
}


