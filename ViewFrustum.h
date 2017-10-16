//
//  ViewFrustum.h
//  dh2323
//
//  Created by Johan Thorell on 2017-06-01.
//  Copyright © 2017 Johan Thorell. All rights reserved.
//

#ifndef ViewFrustum_h
#define ViewFrustum_h

using namespace std;
using glm::vec2;
using glm::vec3;
using glm::ivec2;
using glm::mat3;
using glm::dot;
using glm::length;
using glm::fastNormalize;

class ViewFrustum {
private:
    enum {
        NEARP = 0, FARP, TOP, BOTTOM, LEFT, RIGHT
    };
    
public:
    enum {OUTSIDE, INTERSECT, INSIDE};
    Plane planes[6];
    float nearDist, farDist, nearW, nearH, farW, farH;
    
    void computePlanes(vec3 &p, vec3 &dir, vec3 &up, vec3 &right);
    void setCamInternals(float angle){
        float tang = tan(angle);
        nearH = nearW = 2 * nearDist * tang;
        farH = farW = 2 * farDist * tang;
    }
    bool sphereInFrustum(vec3 &p, float radius);
    void setDist(float nearDist, float farDist){
        this->nearDist = nearDist;
        this->farDist = farDist;
    }
    
};

// PRE: middle point of polygon, and radius of a sphere boundary
bool ViewFrustum::sphereInFrustum(vec3 &p, float radius) {
    
    float distance;
    
    for(int i=0; i < 6; i++) {
        distance = planes[i].signedDistance(p);
        
        if (distance <= -radius){
            //cout<<"false"<<endl;
            return false;
        }
        
        //        else if (distance < radius)
        //            result =  INTERSECT;
    }
    //cout<<"true"<<endl;
    return true;
}

void ViewFrustum::computePlanes(vec3 &p, vec3 &d, vec3 &up, vec3 &right){
    vec3 aux, normal;
    
    vec3 worldForward = d;
    fastNormalize(worldForward);
    
    // X axis of camera with given "up" vector and Z axis
    vec3 worldRight = right;
    
    // the real "up" vector is the cross product of Z and X
    vec3 worldUp = up;

    fastNormalize(worldForward);
    fastNormalize(worldRight);
    fastNormalize(worldUp);
    
    
    vec3 nearCenter = p + worldForward * nearDist;
    vec3 farCenter = p + worldForward * farDist;
    
    
    planes[NEARP].setValues(worldForward, nearCenter);
    planes[FARP].setValues(-worldForward, farCenter);
    
    // Camera position (p) will be on all top, bottom, right, left planes
    aux = (nearCenter + worldRight * nearW/2.f) - p;
    normal = fastNormalize(glm::cross(worldUp,aux));
    planes[RIGHT].setValues(normal, p);
    
    aux = (nearCenter - worldRight * nearW/2.f) - p;
    normal = fastNormalize(glm::cross(aux,worldUp)); // Normal need to point other way than RIGHT
    planes[LEFT].setValues(normal, p);
    
    aux = (nearCenter + worldUp * nearH/2.f) - p;
    normal = fastNormalize(glm::cross(aux,worldRight));
    planes[TOP].setValues(normal, p);
    
    aux = (nearCenter - worldUp * nearH/2.f) - p;
    normal = fastNormalize(glm::cross(worldRight, aux));
    planes[BOTTOM].setValues(normal, p);
}

#endif /* ViewFrustum_h */
