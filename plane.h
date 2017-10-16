//
//  plane.h
//  dh2323
//
//  Created by Johan Thorell on 2017-06-01.
//  Copyright © 2017 Johan Thorell. All rights reserved.
//

#ifndef plane_h
#define plane_h


using glm::vec3;
using glm::dot;

class Plane {
public:
    vec3 normal;
    vec3 point;
    
    void setValues(vec3 n, vec3 p){
        this->normal = n;
        this->point = p;
    }
    
    float signedDistance(vec3 &p){
        float dist = dot(normal, p) + dot(-normal, point);
        return dist;
    }
    
};

#endif /* plane_h */
