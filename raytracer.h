// raytracer.h is for educational purposes 
// not for serious use 
// im learning computer graphics since few weeks 
// so mistakes are expected
// the math logic from articles and community discussion (stackoverflow) 
// all the math functions are from : "LA_MATH.h" my own library 

#ifndef RAYTRACER_
#define RAYTRACER_

// list of includes
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LA_MATH.h"

typedef struct {
    Vec3f origin;
    Vec3f direction;
} Ray;


typedef struct {
    Vec3f center;
    float radius;
} Sphere;

typedef struct {
    Vec3f point;   
    Vec3f normal;  
} Plane;

typedef struct {
    Vec3f origin;
    Vec3f forward;
    Vec3f right;
    Vec3f up;
    float aspect;
    float scale;
} Camera;

typedef struct {
    void *data;  // data (sphere,plane)
    bool (*doIntersect)(Ray, void *data, float *t); // intersection
    Vec3f (*getNormal)(Vec3f point , void* data); // get the normal vector
    
} Hittable;


Ray generateRay(float x , float y ,int w, int h, Camera cam);
Vec3f computeLighting(Vec3f point,Vec3f normal,Vec3f lightpos);
bool sphereIntersect(Ray ray, void *data, float *t);
Vec3f sphereNormal(Vec3f point, void *data);
bool planeIntersect(Ray ray, void *data, float *t);
Vec3f planeNormal(Vec3f point, void *data);

Camera makeCamera(Vec3f origin ,Vec3f lookAt,Vec3f up,float fov,float aspectratio);

Hittable makeSphere(Sphere* sp);
Hittable makePlane(Plane* p);



#endif // RAYTRACER_
#ifdef RAYTRACER_IMPLEMENTATION

Ray generateRay(float x , float y ,int w,int h, Camera cam){
    float px = (2 * ((x + 0.5f) / (float)w) - 1) * cam.aspect * cam.scale;
    float py = (1 - 2 * ((y + 0.5f) / (float)h)) * cam.scale;
    Vec3f dir = vec3f_add(
                    vec3f_add(
                        vec3f_scale(cam.right, px),
                        vec3f_scale(cam.up, py)
                    ), 
                    cam.forward
                );
    return (Ray){ cam.origin, vec3f_normalize(dir) };

}

Vec3f computeLighting(Vec3f point, Vec3f normal, Vec3f lightpos) {
    Vec3f lightDir = vec3f_normalize(vec3f_sub(lightpos, point));
    float diffuse = fmaxf(0.0f, vec3f_dot(normal, lightDir));
    return vec3f_scale((Vec3f){1, 1, 1}, diffuse); // white light
}

bool sphereIntersect(Ray ray, void *data, float *t) {
    Sphere *sp = (Sphere *)data;
    Vec3f oc = vec3f_sub(ray.origin, sp->center);
    float a = vec3f_dot(ray.direction, ray.direction);
    float b = 2.0f * vec3f_dot(oc, ray.direction);
    float c = vec3f_dot(oc, oc) - sp->radius * sp->radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    float discrSQRT = sqrtf(discriminant);
    float t0 = (-b - discrSQRT) / (2 * a);
    float t1 = (-b + discrSQRT) / (2 * a);
    *t = (t0 < t1) ? t0 : t1;
    return true;
}

Vec3f sphereNormal(Vec3f point, void *data) {
    Sphere *sp = (Sphere *)data;
    return vec3f_normalize(vec3f_sub(point, sp->center));
}

bool planeIntersect(Ray ray, void *data, float *t) {
    Plane *pl = (Plane *)data;
    float denom = vec3f_dot(pl->normal, ray.direction);
    if (fabsf(denom) > 1e-6) {
        Vec3f p0l0 = vec3f_sub(pl->point, ray.origin);
        *t = vec3f_dot(p0l0, pl->normal) / denom;
        return *t >= 0;
    }
    return false;
}

Vec3f planeNormal(Vec3f point, void *data) {
    Plane *pl = (Plane *)data;
    return pl->normal; // constant normal
}


Hittable makeSphere(Sphere* sp){
    return (Hittable) {sp,sphereIntersect,sphereNormal};
}
Hittable makePlane(Plane* p){
    return (Hittable){p,planeIntersect,planeNormal};
}

Camera makeCamera(Vec3f origin ,Vec3f lookAt,Vec3f up,float fov,float aspectratio){

    Vec3f forward = vec3f_normalize(vec3f_sub(lookAt, origin));
    Vec3f right = vec3f_normalize(vec3f_cross(forward, up));
    Vec3f camUp = vec3f_cross(right, forward);
    float scale = tanf(fov * 0.5f * LAMATH_PI / 180.0f);
    return (Camera){ origin, forward, right, camUp, aspectratio, scale };
}

#endif // RAYTRACER_IMPLEMENTATION


