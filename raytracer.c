#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LA_MATH.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define RAYTRACER_IMPLEMENTATION
#include "raytracer.h"

#define SCR_WIDTH 800
#define SCR_HEIGHT 600

// template for raytracer.h 
// author : gWall
// date : 10/21/2025
int main() {
    unsigned char *image = malloc(SCR_WIDTH * SCR_HEIGHT * 3);
    if (!image) {
        fprintf(stderr, "Image allocation failed\n");
        return 1;
    }

    // Scene config
    Vec3f camPos = {0, 0, 0};
    Vec3f camTarget = {0, 0, -1};
    Vec3f camUp = {0, 1, 0};
    float fov = 90.0f;
    float aspect = (float)SCR_WIDTH / (float)SCR_HEIGHT;
    Camera cam = makeCamera(camPos, camTarget, camUp, fov, aspect);

    // set the light vector
    Vec3f light = {55, 15, -10};

    // the objects tp render in the scene 
    Sphere sphere = {{0, 0, -5}, 1};
    Plane plane = {{0, -1, 0}, {0, 1, 0}};  // horizontal plane at y = -1
    
    Hittable objects[2];
    objects[0] = makeSphere(&sphere);
    objects[1] = makePlane(&plane);
    int objectCount = 2;
    // render to the screen
    for (int y = 0; y < SCR_HEIGHT; y++) {
        for (int x = 0; x < SCR_WIDTH; x++) {
            Ray ray = generateRay((float)x, (float)y, SCR_WIDTH, SCR_HEIGHT, cam);
            float closestT = INFINITY;
            Hittable *hitObject = NULL;
            Vec3f hitPoint, normal;

            for (int i = 0; i < objectCount; i++) {
                float t;
                if (objects[i].doIntersect(ray, objects[i].data, &t)) {
                    if (t < closestT) {
                        closestT = t;
                        hitObject = &objects[i];
                    }
                }
            }

            int idx = (y * SCR_WIDTH + x) * 3;
            if (hitObject) {
                hitPoint = vec3f_add(ray.origin, vec3f_scale(ray.direction, closestT));
                normal = hitObject->getNormal(hitPoint, hitObject->data);
                Vec3f color = computeLighting(hitPoint, normal, light);
                image[idx + 0] = (unsigned char)(fminf(color.x, 1.0f) * 255);
                image[idx + 1] = (unsigned char)(fminf(color.y, 1.0f) * 255);
                image[idx + 2] = (unsigned char)(fminf(color.z, 1.0f) * 255);
            } else {
                // Background color
                image[idx + 0] = 20;
                image[idx + 1] = 20;
                image[idx + 2] = 20;
            }
        }
    }

    stbi_write_png("output.png", SCR_WIDTH, SCR_HEIGHT, 3, image, SCR_WIDTH * 3);
    free(image);
    printf("Image rendered to output.png\n");
    return 0;
}

