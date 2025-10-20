#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LA_MATH.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define SCR_WIDTH 800
#define SCR_HEIGHT 600

typedef struct {
    Vec3f origin;
    Vec3f direction;
} Ray;

typedef struct {
    Vec3f center;
    float radius;
} Sphere;

// Ray-sphere intersection
bool doIntersect(Ray ray, Sphere sp, float *t) {
    Vec3f oc = vec3f_sub(ray.origin, sp.center);
    float a = vec3f_dot(ray.direction, ray.direction);
    float b = 2.0f * vec3f_dot(oc, ray.direction);
    float c = vec3f_dot(oc, oc) - sp.radius * sp.radius;
    float discriminant = b * b - 4 * a * c;
    if (discriminant < 0) return false;
    float discrSQRT = sqrtf(discriminant);
    float t0 = (-b - discrSQRT) / (2 * a);
    float t1 = (-b + discrSQRT) / (2 * a);
    *t = (t0 < t1) ? t0 : t1;
    return true;
}

int main() {
    unsigned char *image = malloc(SCR_WIDTH * SCR_HEIGHT * 3);
    if (!image) {
        fprintf(stderr, "Image buffer allocation failed\n");
        return 1;
    }

    Vec3f origin = {0, 0, 0};
    Vec3f lookAt = {0, 0, -1};
    Vec3f up = {0, 1, 0};
    float fov = 90.0f;
    float ar = (float)SCR_WIDTH / (float)SCR_HEIGHT;
    float scale = tanf((fov * 0.5f) * LAMATH_PI / 180.0f);

    Vec3f forward = vec3f_normalize(vec3f_sub(lookAt, origin));
    Vec3f right = vec3f_normalize(vec3f_cross(forward, up));
    Vec3f nup = vec3f_cross(right, forward);

    Sphere sphere = {{0, 0, -5}, 1};
    Vec3f light = {5, 5, -10};

    for (int y = 0; y < SCR_HEIGHT; y++) {
        for (int x = 0; x < SCR_WIDTH; x++) {
            float px = (2 * ((x + 0.5f) / (float)SCR_WIDTH) - 1) * ar * scale;
            float py = (1 - 2 * ((y + 0.5f) / (float)SCR_HEIGHT)) * scale;
            Vec3f ray_dir = vec3f_add(vec3f_add(vec3f_scale(right, px), vec3f_scale(nup, py)), forward);
            ray_dir = vec3f_normalize(ray_dir);

            Ray ray = {origin, ray_dir};
            float t;
            int idx = (y * SCR_WIDTH + x) * 3;

            if (doIntersect(ray, sphere, &t)) {
                Vec3f hit_point = vec3f_add(ray.origin, vec3f_scale(ray.direction, t));
                Vec3f normal = vec3f_normalize(vec3f_sub(hit_point, sphere.center));
                Vec3f light_dir = vec3f_normalize(vec3f_sub(light, hit_point));
                float diffuse = fmaxf(0, vec3f_dot(normal, light_dir));
                int intensity = (int)(diffuse * 255);
                image[idx + 0] = intensity;
                image[idx + 1] = intensity;
                image[idx + 2] = intensity;
            } else {
                image[idx + 0] = 20;
                image[idx + 1] = 20;
                image[idx + 2] = 20;
            }
        }
    }

    stbi_write_png("output.png", SCR_WIDTH, SCR_HEIGHT, 3, image, SCR_WIDTH * 3);
    free(image);
    printf("Rendered image saved to output.png\n");
    return 0;
}
