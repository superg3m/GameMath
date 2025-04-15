#pragma once

#include <math.h>
// Vector
// Matrix

#define PI 3.14159265359
#define DEGREES_TO_RAD(degrees) (degrees*(PI/180))
#define RAD_TO_DEGREES(rad) (rad*(180/PI))

float gm_range_map(float x, float s_min, float s_max, float e_min, float e_max);
float gm_move_toward(float current, float target, float delta);


// row major
typedef struct GM_Vec2 {
    union {
        struct {
            float x;
            float y;
        };

        float data[2];
    };
} GM_Vec2;

typedef struct GM_Vec3 {
    union {
        struct {
            float x;
            float y;
            float z;
        };

        float data[3];
    };
} GM_Vec3;

typedef struct GM_Vec4 {
    union {
        struct {
            float x;
            float y;
            float z;
            float w;
        };

        float data[4];
    };
} GM_Vec4;

typedef struct GM_Matix4 {
    union {
        float data[16]; 
        GM_Vec4 v[4];
    };
} GM_Matix4;


// Creation
GM_Matix4 gm_mat4_identity();
GM_Matix4 mat4_translation(GM_Vec3 t);
GM_Matix4 mat4_scale(GM_Vec3 s);

// Convert to radians in the implementation
GM_Matix4 mat4_rotation_x(float degrees);
GM_Matix4 mat4_rotation_y(float degress);
GM_Matix4 mat4_rotation_z(float degress);


GM_Matix4 mat4_perspective(float fov_degrees, float aspect, float near, float far);
GM_Matix4 mat4_orthographic(float left, float right, float bottom, float top, float near, float far);
GM_Matix4 mat4_look_at(GM_Vec3 eye, GM_Vec3 center, GM_Vec3 up);

// Math
GM_Matix4 gm_mat4_mult(GM_Matix4 A, GM_Matix4 B);
//GM_Vec4 mat4_mult_vec4(GM_Matix4 m, GM_Vec4 v);
GM_Matix4 mat4_inverse(GM_Matix4 m);
GM_Matix4 mat4_transpose(GM_Matix4 m);

// float  gm_dot_product()
