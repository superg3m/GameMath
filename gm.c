#define GM_IMPL
#include "gm.h"

float  gm_dot_product() {
}

Matrix4 gm_mat4_rotation_x(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    Matrix4 ret = {
        .data = {
            1, 0,  0, 0,
            0, c, -s, 0,
            0, s,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

Matrix4 gm_mat4_rotation_y(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    Matrix4 ret = {
        .data = {
            c, 0, -s, 0,
            0, 1,  0, 0,
            s, 0,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

Matrix4 gm_mat4_rotation_z(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    Matrix4 ret = {
        .data = {
            c, -s, 0, 0,
            s, c,  0, 0,
            0, 0,  1, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

// Z is negative going away from the viewer here so Right-handed coordinate system OpenGL
Matrix4 gm_mat4_perspective(float fov_degrees, float aspect, float near, float far) {
    float fov_radians = DEGREES_TO_RAD(fov_degrees);
    float p = 1 / (tan(fov_radians) / 2.0f);

    const range = near - far;
    const A = (-far - near) / range; 
    const B = (2 * far * near) / range; 

    Matrix4 ret = {
        .data = {
            p / aspect, 0, 0, 0,
            0, p, 0, 0,
            0, 0, A, B,
            0, 0, 1, 0
        }
    };
    
    return ret;
}