#define GM_IMPL
#include "gm.h"

float gm_range_map(float x, float s_min, float s_max, float e_min, float e_max) {
    x = CLAMP(x, s_min, s_max);
    float s_ratio = (x - s_min) / (s_max - s_min);
    return e_min + (s_ratio * (e_max - e_min));
}

float gm_move_toward(float current, float target, float delta) {
    float diff = target - current;

    if (fabsf(diff) <= delta) {
        return target;
    }

    return current + (diff > 0 ? delta : -delta);
}


float gm_v3_dot_product(GM_Vec3 A, GM_Vec3 B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

float gm_v4_dot_product(GM_Vec4 A, GM_Vec4 B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z) + (A.w * B.w);
}

float gm_v4_magnitude(GM_Vec4 A) {
    return sqrt((A.x*A.x) + (A.x*A.x) + (A.z*A.z));
}

GM_Vec4 gm_v4_normalize(GM_Vec4 A) {
    GM_Vec4 ret;
    const float magnitude = gm_v4_magnitude(A);
    ret.x = A.x / magnitude;
    ret.y = A.y / magnitude;
    ret.z = A.z / magnitude;
    ret.w = A.w / magnitude;

    return ret;
}

GM_Matix4 mat4_identity() {
    GM_Matix4 ret = {
        .data = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_mult(GM_Matix4 A, GM_Matix4 B) {
    GM_Matix4 ret = {0};
    
    for (int i = 0; i < 16; i++) {
        const row = i / 4;
        const column = i % 4;

        const GM_Vec4 a_row = A.v[i / 4];
        const GM_Vec4 b_column = {
            B.data[column + (0 * 4)], 
            B.data[column + (1 * 4)], 
            B.data[column + (2 * 4)], 
            B.data[column + (3 * 4)]
        };

        ret.data[i] = gm_v4_dot_product(a_row, b_column);
    }

    return ret;
}

GM_Matix4 mat4_translation(GM_Vec3 t) {
    GM_Matix4 ret = {
        .data = {
            1, 0, 0, t.x,
            0, 1, 0, t.y,
            0, 0, 1, t.z,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 mat4_scale(GM_Vec3 s) {
    GM_Matix4 ret = {
        .data = {
            s.x, 0, 0, 0,
            0, s.y, 0, 0,
            0, 0, s.z, 0,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_x(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
        .data = {
            1, 0,  0, 0,
            0, c, -s, 0,
            0, s,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_y(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
        .data = {
            c, 0, -s, 0,
            0, 1,  0, 0,
            s, 0,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_z(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
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
GM_Matix4 gm_mat4_perspective(float fov_degrees, float aspect, float near, float far) {
    float fov_radians = DEGREES_TO_RAD(fov_degrees);
    float p = 1 / (tan(fov_radians) / 2.0f);

    const range = near - far;
    const A = (-far - near) / range; 
    const B = (2 * far * near) / range; 

    GM_Matix4 ret = {
        .data = {
            p / aspect, 0, 0, 0,
            0, p, 0, 0,
            0, 0, A, B,
            0, 0, 1, 0
        }
    };
    
    return ret;
}