#pragma once

#include <array>

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_EULER
#define GM_INCLUDE_MATRIX
#define GM_INCLUDE_SHAPES
#define GM_INCLUDE_INTERSECTION
#define GM_INCLUDE_QUATERNION
#define GM_INCLUDE_UTILITY
#define GM_INCLUDE_EASE_FUNCTIONS

#if defined(GM_INCLUDE_TYPES)
    #undef NULLPTR
    #undef PI
    #undef DEGREES_TO_RAD
    #undef RAD_TO_DEGREES
    #undef EPSILON
    #undef NEAR_ZERO
    #undef stringify
    #undef glue
    #undef KiloBytes
    #undef MegaBytes
    #undef GigaBytes
    #undef MIN
    #undef MAX
    #undef CLAMP
    #undef SQUARED
    #undef local_persist
    #undef internal
    #undef FIRST_DIGIT
    #undef GET_BIT
    #undef SET_BIT
    #undef UNSET_BIT
    #undef ArrayCount
    #undef OS_DELIMITER
    #undef CRASH
    #undef UNUSED_FUNCTION

    #include <stdint.h>
    #include <stdio.h>
    #include <stdarg.h>
    #include <stdlib.h>
    #include <stdbool.h>
    #include <math.h>
    
    typedef int8_t  s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef uint8_t  u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    #define NULLPTR 0
    #define PI 3.14159265359f
    #define DEGREES_TO_RAD(degrees) ((degrees)*(PI/180.0f))
    #define RAD_TO_DEGREES(rad) ((rad)*(180.0f/PI))
    #define EPSILON 0.0001f
    #define NEAR_ZERO(x) (fabsf(x) <= EPSILON)

    #define stringify(entry) #entry
    #define glue(a, b) a##b

    #define KiloBytes(value) ((size_t)(value) * 1024L)
    #define MegaBytes(value) ((size_t)KiloBytes(value) * 1024L)
    #define GigaBytes(value) ((size_t)MegaBytes(value) * 1024L)

    #define MIN(a, b) (((a) < (b)) ? (a) : (b))
    #define MAX(a, b) (((a) > (b)) ? (a) : (b))
    #define CLAMP(value, min_value, max_value) (MIN(MAX(value, min_value), max_value))
    #define SQUARED(a) ((a) * (a))

    #define local_persist static
    #define internal static

    #if defined _MSC_VER && !defined _CRT_USE_BUILTIN_OFFSETOF
        #define OFFSET_OF(type, member) (size_t)(&(((type*)0)->member))
    #else
        #define OFFSET_OF(type, member) __builtin_offsetof(type, member)
    #endif
    #define FIRST_DIGIT(number) ((int)number % 10);
    #define GET_BIT(number, bit_to_check) ((number & (1 << bit_to_check)) >> bit_to_check)
    #define SET_BIT(number, bit_to_set) number |= (1 << bit_to_set);
    #define UNSET_BIT(number, bit_to_unset) number &= (~(1 << bit_to_unset));

    #define ArrayCount(array) (sizeof(array) / sizeof(array[0]))

    #if defined(__clang__)
        #define UNUSED_FUNCTION __attribute__((used))
    #elif defined(__GNUC__) || defined(__GNUG__)
        #define UNUSED_FUNCTION __attribute__((used))
    #elif defined(_MSC_VER)
        #define UNUSED_FUNCTION
    #endif
#endif

#if defined(GM_INCLUDE_VECTOR)
    struct GM_Vec2 {
        float x;
        float y;

        GM_Vec2() = default;
        explicit GM_Vec2(float fill);
        explicit GM_Vec2(float x, float y);

        float magnitude();
        float magnitudeSquared();
        GM_Vec2 normalize();
        GM_Vec2 scale(float scale);
        GM_Vec2 scale(GM_Vec2 s);
        GM_Vec2 scale(float scale_x, float scale_y);

        /**
         * @brief The return value tells you:
         * -1: the vectors are 180 degrees from eachother in other words they vectors are pointing in opposite directions
         *  0: the vectors are perpendicular or orthogonal to eachother
         *  1: the vectors are going the same direction
         * 
         * @param a 
         * @param b 
         * @return float 
         */
        static float dot(GM_Vec2 a, GM_Vec2 b);
        static float distanfce(GM_Vec2 a, GM_Vec2 b);
        static float distanfceSquared(GM_Vec2 a, GM_Vec2 b);
        static GM_Vec2 lerp(GM_Vec2 a, GM_Vec2 b, float t);

        GM_Vec2 operator+(const GM_Vec2 &right);
        GM_Vec2& operator+=(const GM_Vec2 &right);

        GM_Vec2 operator-(const GM_Vec2 &right);
        GM_Vec2& operator-=(const GM_Vec2 &right);

        GM_Vec2 operator*(const GM_Vec2 &right);
        GM_Vec2& operator*=(const GM_Vec2 &right);

        GM_Vec2 operator/(const GM_Vec2 &right);
        GM_Vec2& operator/=(const GM_Vec2 &right);

        bool operator==(const GM_Vec2 &right);
        bool operator!=(const GM_Vec2 &right);
    };

    struct GM_Vec3 {
        float x;
        float y;
        float z;

        GM_Vec3() = default;
        explicit GM_Vec3(float fill);
        explicit GM_Vec3(float x, float y, float z);

        float magnitude();
        float magnitudeSquared();
        GM_Vec3 normalize();
        GM_Vec3 scale(float scale);
        GM_Vec3 scale(GM_Vec3 s);
        GM_Vec3 scale(float scale_x, float scale_y, float scale_z);

        /**
         * @brief The return value tells you:
         * -1: the vectors are 180 degrees from eachother in other words they vectors are pointing in opposite directions
         *  0: the vectors are perpendicular or orthogonal to eachother
         *  1: the vectors are going the same direction
         * 
         * @param a 
         * @param b 
         * @return float 
         */
        static float dot(GM_Vec3 a, GM_Vec3 b);
        static float distanfce(GM_Vec3 a, GM_Vec3 b);
        static float distanfceSquared(GM_Vec3 a, GM_Vec3 b);
        static GM_Vec3 lerp(GM_Vec3 a, GM_Vec3 b, float t);
        static GM_Vec3 cross(GM_Vec3 a, GM_Vec3 b);

        GM_Vec3 operator+(const GM_Vec3 &right);
        GM_Vec3& operator+=(const GM_Vec3 &right);

        GM_Vec3 operator-(const GM_Vec3 &right);
        GM_Vec3& operator-=(const GM_Vec3 &right);

        GM_Vec3 operator*(const GM_Vec3 &right);
        GM_Vec3& operator*=(const GM_Vec3 &right);

        GM_Vec3 operator/(const GM_Vec3 &right);
        GM_Vec3& operator/=(const GM_Vec3 &right);

        bool operator==(const GM_Vec3 &right);
        bool operator!=(const GM_Vec3 &right);
    };

    struct GM_Vec4 {
        float x;
        float y;
        float z;
        float w;

        GM_Vec4() = default;
        explicit GM_Vec4(float fill);
        explicit GM_Vec4(float x, float y, float z, float w);

        float magnitude();
        float magnitudeSquared();
        GM_Vec4 normalize();
        GM_Vec4 scale(float scale);
        GM_Vec4 scale(GM_Vec4 s);
        GM_Vec4 scale(float scale_x, float scale_y, float scale_z, float scale_w);

        /**
         * @brief The return value tells you:
         * -1: the vectors are 180 degrees from eachother in other words they vectors are pointing in opposite directions
         *  0: the vectors are perpendicular or orthogonal to eachother
         *  1: the vectors are going the same direction
         * 
         * @param a 
         * @param b 
         * @return float 
         */
        static float dot(GM_Vec4 a, GM_Vec4 b);
        static GM_Vec4 lerp(GM_Vec4 a, GM_Vec4 b, float t);
        static float distanfce(GM_Vec4 a, GM_Vec4 b);
        static float distanfceSquared(GM_Vec4 a, GM_Vec4 b);

        GM_Vec4 operator+(const GM_Vec4 &right);
        GM_Vec4& operator+=(const GM_Vec4 &right);

        GM_Vec4 operator-(const GM_Vec4 &right);
        GM_Vec4& operator-=(const GM_Vec4 &right);

        GM_Vec4 operator*(const GM_Vec4 &right);
        GM_Vec4& operator*=(const GM_Vec4 &right);

        GM_Vec4 operator/(const GM_Vec4 &right);
        GM_Vec4& operator/=(const GM_Vec4 &right);

        bool operator==(const GM_Vec4 &right);
        bool operator!=(const GM_Vec4 &right);
    };
#endif

#if defined(GM_INCLUDE_MATRIX)
    // Matrices are COLUMN-MAJOR (OpenGL convention)
    // data[0-3] = Col 0 (X-axis)
    // data[4-7] = Col 1 (Y-axis)
    // data[8-11] = Col 2 (Z-axis)
    // data[12-15] = Col 3 (Translation)
    struct GM_Matrix4 {
        std::array<GM_Vec4, 4> v;

        GM_Matrix4();
        GM_Matrix4(GM_Vec4 r0, GM_Vec4 r1, GM_Vec4 r2, GM_Vec4 r3);
        GM_Matrix4(float m00, float m01, float m02, float m03,
                float m10, float m11, float m12, float m13,
                float m20, float m21, float m22, float m23,
                float m30, float m31, float m32, float m33);
        GM_Matrix4 transpose();

        static GM_Matrix4 identity();
        static GM_Matrix4 scale(GM_Matrix4 mat, float scale);
        static GM_Matrix4 scale(GM_Matrix4 mat, GM_Vec3 s);
        static GM_Matrix4 scale(GM_Matrix4 mat, float scale_x, float scale_y, float scale_z);
        static GM_Matrix4 rotate(GM_Matrix4 mat, float theta, GM_Vec3 axis);
        static GM_Matrix4 rotate(GM_Matrix4 mat, float theta, float rot_x, float rot_y, float rot_z);
        static GM_Matrix4 translate(GM_Matrix4 mat, GM_Vec3 t);
        static GM_Matrix4 translate(GM_Matrix4 mat, float x, float y, float z);
        static GM_Matrix4 transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t);
        static GM_Matrix4 inverse_transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t);

        static GM_Matrix4 perspective(float fov_degrees, float aspect, float near_plane, float far_plane);
        static GM_Matrix4 orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane);
        static GM_Matrix4 lookat(GM_Vec3 position, GM_Vec3 target, GM_Vec3 world_up);

        static GM_Matrix4 inverse(GM_Matrix4 mat, bool* success);

        GM_Matrix4 operator*(const GM_Matrix4 &right);
        GM_Matrix4& operator*=(const GM_Matrix4 &right);

        bool operator==(const GM_Matrix4 &right);
        bool operator!=(const GM_Matrix4 &right);
    };
#endif

#if defined(GM_INCLUDE_QUATERNION)
    typedef struct GM_Quaternion {
        float w;
        GM_Vec3 v;
    } GM_Quaternion;

    #ifdef __cplusplus
        #define GM_QuaternionLit(w, x, y, z) (GM_Quaternion{w, x, y, z})
    #else
        #define GM_QuaternionLit(w, x, y, z) ((GM_Quaternion){w, x, y, z})
    #endif

    GM_Quaternion gm_quat_create(float theta, GM_Vec3 axis);
    GM_Quaternion gm_quat_inverse(GM_Quaternion quat);
    GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2);
    GM_Quaternion gm_quat_from_euler(GM_Vec3 euler_angles_degrees);
    GM_Quaternion gm_quat_from_angle_axis(float angle, GM_Vec3 axis);
    GM_Matrix4 gm_quat_to_mat4(GM_Quaternion q);
    GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec);
    void gm_quat_to_axis_angle(GM_Quaternion quat, float* theta, GM_Vec3* vec);
    GM_Quaternion gm_quat_sub(GM_Quaternion a, GM_Quaternion b);
    GM_Quaternion gm_quat_add(GM_Quaternion a, GM_Quaternion b);
    GM_Quaternion gm_quat_scale(GM_Quaternion a, float scale);
    GM_Quaternion gm_quat_add_scalar(GM_Quaternion a, float scalar);
    GM_Quaternion gm_quat_sub_scalar(GM_Quaternion a, float scalar);
    GM_Quaternion gm_quat_power(GM_Quaternion q1, float t);
    GM_Quaternion gm_quat_normalize(GM_Quaternion q);
    float gm_quat_dot(GM_Quaternion q, GM_Quaternion r);
    GM_Quaternion gm_quat_slerp(GM_Quaternion q, GM_Quaternion r, float t);
    GM_Quaternion gm_quat_look_at(GM_Vec3 position, GM_Vec3 target, GM_Vec3 up);
#endif

#if defined(GM_INCLUDE_UTILITY)
    float gm_lerp(float a, float b, float t);
    float gm_inverse_lerp(float a, float b, float value);
    GM_Vec3 gm_barycentric(GM_Vec3 a, GM_Vec3 b, GM_Vec3 c, float u, float v);

    float gm_remap(float x, float s_min, float s_max, float e_min, float e_max);
    float gm_move_toward(float current, float target, float delta);

    float gm_smoothstep(float edge0, float edge1, float x);
    float gm_smootherstep(float edge0, float edge1, float x);

    bool gm_near_equal(float a, float b);
#endif

#if defined(GM_INCLUDE_EASE_FUNCTIONS)
    // Date: May 18, 2025
    // NOTE(Jovanni): Visualize these at: https://easinfgs.net/
    float gm_ease_in_sinfe(float t);
    float gm_ease_out_sinfe(float t);
    float gm_ease_in_out_sinfe(float t);
    float gm_ease_in_quad(float t);
    float gm_ease_out_quad(float t);
    float gm_ease_in_cubic(float t);
    float gm_ease_out_cubic(float t);
    float gm_ease_in_out_cubic(float t);
    float gm_ease_in_quart(float t);
    float gm_ease_out_quart(float t);
    float gm_ease_in_out_quart(float t);
    float gm_ease_in_quint(float t);
    float gm_ease_out_quint(float t);
    float gm_ease_in_out_quint(float t);
    float gm_ease_in_expo(float t);
    float gm_ease_out_expo(float t);
    float gm_ease_in_out_expo(float t);
    float gm_ease_in_circ(float t);
    float gm_ease_out_circ(float t);
    float gm_ease_in_out_circ(float t);
    float gm_ease_in_back(float t);
    float gm_ease_out_back(float t);
    float gm_ease_in_out_back(float t);
    float gm_ease_in_elastic(float t);
    float gm_ease_out_elastic(float t);
    float gm_ease_in_out_elastic(float t);
    float gm_ease_in_bounce(float t);
    float gm_ease_out_bounce(float t);
    float gm_ease_in_out_bounce(float t);
#endif

#if defined(GM_INCLUDE_EULER)
    GM_Vec2 gm_euler_to_vec2(float yaw, float pitch);
    GM_Vec3 gm_euler_to_vec3(float yaw, float pitch);
#endif

//
// ===================================================== GM_IMPL =====================================================
//

#if defined(GM_IMPL_EULER)
    GM_Vec2 gm_euler_to_vec2(float yaw, float pitch) {
        GM_Vec2 ret = GM_Vec2Lit(0, 0);
        ret.x = cosf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));
        ret.y = sinf(DEGREES_TO_RAD(pitch));

        return ret;
    }

    GM_Vec3 gm_euler_to_vec3(float yaw, float pitch) {
        GM_Vec3 ret = GM_Vec3Lit(0, 0, 0);
        ret.x = cosf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));
        ret.y = sinf(DEGREES_TO_RAD(pitch));
        ret.z = sinf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));

        return ret;
    }
#endif

#if defined(GM_IMPL_MATRIX)
    GM_Matrix4 gm_mat4_identity() {
        GM_Matrix4 ret = {
            .data = {
                1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f
            }
        };
        return ret;
    }



    GM_Matrix4 gm_mat4_perspective(float fov_degrees, float aspect, float near_plane, float far_plane) {
        float fov_radians = DEGREES_TO_RAD(fov_degrees);

        const float t = tanf(fov_radians / 2) * near_plane;
        const float b = -t;
        const float r = t * aspect;
        const float l = -t * aspect;

        const float p = (2.0f * near_plane);

        const float A = p / (r - l);
        const float B = p / (t - b);
        const float C = -((far_plane + near_plane) / (far_plane - near_plane));
        const float D = -((p * far_plane) / (far_plane - near_plane));

        GM_Matrix4 ret = {
            .data = {
                A,  0,  0,  0,
                0,  B,  0,  0,
                0,  0,  C,  D,
                0,  0, -1,  0
            }
        };

        return ret;
    }

    // Found at: https://en.wikipedia.org/wiki/Orthographic_projection
    GM_Matrix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane) {
        const float A = 2.0f / (right - left);
        const float B = 2.0f / (top - bottom);
        const float C = -2.0f / (far_plane - near_plane);
        const float D = -(right + left) / (right - left);
        const float E = -(top + bottom) / (top - bottom);
        const float F = -(far_plane + near_plane) / (far_plane - near_plane);

        GM_Matrix4 ret = {
            .data = {
                A,  0,  0,  D,
                0,  B,  0,  E,
                0,  0,  C,  F,
                0,  0,  0,  1
            }
        };

        return ret;
    }

    // Found at: https://www.khronos.org/opengl/wiki/GluLookAt_code
    GM_Matrix4 gm_mat4_look_at(GM_Vec3 position, GM_Vec3 target, GM_Vec3 world_up) {
        GM_Vec3 forward = gm_vec3_normalize(gm_vec3_sub(position, target));
        GM_Vec3 right   = gm_vec3_normalize(gm_vec3_cross(world_up, forward));
        GM_Vec3 up      = gm_vec3_cross(forward, right);

        GM_Matrix4 rotation = {
            .data = {
                right.x,   right.y,   right.z,   0,
                up.x,      up.y,      up.z,      0,
                forward.x, forward.y, forward.z, 0,
                0.0f,      0.0f,      0.0f,      1.0f
            }
        };
        
        GM_Matrix4 translation = gm_mat4_translate_xyz(gm_mat4_identity(), -position.x, -position.y, -position.z);

        return gm_mat4_mult(rotation, translation);
    }

    #define M_ELEM(mat, row, col) (mat.data[((row) * 4) + (col)])

    float gm_mat3_determinant_helper(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
        return a * (e * i - f * h) -
               b * (d * i - f * g) +
               c * (d * h - e * g);
    }


    GM_Matrix4 gm_mat4_inverse(GM_Matrix4 m, bool* success) {
        if (success) {
            *success = false;
        }

        float m00 = M_ELEM(m, 0, 0), m01 = M_ELEM(m, 0, 1), m02 = M_ELEM(m, 0, 2), m03 = M_ELEM(m, 0, 3);
        float m10 = M_ELEM(m, 1, 0), m11 = M_ELEM(m, 1, 1), m12 = M_ELEM(m, 1, 2), m13 = M_ELEM(m, 1, 3);
        float m20 = M_ELEM(m, 2, 0), m21 = M_ELEM(m, 2, 1), m22 = M_ELEM(m, 2, 2), m23 = M_ELEM(m, 2, 3);
        float m30 = M_ELEM(m, 3, 0), m31 = M_ELEM(m, 3, 1), m32 = M_ELEM(m, 3, 2), m33 = M_ELEM(m, 3, 3);

        float c00 = gm_mat3_determinant_helper(m11, m12, m13, m21, m22, m23, m31, m32, m33);
        float c01 = gm_mat3_determinant_helper(m10, m12, m13, m20, m22, m23, m30, m32, m33);
        float c02 = gm_mat3_determinant_helper(m10, m11, m13, m20, m21, m23, m30, m31, m33);
        float c03 = gm_mat3_determinant_helper(m10, m11, m12, m20, m21, m22, m30, m31, m32);

        float det = m00 * c00 - m01 * c01 + m02 * c02 - m03 * c03;
        if (NEAR_ZERO(det)) {
            return gm_mat4_identity();
        }

        float invDet = 1.0f / det;

        GM_Matrix4 inv;

        // Row 0
        M_ELEM(inv, 0, 0) = invDet * c00;
        M_ELEM(inv, 0, 1) = invDet * (-gm_mat3_determinant_helper(m01, m02, m03, m21, m22, m23, m31, m32, m33));
        M_ELEM(inv, 0, 2) = invDet * gm_mat3_determinant_helper(m01, m02, m03, m11, m12, m13, m31, m32, m33);
        M_ELEM(inv, 0, 3) = invDet * (-gm_mat3_determinant_helper(m01, m02, m03, m11, m12, m13, m21, m22, m23));

        // Row 1
        M_ELEM(inv, 1, 0) = invDet * (-c01);
        M_ELEM(inv, 1, 1) = invDet * gm_mat3_determinant_helper(m00, m02, m03, m20, m22, m23, m30, m32, m33);
        M_ELEM(inv, 1, 2) = invDet * (-gm_mat3_determinant_helper(m00, m02, m03, m10, m12, m13, m30, m32, m33));
        M_ELEM(inv, 1, 3) = invDet * gm_mat3_determinant_helper(m00, m02, m03, m10, m12, m13, m20, m22, m23);

        // Row 2
        M_ELEM(inv, 2, 0) = invDet * c02;
        M_ELEM(inv, 2, 1) = invDet * (-gm_mat3_determinant_helper(m00, m01, m03, m20, m21, m23, m30, m31, m33));
        M_ELEM(inv, 2, 2) = invDet * gm_mat3_determinant_helper(m00, m01, m03, m10, m11, m13, m30, m31, m33);
        M_ELEM(inv, 2, 3) = invDet * (-gm_mat3_determinant_helper(m00, m01, m03, m10, m11, m13, m20, m21, m23));

        // Row 3
        M_ELEM(inv, 3, 0) = invDet * (-c03);
        M_ELEM(inv, 3, 1) = invDet * gm_mat3_determinant_helper(m00, m01, m02, m20, m21, m22, m30, m31, m32);
        M_ELEM(inv, 3, 2) = invDet * (-gm_mat3_determinant_helper(m00, m01, m02, m10, m11, m12, m30, m31, m32));
        M_ELEM(inv, 3, 3) = invDet * gm_mat3_determinant_helper(m00, m01, m02, m10, m11, m12, m20, m21, m22);

        if (success) {
            *success = true;
        }
        return inv;
    }

    #undef M_ELEM // Cleanup the macro

    
    GM_Matrix4 gm_mat4_transpose(GM_Matrix4 m) {
        GM_Matrix4 ret = {0};

        ret.v[0].x = m.v[0].x;
        ret.v[0].y = m.v[1].x;
        ret.v[0].z = m.v[2].x;
        ret.v[0].w = m.v[3].x;

        ret.v[1].x = m.v[0].y;
        ret.v[1].y = m.v[1].y;
        ret.v[1].z = m.v[2].y;
        ret.v[1].w = m.v[3].y;

        ret.v[2].x = m.v[0].z;
        ret.v[2].y = m.v[1].z;
        ret.v[2].z = m.v[2].z;
        ret.v[2].w = m.v[3].z;

        ret.v[3].x = m.v[0].w;
        ret.v[3].y = m.v[1].w;
        ret.v[3].z = m.v[2].w;
        ret.v[3].w = m.v[3].w;      

        return ret;
    }

    
    bool gm_mat4_equal(GM_Matrix4 A, GM_Matrix4 B) {
        for (int i = 0; i < 16; i++) {
            if (!NEAR_ZERO(A.data[i] - B.data[i])) {
                return false;
            }
        }

        return true;
    }


#endif

#if defined(GM_IMPL_SHAPES)
    GM_Rectanfgle2D gm_rectanfgle2d_create(float x, float y, float width, float height) {
        GM_Rectanfgle2D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    GM_RectanfgleReference2D gm_rectanfgle_reference2d_create(GM_Vec2* center, float width, float height) {
        GM_RectanfgleReference2D ret;
        ret.center = center;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    GM_Rectanfgle3D gm_rectanfgle3d_create(float x, float y, float z, float width, float height, float length) {
        GM_Rectanfgle3D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.center.z = z;
        ret.width = width;
        ret.height = height;
        ret.length = length;

        return ret;
    }

    GM_RectanfgleReference3D gm_rectanfgle_reference3d_create(GM_Vec3* center, float width, float height, float length) {
        GM_RectanfgleReference3D ret;
        ret.center = center;
        ret.width = width;
        ret.height = height;
        ret.length = length;

        return ret;
    }

    GM_Circle2D gm_circle2d_create(float x, float y, float radius) {
        GM_Circle2D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.radius = radius;

        return ret;
    }

    GM_CircleReference2D gm_circle_reference2d_create(GM_Vec2* center, float radius) {
        GM_CircleReference2D ret;
        ret.center = center;
        ret.radius = radius;

        return ret;
    }

    GM_Circle3D gm_circle3d_create(float x, float y, float z, float radius) {
        GM_Circle3D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.center.z = z;
        ret.radius = radius;

        return ret;
    }

    GM_CircleReference3D gm_circle_reference3d_create(GM_Vec3* center, float radius) {
        GM_CircleReference3D ret;
        ret.center = center;
        ret.radius = radius;

        return ret;
    }

    /*
    bool gm_rectanfgle_check_aabb_collision(GM_RectanfgleReference2D rect1, GM_RectanfgleReference2D rect2) {
        if (rect1.position->x < rect2.position->x + rect2.width && rect1.position->x + rect1.width > rect2.position->x &&
            rect1.position->y < rect2.position->y + rect2.height && rect1.position->y + rect1.height > rect2.position->y) {
            return true;
        }
        return false;
    }
    */

#endif

#if defined(GM_IMPL_INTERSECTION)
    // Found at: https://imois.in/posts/line-intersections-with-cross-products/
    bool gm_intersection2d_line_line(GM_Vec2 a, GM_Vec2 b, GM_Vec2 c, GM_Vec2 d, GM_Vec2* intersection) {
        GM_Vec3 line1 = gm_vec3_cross(GM_Vec3Lit(a.x, a.y, 1), GM_Vec3Lit(b.x, b.y, 1));
        GM_Vec3 line2 = gm_vec3_cross(GM_Vec3Lit(c.x, c.y, 1), GM_Vec3Lit(d.x, d.y, 1));
        GM_Vec3 solution = gm_vec3_cross(line1, line2);

        if (solution.z == 0) {
            return false;
        } else {
            GM_Vec2 temp = GM_Vec2Lit(solution.x / solution.z, solution.y / solution.z);

            if (intersection) {
                *intersection = temp;
            }

            return (
                (temp.x < a.x) != (temp.x < b.x) && 
                (temp.y < a.y) != (temp.y < b.y) && 
                (temp.x < c.x) != (temp.x < d.x) && 
                (temp.y < c.y) != (temp.y < d.y)
            );
        }
    }


    bool gm_intersection2d_line_aabb(GM_Vec2 p0, GM_Vec2 p1, GM_RectanfgleReference2D aabb, GM_Vec2* inPoint, GM_Vec2* outPoint) {
        float t0 = 0.0f;
        float t1 = 1.0f;
        float dx = p1.x - p0.x;
        float dy = p1.y - p0.y;

        #define CLIP(p, q)                   \
        do {                                 \
            if ((p) == 0.0f && (q) < 0.0f) { \
                return false;                \
            }                                \
            if ((p) < 0.0f) {                \
                float r = (q) / (p);         \
                if (r > t1) return false;    \
                if (r > t0) t0 = r;          \
            } else if ((p) > 0.0f) {         \
                float r = (q) / (p);         \
                if (r < t0) return false;    \
                if (r < t1) t1 = r;          \
            }                                \
        } while (0)                          \

        float half_width  = (aabb.width  / 2.0f);
        float half_height = (aabb.height / 2.0f);

        GM_Vec2 aabb_min = GM_Vec2Lit(aabb.center->x - half_width, aabb.center->y - half_height);
        GM_Vec2 aabb_max = GM_Vec2Lit(aabb.center->x + half_width, aabb.center->y + half_height);

        // X-axis
        CLIP(-dx, p0.x - aabb_min.x); // Left
        CLIP( dx, aabb_max.x - p0.x); // Right

        // Y-axis
        CLIP(-dy, p0.y - aabb_min.y); // Bottom
        CLIP( dy, aabb_max.y - p0.y); // Top

        #undef CLIP

        if (t1 < t0) {
            return false;
        }

        if (inPoint) {
            inPoint->x = p0.x + t0 * dx;
            inPoint->y = p0.y + t0 * dy;
        }

        if (outPoint) {
            outPoint->x = p0.x + t1 * dx;
            outPoint->y = p0.y + t1 * dy;
        }

        return true;
    }

    bool gm_intersection3d_line_aabb(GM_Vec3 p0, GM_Vec3 p1, GM_RectanfgleReference3D aabb, GM_Vec3* inPoint, GM_Vec3* outPoint) {
        float t0 = 0.0f;
        float t1 = 1.0f;

        float dx = p1.x - p0.x;
        float dy = p1.y - p0.y;
        float dz = p1.z - p0.z;

        #define CLIP(p, q)                   \
        do {                                 \
            if ((p) == 0.0f && (q) < 0.0f) { \
                return false;                \
            }                                \
            if ((p) < 0.0f) {                \
                float r = (q) / (p);         \
                if (r > t1) return false;    \
                if (r > t0) t0 = r;          \
            } else if ((p) > 0.0f) {         \
                float r = (q) / (p);         \
                if (r < t0) return false;    \
                if (r < t1) t1 = r;          \
            }                                \
        } while (0)                          \

        float half_width  = (aabb.width  / 2.0f);
        float half_height = (aabb.height / 2.0f);
        float half_length = (aabb.length / 2.0f);

        GM_Vec3 aabb_min = GM_Vec3Lit(aabb.center->x - half_width, aabb.center->y - half_height, aabb.center->z - half_length);
        GM_Vec3 aabb_max = GM_Vec3Lit(aabb.center->x + half_width, aabb.center->y + half_height, aabb.center->z + half_length);

        // X-axis
        CLIP(-dx, p0.x - aabb_min.x); // Left
        CLIP( dx, aabb_max.x - p0.x); // Right

        // Y-axis
        CLIP(-dy, p0.y - aabb_min.y); // Bottom
        CLIP( dy, aabb_max.y - p0.y); // Top

        // Z-axis
        CLIP(-dz, p0.z - aabb_min.z); // Near
        CLIP( dz, aabb_max.z - p0.z); // Far

        #undef CLIP

        if (t1 < t0) {
            return false;
        }

        if (inPoint) {
            inPoint->x = p0.x + t0 * dx;
            inPoint->y = p0.y + t0 * dy;
            inPoint->z = p0.z + t0 * dz;
        }

        if (outPoint) {
            outPoint->x = p0.x + t1 * dx;
            outPoint->y = p0.y + t1 * dy;
            outPoint->z = p0.z + t1 * dz;
        }

        return true;
    }
#endif

#if defined(GM_IMPL_QUATERNION)
    GM_Quaternion gm_quat_create(float theta, GM_Vec3 axis) {
        GM_Quaternion ret = {0};

        float radians = DEGREES_TO_RAD(theta);
        ret.w = cosf(radians / 2.0f);

        if (NEAR_ZERO(ret.w)) {
            ret.w = 0.0f;
        }

        GM_Vec3 norm_axis = gm_vec3_normalize(axis);
        ret.v = gm_vec3_scale(norm_axis, sinf(radians / 2.0f));

        return ret;
    }


    GM_Quaternion gm_quat_inverse(GM_Quaternion quat) {
        GM_Quaternion ret = {0};

        float mag_sq = (quat.w * quat.w) + gm_vec3_dot(quat.v, quat.v);
        if (mag_sq == 0.0f) { 
            return (GM_Quaternion){0,0,0,0};
        }

        ret.w = quat.w / mag_sq;
        ret.v = gm_vec3_scale(quat.v, -1.0f / mag_sq);

        return ret;
    }

    GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2) {
        GM_Quaternion ret;
        ret.w = q1.w * q2.w - gm_vec3_dot(q1.v, q2.v);
        ret.v = gm_vec3_add(gm_vec3_add(gm_vec3_scale(q1.v, q2.w), gm_vec3_scale(q2.v, q1.w)), gm_vec3_cross(q1.v, q2.v));
        return ret;
    }

    GM_Quaternion gm_quat_from_euler(GM_Vec3 euler_angles_degrees) {
        float roll_rad_half = DEGREES_TO_RAD(euler_angles_degrees.x) * 0.5ff;
        float pitch_rad_half = DEGREES_TO_RAD(euler_angles_degrees.y) * 0.5ff;
        float yaw_rad_half = DEGREES_TO_RAD(euler_angles_degrees.z) * 0.5ff;

        float cx = cosf(roll_rad_half);
        float sx = sinf(roll_rad_half);
        float cy = cosf(pitch_rad_half);
        float sy = sinf(pitch_rad_half);
        float cz = cosf(yaw_rad_half);
        float sz = sinf(yaw_rad_half);

        GM_Quaternion q;

        q.w = cx * cy * cz + sx * sy * sz;
        q.v.x = sx * cy * cz - cx * sy * sz;
        q.v.y = cx * sy * cz + sx * cy * sz;
        q.v.z = cx * cy * sz - sx * sy * cz;

        return q;
    }

    GM_Quaternion gm_quat_from_angle_axis(float angle, GM_Vec3 axis) {
        float half_angle = DEGREES_TO_RAD(angle) * 0.5ff;
        float sinf_half = sinf(half_angle);
        float cosf_half = cosf(half_angle);

        GM_Quaternion q;
        axis = gm_vec3_normalize(axis);

        q.w     = cosf_half;
        q.v.x   = axis.x * sinf_half;
        q.v.y   = axis.y * sinf_half;
        q.v.z   = axis.z * sinf_half;

        return q;
    }

    GM_Matrix4 gm_quat_to_mat4(GM_Quaternion q) {
        GM_Matrix4 result = gm_mat4_identity();

        float x2 = q.v.x * q.v.x;
        float y2 = q.v.y * q.v.y;
        float z2 = q.v.z * q.v.z;

        float xy = q.v.x * q.v.y;
        float xz = q.v.x * q.v.z;
        float yz = q.v.y * q.v.z;
        float xw = q.v.x * q.w;
        float yw = q.v.y * q.w;
        float zw = q.v.z * q.w;

        result.data[0] = 1.0f - 2.0f * (y2 + z2);  // m00
        result.data[1] = 2.0f * (xy - zw);         // m01
        result.data[2] = 2.0f * (xz + yw);         // m02
        result.data[3] = 0.0f;                     // m03

        result.data[4] = 2.0f * (xy + zw);         // m10
        result.data[5] = 1.0f - 2.0f * (x2 + z2);  // m11
        result.data[6] = 2.0f * (yz - xw);         // m12
        result.data[7] = 0.0f;                     // m13

        result.data[8] = 2.0f * (xz - yw);         // m20
        result.data[9] = 2.0f * (yz + xw);         // m21
        result.data[10] = 1.0f - 2.0f * (x2 + y2); // m22
        result.data[11] = 0.0f;                    // m23

        result.data[12] = 0.0f;                    // m30
        result.data[13] = 0.0f;                    // m31
        result.data[14] = 0.0f;                    // m32
        result.data[15] = 1.0f;                    // m33

        return result;
    }

    // Rotate a vector with this quaternion (q * p * q_inverse)
    GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec) {
        GM_Quaternion p;
        p.w = 0.0f;
        p.v = vec;

        GM_Quaternion temp = gm_quat_mult(quat, p);
        return gm_quat_mult(temp, gm_quat_inverse(quat)).v;
    }

    void gm_quat_to_axis_angle(GM_Quaternion quat, float* theta, GM_Vec3* vec) {
        quat = gm_quat_normalize(quat);
        float sinf_half_theta = gm_vec3_magnitude(quat.v);

        if (vec) {
            if (sinf_half_theta < EPSILON) {
                *vec = GM_Vec3Lit(1, 0, 0);
            } else {
                *vec = gm_vec3_scale(quat.v, 1.0f / sinf_half_theta);
            }
        }

        if (theta) {
            // Clamp w to [-1, 1] to avoid NaNs due to precision issues
            float w = CLAMP(quat.w, -1.0f, 1.0f);
            *theta = 2.0f * acosf(w);
            *theta = RAD_TO_DEGREES(*theta);
        }
    }

    GM_Quaternion gm_quat_scale(GM_Quaternion q, float scale) {
        GM_Quaternion ret;

        ret.w   = q.w   * scale;
        ret.v.x = q.v.x * scale;
        ret.v.y = q.v.y * scale;
        ret.v.z = q.v.z * scale;

        return ret;
    }

    GM_Quaternion gm_quat_add_scalar(GM_Quaternion a, float scalar) {
        GM_Quaternion ret;

        ret.w   = a.w   + scalar;
        ret.v.x = a.v.x + scalar;
        ret.v.y = a.v.y + scalar;
        ret.v.z = a.v.z + scalar;

        return ret;
    }

    GM_Quaternion gm_quat_sub_scalar(GM_Quaternion a, float scalar) {
        GM_Quaternion ret;

        ret.w   = a.w   - scalar;
        ret.v.x = a.v.x - scalar;
        ret.v.y = a.v.y - scalar;
        ret.v.z = a.v.z - scalar;

        return ret;
    }

    GM_Quaternion gm_quat_add(GM_Quaternion a, GM_Quaternion b) {
        GM_Quaternion ret;

        ret.w   = a.w   + b.w;
        ret.v.x = a.v.x + b.v.x;
        ret.v.y = a.v.y + b.v.y;
        ret.v.z = a.v.z + b.v.z;

        return ret;
    }

    GM_Quaternion gm_quat_sub(GM_Quaternion a, GM_Quaternion b) {
        GM_Quaternion ret;

        ret.w   = a.w   - b.w;
        ret.v.x = a.v.x - b.v.x;
        ret.v.y = a.v.y - b.v.y;
        ret.v.z = a.v.z - b.v.z;

        return ret;
    }

    GM_Quaternion gm_quat_power(GM_Quaternion q1, float t) {
        float a;
        GM_Vec3 n;

        gm_quat_to_axis_angle(q1, &a, &n);
        return gm_quat_create(a * t, n);
    }

    GM_Quaternion gm_quat_normalize(GM_Quaternion q) {
        GM_Vec4 temp = gm_vec4_normalize(GM_Vec4Lit(q.w, q.v.x, q.v.y, q.v.z));
        
        GM_Quaternion ret;
        ret.w = temp.x;
        ret.v.x = temp.y;
        ret.v.y = temp.z;
        ret.v.z = temp.w;

        return ret;
    }

    float gm_quat_dot(GM_Quaternion a, GM_Quaternion b) {
        float dot = a.w   * b.w   +
                    a.v.x * b.v.x +
                    a.v.y * b.v.y +
                    a.v.z * b.v.z;

        return dot;
    }

    GM_Quaternion gm_quat_slerp(GM_Quaternion q, GM_Quaternion r, float t) {
        q = gm_quat_normalize(q);
        r = gm_quat_normalize(r);
        float dot = gm_quat_dot(q, r);

        if (dot < 0.0f) {
            r = gm_quat_scale(r, -1.0f);
            dot = -dot;
        }

        if (dot > 0.9995f) {
            GM_Quaternion lerp = gm_quat_add(q, gm_quat_scale(gm_quat_sub(r, q), t));
            return gm_quat_normalize(lerp);
        }

        float theta_0 = acosf(dot);
        float theta = theta_0 * t;

        GM_Quaternion q3 = gm_quat_sub(r, gm_quat_scale(q, dot));
        q3 = gm_quat_normalize(q3);

        GM_Quaternion term1 = gm_quat_scale(q, cosf(theta));
        GM_Quaternion term2 = gm_quat_scale(q3, sinf(theta));
        return gm_quat_add(term1, term2);
    }

    GM_Quaternion gm_quat_look_at(GM_Vec3 position, GM_Vec3 target, GM_Vec3 up) {
        GM_Vec3 direction = gm_vec3_normalize(gm_vec3_sub(target, position));
        GM_Vec3 forward = {0, 0, 1};

        if (NEAR_ZERO(gm_vec3_magnitude_squared(gm_vec3_sub(direction, forward)))) {
            return GM_QuaternionLit(1, 0, 0, 0);
        }

        GM_Vec3 axis = gm_vec3_cross(forward, direction);
        float angle = acosf(gm_vec3_dot(forward, direction));

        if (NEAR_ZERO(gm_vec3_magnitude_squared(axis))) {
            axis = gm_vec3_normalize(up);
        } else {
            axis = gm_vec3_normalize(axis);
        }

        return gm_quat_from_angle_axis(angle, axis);
    }
#endif

#if defined(GM_IMPL_UTILITY)
    float gm_lerp(float a, float b, float t) {
        return a + ((b - a) * t);
    }

    float gm_inverse_lerp(float a, float b, float value) {
        if (NEAR_ZERO(a - b)) {
            return 0.0f; // Avoid division by zero
        }

        return (value - a) / (b - a);
    }

    float gm_remap(float x, float s_min, float s_max, float e_min, float e_max) {
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

    GM_Vec3 gm_barycentric(GM_Vec3 a, GM_Vec3 b, GM_Vec3 c, float u, float v);
    float gm_smoothstep(float edge0, float edge1, float x);
    float gm_smootherstep(float edge0, float edge1, float x);

    bool gm_near_equal(float a, float b) {
        return fabsf(a - b) < EPSILON;
    }
#endif

#if defined(GM_IMPL_EASE_FUNCTIONS)
    float gm_ease_in_sinfe(float t) {
        return 1.0f - cosf((t * PI) / 2.0f);
    }

    float gm_ease_out_sinfe(float t) {
        return sinf((t * PI) / 2.0f);
    }

    float gm_ease_in_out_sinfe(float t) {
        return -(cosf(PI * t) - 1.0f) / 2.0f;
    }

    float gm_ease_in_quad(float t) {
        return t * t;
    }

    float gm_ease_out_quad(float t) {
        return 1.0f - (1.0f - t) * (1.0f - t);
    }

    float gm_ease_in_cubic(float t) {
        return t * t * t;
    }

    float gm_ease_out_cubic(float t) {
        return 1.0f - powf(1.0f - t, 3.0f);
    }

    float gm_ease_in_out_cubic(float t) {
        return t < 0.5ff ? 4.0f * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 3.0f) / 2.0f;
    }

    float gm_ease_in_quart(float t) {
        return t * t * t * t;
    }

    float gm_ease_out_quart(float t) {
        return 1.0f - powf(1.0f - t, 4.0f);
    }

    float gm_ease_in_out_quart(float t) {
        return t < 0.5ff ? 8.0f * t * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 4.0f) / 2.0f;
    }

    float gm_ease_in_quint(float t) {
        return t * t * t * t * t;
    }

    float gm_ease_out_quint(float t) {
        return 1.0f + powf(t - 1.0f, 5.0f);
    }

    float gm_ease_in_out_quint(float t) {
        return t < 0.5ff ? 16.0f * t * t * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 5.0f) / 2.0f;
    }

    float gm_ease_in_expo(float t) {
        return t == 0.0f ? 0.0f : powf(2.0f, 10.0f * t - 10.0f);
    }

    float gm_ease_out_expo(float t) {
        return t == 1.0f ? 1.0f : 1.0f - powf(2.0f, -10.0f * t);
    }

    float gm_ease_in_out_expo(float t) {
        if (t == 0.0f) return 0.0f;
        if (t == 1.0f) return 1.0f;
        if (t < 0.5ff) return powf(2.0f, 20.0f * t - 10.0f) / 2.0f;
        return (2.0f - powf(2.0f, -20.0f * t + 10.0f)) / 2.0f;
    }

    float gm_ease_in_circ(float t) {
        return 1.0f - sqrtf(1.0f - (t * t));
    }

    float gm_ease_out_circ(float t) {
        return sqrtf(1.0f - powf(t - 1.0f, 2.0f));
    }

    float gm_ease_in_out_circ(float t) {
        if (t < 0.5ff) return (1.0f - sqrtf(1.0f - powf(2.0f * t, 2.0f))) / 2.0f;
        return (sqrtf(1.0f - powf(-2.0f * t + 2.0f, 2.0f)) + 1.0f) / 2.0f;
    }

    float gm_ease_in_back(float t) {
        float c1 = 1.70158f;
        float c3 = c1 + 1.0f;
        return c3 * t * t * t - c1 * t * t;
    }

    float gm_ease_out_back(float t) {
        float c1 = 1.70158f;
        float c3 = c1 + 1.0f;
        return 1.0f + c3 * powf(t - 1.0f, 3.0f) + c1 * powf(t - 1.0f, 2.0f);
    }

    float gm_ease_in_out_back(float t) {
        float c1 = 1.70158f;
        float c2 = c1 * 1.5f25f;
        if (t < 0.5ff) return (powf(2.0f * t, 2.0f) * ((c2 + 1.0f) * 2.0f * t - c2)) / 2.0f;
        return (powf(2.0f * t - 2.0f, 2.0f) * ((c2 + 1.0f) * (t * 2.0f - 2.0f) + c2) + 2.0f) / 2.0f;
    }

    float gm_ease_in_elastic(float t) {
        float c4 = (2.0f * PI) / 3.0f;
        if (t == 0.0f) return 0.0f;
        if (t == 1.0f) return 1.0f;
        return -powf(2.0f, 10.0f * t - 10.0f) * sinf((t * 10.0f - 10.75f) * c4);
    }

    float gm_ease_out_elastic(float t) {
        float c4 = (2.0f * PI) / 3.0f;
        if (t == 0.0f) return 0.0f;
        if (t == 1.0f) return 1.0f;
        return powf(2.0f, -10.0f * t) * sinf((t * 10.0f - 0.75f) * c4) + 1.0f;
    }

    float gm_ease_in_out_elastic(float t) {
        float c5 = (2.0f * PI) / 4.5ff;
        if (t == 0.0f) return 0.0f;
        if (t == 1.0f) return 1.0f;
        if (t < 0.5ff) return -(powf(2.0f, 20.0f * t - 10.0f) * sinf((20.0f * t - 11.125f) * c5)) / 2.0f;
        return (powf(2.0f, -20.0f * t + 10.0f) * sinf((20.0f * t - 11.125f) * c5)) / 2.0f + 1.0f;
    }

    float gm_ease_in_bounce(float t) {
        return 1.0f - gm_ease_out_bounce(1.0f - t);
    }

    float gm_ease_out_bounce(float t) {
        float n1 = 7.5f625f;
        float d1 = 2.75f;
        if (t < 1.0f / d1) {
            return n1 * t * t;
        } else if (t < 2.0f / d1) {
            return n1 * (t -= 1.5ff / d1) * t + 0.75f;
        } else if (t < 2.5ff / d1) {
            return n1 * (t -= 2.25f / d1) * t + 0.9375f;
        } else {
            return n1 * (t -= 2.625f / d1) * t + 0.984375f;
        }
    }

    float gm_ease_in_out_bounce(float t) {
        return t < 0.5ff
            ? (1.0f - gm_ease_out_bounce(1.0f - 2.0f * t)) / 2.0f
            : (1.0f + gm_ease_out_bounce(2.0f * t - 1.0f)) / 2.0f;
    }
#endif