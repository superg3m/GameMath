#pragma once

#include <array>
#include <limits>

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_AABB
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
    #undef stringify
    #undef glue
    #undef KiloBytes
    #undef MegaBytes
    #undef GigaBytes
    #undef MIN
    #undef MAX
    #undef CLAMP
    #undef SQUARED
    #undef FIRST_DIGIT
    #undef GET_BIT
    #undef SET_BIT
    #undef UNSET_BIT
    #undef ArrayCount
    #undef PLATFORM_WINDOWS
    #undef NOMINMAX
    #undef WIN32_LEAN_AND_MEAN
    #undef WIN32_LEAN_AND_MEAN
    #undef PLATFORM_APPLE
    #undef PLATFORM_LINUX
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

        GM_Vec2();
        explicit GM_Vec2(float fill);
        explicit GM_Vec2(float x, float y);

        float magnitude();
        float magnitudeSquared();
        GM_Vec2 normalize();
        GM_Vec2 scale(float scale) const;
        GM_Vec2 scale(GM_Vec2 s) const;
        GM_Vec2 scale(float scale_x, float scale_y) const;

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
        static float distance(GM_Vec2 a, GM_Vec2 b);
        static float distanceSquared(GM_Vec2 a, GM_Vec2 b);
        static GM_Vec2 lerp(GM_Vec2 a, GM_Vec2 b, float t);
        static GM_Vec2 euler(float yaw, float pitch);

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

    typedef struct GM_Vec4 GM_Vec4;
    struct GM_Vec3 {
        union {
            struct {
                float x;
                float y;
                float z;
            };

            struct {
                float r;
                float g;
                float b;
            };
        };

        GM_Vec3();
        explicit GM_Vec3(float fill);
        explicit GM_Vec3(float x, float y, float z);
        explicit GM_Vec3(GM_Vec2 v, float z);
        explicit GM_Vec3(GM_Vec4 v);

        float magnitude();
        float magnitudeSquared();
        GM_Vec3 normalize();
        GM_Vec3 scale(float scale) const;
        GM_Vec3 scale(GM_Vec3 s) const;
        GM_Vec3 scale(float scale_x, float scale_y, float scale_z) const;

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
        static float distance(GM_Vec3 a, GM_Vec3 b);
        static float distanceSquared(GM_Vec3 a, GM_Vec3 b);
        static GM_Vec3 lerp(GM_Vec3 a, GM_Vec3 b, float t);
        static GM_Vec3 cross(GM_Vec3 a, GM_Vec3 b);
        static GM_Vec3 euler(float yaw, float pitch);

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
        union {
            struct {
                float x;
                float y;
                float z;
                float w;
            };

            struct {
                float r;
                float g;
                float b;
                float a;
            };
        };

        GM_Vec4();
        explicit GM_Vec4(float fill);
        explicit GM_Vec4(float x, float y, float z, float w);
        explicit GM_Vec4(GM_Vec3 v, float w);

        float magnitude();
        float magnitudeSquared();
        GM_Vec4 normalize();
        GM_Vec4 scale(float scale) const;
        GM_Vec4 scale(GM_Vec4 s) const;
        GM_Vec4 scale(float scale_x, float scale_y, float scale_z, float scale_w) const;

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
        static float distance(GM_Vec4 a, GM_Vec4 b);
        static float distanceSquared(GM_Vec4 a, GM_Vec4 b);

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

#if defined(GM_INCLUDE_AABB)
    struct GM_AABB {
        GM_Vec3 min;
        GM_Vec3 max;

        GM_AABB();
        GM_AABB(GM_Vec3 min, GM_Vec3 max);
        GM_AABB(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z);

        GM_Vec3 getCenter();
        GM_Vec3 getExtents();

        static GM_AABB fromCenterExtents(GM_Vec3 center, GM_Vec3 extents);
        static bool intersection(GM_AABB aabb, GM_Vec3 p0, GM_Vec3 p1);
    };
#endif

#if defined(GM_INCLUDE_MATRIX)
    // Matrices are COLUMN-MAJOR (OpenGL convention)
    // data[0-3] = Col 0 (X-axis)
    // data[4-7] = Col 1 (Y-axis)
    // data[8-11] = Col 2 (Z-axis)
    // data[12-15] = Col 3 (Translation)
    typedef struct GM_Quaternion GM_Quaternion;

    struct GM_Matrix4 {
        std::array<GM_Vec4, 4> v;

        GM_Matrix4();
        GM_Matrix4(GM_Vec4 r0, GM_Vec4 r1, GM_Vec4 r2, GM_Vec4 r3);
        GM_Matrix4(float m00, float m01, float m02, float m03,
                float m10, float m11, float m12, float m13,
                float m20, float m21, float m22, float m23,
                float m30, float m31, float m32, float m33);
        GM_Matrix4 transpose();
        GM_Matrix4 inverse(bool &success);

        static GM_Matrix4 identity();
        static GM_Matrix4 fromColumnMajor(const float mat[16]);
        static GM_Matrix4 scale(GM_Matrix4 mat, float scale);
        static GM_Matrix4 scale(GM_Matrix4 mat, GM_Vec3 s);
        static GM_Matrix4 scale(GM_Matrix4 mat, float scale_x, float scale_y, float scale_z);
        static GM_Matrix4 rotate(GM_Matrix4 mat, float theta, GM_Vec3 axis);
        static GM_Matrix4 rotate(GM_Matrix4 mat, float theta, float rot_x, float rot_y, float rot_z);
        static GM_Matrix4 rotate(GM_Matrix4 mat, GM_Quaternion quat);
        static GM_Matrix4 translate(GM_Matrix4 mat, GM_Vec3 t);
        static GM_Matrix4 translate(GM_Matrix4 mat, float x, float y, float z);
        static GM_Matrix4 transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t);
        static GM_Matrix4 inverse_transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t);
        static void decompose(GM_Matrix4 mat, GM_Vec3* out_position, GM_Quaternion* out_orientation, GM_Vec3* out_scale);

        static GM_Matrix4 perspective(float fov_degrees, float aspect, float near_plane, float far_plane);
        static GM_Matrix4 orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane);
        static GM_Matrix4 lookat(GM_Vec3 position, GM_Vec3 target, GM_Vec3 world_up);

        GM_Vec4 operator*(const GM_Vec4 &right);
        GM_Matrix4 operator*(const GM_Matrix4 &right);
        GM_Matrix4& operator*=(const GM_Matrix4 &right);

        bool operator==(const GM_Matrix4 &right);
        bool operator!=(const GM_Matrix4 &right);
    };
#endif

#if defined(GM_INCLUDE_QUATERNION)
    struct GM_Quaternion {
        float w;
        GM_Vec3 v;

        GM_Quaternion();
        GM_Quaternion(float theta, GM_Vec3 axis);
        GM_Quaternion(float theta, float x, float y, float z);

        GM_Quaternion inverse();
        GM_Quaternion normalize();
        GM_Quaternion scale(float scale);
        GM_Matrix4 toMatrix4();
        void toAngleAxis(float &theta, GM_Vec3 &vec);

        static GM_Quaternion identity();
        static GM_Quaternion literal(float w, GM_Vec3 axis);
        static GM_Quaternion literal(float w, float x, float y, float z);
        static GM_Quaternion fromEuler(GM_Vec3 euler_angles_degrees);
        static GM_Quaternion fromAngleAxis(float w, GM_Vec3 axis);
        static GM_Quaternion fromRotationMatrix(const float m[16]);
        static GM_Quaternion fromRotationMatrix(GM_Matrix4 mat);
        static GM_Quaternion slerp(GM_Quaternion q, GM_Quaternion r, float t);
        static float dot(GM_Quaternion a, GM_Quaternion b);

        GM_Quaternion operator+(const GM_Quaternion &right);
        GM_Quaternion& operator+=(const GM_Quaternion &right);
        
        GM_Quaternion operator-(const GM_Quaternion &right);
        GM_Quaternion& operator-=(const GM_Quaternion &right);

        GM_Quaternion operator*(const GM_Quaternion &right);
        GM_Quaternion& operator*=(const GM_Quaternion &right);

        GM_Vec3 operator*(const GM_Vec3 &right);

        bool operator==(const GM_Quaternion &right);
        bool operator!=(const GM_Quaternion &right);
    };
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