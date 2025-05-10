#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#if defined(CKG_INCLUDE_TYPES)
    #include <stdint.h>
    #include <stdio.h>
    #include <stdarg.h>
    #include <stdlib.h>
    #include <stdbool.h>

    typedef int8_t  s8;
    typedef int16_t s16;
    typedef int32_t s32;
    typedef int64_t s64;

    typedef uint8_t  u8;
    typedef uint16_t u16;
    typedef uint32_t u32;
    typedef uint64_t u64;

    #define NULLPTR 0
    #define PI 3.14159265359
    #define DEGREES_TO_RAD(degrees) ((degrees)*(PI/180))
    #define RAD_TO_DEGREES(rad) ((rad)*(180/PI))

    #define stringify(entry) #entry
    #define glue(a, b) a##b

    #define KiloBytes(value) ((size_t)(value) * 1024L)
    #define MegaBytes(value) ((size_t)KiloBytes(value) * 1024L)
    #define GigaBytes(value) ((size_t)MegaBytes(value) * 1024L)

    #define MIN(a, b) (((a) < (b)) ? (a) : (b))
    #define MAX(a, b) (((a) > (b)) ? (a) : (b))
    #define CLAMP(value, min_value, max_value) (MIN(MAX(value, min_value), max_value))

    #define local_persist static
    #define internal static

    // Date: July 12, 2024
    // TODO(Jovanni): Test this to make sure its actually works but it makes sense to me
    #define OFFSET_OF(type, member) sizeof((size_t)(&(((type*)0)->member)))
    #define FIRST_DIGIT(number) ((int)number % 10);
    #define GET_BIT(number, bit_to_check) ((number & (1 << bit_to_check)) >> bit_to_check)
    #define SET_BIT(number, bit_to_set) number |= (1 << bit_to_set);
    #define UNSET_BIT(number, bit_to_unset) number &= (~(1 << bit_to_unset));

    #define ArrayCount(array) (sizeof(array) / sizeof(array[0]))

    #define PLATFORM_MAX_PATH 256

    #if defined(_WIN32)
        #define NOMINMAX
        #define WIN32_LEAN_AND_MEAN
        #include <windows.h>
        #define PLATFORM_WINDOWS
        #define OS_DELIMITER '\\'
        #define CRASH __debugbreak()
    #elif defined(__APPLE__)
        #define PLATFORM_APPLE
        #define OS_DELIMITER '/'
        #define CRASH __builtin_trap()
    #elif defined(__linux__) || defined(__unix__) || defined(__POSIX__)
        #define PLATFORM_LINUX
        #define OS_DELIMITER '/'
        #define CRASH __builtin_trap()
    #else
        #error "Unknown Platform???"
    #endif

    #if defined(__clang__)
        #define UNUSED_FUNCTION __attribute__((used))
        #define WRITE_FENCE() __asm__ volatile("" ::: "memory"); __asm__ volatile("sfence" ::: "memory")
        #define READ_FENCE() __asm__ volatile("" ::: "memory");
    #elif defined(__GNUC__) || defined(__GNUG__)
        #define UNUSED_FUNCTION __attribute__((used))
        #define WRITE_FENCE() __asm__ volatile("" ::: "memory"); __asm__ volatile("sfence" ::: "memory")
        #define READ_FENCE() __asm__ volatile("" ::: "memory");
    #elif defined(_MSC_VER)
        #define UNUSED_FUNCTION
        #define WRITE_FENCE() _WriteBarrier(); _mm_sfence()
        #define READ_FENCE() _ReadBarrier()
    #endif

    CKG_API void ckg_stack_trace_dump();
#endif

#include <math.h>

float gm_range_map(float x, float s_min, float s_max, float e_min, float e_max);
float gm_move_toward(float current, float target, float delta);

typedef struct CKIT_Vector2 {
    union {
        struct {
            double x;
            double y;
        };
        double v[2];
    };
} CKIT_Vector2;

typedef struct CKIT_Vector3 {
    union {
        struct {
            double x;
            double y;
            double z;
        };
        struct {
            double r;
            double g;
            double b;
        };
        double v[3];
    };
} CKIT_Vector3;

typedef struct CKIT_Vector4 {
    union {
        struct {
            double x;
            double y;
            double z;
            double w;
        };
        struct {
            double r;
            double g;
            double b;
            double a;
        };
        double v[4];
    };
} CKIT_Vector4;


typedef struct GM_Rectangle2D {
    CKIT_Vector2 position;
    u32 width;
    u32 height;
} CKIT_Rectangle2D;

typedef struct GM_Rectangle3D {
    CKIT_Vector3 position;
    u32 length;
    u32 width;
    u32 height;
} CKIT_Rectangle3D;

typedef struct GM_Circle2D {
    CKIT_Vector2 position;
    u32 radius;
} CKIT_Circle2D;

typedef struct GM_Circle3D {
    CKIT_Vector3 position;
    u32 radius;
} CKIT_Circle3D;

GM_API double gm_lerp(double a, double b, double t);
GM_API GM_Vector2 gm_vector2_lerp(GM_Vector2 a, GM_Vector2 b, double t);
GM_API GM_Vector3 gm_vector3_lerp(GM_Vector3 a, GM_Vector3 b, double t);
GM_API GM_Vector4 gm_vector4_lerp(GM_Vector4 a, GM_Vector4 b, double t);
GM_API GM_Vector2 gm_vector2_spline_point(GM_Vector2* spline_points, u32 spline_points_count, double t);


// row major
typedef struct GM_Vec2 {
    union {
        struct {
            double x;
            double y;
        };
        double v[2];
    };
} GM_Vec2;

typedef struct GM_Vec3 {
    union {
        struct {
            double x;
            double y;
            double z;
        };
        struct {
            double r;
            double g;
            double b;
        };
        double v[3];
    };
} GM_Vec3;

typedef struct GM_Vec4 {
    union {
        struct {
            double x;
            double y;
            double z;
            double w;
        };
        struct {
            double r;
            double g;
            double b;
            double a;
        };
        double v[4];
    };
} GM_Vec4;

GM_API CKIT_Vector2 ckit_vector2_lerp(CKIT_Vector2 a, CKIT_Vector2 b, double t);
GM_API CKIT_Vector3 ckit_vector3_lerp(CKIT_Vector3 a, CKIT_Vector3 b, double t);
GM_API CKIT_Vector4 ckit_vector4_lerp(CKIT_Vector4 a, CKIT_Vector4 b, double t);

GM_API CKIT_Vector2 ckit_vector2_spline_point(CKIT_Vector2* spline_points, u32 spline_points_count, double t);


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
