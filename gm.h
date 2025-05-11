#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_MATRIX
#define GM_INCLUDE_UTILITY
#define GM_INCLUDE_SMOOTHING_FUCNTIONS

#if defined(GM_INCLUDE_TYPES)
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
#endif

#include <math.h>

typedef struct CKIT_Vector2 {
    union {
        struct {
            float x;
            float y;
        };
        float v[2];
    };
} CKIT_Vector2;

typedef struct CKIT_Vector3 {
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
        float v[3];
    };
} CKIT_Vector3;

typedef struct CKIT_Vector4 {
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
        float v[4];
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


// row major
typedef struct GM_Vec2 {
    union {
        struct {
            float x;
            float y;
        };
        float v[2];
    };
} GM_Vec2;

typedef struct GM_Vec3 {
    union {
        struct {
            float x;
            float y;
            float z;
        };
        float v[3];
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
        float v[4];
    };
} GM_Vec4;


typedef struct GM_RGB {
    union {
        union {
            struct {
                u8 r;
                u8 g;
                u8 b;
            };
            u8 c[3];
        } rgb;
    };
} GM_RGB;

typedef struct GM_RGBA {
    union {
        union {
            struct {
                u8 a;
                u8 r;
                u8 g;
                u8 b;
            };
            u32 hex;
            u8 c[4];
        } argb;

        union {
            struct {
                u8 r;
                u8 g;
                u8 b;
                u8 a;
            };
            u32 hex;
            u8 c[4];
        } rgba;
    };
} GM_RGBA;

typedef struct GM_AABB {
    GM_Vec3 min;
    GM_Vec3 max;
} GM_AABB;

GM_API bool gm_aabb_point_colliding(GM_Vec3 point, GM_AABB aabb);
GM_API bool gm_aabb_aabb_colliding(GM_AABB a, GM_AABB b);

GM_API float gm_lerp(float a, float b, float t);
GM_API GM_Vec2 gm_vector2_lerp(GM_Vec2 a, GM_Vec2 b, float t);
GM_API GM_Vec3 gm_vector3_lerp(GM_Vec3 a, GM_Vec3 b, float t);
GM_API GM_Vec4 gm_vector4_lerp(GM_Vec4 a, GM_Vec4 b, float t);
GM_API GM_Vec2 gm_vector2_spline_point(GM_Vec2* spline_points, u32 spline_points_count, float t);


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

float gm_range_mapper(float x, float s_min, float s_max, float e_min, float e_max);
float gm_move_toward(float current, float target, float delta);