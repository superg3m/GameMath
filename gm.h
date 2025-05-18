#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_COLOR
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_COLLISION
#define GM_INCLUDE_MATRIX
#define GM_INCLUDE_SHAPES
#define GM_INCLUDE_EASE_FUNCTIONS

#if defined(GM_INCLUDE_TYPES)
    #undef NULLPTR
    #undef PI
    #undef stringify
    #undef glue
    #undef KiloBytes
    #undef MegaBytes
    #undef GigaBytes
    #undef MIN
    #undef MAX
    #undef CLAMP
    #undef local_persist
    #undef internal
    #undef OFFSET_OF
    #undef FIRST_DIGIT
    #undef GET_BIT
    #undef SET_BIT
    #undef UNSET_BIT
    #undef ArrayCount
    #undef PLATFORM_MAX_PATH
    #undef PLATFORM_WINDOWS
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

    #define PLATFORM_MAX_PATH 256

    #if defined(_WIN32)
        #define NOMINMAX
        #define WIN32_LEAN_AND_MEAN
        #include <windows.h>
        #define PLATFORM_WINDOWS
        #define OS_DELIMITER '\\'
        #define CRASH __debugbreak()
    #elif defined(__APPLE__)
        #include <dlfcn.h>
        #define PLATFORM_APPLE
        #define OS_DELIMITER '/'
        #define CRASH __builtin_trap()
    #elif defined(__linux__) || defined(__unix__) || defined(__POSIX__)
        #include <dlfcn.h>
        #define PLATFORM_LINUX
        #define OS_DELIMITER '/'
        #define CRASH __builtin_trap()
    #else
        #error "Unknown Platform???"
    #endif

    #if defined(__clang__)
        #define UNUSED_FUNCTION __attribute__((used))
    #elif defined(__GNUC__) || defined(__GNUG__)
        #define UNUSED_FUNCTION __attribute__((used))
    #elif defined(_MSC_VER)
        #define UNUSED_FUNCTION
    #endif
#endif

#include <math.h>

#if defined(GM_INCLUDE_COLOR)
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

    GM_API u32 ckit_color_to_u32(GM_RGBA color);
    GM_API GM_RGBA ckit_color_from_u32(u32 color);
    GM_API GM_RGBA ckit_color_alpha_blend(GM_RGBA front_color, GM_RGBA back_color);
    GM_API GM_RGBA ckit_color_u32_alpha_blend(u32 front_color_u32, u32 back_color_u32);
    /**
	 * @brief value from 0.0 to 1.0
	 * 
	 * @param color 
	 * @param value 
	 * @return CKIT_Color 
	 */
	GM_API GM_RGBA ckit_color_multiply(GM_RGBA color, float value);

    #define GM_COLOR_BLACK ((GM_RGBA){0, 0, 0, 255})
    #define GM_COLOR_RED ((GM_RGBA){255, 0, 0, 255})
    #define GM_COLOR_BLUE ((GM_RGBA){0, 0, 255, 255})
    #define GM_COLOR_GREEN ((GM_RGBA){0, 255, 0, 255})
    #define GM_COLOR_WHITE ((GM_RGBA){255, 255, 255, 255})
    #define GM_COLOR_PINK ((GM_RGBA){255, 105, 180, 255})
    #define GM_COLOR_LIME ((GM_RGBA){0, 255, 128, 255})
    #define GM_COLOR_CYAN ((GM_RGBA){0, 255, 255, 255})
    #define GM_COLOR_PURPLE ((GM_RGBA){128, 0, 128, 255})
    #define GM_COLOR_YELLOW ((GM_RGBA){255, 255, 0, 255})
#endif

#if defined(GM_INCLUDE_COLLISION)
    typedef struct GM_AABB {
        GM_Vec3 min;
        GM_Vec3 max;
    } GM_AABB;

    // SAT collision

    GM_API bool gm_aabb_point_colliding(GM_Vec3 point, GM_AABB aabb);
    GM_API bool gm_aabb_aabb_colliding(GM_AABB a, GM_AABB b);
#endif

// https://www.youtube.com/watch?v=QS3677PlIos
// https://www.youtube.com/watch?v=YJB1QnEmlTs
#if defined(GM_INCLUDE_INTERPOLATION)
    GM_API float gm_lerp(float a, float b, float t);
    GM_API Quaternion gm_slerp(Quaternion a, Quaternion b, float t);
    GM_API Vector3 gm_barycentric(Vector3 a, Vector3 b, Vector3 c, float u, float v);
    GM_API float gm_inverse_lerp(float a, float b, float value);
    GM_API float gm_remap(float x, float s_min, float s_max, float e_min, float e_max);
    GM_API float gm_move_toward(float current, float target, float delta);
    GM_API float gm_smoothstep(float edge0, float edge1, float x);
    GM_API float gm_smootherstep(float edge0, float edge1, float x);
#endif

#if defined(GM_INCLUDE_VECTOR)
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


    GM_API GM_Vec2 gm_vector2_lerp(GM_Vec2 a, GM_Vec2 b, float t);
    GM_API GM_Vec3 gm_vector3_lerp(GM_Vec3 a, GM_Vec3 b, float t);
    GM_API GM_Vec4 gm_vector4_lerp(GM_Vec4 a, GM_Vec4 b, float t);
    GM_API GM_Vec2 gm_vector2_spline_point(GM_Vec2* spline_points, u32 spline_points_count, float t);
    GM_API float gm_v3_dot_product(GM_Vec3 A, GM_Vec3 B);
    GM_API float gm_v4_dot_product(GM_Vec4 A, GM_Vec4 B);
#endif

#if defined(CKG_INCLUDE_SHAPES)
    typedef struct GM_Rectangle2D {
        GM_Vec2 position;
        u32 width;
        u32 height;
    } CKIT_Rectangle2D;

    typedef struct GM_Rectangle3D {
        GM_Vec3 position;
        u32 length;
        u32 width;
        u32 height;
    } CKIT_Rectangle3D;

    typedef struct GM_Circle2D {
        GM_Vec2 position;
        u32 radius;
    } CKIT_Circle2D;

    typedef struct GM_Circle3D {
        GM_Vec3 position;
        u32 radius;
    } CKIT_Circle3D;
#endif

#if defined(CKG_INCLUDE_MATRIX)
    typedef struct GM_Matix4 {
        union {
            float data[16]; 
            GM_Vec4 v[4];
        };
    } GM_Matix4;

    GM_API GM_Matix4 gm_mat4_identity();
    GM_API GM_Matix4 gm_mat4_translation(GM_Vec3 t);
    GM_API GM_Matix4 gm_mat4_scale(GM_Vec3 s);
    GM_API GM_Matix4 gm_mat4_scale_xyz(float x, float y, float z);

    GM_API GM_Matix4 gm_mat4_rotation_x(float degrees);
    GM_API GM_Matix4 gm_mat4_rotation_y(float degress);
    GM_API GM_Matix4 gm_mat4_rotation_z(float degress);

    GM_API GM_Matix4 gm_mat4_perspective(float fov_degrees, float aspect, float near, float far);
    GM_API GM_Matix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near, float far);
    GM_API GM_Matix4 gm_mat4_look_at(GM_Vec3 eye, GM_Vec3 center, GM_Vec3 up);

    GM_API GM_Matix4 gm_mat4_mult(GM_Matix4 A, GM_Matix4 B);
    GM_API GM_Matix4 gm_mat4_inverse(GM_Matix4 m);
    GM_API GM_Matix4 gm_mat4_transpose(GM_Matix4 m);
#endif



#if defined(CKG_INCLUDE_EASE_FUNCTIONS)
    // Date: May 18, 2025
    // NOTE(Jovanni): Visualize these at: https://easings.net/
    GM_API float gm_ease_in_sine(float t);
    GM_API float gm_ease_out_sine(float t);
    GM_API float gm_ease_in_out_sine(float t);
    GM_API float gm_ease_in_quad(float t);
    GM_API float gm_ease_out_quad(float t);
    GM_API float gm_ease_in_cubic(float t);
    GM_API float gm_ease_out_cubic(float t);
    GM_API float gm_ease_in_out_cubic(float t);
    GM_API float gm_ease_in_quart(float t);
    GM_API float gm_ease_out_quart(float t);
    GM_API float gm_ease_in_out_quart(float t);
    GM_API float gm_ease_in_quint(float t);
    GM_API float gm_ease_out_quint(float t);
    GM_API float gm_ease_in_out_quint(float t);
    GM_API float gm_ease_in_expo(float t);
    GM_API float gm_ease_out_expo(float t);
    GM_API float gm_ease_in_out_expo(float t);
    GM_API float gm_ease_in_circ(float t);
    GM_API float gm_ease_out_circ(float t);
    GM_API float gm_ease_in_out_circ(float t);
    GM_API float gm_ease_in_back(float t);
    GM_API float gm_ease_out_back(float t);
    GM_API float gm_ease_in_out_back(float t);
    GM_API float gm_ease_in_elastic(float t);
    GM_API float gm_ease_out_elastic(float t);
    GM_API float gm_ease_in_out_elastic(float t);
    GM_API float gm_ease_in_bounce(float t);
    GM_API float gm_ease_out_bounce(float t);
    GM_API float gm_ease_in_out_bounce(float t);
#endif