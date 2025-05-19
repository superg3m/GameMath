#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_COLOR
#define GM_INCLUDE_INTERPOLATION
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_MATRIX
#define CKG_INCLUDE_QUATERNION
#define GM_INCLUDE_COLLISION
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

    GM_API GM_Vec2 gm_v2_create(float x, float y);
    GM_API GM_Vec3 gm_v3_create(float x, float y, float z);
    GM_API GM_Vec4 gm_v4_create(float x, float y, float z, float w);

    GM_API GM_Vec2 gm_v2_scale(GM_Vec2 v, float scale);
    GM_API GM_Vec3 gm_v3_scale(GM_Vec3 v, float scale);
    GM_API GM_Vec4 gm_v4_scale(GM_Vec4 v, float scale);

    GM_API GM_Vec2 gm_v2_lerp(GM_Vec2 A, GM_Vec2 B, float t);
    GM_API GM_Vec3 gm_v3_lerp(GM_Vec3 A, GM_Vec3 B, float t);
    GM_API GM_Vec4 gm_v4_lerp(GM_Vec4 A, GM_Vec4 B, float t);

    GM_API GM_Vec2 gm_v2_spline_point(GM_Vec2* spline_points, u32 spline_points_count, float t);

    GM_API GM_Vec2 gm_v2_projection(GM_Vec2 A, GM_Vec2 B);
    GM_API GM_Vec3 gm_v3_projection(GM_Vec3 A, GM_Vec3 B);

    GM_API GM_Vec2 gm_v2_normalize(GM_Vec2 A);
    GM_API GM_Vec3 gm_v3_normalize(GM_Vec3 A);
    GM_API GM_Vec4 gm_v4_normalize(GM_Vec4 A);

    GM_API float gm_v2_magnitude(GM_Vec2 A);
    GM_API float gm_v3_magnitude(GM_Vec3 A);
    GM_API float gm_v4_magnitude(GM_Vec4 A);

    GM_API float gm_v2_dot(GM_Vec2 A, GM_Vec2 B);
    GM_API float gm_v3_dot(GM_Vec3 A, GM_Vec3 B);
    GM_API float gm_v4_dot(GM_Vec4 A, GM_Vec4 B);

    GM_API float gm_v3_cross(GM_Vec4 A, GM_Vec4 B);
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
    GM_API GM_Matix4 gm_mat4_inverse(GM_Matix4 m, bool* success);
    GM_API GM_Matix4 gm_mat4_transpose(GM_Matix4 m);
#endif

#if defined(CKG_INCLUDE_QUATERNION)
    typedef struct GM_Quaternion {
        float w;
        GM_Vec3 v;
    } GM_Quaternion;

    GM_Quaternion gm_quat_create(GM_Vec3 axis, float theta);
    GM_Quaternion gm_quat_inverse(GM_Quaternion quat);
    GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2);
    GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec);
#endif

// https://www.youtube.com/watch?v=QS3677PlIos
// https://www.youtube.com/watch?v=YJB1QnEmlTs
#if defined(GM_INCLUDE_INTERPOLATION)
    GM_API float gm_lerp(float a, float b, float t);
    GM_API float gm_inverse_lerp(float a, float b, float value);
    GM_API GM_Vec3 gm_barycentric(GM_Vec3 a, GM_Vec3 b, GM_Vec3 c, float u, float v);

    GM_API float gm_remap(float x, float s_min, float s_max, float e_min, float e_max);
    GM_API float gm_move_toward(float current, float target, float delta);

    GM_API float gm_smoothstep(float edge0, float edge1, float x);
    GM_API float gm_smootherstep(float edge0, float edge1, float x);
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

    GM_Rectangle2D gm_rectangle2d_create(s32 x, s32 y, u32 width, u32 height) {
        GM_Rectangle2D ret = {0};
        ret.position.x = x;
        ret.position.y = y;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    GM_Rectangle3D gm_rectangle3d_create(s32 x, s32 y, s32 z, u32 length, u32 width, u32 height) {
        GM_Rectangle3D ret = {0};
        ret.position.x = x;
        ret.position.y = y;
        ret.position.z = z;
        ret.length = length;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    bool gm_rectangle_check_aabb_collision(GM_Rectangle2D rect1, GM_Rectangle2D rect2) {
        if (rect1.position.x < rect2.position.x + rect2.width && rect1.position.x + rect1.width > rect2.position.x &&
            rect1.position.y < rect2.position.y + rect2.height && rect1.position.y + rect1.height > rect2.position.y) {
            return true;
        }

        return false;
    }

    CKIT_Rectangle2D gm_rectangle_get_aabb_collision(GM_Rectangle2D rect1, GM_Rectangle2D rect2) {
        GM_Rectangle2D ret = {0};
        (void)rect1;
        (void)rect2;
        
        return ret;
    }

    GM_Circle2D gm_circle2d_create(s32 x, s32 y, u32 radius) {
        GM_Circle2D ret = {0};
        ret.position.x = x;
        ret.position.y = y;
        ret.radius = radius;

        return ret;
    }

    GM_Circle3D gm_circle3d_create(s32 x, s32 y, s32 z, u32 radius) {
        GM_Circle3D ret = {0};
        ret.position.x = x;
        ret.position.y = y;
        ret.position.z = z;
        ret.radius = radius;

        return ret;
    }
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

//
// ===================================================== CKIT_IMPL =====================================================
//

#if defined(CKG_IMPL_VECTOR) 
    float gm_v2_dot(GM_Vec2 A, GM_Vec2 B) {
        return (A.x * B.x) + (A.y * B.y);
    }

    float gm_v3_dot(GM_Vec3 A, GM_Vec3 B) {
        return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
    }

    float gm_v4_dot(GM_Vec4 A, GM_Vec4 B) {
        return (A.x * B.x) + (A.y * B.y) + (A.z * B.z) + (A.w * B.w);
    }

    float gm_v3_cross(GM_Vec3 A, GM_Vec3 B) {
        GM_Vec3 ret;

        ret.x = A.y * B.z - A.z * B.y;
        ret.y = A.z * B.x - A.x * B.z;
        ret.z = A.x * B.y - A.y * B.x;

        return ret;
    }

    float gm_v2_magnitude(GM_Vec2 A) {
        return sqrt((A.x*A.x) + (A.x*A.x));
    }

    float gm_v3_magnitude(GM_Vec3 A) {
        return sqrt((A.x*A.x) + (A.x*A.x) + (A.z*A.z));
    }

    float gm_v4_magnitude(GM_Vec4 A) {
        return sqrt((A.x*A.x) + (A.x*A.x) + (A.z*A.z) + (A.w*A.w));
    }

    GM_Vec2 gm_v2_normalize(GM_Vec2 A) {
        GM_Vec2 ret;
        const float magnitude = gm_v4_magnitude(A);
        ret.x = A.x / magnitude;
        ret.y = A.y / magnitude;

        return ret;
    }

    GM_Vec3 gm_v3_normalize(GM_Vec3 A) {
        GM_Vec3 ret;
        const float magnitude = gm_v4_magnitude(A);
        ret.x = A.x / magnitude;
        ret.y = A.y / magnitude;
        ret.z = A.z / magnitude;

        return ret;
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

    GM_Vec2 gm_v2_create(float x, float y) {
        GM_Vec2 ret = {0};
        ret.x = x;
        ret.y = y;

        return ret;
    }

    GM_Vec3 gm_v3_create(float x, float y, float z) {
        GM_Vec3 ret = {0};
        ret.x = x;
        ret.y = y;
        ret.z = z; 

        return ret;
    }

    GM_Vec4 gm_v4_create(float x, float y, float z, float w) {
        GM_Vec4 ret = {0};
        ret.x = x;
        ret.y = y;
        ret.z = z; 
        ret.w = w; 

        return ret;
    }

    GM_Vec2 gm_v2_scale(GM_Vec2 v, float scale); {
        GM_Vec2 ret = {0};
        ret.x = v.x * scale;
        ret.y = v.y * scale;

        return ret;
    }

    GM_Vec3 gm_v3_scale(GM_Vec3 v, float scale) {
        GM_Vec3 ret = {0};
        ret.x = v.x * scale;
        ret.y = v.y * scale;
        ret.z = v.z * scale; 

        return ret;
    }

    GM_Vec4 gm_v4_scale(GM_Vec4 v, float scale) {
        GM_Vec4 ret = {0};
        ret.x = v.x * scale;
        ret.y = v.y * scale;
        ret.z = v.z * scale; 
        ret.w = v.w * scale; 

        return ret;
    }
#endif

#if defined(CKG_IMPL_INTERPOLATION) 
    float gm_lerp(float a, float b, float t) {
        return a + ((b - a) * t);
    }

    float gm_inverse_lerp(float a, float b, float value) {
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

    GM_Quaternion gm_slerp(GM_Quaternion a, GM_Quaternion b, float t);
    GM_Vec3 gm_barycentric(GM_Vec3 a, GM_Vec3 b, GM_Vec3 c, float u, float v);
    float gm_smoothstep(float edge0, float edge1, float x);
    float gm_smootherstep(float edge0, float edge1, float x);
#endif

#if defined(CKG_IMPL_EASE_FUNCTION)
    float gm_ease_in_sine(float t) {
        return 1 - cos((t * PI) / 2);
    }

    float gm_ease_out_sine(float t) {
            return sin((t * PI) / 2);
    }

    float gm_ease_in_out_sine(float t) {
        return -(cos(PI * t) - 1) / 2;
    }

    float gm_ease_in_quad(float t) {
        return t * t;
    }

    float gm_ease_out_quad(float t) {
        return 1 - (1 - t) * (1 - t);
    }

    float gm_ease_in_cubic(float t) {
        return t * t * t;
    }

    float gm_ease_out_cubic(float t) {
        return 1 - pow(1 - t, 3);
    }

    float gm_ease_in_out_cubic(float t) {
        return t < 0.5 ? 4 * t * t * t : 1 - pow(-2 * t + 2, 3) / 2;
    }

    float gm_ease_in_quart(float t) {
        return t * t * t * t;
    }

    float gm_ease_out_quart(float t) {
        return 1 - pow(1 - t, 4);
    }

    float gm_ease_in_out_quart(float t) {
        return t < 0.5 ? 8 * t * t * t * t : 1 - pow(-2 * t + 2, 4) / 2;
    }

    float gm_ease_in_quint(float t) {
        return t * t * t * t * t;
    }

    float gm_ease_out_quint(float t) {
        return 1 + pow(t - 1, 5);
    }

    float gm_ease_in_out_quint(float t) {
        return t < 0.5 ? 16 * t * t * t * t * t : 1 - pow(-2 * t + 2, 5) / 2;
    }

    float gm_ease_in_expo(float t) {
        return t == 0 ? 0 : pow(2, 10 * t - 10);
    }

    float gm_ease_out_expo(float t) {
        return t == 1 ? 1 : 1 - pow(2, -10 * t);
    }

    float gm_ease_in_out_expo(float t) {
        if (t == 0) return 0;
        if (t == 1) return 1;
        if (t < 0.5) return pow(2, 20 * t - 10) / 2;
        return (2 - pow(2, -20 * t + 10)) / 2;
    }

    float gm_ease_in_circ(float t) {
        return 1 - sqrt(1 - (t * t));
    }

    float gm_ease_out_circ(float t) {
        return sqrt(1 - pow(t - 1, 2));
    }

    float gm_ease_in_out_circ(float t) {
        if (t < 0.5) return (1 - sqrt(1 - pow(2 * t, 2))) / 2;
        return (sqrt(1 - pow(-2 * t + 2, 2)) + 1) / 2;
    }

    float gm_ease_in_back(float t) {
        float c1 = 1.70158;
        float c3 = c1 + 1;
        return c3 * t * t * t - c1 * t * t;
    }

    float gm_ease_out_back(float t) {
        float c1 = 1.70158;
        float c3 = c1 + 1;
        return 1 + c3 * pow(t - 1, 3) + c1 * pow(t - 1, 2);
    }

    float gm_ease_in_out_back(float t) {
        float c1 = 1.70158;
        float c2 = c1 * 1.525;
        if (t < 0.5) return (pow(2 * t, 2) * ((c2 + 1) * 2 * t - c2)) / 2;
        return (pow(2 * t - 2, 2) * ((c2 + 1) * (t * 2 - 2) + c2) + 2) / 2;
    }

    float gm_ease_in_elastic(float t) {
        float c4 = (2 * PI) / 3;
        if (t == 0) return 0;
        if (t == 1) return 1;
        return -pow(2, 10 * t - 10) * sin((t * 10 - 10.75) * c4);
    }

    float gm_ease_out_elastic(float t) {
        float c4 = (2 * PI) / 3;
        if (t == 0) return 0;
        if (t == 1) return 1;
        return pow(2, -10 * t) * sin((t * 10 - 0.75) * c4) + 1;
    }

    float gm_ease_in_out_elastic(float t) {
        float c5 = (2 * PI) / 4.5;
        if (t == 0) return 0;
        if (t == 1) return 1;
        if (t < 0.5) return -(pow(2, 20 * t - 10) * sin((20 * t - 11.125) * c5)) / 2;
        return (pow(2, -20 * t + 10) * sin((20 * t - 11.125) * c5)) / 2 + 1;
    }

    float gm_ease_in_bounce(float t) {
        return 1 - gm_ease_out_bounce(1 - t);
    }

    float gm_ease_out_bounce(float t) {
        float n1 = 7.5625;
        float d1 = 2.75;
        if (t < 1 / d1) {
            return n1 * t * t;
        } else if (t < 2 / d1) {
            return n1 * (t -= 1.5 / d1) * t + 0.75;
        } else if (t < 2.5 / d1) {
            return n1 * (t -= 2.25 / d1) * t + 0.9375;
        } else {
            return n1 * (t -= 2.625 / d1) * t + 0.984375;
        }
    }

    float gm_ease_in_out_bounce(float t) {
        return t < 0.5
            ? (1 - gm_ease_out_bounce(1 - 2 * t)) / 2
            : (1 + gm_ease_out_bounce(2 * t - 1)) / 2;
    }
#endif

#if defined(CKG_IMPL_QUATERNION)
    GM_Quaternion gm_quat_create(GM_Vec3 axis, float theta) {
        GM_Quaternion ret = {0};

        float radians = DEGREES_TO_RAD(theta);
        ret.w = cos(radians / 2);
        ret.v = gm_v3_scale(axis, sin(radians / 2));

        return ret;
    }

    GM_Quaternion gm_quat_inverse(GM_Quaternion quat) {
        GM_Quaternion ret = {0};

        ret.w = quat.w;
        ret.v = gm_v3_scale(quat.v, -1);

        return ret;
    }

    GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2) {
        GM_Quaternion ret;

        ret.w = (q1.w * q2.w) + gm_v3_dot(q1.v, q2.v);
        ret.v = gm_v3_scale(q1.v, q2.w) + gm_v3_scale(q2.v, q1.w) + gm_v3_cross(q1.v, q2.v);

        return ret;
    }

    // Rotate a vector with this quaternion.
    GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec) {
        GM_Quaternion p;
        p.w = 0;
        p.v = vec;

        return gm_quat_mult(gm_quat_mult(quat, p),  gm_quat_inverse(quat)).v;
    }
#endif