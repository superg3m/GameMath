#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#if defined(GM_IMPL)
    #define GM_IMPL_COLOR
    #define GM_IMPL_VECTOR
    #define GM_IMPL_MATRIX
    #define GM_IMPL_QUATERNION
    #define GM_IMPL_UTILITY
    #define GM_IMPL_EASE_FUNCTIONS
    #define GM_IMPL_COLLISION
    #define GM_IMPL_SHAPES
    #define GM_IMPL_PHYSICS
#endif

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_COLOR
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_MATRIX
#define GM_INCLUDE_QUATERNION
#define GM_INCLUDE_UTILITY
#define GM_INCLUDE_EASE_FUNCTIONS
#define GM_INCLUDE_COLLISION
#define GM_INCLUDE_SHAPES
#define GM_INCLUDE_PHYSICS

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
    #undef SQUARE
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
    #define SQUARE(a) ((a) * (a))

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

    // NOTE: Union for hex in RGBA is potentially problematic due to endianness.
    // It's safer to use explicit conversion functions (like gm_rgba_to_u32).
    typedef struct GM_RGBA {
        union {
            struct {
                u8 r;
                u8 g;
                u8 b;
                u8 a;
            };
            u8 c[4];
            u32 hex_rgba; // Explicitly name for rgba order for hex
        };
    } GM_RGBA;

    GM_API u32 gm_rgba_to_u32(GM_RGBA color);
    GM_API GM_RGBA gm_rgba_from_u32(u32 color);
    GM_API GM_RGBA gm_rgba_alpha_blend(GM_RGBA front_color, GM_RGBA back_color);
    GM_API GM_RGBA gm_rgba_u32_alpha_blend(u32 front_color_u32, u32 back_color_u32);
    /**
     * @brief value from 0.0 to 1.0
     *
     * @param color
     * @param value
     * @return GM_RGBA
     */
    GM_API GM_RGBA gm_rgba_multiply(GM_RGBA color, float value);

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

#if defined(GM_INCLUDE_VECTOR)
    // Vectors are typically treated as column vectors for post-multiplication (M * v)
    // Their memory layout here is simple sequential (x, y, z, w)
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

    GM_API GM_Vec2 gm_vec2_create(float x, float y);
    GM_API GM_Vec3 gm_vec3_create(float x, float y, float z);
    GM_API GM_Vec4 gm_vec4_create(float x, float y, float z, float w);

    GM_API GM_Vec2 gm_vec2_add(GM_Vec2 A, GM_Vec2 B);
    GM_API GM_Vec3 gm_vec3_add(GM_Vec3 A, GM_Vec3 B);
    GM_API GM_Vec4 gm_vec4_add(GM_Vec4 A, GM_Vec4 B);

    GM_API GM_Vec2 gm_vec2_sub(GM_Vec2 A, GM_Vec2 B);
    GM_API GM_Vec3 gm_vec3_sub(GM_Vec3 A, GM_Vec3 B);
    GM_API GM_Vec4 gm_vec4_sub(GM_Vec4 A, GM_Vec4 B);

    GM_API GM_Vec2 gm_vec2_scale(GM_Vec2 v, float scale);
    GM_API GM_Vec3 gm_vec3_scale(GM_Vec3 v, float scale);
    GM_API GM_Vec4 gm_vec4_scale(GM_Vec4 v, float scale);

    GM_API GM_Vec2 gm_vec2_lerp(GM_Vec2 A, GM_Vec2 B, float t);
    GM_API GM_Vec3 gm_vec3_lerp(GM_Vec3 A, GM_Vec3 B, float t);
    GM_API GM_Vec4 gm_vec4_lerp(GM_Vec4 A, GM_Vec4 B, float t);

    GM_API GM_Vec2 gm_vec2_projection(GM_Vec2 A, GM_Vec2 B);
    GM_API GM_Vec3 gm_vec3_projection(GM_Vec3 A, GM_Vec3 B);

    GM_API GM_Vec2 gm_vec2_normalize(GM_Vec2 A);
    GM_API GM_Vec3 gm_vec3_normalize(GM_Vec3 A);
    GM_API GM_Vec4 gm_vec4_normalize(GM_Vec4 A);

    GM_API float gm_vec2_distance(GM_Vec2 A, GM_Vec2 B);
    GM_API float gm_vec3_distance(GM_Vec3 A, GM_Vec3 B);
    GM_API float gm_vec4_distance(GM_Vec4 A, GM_Vec4 B);

    GM_API float gm_vec2_magnitude(GM_Vec2 A);
    GM_API float gm_vec3_magnitude(GM_Vec3 A);
    GM_API float gm_vec4_magnitude(GM_Vec4 A);

    GM_API float gm_vec2_dot(GM_Vec2 A, GM_Vec2 B);
    GM_API float gm_vec3_dot(GM_Vec3 A, GM_Vec3 B);
    GM_API float gm_vec4_dot(GM_Vec4 A, GM_Vec4 B);

    GM_API GM_Vec3 gm_vec3_cross(GM_Vec3 A, GM_Vec3 B);
    GM_API GM_Vec2 gm_vec2_spline_point(GM_Vec2* spline_points, u32 spline_points_count, float t); // Declaration
#endif

#if defined(GM_INCLUDE_MATRIX)
    // Matrices are COLUMN-MAJOR (OpenGL convention)
    // data[0-3] = Col 0 (X-axis)
    // data[4-7] = Col 1 (Y-axis)
    // data[8-11] = Col 2 (Z-axis)
    // data[12-15] = Col 3 (Translation)
    typedef struct GM_Matrix4 {
        union {
            float data[16];
            GM_Vec4 v[4];
        };
    } GM_Matrix4;

    GM_API GM_Matrix4 gm_mat4_identity();
    GM_API GM_Matrix4 gm_mat4_translate(GM_Matrix4 mat, GM_Vec3 t);
    GM_API GM_Matrix4 gm_mat4_translate_xyz(GM_Matrix4 mat, float x, float y, float z);
    GM_API GM_Matrix4 gm_mat4_scale(GM_Matrix4 mat, GM_Vec3 s);
    GM_API GM_Matrix4 gm_mat4_scale_xyz(GM_Matrix4 mat, float x, float y, float z);

    GM_API GM_Matrix4 gm_mat4_rotate(GM_Matrix4 mat, float degrees, GM_Vec3 axis);
    GM_API GM_Matrix4 gm_mat4_rotate_xyz(GM_Matrix4 mat, float degrees, float x, float y, float z);

    GM_API GM_Matrix4 gm_mat4_perspective(float fov_degrees, float aspect, float near_plane, float far_plane);
    GM_API GM_Matrix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane);
    GM_API GM_Matrix4 gm_mat4_look_at(GM_Vec3 camera_position, GM_Vec3 target_position, GM_Vec3 world_up);

    GM_API GM_Matrix4 gm_mat4_mult(GM_Matrix4 A, GM_Matrix4 B);
    GM_API GM_Matrix4 gm_mat4_inverse(GM_Matrix4 m, bool* success);
    GM_API GM_Matrix4 gm_mat4_transpose(GM_Matrix4 m);
#endif

#if defined(GM_INCLUDE_QUATERNION)
    typedef struct GM_Quaternion {
        float w;
        GM_Vec3 v;
    } GM_Quaternion;

    GM_API GM_Quaternion gm_quat_create(GM_Vec3 axis, float theta);
    GM_API GM_Quaternion gm_quat_inverse(GM_Quaternion quat);
    GM_API GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2);
    GM_API GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec);
    GM_API GM_Quaternion gm_slerp(GM_Quaternion a, GM_Quaternion b, float t);
#endif

#if defined(GM_INCLUDE_UTILITY)
    GM_API float gm_lerp(float a, float b, float t);
    GM_API float gm_inverse_lerp(float a, float b, float value);
    GM_API GM_Vec3 gm_barycentric(GM_Vec3 a, GM_Vec3 b, GM_Vec3 c, float u, float v);

    GM_API float gm_remap(float x, float s_min, float s_max, float e_min, float e_max);
    GM_API float gm_move_toward(float current, float target, float delta);

    GM_API float gm_smoothstep(float edge0, float edge1, float x);
    GM_API float gm_smootherstep(float edge0, float edge1, float x);

    GM_API bool gm_near_equal(float a, float b);
#endif

#if defined(GM_INCLUDE_EASE_FUNCTIONS)
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

#if defined(GM_INCLUDE_SHAPES)
    typedef struct GM_Rectangle2D {
        GM_Vec2 position;
        u32 width;
        u32 height;
    } GM_Rectangle2D;

    typedef struct GM_Rectangle3D {
        GM_Vec3 position;
        u32 length;
        u32 width;
        u32 height;
    } GM_Rectangle3D;

    typedef struct GM_Circle2D {
        GM_Vec2 position;
        u32 radius;
    } GM_Circle2D;

    typedef struct GM_Circle3D {
        GM_Vec3 position;
        u32 radius;
    } GM_Circle3D;

    GM_API GM_Rectangle2D gm_rectangle2d_create(float x, float y, u32 width, u32 height);
    GM_API GM_Rectangle3D gm_rectangle3d_create(float x, float y, float z, u32 length, u32 width, u32 height);
    GM_API bool gm_rectangle_check_aabb_collision(GM_Rectangle2D rect1, GM_Rectangle2D rect2);
    GM_API GM_Circle2D gm_circle2d_create(float x, float y, u32 radius);
    GM_API GM_Circle3D gm_circle3d_create(float x, float y, float z, u32 radius);
#endif

#if defined(GM_INCLUDE_COLLISION)

    typedef struct GM_CollisionInfo2D {
        GM_Vec2 normal;
        float depth;
    } GM_CollisionInfo2D;

    typedef enum GM_Collider2DType { 
        GM_COLLIDER_AABB,
        GM_COLLIDER_POINT,
        GM_COLLIDER_CIRCLE,
    } GM_Collider2DType;

    typedef struct GM_Collider2D {
        GM_Collider2DType type;
        union {
            GM_Rectangle2D aabb;
            GM_Circle2D circle;
        };

        u64 collision_mask; // Checking other layers to see if should actually collide or ignore based on their layer_mask
        u64 layer_mask; // is used by other colliders to see if it should actually collide or ignore

        // collision_enter(GM_Collider2D* other, GM_Collider2D* other);
        // collision_persisting(GM_Collider2D* other, GM_Collider2D* other, float time);
        // collision_function_leave(GM_Collider2D* other, GM_Collider2D* other);
    } GM_Collider2D;

    GM_API bool gm_collision2d_circles(GM_Circle2D c1, GM_Circle2D c2, GM_CollisionInfo2D* collision_info);

    // Date: May 28, 2025
    // TODO(Jovanni): SAT collision

    GM_API GM_Collider2D gm_collider2d_circle_create(GM_Circle2D circle);
    GM_API GM_Collider2D gm_collider2d_aabb_create(GM_Rectangle2D aabb);

    // 0 - 63
    GM_API GM_Collider2D gm_collider2d_set_layer_bit(u8 layer_bit);
    GM_API GM_Collider2D gm_collider2d_unset_layer_bit(u8 layer_bit);

    // 0 - 63
    GM_API GM_Collider2D gm_collider2d_set_mask_bit(u8 mask_bit_index);
    GM_API GM_Collider2D gm_collider2d_unset_mask_bit(u8 mask_bit_index);

    // GM_API bool gm_aabb_point_colliding(GM_Vec3 point, GM_AABB aabb);
    // GM_API bool gm_aabb_aabb_colliding(GM_AABB a, GM_AABB b);
#endif

#if defined(GM_INCLUDE_PHYSICS)
    typedef struct GM_RigidBody2D {
        GM_Vec2 position;
        GM_Vec2 velocity;
        GM_Vec2 acceleration;
        GM_Vec2 force;
        float mass; // In kilograms (KG)
    } GM_RigidBody2D;

    typedef struct GM_PhysicsObject2D {
        GM_RigidBody2D rb;
        GM_Collider2D collider;
    } GM_PhysicsObject2D;

    GM_API GM_RigidBody2D gm_physics2d_rb_create(GM_Vec2 position, float mass);
    GM_API GM_RigidBody2D gm_physics2d_rb_create_xy(float x, float y, float mass);

    GM_API GM_PhysicsObject2D gm_physics2d_object_create(GM_Vec2 position, float mass, GM_Collider2D collider);
    GM_API GM_PhysicsObject2D gm_physics2d_object_create_xy(float x, float y, float mass, GM_Collider2D collider);

    GM_API void gm_physics2d_add_velocity(GM_PhysicsObject2D* object, GM_Vec2 velocity);
    GM_API void gm_physics2d_add_velocity_xy(GM_PhysicsObject2D* object, float velocity_x, float velocity_y);

    GM_API void gm_physics2d_add_force(GM_PhysicsObject2D* object, GM_Vec2 force);
    GM_API void gm_physics2d_add_force_xy(GM_PhysicsObject2D* object, float force_x, float force_y);

    GM_API void gm_physics2d_simulate(GM_PhysicsObject2D* objects, int rb_count, float dt);
#endif

//
// ===================================================== GM_IMPL =====================================================
//

#if defined(GM_IMPL_COLOR)
    GM_RGBA gm_rgba_from_u32(u32 color) {
        GM_RGBA ret = {0};
        ret.r = (u8)((color >> 24) & 0xFF);
        ret.g = (u8)((color >> 16) & 0xFF);
        ret.b = (u8)((color >> 8) & 0xFF);
        ret.a = (u8)((color >> 0) & 0xFF);

        return ret;
    }

    GM_RGBA gm_rgba_multiply(GM_RGBA color, float value) {
        GM_RGBA ret = {0};
        ret.r = (u8)CLAMP(color.r * value, 0, 255);
        ret.g = (u8)CLAMP(color.g * value, 0, 255);
        ret.b = (u8)CLAMP(color.b * value, 0, 255);
        ret.a = (u8)CLAMP(color.a * value, 0, 255);

        return ret;
    }

    u32 gm_rgba_u32_multiply(u32 color, float value) {
        u8 r = (u8)(((color >> 24) & 0xFF) * value);
        u8 g = (u8)(((color >> 16) & 0xFF) * value);
        u8 b = (u8)(((color >> 8) & 0xFF) * value);
        u8 a = (u8)(((color >> 0) & 0xFF) * value);

        u32 result = (u32)(r << 24 | g << 16 | b << 8 | a);
        return result;
    }

    GM_RGBA gm_rgba_alpha_blend(GM_RGBA front_color, GM_RGBA back_color) {
        GM_RGBA ret = {0};

        float normalized_front_alpha = (float)front_color.a / 255.0f; // Alpha of the front color
        float normalized_back_alpha = (float)back_color.a / 255.0f;  // Alpha of the back color

        ret.a = (u8)CLAMP((front_color.a + back_color.a * (1.0f - normalized_front_alpha)), 0, 255);
        ret.r = (u8)CLAMP((front_color.r * normalized_front_alpha + back_color.r * normalized_back_alpha * (1.0f - normalized_front_alpha)), 0, 255);
        ret.g = (u8)CLAMP((front_color.g * normalized_front_alpha + back_color.g * normalized_back_alpha * (1.0f - normalized_front_alpha)), 0, 255);
        ret.b = (u8)CLAMP((front_color.b * normalized_front_alpha + back_color.b * normalized_back_alpha * (1.0f - normalized_front_alpha)), 0, 255);

        return ret;
    }

    GM_RGBA gm_rgba_u32_alpha_blend(u32 front_color_u32, u32 back_color_u32) {
        GM_RGBA front_color = gm_rgba_from_u32(front_color_u32);
        GM_RGBA back_color = gm_rgba_from_u32(back_color_u32);
        
        return gm_rgba_alpha_blend(front_color, back_color);
    }
#endif

#if defined(GM_IMPL_VECTOR)
    float gm_vec2_dot(GM_Vec2 A, GM_Vec2 B) {
        return (A.x * B.x) + (A.y * B.y);
    }

    float gm_vec3_dot(GM_Vec3 A, GM_Vec3 B) {
        return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
    }

    float gm_vec4_dot(GM_Vec4 A, GM_Vec4 B) {
        return (A.x * B.x) + (A.y * B.y) + (A.z * B.z) + (A.w * B.w);
    }

    GM_Vec3 gm_vec3_cross(GM_Vec3 A, GM_Vec3 B) {
        GM_Vec3 ret;
        ret.x = A.y * B.z - A.z * B.y;
        ret.y = A.z * B.x - A.x * B.z;
        ret.z = A.x * B.y - A.y * B.x;
        return ret;
    }

    float gm_vec2_magnitude(GM_Vec2 A) {
        return sqrtf((A.x*A.x) + (A.y*A.y));
    }

    float gm_vec3_magnitude(GM_Vec3 A) {
        return sqrtf((A.x*A.x) + (A.y*A.y) + (A.z*A.z));
    }

    float gm_vec4_magnitude(GM_Vec4 A) {
        return sqrtf((A.x*A.x) + (A.y*A.y) + (A.z*A.z) + (A.w*A.w));
    }

    GM_Vec2 gm_vec2_normalize(GM_Vec2 A) {
        GM_Vec2 ret;
        const float magnitude = gm_vec2_magnitude(A);
        if (magnitude == 0) return (GM_Vec2){0,0};
        ret.x = A.x / magnitude;
        ret.y = A.y / magnitude;

        return ret;
    }

    GM_Vec3 gm_vec3_normalize(GM_Vec3 A) {
        GM_Vec3 ret;
        const float magnitude = gm_vec3_magnitude(A);
        if (magnitude == 0) return (GM_Vec3){0,0,0};
        ret.x = A.x / magnitude;
        ret.y = A.y / magnitude;
        ret.z = A.z / magnitude;

        return ret;
    }

    GM_Vec4 gm_vec4_normalize(GM_Vec4 A) {
        GM_Vec4 ret;
        const float magnitude = gm_vec4_magnitude(A);
        if (magnitude == 0) return (GM_Vec4){0,0,0,0};
        ret.x = A.x / magnitude;
        ret.y = A.y / magnitude;
        ret.z = A.z / magnitude;
        ret.w = A.w / magnitude;

        return ret;
    }

    float gm_vec2_distance(GM_Vec2 A, GM_Vec2 B) {
        return sqrtf(SQUARE(B.x - A.x) + SQUARE(B.y - A.y));
    }

    float gm_vec3_distance(GM_Vec3 A, GM_Vec3 B) {
        return sqrtf(SQUARE(B.x - A.x) + SQUARE(B.y - A.y) + SQUARE(B.z - A.z));
    }

    float gm_vec4_distance(GM_Vec4 A, GM_Vec4 B) {
        return sqrtf(SQUARE(B.x - A.x) + SQUARE(B.y - A.y) + SQUARE(B.z - A.z) + SQUARE(B.w - A.w));
    }
    
    GM_Vec2 gm_vec2_create(float x, float y) {
        GM_Vec2 ret = { .x = x, .y = y };
        return ret;
    }

    GM_Vec3 gm_vec3_create(float x, float y, float z) {
        GM_Vec3 ret = { .x = x, .y = y, .z = z };
        return ret;
    }

    GM_Vec4 gm_vec4_create(float x, float y, float z, float w) {
        GM_Vec4 ret = { .x = x, .y = y, .z = z, .w = w };
        return ret;
    }

    GM_Vec2 gm_vec2_add(GM_Vec2 A, GM_Vec2 B) {
        GM_Vec2 ret = { .x = A.x + B.x, .y = A.y + B.y };
        return ret;
    }

    GM_Vec3 gm_vec3_add(GM_Vec3 A, GM_Vec3 B) {
        GM_Vec3 ret = { .x = A.x + B.x, .y = A.y + B.y, .z = A.z + B.z };
        return ret;
    }

    GM_Vec4 gm_vec4_add(GM_Vec4 A, GM_Vec4 B) {
        GM_Vec4 ret = { .x = A.x + B.x, .y = A.y + B.y, .z = A.z + B.z, .w = A.w + B.w };
        return ret;
    }

    GM_Vec2 gm_vec2_sub(GM_Vec2 A, GM_Vec2 B) {
        GM_Vec2 ret = { .x = A.x - B.x, .y = A.y - B.y };
        return ret;
    }

    GM_Vec3 gm_vec3_sub(GM_Vec3 A, GM_Vec3 B) {
        GM_Vec3 ret = { .x = A.x - B.x, .y = A.y - B.y, .z = A.z - B.z };
        return ret;
    }

    GM_Vec4 gm_vec4_sub(GM_Vec4 A, GM_Vec4 B) {
        GM_Vec4 ret = { .x = A.x - B.x, .y = A.y - B.y, .z = A.z - B.z, .w = A.w - B.w };
        return ret;
    }

    GM_Vec2 gm_vec2_scale(GM_Vec2 v, float scale) {
        GM_Vec2 ret = { .x = v.x * scale, .y = v.y * scale };
        return ret;
    }

    GM_Vec3 gm_vec3_scale(GM_Vec3 v, float scale) {
        GM_Vec3 ret = { .x = v.x * scale, .y = v.y * scale, .z = v.z * scale };
        return ret;
    }

    GM_Vec4 gm_vec4_scale(GM_Vec4 v, float scale) {
        GM_Vec4 ret = { .x = v.x * scale, .y = v.y * scale, .z = v.z * scale, .w = v.w * scale };
        return ret;
    }

    GM_Vec2 gm_vec2_lerp(GM_Vec2 A, GM_Vec2 B, float t) {
        return gm_vec2_add(A, gm_vec2_scale(gm_vec2_add(B, gm_vec2_scale(A, -1.0f)), t));
    }

    GM_Vec3 gm_vec3_lerp(GM_Vec3 A, GM_Vec3 B, float t) {
        return gm_vec3_add(A, gm_vec3_scale(gm_vec3_add(B, gm_vec3_scale(A, -1.0f)), t));
    }

    GM_Vec4 gm_vec4_lerp(GM_Vec4 A, GM_Vec4 B, float t) {
        return gm_vec4_add(A, gm_vec4_scale(gm_vec4_add(B, gm_vec4_scale(A, -1.0f)), t));
    }

    GM_Vec2 gm_vec2_projection(GM_Vec2 A, GM_Vec2 B) {
        float dot_product = gm_vec2_dot(A, B);
        float B_mag_sq = gm_vec2_dot(B, B);
        if (B_mag_sq == 0) return (GM_Vec2){0,0};
        return gm_vec2_scale(B, dot_product / B_mag_sq);
    }

    GM_Vec3 gm_vec3_projection(GM_Vec3 A, GM_Vec3 B) {
        float dot_product = gm_vec3_dot(A, B);
        float B_mag_sq = gm_vec3_dot(B, B);
        if (B_mag_sq == 0) return (GM_Vec3){0,0,0};

        return gm_vec3_scale(B, dot_product / B_mag_sq);
    }

    GM_Vec2 gm_vec2_spline_point(GM_Vec2* spline_points, u32 spline_points_count, float t);
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

    GM_Matrix4 gm_mat4_mult(GM_Matrix4 A, GM_Matrix4 B) {
        GM_Matrix4 C = {0};

        for (int i = 0; i < 4; i++) {
            C.v[i].x += A.v[i].x * B.data[0];
            C.v[i].x += A.v[i].y * B.data[4];
            C.v[i].x += A.v[i].z * B.data[8];
            C.v[i].x += A.v[i].w * B.data[12];

            C.v[i].y += A.v[i].x * B.data[1];
            C.v[i].y += A.v[i].y * B.data[5];
            C.v[i].y += A.v[i].z * B.data[9];
            C.v[i].y += A.v[i].w * B.data[13];

            C.v[i].z += A.v[i].x * B.data[2];
            C.v[i].z += A.v[i].y * B.data[6];
            C.v[i].z += A.v[i].z * B.data[10];
            C.v[i].z += A.v[i].w * B.data[14];
            
            C.v[i].w += A.v[i].x * B.data[3];
            C.v[i].w += A.v[i].y * B.data[7];
            C.v[i].w += A.v[i].z * B.data[11];
            C.v[i].w += A.v[i].w * B.data[15];
        }
        
        return C;
    }

    GM_Matrix4 gm_mat4_translate(GM_Matrix4 mat, GM_Vec3 t) {
        GM_Matrix4 translate_matrix = {
            .data = {
                1.0f, 0.0f, 0.0f, t.x,
                0.0f, 1.0f, 0.0f, t.y,
                0.0f, 0.0f, 1.0f, t.z,
                0.0f, 0.0f, 0.0f, 1.0f
            }
        };

        return gm_mat4_mult(translate_matrix, mat);
    }

    GM_Matrix4 gm_mat4_translate_xyz(GM_Matrix4 mat, float x, float y, float z) {
        return gm_mat4_translate(mat, (GM_Vec3){.x=x, .y=y, .z=z});
    }

    GM_Matrix4 gm_mat4_scale(GM_Matrix4 mat, GM_Vec3 s) {
        GM_Matrix4 scale_matrix = {
            .data = {
                s.x,  0.0f, 0.0f, 0.0f,
                0.0f, s.y,  0.0f, 0.0f,
                0.0f, 0.0f, s.z,  0.0f,
                0.0f, 0.0f, 0.0f, 1.0f 
            }
        };
        return gm_mat4_mult(scale_matrix, mat);
    }

    GM_Matrix4 gm_mat4_scale_xyz(GM_Matrix4 mat, float x, float y, float z) {
        return gm_mat4_scale(mat, (GM_Vec3){.x=x, .y=y, .z=z});
    }

    GM_Matrix4 gm_mat4_rotate(GM_Matrix4 mat, float degrees, GM_Vec3 axis) {
        float rad = DEGREES_TO_RAD(degrees);
        float c = cosf(rad);
        float s = sinf(rad);
        float t = 1.0f - c;

        GM_Vec3 norm_axis = gm_vec3_normalize(axis);
        float x = norm_axis.x;
        float y = norm_axis.y;
        float z = norm_axis.z;

        GM_Matrix4 rot = {
            .data = {
                t * x * x + c,     t * x * y - z * s, t * x * z + y * s, 0.0f,
                t * x * y + z * s, t * y * y + c,     t * y * z - x * s, 0.0f,
                t * x * z - y * s, t * y * z + x * s, t * z * z + c,     0.0f,
                0.0f,              0.0f,              0.0f,              1.0f
            }
        };

        return gm_mat4_mult(rot, mat);
    }


    GM_API GM_Matrix4 gm_mat4_rotate_xyz(GM_Matrix4 mat, float degrees, float x1, float y1, float z1) {
        return gm_mat4_rotate(mat, degrees, (GM_Vec3){.x=x1, .y=y1, .z=z1});
    }

    GM_API GM_Matrix4 gm_mat4_perspective(float fov_degrees, float aspect, float near_plane, float far_plane) {
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

    GM_Matrix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane);

    GM_Matrix4 gm_mat4_look_at(GM_Vec3 camera_position, GM_Vec3 target_position, GM_Vec3 world_up) {
        GM_Vec3 forward = gm_vec3_normalize(gm_vec3_sub(camera_position, target_position));
        GM_Vec3 right   = gm_vec3_normalize(gm_vec3_cross(world_up, forward));
        GM_Vec3 up      = gm_vec3_cross(forward, right);

        float dot_right   = -gm_vec3_dot(right, camera_position);
        float dot_up      = -gm_vec3_dot(up, camera_position);
        float dot_forward = -gm_vec3_dot(forward, camera_position);

        GM_Matrix4 ret = {
            .data = {
                right.x,   right.y,   right.z,   dot_right,
                up.x,      up.y,      up.z,      dot_up,
                forward.x, forward.y, forward.z, dot_forward,
                0.0f,      0.0f,      0.0f,      1.0f
            }
        };

        return ret;
    }

    GM_Matrix4 gm_mat4_inverse(GM_Matrix4 m, bool* success) {
        if (success) {
            *success = false;
        }
        
        return gm_mat4_identity();
    }

    GM_API GM_Matrix4 gm_mat4_transpose(GM_Matrix4 m) {
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
#endif

#if defined(GM_IMPL_QUATERNION)
    GM_Quaternion gm_quat_create(GM_Vec3 axis, float theta) {
        GM_Quaternion ret = {0};

        float radians = DEGREES_TO_RAD(theta);
        float half_angle = radians / 2.0f;
        ret.w = cosf(half_angle);
        GM_Vec3 norm_axis = gm_vec3_normalize(axis);
        ret.v = gm_vec3_scale(norm_axis, sinf(half_angle));

        return ret;
    }

    GM_Quaternion gm_quat_inverse(GM_Quaternion quat) {
        GM_Quaternion ret = {0};

        float mag_sq = (quat.w * quat.w) + gm_vec3_dot(quat.v, quat.v);
        if (mag_sq == 0.0f) return (GM_Quaternion){0,0,0,0};

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

    // Rotate a vector with this quaternion (q * p * q_inverse)
    GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec) {
        GM_Quaternion p;
        p.w = 0.0f;
        p.v = vec;

        // Efficient vector rotation:
        // Quat * Vec * Quat_Inverse
        // vprime = q * v * q_inv
        // Where v is a pure quaternion (0, vec.x, vec.y, vec.z)
        GM_Quaternion temp = gm_quat_mult(quat, p);
        return gm_quat_mult(temp, gm_quat_inverse(quat)).v;
    }

    GM_Quaternion gm_slerp(GM_Quaternion a, GM_Quaternion b, float t);
#endif

#if defined(GM_IMPL_UTILITY)
    float gm_lerp(float a, float b, float t) {
        return a + ((b - a) * t);
    }

    float gm_inverse_lerp(float a, float b, float value) {
        if (fabsf(a - b) < 0.00001f) return 0.0f; // Avoid division by zero
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
    float gm_ease_in_sine(float t) {
        return 1.0f - cosf((t * PI) / 2.0f);
    }

    float gm_ease_out_sine(float t) {
        return sinf((t * PI) / 2.0f);
    }

    float gm_ease_in_out_sine(float t) {
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
        return t < 0.5f ? 4.0f * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 3.0f) / 2.0f;
    }

    float gm_ease_in_quart(float t) {
        return t * t * t * t;
    }

    float gm_ease_out_quart(float t) {
        return 1.0f - powf(1.0f - t, 4.0f);
    }

    float gm_ease_in_out_quart(float t) {
        return t < 0.5f ? 8.0f * t * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 4.0f) / 2.0f;
    }

    float gm_ease_in_quint(float t) {
        return t * t * t * t * t;
    }

    float gm_ease_out_quint(float t) {
        return 1.0f + powf(t - 1.0f, 5.0f);
    }

    float gm_ease_in_out_quint(float t) {
        return t < 0.5f ? 16.0f * t * t * t * t * t : 1.0f - powf(-2.0f * t + 2.0f, 5.0f) / 2.0f;
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
        if (t < 0.5f) return powf(2.0f, 20.0f * t - 10.0f) / 2.0f;
        return (2.0f - powf(2.0f, -20.0f * t + 10.0f)) / 2.0f;
    }

    float gm_ease_in_circ(float t) {
        return 1.0f - sqrtf(1.0f - (t * t));
    }

    float gm_ease_out_circ(float t) {
        return sqrtf(1.0f - powf(t - 1.0f, 2.0f));
    }

    float gm_ease_in_out_circ(float t) {
        if (t < 0.5f) return (1.0f - sqrtf(1.0f - powf(2.0f * t, 2.0f))) / 2.0f;
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
        float c2 = c1 * 1.525f;
        if (t < 0.5f) return (powf(2.0f * t, 2.0f) * ((c2 + 1.0f) * 2.0f * t - c2)) / 2.0f;
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
        float c5 = (2.0f * PI) / 4.5f;
        if (t == 0.0f) return 0.0f;
        if (t == 1.0f) return 1.0f;
        if (t < 0.5f) return -(powf(2.0f, 20.0f * t - 10.0f) * sinf((20.0f * t - 11.125f) * c5)) / 2.0f;
        return (powf(2.0f, -20.0f * t + 10.0f) * sinf((20.0f * t - 11.125f) * c5)) / 2.0f + 1.0f;
    }

    float gm_ease_in_bounce(float t) {
        return 1.0f - gm_ease_out_bounce(1.0f - t);
    }

    float gm_ease_out_bounce(float t) {
        float n1 = 7.5625f;
        float d1 = 2.75f;
        if (t < 1.0f / d1) {
            return n1 * t * t;
        } else if (t < 2.0f / d1) {
            return n1 * (t -= 1.5f / d1) * t + 0.75f;
        } else if (t < 2.5f / d1) {
            return n1 * (t -= 2.25f / d1) * t + 0.9375f;
        } else {
            return n1 * (t -= 2.625f / d1) * t + 0.984375f;
        }
    }

    float gm_ease_in_out_bounce(float t) {
        return t < 0.5f
            ? (1.0f - gm_ease_out_bounce(1.0f - 2.0f * t)) / 2.0f
            : (1.0f + gm_ease_out_bounce(2.0f * t - 1.0f)) / 2.0f;
    }
#endif

#if defined(GM_IMPL_COLLISION)
    // SAT collision

    bool gm_collision2d_circles(GM_Circle2D c1, GM_Circle2D c2, GM_CollisionInfo2D* collision_info) {
        float distance = gm_vec2_distance(c1.position, c2.position);

        float total_radius = (float)c1.radius + (float)c2.radius;
        if (distance >= total_radius) {
            return false;
        }

        if (collision_info) {
            collision_info->normal = gm_vec2_sub(c2.position, c1.position);
            collision_info->depth = total_radius - distance;
        }

        return true;
    }

    GM_Collider2D gm_collider2d_circle_create(GM_Circle2D circle) {
        GM_Collider2D ret = {0};
        ret.type = GM_COLLIDER_CIRCLE;
        ret.circle = circle;

        return ret;
    }

    GM_Collider2D gm_collider2d_aabb_create(GM_Rectangle2D aabb) {
        GM_Collider2D ret = {0};
        ret.type = GM_COLLIDER_AABB;
        ret.aabb = aabb;

        return ret;
    }

    /*
    bool gm_aabb_point_colliding(GM_Vec3 point, GM_AABB aabb) {
        return (
            point.x >= aabb.min.x &&
            point.x <= aabb.max.x &&
            point.y >= aabb.min.y &&
            point.y <= aabb.max.y &&
            point.z >= aabb.min.z &&
            point.z <= aabb.max.z
        );
    }

    bool gm_aabb_aabb_colliding(GM_AABB a, GM_AABB b) {
        return (
            a.min.x <= b.max.x &&
            a.max.x >= b.min.x &&
            a.min.y <= b.max.y &&
            a.max.y >= b.min.y &&
            a.min.z <= b.max.z &&
            a.max.z >= b.min.z
        );
    }
    */
#endif

#if defined(GM_IMPL_SHAPES)
    GM_Rectangle2D gm_rectangle2d_create(float x, float y, u32 width, u32 height) {
        GM_Rectangle2D ret = { .position.x = x, .position.y = y, .width = width, .height = height };
        return ret;
    }

    GM_Rectangle3D gm_rectangle3d_create(float x, float y, float z, u32 length, u32 width, u32 height) {
        GM_Rectangle3D ret = { .position.x = x, .position.y = y, .position.z = z, .length = length, .width = width, .height = height };
        return ret;
    }

    bool gm_rectangle_check_aabb_collision(GM_Rectangle2D rect1, GM_Rectangle2D rect2) {
        if (rect1.position.x < rect2.position.x + rect2.width && rect1.position.x + rect1.width > rect2.position.x &&
            rect1.position.y < rect2.position.y + rect2.height && rect1.position.y + rect1.height > rect2.position.y) {
            return true;
        }
        return false;
    }

    GM_Circle2D gm_circle2d_create(float x, float y, u32 radius) {
        GM_Circle2D ret = { .position.x = x, .position.y = y, .radius = radius };
        return ret;
    }

    GM_Circle3D gm_circle3d_create(float x, float y, float z, u32 radius) {
        GM_Circle3D ret = { .position.x = x, .position.y = y, .position.z = z, .radius = radius };
        return ret;
    }
#endif

#if defined(GM_IMPL_PHYSICS)
    internal GM_Vec2 gravity_vector = {0.0f, 9.81f};

    internal void gm_physics2d_resolve_collisions(GM_PhysicsObject2D* objects, int object_count) {
        for (int i = 0; i < object_count - 1; i++) {
            GM_PhysicsObject2D* object_a = &objects[i];
            for (int j = i + 1; j < object_count; j++) {
                GM_PhysicsObject2D* object_b = &objects[j];

                if (object_a->collider.type == GM_COLLIDER_CIRCLE && object_b->collider.type == GM_COLLIDER_CIRCLE) {
                    GM_CollisionInfo2D collision_info = {0};
                    if (gm_collision2d_circles(object_a->collider.circle, object_b->collider.circle, &collision_info)) {
                        object_a->rb.position = gm_vec2_add(object_a->rb.position, gm_vec2_scale(collision_info.normal, collision_info.depth));
                    }
                }
            }
        }
    }

    GM_RigidBody2D gm_physics2d_rb_create(GM_Vec2 position, float mass) {
        GM_RigidBody2D ret;
        ret.position = position;
        ret.velocity = gm_vec2_create(0, 0);
        ret.acceleration = gm_vec2_create(0, 0);
        ret.mass = mass;

        return ret;
    }

    GM_RigidBody2D gm_physics2d_rb_create_xy(float x, float y, float mass) {
        GM_RigidBody2D ret;
        ret.position.x = x;
        ret.position.y = y;
        ret.velocity = gm_vec2_create(0, 0);
        ret.acceleration = gm_vec2_create(0, 0);
        ret.mass = mass;

        return ret;
    }

    GM_PhysicsObject2D gm_physics2d_object_create(GM_Vec2 position, float mass, GM_Collider2D collider) {
        GM_PhysicsObject2D ret = {0};
        ret.rb = gm_physics2d_rb_create(position, mass);
        ret.collider = collider;
        if (collider.type == GM_COLLIDER_CIRCLE) {
            ret.collider.circle.position = position;
        } else if (collider.type == GM_COLLIDER_AABB) {
            ret.collider.aabb.position = position;
        }

        return ret;
    }

    GM_PhysicsObject2D gm_physics2d_object_create_xy(float x, float y, float mass, GM_Collider2D collider) {
        GM_PhysicsObject2D ret = {0};
        ret.rb = gm_physics2d_rb_create_xy(x, y, mass);
        ret.collider = collider;
        if (collider.type == GM_COLLIDER_CIRCLE) {
            ret.collider.circle.position.x = x;
            ret.collider.circle.position.y = y;
        } else if (collider.type == GM_COLLIDER_AABB) {
            ret.collider.aabb.position.x = x;
            ret.collider.aabb.position.y = y;
        }

        return ret;
    }

    void gm_physics2d_rb_apply_velocity(GM_PhysicsObject2D* obj, GM_Vec2 velocity) {
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, velocity);
    }

    void gm_physics2d_apply_velocity_xy(GM_PhysicsObject2D* obj, float velocity_x, float velocity_y) {
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, gm_vec2_create(velocity_x, velocity_y));
    }

    void gm_physics2d_apply_force(GM_PhysicsObject2D* obj, GM_Vec2 force) {
        obj->rb.force = gm_vec2_add(obj->rb.force, force);
    }

    void gm_physics2d_apply_force_xy(GM_PhysicsObject2D* obj, float force_x, float force_y) {
        obj->rb.force = gm_vec2_add(obj->rb.force, gm_vec2_create(force_x, force_y));
    }

    void gm_physics2d_update(GM_PhysicsObject2D* obj, float dt) {
        obj->rb.acceleration = gm_vec2_scale(obj->rb.force, 1 / obj->rb.mass);
        obj->rb.position = gm_vec2_add(obj->rb.position, gm_vec2_scale(obj->rb.velocity, dt));
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, gm_vec2_scale(obj->rb.acceleration, dt));
    }

    void gm_physics2d_simulate(GM_PhysicsObject2D* objects, int object_count, float dt) {
        for (int i = 0; i < object_count; i++) {
            GM_PhysicsObject2D* obj = &objects[i];
            gm_physics2d_apply_force(obj, gm_vec2_scale(gravity_vector, obj->rb.mass));

            gm_physics2d_resolve_collisions(objects, object_count);

            gm_physics2d_update(obj, dt);
        }
    }
#endif