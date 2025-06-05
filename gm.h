#pragma once

#ifdef __cplusplus
    #define GM_API extern "C"
#else
    #define GM_API
#endif

#if defined(GM_IMPL)
    #define GM_IMPL_COLOR
    #define GM_IMPL_VECTOR
    #define GM_IMPL_EULER
    #define GM_IMPL_MATRIX
    #define GM_IMPL_SHAPES
    #define GM_IMPL_INTERSECTION
    #define GM_IMPL_QUATERNION
    #define GM_IMPL_UTILITY
    #define GM_IMPL_EASE_FUNCTIONS
    #define GM_IMPL_COLLISION
    #define GM_IMPL_PHYSICS
#endif

#define GM_INCLUDE_TYPES
#define GM_INCLUDE_COLOR
#define GM_INCLUDE_VECTOR
#define GM_INCLUDE_EULER
#define GM_INCLUDE_MATRIX
#define GM_INCLUDE_SHAPES
#define GM_INCLUDE_INTERSECTION
#define GM_INCLUDE_QUATERNION
#define GM_INCLUDE_UTILITY
#define GM_INCLUDE_EASE_FUNCTIONS
#define GM_INCLUDE_COLLISION
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
    typedef struct GM_RGBA {
       union {
            struct {
                u8 r;
                u8 g;
                u8 b;
                u8 a;
            };

            u8 c[4];
            u32 hex;
        };
    } GM_RGBA;

    #ifdef __cplusplus
        #define GM_RGBALit(r, g, b, a) (GM_RGBA{r, g, b, a})
    #else
        #define GM_RGBALit(r, g, b, a) ((GM_RGBA){r, g, b, a})
    #endif

    GM_API GM_RGBA gm_rgba_from_u32(u32 color);
    GM_API GM_RGBA gm_argb_from_u32(u32 color);
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

    #define GM_RGBA_BLACK GM_RGBALit(0, 0, 0, 255)
    #define GM_RGBA_RED GM_RGBALit(255, 0, 0, 255)
    #define GM_RGBA_BLUE GM_RGBALit(0, 0, 255, 255)
    #define GM_RGBA_GREEN GM_RGBALit(0, 255, 0, 255)
    #define GM_RGBA_WHITE GM_RGBALit(255, 255, 255, 255)
    #define GM_RGBA_PINK GM_RGBALit(255, 105, 180, 255)
    #define GM_RGBA_LIME GM_RGBALit(0, 255, 128, 255)
    #define GM_RGBA_CYAN GM_RGBALit(0, 255, 255, 255)
    #define GM_RGBA_PURPLE GM_RGBALit(128, 0, 128, 255)
    #define GM_RGBA_YELLOW GM_RGBALit(255, 255, 0, 255)
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

    #ifdef __cplusplus
        #define GM_Vec2Lit(x, y) GM_Vec2{x, y}
        #define GM_Vec3Lit(x, y, z) GM_Vec3{x, y, z}
        #define GM_Vec4Lit(x, y, z) GM_Vec4{x, y, z, w}
    #else
        #define GM_Vec2Lit(x, y) (GM_Vec2){x, y}
        #define GM_Vec3Lit(x, y, z) (GM_Vec3){x, y, z}
        #define GM_Vec4Lit(x, y, z, w) (GM_Vec4){x, y, z, w}
    #endif

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
    GM_API GM_Vec2 gm_vec2_scale_non_uniform(GM_Vec2 v, GM_Vec2 s);
    GM_API GM_Vec3 gm_vec3_scale_non_uniform(GM_Vec3 v, GM_Vec3 s);
    GM_API GM_Vec4 gm_vec4_scale_non_uniform(GM_Vec4 v, GM_Vec4 s);
    GM_API GM_Vec2 gm_vec2_scale_xyz(GM_Vec2 v, float scale_x, float scale_y);
    GM_API GM_Vec3 gm_vec3_scale_xyz(GM_Vec3 v, float scale_x, float scale_y, float sacle_z);
    GM_API GM_Vec4 gm_vec4_scale_xyz(GM_Vec4 v, float scale_x, float scale_y, float sacle_z, float scale_w);

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
    GM_API float gm_vec2_magnitude_squared(GM_Vec2 A);

    GM_API float gm_vec3_magnitude(GM_Vec3 A);
    GM_API float gm_vec3_magnitude_squared(GM_Vec3 A);

    GM_API float gm_vec4_magnitude(GM_Vec4 A);
    GM_API float gm_vec4_magnitude_squared(GM_Vec4 A);

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

    /**
     * @brief ensure that left, right, top, and bottom account for aspect ratio by multiply it
     * 
     * @param left 
     * @param right 
     * @param bottom 
     * @param top 
     * @param near_plane 
     * @param far_plane 
     * @return GM_API 
     */
    GM_API GM_Matrix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane);
    GM_API GM_Matrix4 gm_mat4_look_at(GM_Vec3 camera_position, GM_Vec3 target_position, GM_Vec3 world_up);

    GM_API GM_Matrix4 gm_mat4_mult(GM_Matrix4 A, GM_Matrix4 B);
    GM_API GM_Matrix4 gm_mat4_inverse(GM_Matrix4 m, bool* success);
    GM_API GM_Matrix4 gm_mat4_transpose(GM_Matrix4 m);

    GM_API bool gm_mat4_equal(GM_Matrix4 A, GM_Matrix4 B);
    GM_API GM_Matrix4 gm_mat4_transform(GM_Vec3 scale, float theta, GM_Vec3 rotation_axis, GM_Vec3 translation);
    GM_API GM_Matrix4 gm_mat4_inverse_transform(GM_Vec3 scale, float theta, GM_Vec3 rotation_axis, GM_Vec3 translation);
#endif

#if defined(GM_INCLUDE_SHAPES)
    typedef struct GM_Rectangle2D {
        GM_Vec2 center;   // Center in 3D space (x, y, z)
        float width;  // along the X-axis
        float height; // along the Y-axis
    } GM_Rectangle2D;

    typedef struct GM_RectangleReference2D {
        GM_Vec2* center;   // Center in 3D space (x, y, z)
        float width;  // along the X-axis
        float height; // along the Y-axis
    } GM_RectangleReference2D;

    typedef struct GM_Rectangle3D {
        GM_Vec3 center;       // Center in 3D space (x, y, z)
        float width;  // along the X-axis
        float height; // along the Y-axis
        float length; // along the Z-axis
    } GM_Rectangle3D;

    typedef struct GM_RectangleReference3D {
        GM_Vec3* center;   // Center in 3D space (x, y, z)
        float width;  // along the X-axis
        float height; // along the Y-axis
        float length; // along the Z-axis
    } GM_RectangleReference3D;

    typedef struct GM_Circle2D {
        GM_Vec2 center;
        float radius;
    } GM_Circle2D;

    typedef struct GM_CircleReference2D {
        GM_Vec2* center;
        float radius;
    } GM_CircleReference2D;

    typedef struct GM_Circle3D {
        GM_Vec3 center;
        float radius;
    } GM_Circle3D;

    typedef struct GM_CircleReference3D {
        GM_Vec3* center;
        float radius;
    } GM_CircleReference3D;

    GM_API GM_Rectangle2D gm_rectangle2d_create(float x, float y, float width, float height);
    GM_API GM_RectangleReference2D gm_rectangle_reference2d_create(GM_Vec2* center, float width, float height);
    GM_API GM_Rectangle3D gm_rectangle3d_create(float x, float y, float z, float width, float height, float length);
    GM_API GM_RectangleReference3D gm_rectangle_reference3d_create(GM_Vec3* center, float width, float height, float length);
    GM_API bool gm_rectangle_check_aabb_collision(GM_RectangleReference2D rect1, GM_RectangleReference2D rect2);
    GM_API GM_Circle2D gm_circle2d_create(float x, float y, float radius);
    GM_API GM_CircleReference2D gm_circle_reference2d_create(GM_Vec2* center, float radius);
    GM_API GM_Circle3D gm_circle3d_create(float x, float y, float z, float radius);
    GM_API GM_CircleReference3D gm_circle_reference3d_create(GM_Vec3* center, float radius);
#endif

#if defined(GM_INCLUDE_INTERSECTION)
    // Found at: https://imois.in/posts/line-intersections-with-cross-products/
    GM_API bool gm_intersection2d_line_line(GM_Vec2 a, GM_Vec2 b, GM_Vec2 c, GM_Vec2 d, GM_Vec2* intersection);
    GM_API bool gm_intersection2d_line_aabb(GM_Vec2 p0, GM_Vec2 p1, GM_RectangleReference2D aabb, GM_Vec2* inPoint, GM_Vec2* outPoint);
    GM_API bool gm_intersection3d_line_aabb(GM_Vec3 p0, GM_Vec3 p1, GM_RectangleReference3D aabb, GM_Vec3* inPoint, GM_Vec3* outPoint);
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

    GM_API GM_Quaternion gm_quat_create(float theta, GM_Vec3 axis);
    GM_API GM_Quaternion gm_quat_inverse(GM_Quaternion quat);
    GM_API GM_Quaternion gm_quat_mult(GM_Quaternion q1, GM_Quaternion q2);
    GM_API GM_Quaternion gm_quat_from_euler(GM_Vec3 euler_angles_degrees);
    GM_API GM_Matrix4 gm_quat_to_mat4(GM_Quaternion q);
    GM_API GM_Vec3 gm_quat_vector_mult(GM_Quaternion quat, GM_Vec3 vec);
    GM_API void gm_quat_to_axis_angle(GM_Quaternion quat, float* theta, GM_Vec3* vec);
    GM_API GM_Quaternion gm_quat_sub(GM_Quaternion a, GM_Quaternion b);
    GM_API GM_Quaternion gm_quat_add(GM_Quaternion a, GM_Quaternion b);
    GM_API GM_Quaternion gm_quat_scale(GM_Quaternion a, float scale);
   GM_API  GM_Quaternion gm_quat_add_scalar(GM_Quaternion a, float scalar);
   GM_API  GM_Quaternion gm_quat_sub_scalar(GM_Quaternion a, float scalar);
    GM_API GM_Quaternion gm_quat_power(GM_Quaternion q1, float t);
    GM_API GM_Quaternion gm_quat_normalize(GM_Quaternion q);
    GM_API float gm_quat_dot(GM_Quaternion q, GM_Quaternion r);
    GM_API GM_Quaternion gm_quat_slerp(GM_Quaternion q, GM_Quaternion r, float t);
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

#if defined(GM_INCLUDE_EULER)
    GM_API GM_Vec2 gm_euler_to_vec2(float yaw, float pitch);
    GM_API GM_Vec3 gm_euler_to_vec3(float yaw, float pitch);
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
            GM_RectangleReference2D aabb;
            GM_CircleReference2D circle;
        };

        u64 collision_mask; // Checking other layers to see if should actually collide or ignore based on their layer_mask
        u64 layer_mask; // is used by other colliders to see if it should actually collide or ignore

        // collision_enter(GM_Collider2D* other, GM_Collider2D* other);
        // collision_persisting(GM_Collider2D* other, GM_Collider2D* other, float time);
        // collision_function_leave(GM_Collider2D* other, GM_Collider2D* other);
    } GM_Collider2D;

    GM_API bool gm_collision2d_circles(GM_CircleReference2D c1, GM_CircleReference2D c2, GM_CollisionInfo2D* collision_info);

    // Date: May 28, 2025
    // TODO(Jovanni): SAT collision

    GM_API GM_Collider2D gm_collider2d_circle_create(GM_CircleReference2D circle);
    GM_API GM_Collider2D gm_collider2d_aabb_create(GM_RectangleReference2D aabb);

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
        GM_Vec2* position; // reference to the position
        GM_Vec2 velocity;
        GM_Vec2 acceleration;
        GM_Vec2 force;
        float mass; // In kilograms (KG)
    } GM_RigidBody2D;

    typedef struct GM_PhysicsObject2D {
        GM_RigidBody2D rb;
        GM_Collider2D collider;
    } GM_PhysicsObject2D;

    GM_API void gm_physics2d_resolve_collisions(GM_PhysicsObject2D* objects, int object_count);

    GM_API GM_RigidBody2D gm_physics2d_rb_create(GM_Vec2* position, float mass);

    GM_API GM_PhysicsObject2D gm_physics2d_object_create(GM_Vec2* position, float mass, GM_Collider2D collider);

    GM_API void gm_physics2d_apply_velocity(GM_PhysicsObject2D* object, GM_Vec2 velocity);
    GM_API void gm_physics2d_apply_velocity_xy(GM_PhysicsObject2D* object, float velocity_x, float velocity_y);

    GM_API void gm_physics2d_apply_force(GM_PhysicsObject2D* object, GM_Vec2 force);
    GM_API void gm_physics2d_apply_force_xy(GM_PhysicsObject2D* object, float force_x, float force_y);

    GM_API void gm_physics2d_update(GM_PhysicsObject2D* obj, float dt);

    /**
     * @brief The gravity force vector points towards position_a or mass_a
     * 
     * @param G 
     * @param position_a
     * @param position_a 
     * @param mass_a 
     * @param position_b 
     * @param mass_b 
     * @return GM_API 
     */
    GM_API GM_Vec2 gm_physics2d_gravity_force(const double G, float min_distance, GM_Vec2 position_a, float mass_a, GM_Vec2 position_b, float mass_b);
#endif

//
// ===================================================== GM_IMPL =====================================================
//

#if defined(GM_IMPL_COLOR)
    GM_RGBA gm_rgba_from_u32(u32 color) {
        GM_RGBA ret = {0};
        ret.r = (u8)((color >> 0) & 0xFF);
        ret.g = (u8)((color >> 8) & 0xFF);
        ret.b = (u8)((color >> 16) & 0xFF);
        ret.a = (u8)((color >> 24) & 0xFF);

        return ret;
    }

    GM_RGBA gm_argb_from_u32(u32 color) {
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

    float gm_vec2_magnitude_squared(GM_Vec2 A) {
        return (A.x*A.x) + (A.y*A.y);
    }

    float gm_vec3_magnitude(GM_Vec3 A) {
        return sqrtf((A.x*A.x) + (A.y*A.y) + (A.z*A.z));
    }

    float gm_vec3_magnitude_squared(GM_Vec3 A) {
        return (A.x*A.x) + (A.y*A.y) + (A.z*A.z);
    }

    float gm_vec4_magnitude(GM_Vec4 A) {
        return sqrtf((A.x*A.x) + (A.y*A.y) + (A.z*A.z) + (A.w*A.w));
    }

    float gm_vec4_magnitude_squared(GM_Vec4 A) {
        return (A.x*A.x) + (A.y*A.y) + (A.z*A.z) + (A.w*A.w);
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

    GM_Vec2 gm_vec2_scale_non_uniform(GM_Vec2 v, GM_Vec2 s) {
        GM_Vec2 ret = { .x = v.x * s.x, .y = v.y * s.y };
        return ret;
    }

    GM_Vec2 gm_vec2_scale_xyz(GM_Vec2 v, float scale_x, float scale_y) {
        GM_Vec2 ret = { .x = v.x * scale_x, .y = v.y * scale_y };
        return ret;
    }

    GM_Vec3 gm_vec3_scale(GM_Vec3 v, float scale) {
        GM_Vec3 ret = { .x = v.x * scale, .y = v.y * scale, .z = v.z * scale };
        return ret;
    }

    GM_Vec3 gm_vec3_scale_non_uniform(GM_Vec3 v, GM_Vec3 s) {
        GM_Vec3 ret = { .x = v.x * s.x, .y = v.y * s.y, .z = v.z * s.z };
        return ret;
    }

    GM_Vec3 gm_vec3_scale_xyz(GM_Vec3 v, float scale_x, float scale_y, float sacle_z) {
        GM_Vec3 ret = { .x = v.x * scale_x, .y = v.y * scale_y, .z = v.z * sacle_z };
        return ret;
    }

    GM_Vec4 gm_vec4_scale(GM_Vec4 v, float scale) {
        GM_Vec4 ret = { .x = v.x * scale, .y = v.y * scale, .z = v.z * scale, .w = v.w * scale };
        return ret;
    }

    GM_Vec4 gm_vec4_scale_non_uniform(GM_Vec4 v, GM_Vec4 s) {
        GM_Vec4 ret = { .x = v.x * s.x, .y = v.y * s.y, .z = v.z * s.z, .w = v.w * s.w };
        return ret;
    }

    GM_Vec4 gm_vec4_scale_xyz(GM_Vec4 v, float scale_x, float scale_y, float sacle_z, float scale_w) {
        GM_Vec4 ret = { .x = v.x * scale_x, .y = v.y * scale_y, .z = v.z * sacle_z, .w = v.w * scale_w };
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
        float c = fabsf(cosf(rad));
        float s = fabsf(sinf(rad));
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

    // Found at: https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/orthographic-projection-matrix.html
    GM_Matrix4 gm_mat4_orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane) {
        const float A = 2 / (right - left);
        const float B = 2 / (top - bottom);
        const float C = -2 / (far_plane - near_plane);
        const float D = -((right + left) / (right - left));
        const float E = -((top + bottom) / (top - bottom));
        const float F = -((far_plane + near_plane) / (far_plane - near_plane));

        GM_Matrix4 ret = {
            .data = {
                A,  0,  0,  0,
                0,  B,  0,  0,
                0,  0,  C,  0,
                D,  E,  F,  1
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

    
    bool gm_mat4_equal(GM_Matrix4 A, GM_Matrix4 B) {
        for (int i = 0; i < 16; i++) {
            if (!NEAR_ZERO(A.data[i] - B.data[i])) {
                return false;
            }
        }

        return true;
    }

    GM_Matrix4 gm_mat4_transform(GM_Vec3 scale, float theta, GM_Vec3 rotation_axis, GM_Vec3 translation) {
        GM_Matrix4 scale_matrix = gm_mat4_scale(gm_mat4_identity(), scale);
        GM_Matrix4 rotation_matrix = gm_mat4_rotate(gm_mat4_identity(), theta, rotation_axis);
        GM_Matrix4 translation_matrix = gm_mat4_translate(gm_mat4_identity(), translation);

        return gm_mat4_mult(translation_matrix, gm_mat4_mult(rotation_matrix, scale_matrix));
    }

    GM_Matrix4 gm_mat4_inverse_transform(GM_Vec3 scale, float theta, GM_Vec3 rotation_axis, GM_Vec3 translation) {
        GM_Matrix4 inverse_scale_matrix = gm_mat4_scale(gm_mat4_identity(), gm_vec3_scale_xyz(scale, 1 / scale.x, 1 / scale.y, 1 / scale.z));
        GM_Matrix4 inverse_rotation_matrix = gm_mat4_transpose(gm_mat4_rotate(gm_mat4_identity(), theta, rotation_axis));
        GM_Matrix4 inverse_translation_matrix = gm_mat4_translate(gm_mat4_identity(), gm_vec3_scale(translation, -1));

        return gm_mat4_mult(inverse_scale_matrix, gm_mat4_mult(inverse_rotation_matrix, inverse_translation_matrix));
    }
#endif

#if defined(GM_IMPL_SHAPES)
    GM_Rectangle2D gm_rectangle2d_create(float x, float y, float width, float height) {
        GM_Rectangle2D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    GM_RectangleReference2D gm_rectangle_reference2d_create(GM_Vec2* center, float width, float height) {
        GM_RectangleReference2D ret;
        ret.center = center;
        ret.width = width;
        ret.height = height;

        return ret;
    }

    GM_Rectangle3D gm_rectangle3d_create(float x, float y, float z, float width, float height, float length) {
        GM_Rectangle3D ret;
        ret.center.x = x;
        ret.center.y = y;
        ret.center.z = z;
        ret.width = width;
        ret.height = height;
        ret.length = length;

        return ret;
    }

    GM_RectangleReference3D gm_rectangle_reference3d_create(GM_Vec3* center, float width, float height, float length) {
        GM_RectangleReference3D ret;
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
    bool gm_rectangle_check_aabb_collision(GM_RectangleReference2D rect1, GM_RectangleReference2D rect2) {
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


    bool gm_intersection2d_line_aabb(GM_Vec2 p0, GM_Vec2 p1, GM_RectangleReference2D aabb, GM_Vec2* inPoint, GM_Vec2* outPoint) {
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

    bool gm_intersection3d_line_aabb(GM_Vec3 p0, GM_Vec3 p1, GM_RectangleReference3D aabb, GM_Vec3* inPoint, GM_Vec3* outPoint) {
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
        float roll_rad_half = DEGREES_TO_RAD(euler_angles_degrees.x) * 0.5f;
        float pitch_rad_half = DEGREES_TO_RAD(euler_angles_degrees.y) * 0.5f;
        float yaw_rad_half = DEGREES_TO_RAD(euler_angles_degrees.z) * 0.5f;

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
        if (NEAR_ZERO(gm_vec3_magnitude_squared(quat.v)))
            if (vec) {
                *vec = GM_Vec3Lit(1, 0, 0);
            }
        else {
            if (vec) {
                *vec = gm_vec3_normalize(quat.v);
            }
        }

        if (theta) {
            *theta = acosf(quat.w) * 2.0f;
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

    bool gm_collision2d_circles(GM_CircleReference2D c1, GM_CircleReference2D c2, GM_CollisionInfo2D* collision_info) {
        float distance = gm_vec2_distance(*c1.center, *c2.center);
        float total_radius = c1.radius + c2.radius;
        if (distance >= total_radius) {
            return false;
        }

        if (collision_info) {
            collision_info->normal = gm_vec2_normalize(gm_vec2_sub(*c1.center, *c2.center));
            collision_info->depth = total_radius - distance;
        }

        return true;
    }

    GM_Collider2D gm_collider2d_circle_create(GM_CircleReference2D circle) {
        GM_Collider2D ret = {0};
        ret.type = GM_COLLIDER_CIRCLE;
        ret.circle = circle;

        return ret;
    }

    GM_Collider2D gm_collider2d_aabb_create(GM_RectangleReference2D aabb) {
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

#if defined(GM_IMPL_PHYSICS)
    void gm_physics2d_resolve_collisions(GM_PhysicsObject2D* objects, int object_count) {
        for (int i = 0; i < object_count - 1; i++) {
            GM_PhysicsObject2D* object_a = &objects[i];
            for (int j = i + 1; j < object_count; j++) {
                GM_PhysicsObject2D* object_b = &objects[j];

                if (object_a->collider.type == GM_COLLIDER_CIRCLE && object_b->collider.type == GM_COLLIDER_CIRCLE) {
                    GM_CollisionInfo2D collision_info = {0};
                    if (gm_collision2d_circles(object_a->collider.circle, object_b->collider.circle, &collision_info)) {
                        if (collision_info.depth / object_a->collider.circle.radius >= 0.70) {
                            *object_a->rb.position = gm_vec2_add(*object_a->rb.position, gm_vec2_scale(collision_info.normal, collision_info.depth / 2.0f));
                            *object_b->rb.position = gm_vec2_add(*object_b->rb.position, gm_vec2_scale(collision_info.normal, -collision_info.depth / 2.0f));
                        } else {
                            gm_physics2d_apply_velocity(object_a, gm_vec2_scale(collision_info.normal, collision_info.depth / 2.0f));
                            gm_physics2d_apply_velocity(object_b, gm_vec2_scale(collision_info.normal, -collision_info.depth / 2.0f));
                        }
                    }
                }
            }
        }
    }

    GM_RigidBody2D gm_physics2d_rb_create(GM_Vec2* position, float mass) {
        GM_RigidBody2D ret;
        ret.position = position;
        ret.velocity = GM_Vec2Lit(0, 0);
        ret.acceleration = GM_Vec2Lit(0, 0);
        ret.mass = mass;

        return ret;
    }

    GM_PhysicsObject2D gm_physics2d_object_create(GM_Vec2* position, float mass, GM_Collider2D collider) {
        GM_PhysicsObject2D ret = {0};
        ret.rb = gm_physics2d_rb_create(position, mass);
        ret.collider = collider;
        if (collider.type == GM_COLLIDER_CIRCLE) {
            ret.collider.circle.center = position;
        } else if (collider.type == GM_COLLIDER_AABB) {
            ret.collider.aabb.center = position;
        }

        return ret;
    }

    void gm_physics2d_apply_velocity(GM_PhysicsObject2D* obj, GM_Vec2 velocity) {
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, velocity);
    }

    void gm_physics2d_apply_velocity_xy(GM_PhysicsObject2D* obj, float velocity_x, float velocity_y) {
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, GM_Vec2Lit(velocity_x, velocity_y));
    }

    void gm_physics2d_apply_force(GM_PhysicsObject2D* obj, GM_Vec2 force) {
        obj->rb.force = gm_vec2_add(obj->rb.force, force);
    }

    void gm_physics2d_apply_force_xy(GM_PhysicsObject2D* obj, float force_x, float force_y) {
        obj->rb.force = gm_vec2_add(obj->rb.force, GM_Vec2Lit(force_x, force_y));
    }

    void gm_physics2d_update(GM_PhysicsObject2D* obj, float dt) {
        obj->rb.acceleration = gm_vec2_scale(obj->rb.force, 1 / obj->rb.mass);
        *obj->rb.position = gm_vec2_add(*obj->rb.position, gm_vec2_scale(obj->rb.velocity, dt));
        obj->rb.velocity = gm_vec2_add(obj->rb.velocity, gm_vec2_scale(obj->rb.acceleration, dt));

        if (obj->collider.type == GM_COLLIDER_CIRCLE) {
            obj->collider.circle.center = obj->rb.position;
        } else if (obj->collider.type == GM_COLLIDER_AABB) {
            obj->collider.aabb.center = obj->rb.position;
        }
    }

    GM_Vec2 gm_physics2d_gravity_force(const double G, float min_distance, GM_Vec2 position_a, float mass_a, GM_Vec2 position_b, float mass_b) {
        GM_Vec2 AB = gm_vec2_normalize(gm_vec2_sub(position_b, position_a));
        float distance = gm_vec2_distance(position_a, position_b);

        if (distance < min_distance) {
            distance = min_distance;
        }

        float force_magnitude = (float)(-G * ((mass_a * mass_b) / (distance * distance)));
        GM_Vec2 gravity_force = gm_vec2_scale(AB, force_magnitude);

        return gravity_force;
    }
#endif