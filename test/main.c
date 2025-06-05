#include <gm.h>

// --- Global Test Counters ---
static int g_tests_run = 0;
static int g_tests_passed = 0;
static const char* g_current_test_name = NULL;
static const char* g_current_suite_name = NULL;

// --- Helper Functions for Assertions ---
bool float_equals_near(float a, float b) {
    return fabsf(a - b) <= EPSILON;
}

bool vec2_equals_near(GM_Vec2 a, GM_Vec2 b) {
    return float_equals_near(a.x, b.x) && float_equals_near(a.y, b.y);
}

bool vec3_equals_near(GM_Vec3 a, GM_Vec3 b) {
    return float_equals_near(a.x, b.x) && float_equals_near(a.y, b.y) && float_equals_near(a.z, b.z);
}

bool vec4_equals_near(GM_Vec4 a, GM_Vec4 b) {
    return float_equals_near(a.x, b.x) && float_equals_near(a.y, b.y) && float_equals_near(a.z, b.z) && float_equals_near(a.w, b.w);
}

bool rgba_equals(GM_RGBA a, GM_RGBA b) {
    return a.r == b.r && a.g == b.g && a.b == b.b && a.a == b.a;
}

bool quat_equals_near(GM_Quaternion a, GM_Quaternion b) {
    if (gm_quat_dot(a, b) < 0.0f) {
        b.w = -b.w;
        b.v.x = -b.v.x;
        b.v.y = -b.v.y;
        b.v.z = -b.v.z;
    }

    return fabsf(a.w - b.w) < EPSILON &&
           fabsf(a.v.x - b.v.x) < EPSILON &&
           fabsf(a.v.y - b.v.y) < EPSILON &&
           fabsf(a.v.z - b.v.z) < EPSILON;
}

// --- Assertion Macros ---
#define _FAIL_TEST() \
    do { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__); \
        return; /* Exit the current test function */ \
    } while (0)

#define ASSERT_TRUE(condition) \
    do { if (!(condition)) { _FAIL_TEST(); } } while (0)

#define ASSERT_FALSE(condition) \
    do { if (condition) { _FAIL_TEST(); } } while (0)

#define ASSERT_NEAR(actual, expected) \
    do { if (!float_equals_near(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected %f, got %f\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected, actual); \
        return; \
    } } while (0)

#define ASSERT_VEC2_NEAR(actual, expected) \
    do { if (!vec2_equals_near(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f), got (%.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, actual.x, actual.y); \
        return; \
    } } while (0)

#define ASSERT_VEC3_NEAR(actual, expected) \
    do { if (!vec3_equals_near(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f, %.4f), got (%.4f, %.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, expected.z, actual.x, actual.y, actual.z); \
        return; \
    } } while (0)

#define ASSERT_VEC4_NEAR(actual, expected) \
    do { if (!vec4_equals_near(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f, %.4f, %.4f), got (%.4f, %.4f, %.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, expected.z, expected.w, actual.x, actual.y, actual.z, actual.w); \
        return; \
    } } while (0)

#define ASSERT_RGBA_EQUAL(actual, expected) \
    do { if (!rgba_equals(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected RGBA(%u, %u, %u, %u), got RGBA(%u, %u, %u, %u)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, (unsigned int)expected.r, (unsigned int)expected.g, (unsigned int)expected.b, (unsigned int)expected.a, (unsigned int)actual.r, (unsigned int)actual.g, (unsigned int)actual.b, (unsigned int)actual.a); \
        return; \
    } } while (0)

#define ASSERT_MATRIX4_NEAR(actual, expected) \
    do { if (!gm_mat4_equal(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Matrix mismatch.\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__); \
        return; \
    } } while (0)

#define ASSERT_QUAT_NEAR(actual, expected) \
    do { if (!quat_equals_near(actual, expected)) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (w=%.4f, x=%.4f, y=%.4f, z=%.4f), got (w=%.4f, x=%.4f, y=%.4f, z=%.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.w, expected.v.x, expected.v.y, expected.v.z, actual.w, actual.v.x, actual.v.y, actual.v.z); \
        return; \
    } } while (0)


// --- Test Runner Functions ---
void test_start(const char* name) {
    g_current_test_name = name;
    g_tests_run++;
}

void test_end() {
    printf("    \033[32mPASS\033[0m: %s.%s\n", g_current_suite_name, g_current_test_name);
    g_tests_passed++;
}

void test_suite_start(const char* name) {
    g_current_suite_name = name;
    printf("\n--- Running Test Suite: %s ---\n", name);
}

void test_suite_end() {
    printf("--- Finished Test Suite: %s ---\n", g_current_suite_name);
}

// --- Individual Test Functions (Declarations) ---
// COLOR tests
void test_rgba_functions();

// VECTOR tests
void test_vec2_functions();
void test_vec3_functions();
void test_vec4_functions();

// MATRIX tests
void test_matrix4_identity();
void test_matrix4_multiply();
void test_matrix4_transpose();
void test_matrix4_inverse();
void test_matrix4_transform(); // And transform_inverse
void test_matrix4_translation_scale_rotation(); // Translate, Scale, Rotate

// QUATERNION tests
void test_quaternion_functions();

// UTILITY tests
void test_utility_functions();

// EULER tests
void test_euler_functions();

void test_intersection2d_line_aabb();
void test_intersection3d_line_aabb();

// --- Main Test Runner ---
int main() {
    printf("Starting all GM Math Unit Tests...\n");

    // Call all test suites
    test_rgba_functions();
    test_vec2_functions();
    test_vec3_functions();
    test_vec4_functions();
    test_matrix4_identity();
    test_matrix4_multiply();
    test_matrix4_transpose();
    test_matrix4_inverse();
    test_matrix4_translation_scale_rotation();
    test_matrix4_transform(); // and transform_inverse
    test_intersection2d_line_aabb();
    test_intersection3d_line_aabb();
    test_quaternion_functions();
    test_utility_functions();
    test_euler_functions();


    printf("\n--- All Tests Finished ---\n");
    printf("Total tests run: %d\n", g_tests_run);
    printf("Tests passed:    \033[32m%d\033[0m\n", g_tests_passed);
    printf("Tests failed:    \033[31m%d\033[0m\n", g_tests_run - g_tests_passed);

    return (g_tests_run == g_tests_passed) ? 0 : 1; // Return 0 on success, 1 on failure
}


// --- TEST FUNCTION IMPLEMENTATIONS ---

// COLOR tests
void test_rgba_functions() {
    test_suite_start("GM_RGBA Functions");

    test_start("gm_rgba_to_u32 and gm_rgba_from_u32 Roundtrip");
    GM_RGBA original = GM_RGBALit(123, 200, 50, 255);
    u32 u32_color = original.hex;
    GM_RGBA converted = gm_rgba_from_u32(u32_color);
    ASSERT_RGBA_EQUAL(converted, original);
    test_end();

    test_start("gm_rgba_alpha_blend - Opaque on Opaque");
    GM_RGBA front = GM_RGBALit(255, 0, 0, 255); // Red
    GM_RGBA back = GM_RGBALit(0, 0, 255, 255);  // Blue
    GM_RGBA blended = gm_rgba_alpha_blend(front, back);
    ASSERT_RGBA_EQUAL(blended, GM_RGBALit(255, 0, 0, 255)); // Should be red (front is opaque)
    test_end();

    test_start("gm_rgba_alpha_blend - Semi-transparent on Opaque");
    front = GM_RGBALit(255, 0, 0, 128); // 50% Red
    back = GM_RGBALit(0, 0, 255, 255);  // Blue
    blended = gm_rgba_alpha_blend(front, back);
    ASSERT_RGBA_EQUAL(blended, GM_RGBALit(128, 0, 126, 255));
    test_end();
    
    test_start("gm_rgba_multiply");
    GM_RGBA color_orig = GM_RGBALit(100, 200, 50, 255);
    GM_RGBA multiplied = gm_rgba_multiply(color_orig, 0.5f);
    ASSERT_RGBA_EQUAL(multiplied, GM_RGBALit(50, 100, 25, 127));
    multiplied = gm_rgba_multiply(color_orig, 2.0f);
    ASSERT_RGBA_EQUAL(multiplied, GM_RGBALit(200, 255, 100, 255));
    test_end();

    test_suite_end();
}

// VECTOR tests
void test_vec2_functions() {
    test_suite_start("GM_Vec2 Functions");

    test_start("gm_vec2_create");
    GM_Vec2 v = gm_vec2_create(1.0f, 2.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    test_end();

    test_start("gm_vec2_add");
    GM_Vec2 v1 = GM_Vec2Lit(1, 2);
    GM_Vec2 v2 = GM_Vec2Lit(3, 4);
    GM_Vec2 sum = gm_vec2_add(v1, v2);
    ASSERT_VEC2_NEAR(sum, GM_Vec2Lit(4, 6));
    test_end();

    test_start("gm_vec2_sub");
    GM_Vec2 diff = gm_vec2_sub(v2, v1);
    ASSERT_VEC2_NEAR(diff, GM_Vec2Lit(2, 2));
    test_end();

    test_start("gm_vec2_scale");
    GM_Vec2 scaled = gm_vec2_scale(v1, 5.0f);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2Lit(5, 10));
    test_end();

    test_start("gm_vec2_scale_non_uniform");
    GM_Vec2 s_nu = GM_Vec2Lit(2.0f, 3.0f);
    scaled = gm_vec2_scale_non_uniform(v1, s_nu);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2Lit(2.0f, 6.0f));
    test_end();

    test_start("gm_vec2_scale_xyz"); // (Note: Should be _xy for Vec2, but testing as declared)
    scaled = gm_vec2_scale_xyz(v1, 2.0f, 3.0f);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2Lit(2.0f, 6.0f));
    test_end();

    test_start("gm_vec2_lerp");
    GM_Vec2 a = GM_Vec2Lit(0, 0);
    GM_Vec2 b = GM_Vec2Lit(10, 10);
    GM_Vec2 lerp_res = gm_vec2_lerp(a, b, 0.5f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2Lit(5, 5));
    lerp_res = gm_vec2_lerp(a, b, 0.0f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2Lit(0, 0));
    lerp_res = gm_vec2_lerp(a, b, 1.0f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2Lit(10, 10));
    test_end();

    test_start("gm_vec2_projection");
    GM_Vec2 proj_a = GM_Vec2Lit(1, 1);
    GM_Vec2 proj_b = GM_Vec2Lit(5, 0);
    GM_Vec2 proj_res = gm_vec2_projection(proj_a, proj_b);
    ASSERT_VEC2_NEAR(proj_res, GM_Vec2Lit(1, 0));
    test_end();

    test_start("gm_vec2_normalize");
    GM_Vec2 norm_v = GM_Vec2Lit(3, 4); // Magnitude 5
    GM_Vec2 normalized = gm_vec2_normalize(norm_v);
    ASSERT_VEC2_NEAR(normalized, GM_Vec2Lit(0.6f, 0.8f));
    test_end();

    test_start("gm_vec2_distance");
    GM_Vec2 d1 = GM_Vec2Lit(0, 0);
    GM_Vec2 d2 = GM_Vec2Lit(3, 4);
    ASSERT_NEAR(gm_vec2_distance(d1, d2), 5.0f);
    test_end();

    test_start("gm_vec2_magnitude and gm_vec2_magnitude_squared");
    GM_Vec2 mag_v = GM_Vec2Lit(3, -4);
    ASSERT_NEAR(gm_vec2_magnitude(mag_v), 5.0f);
    ASSERT_NEAR(gm_vec2_magnitude_squared(mag_v), 25.0f);
    test_end();

    test_start("gm_vec2_dot");
    GM_Vec2 dot_v1 = GM_Vec2Lit(1, 2);
    GM_Vec2 dot_v2 = GM_Vec2Lit(3, 4);
    ASSERT_NEAR(gm_vec2_dot(dot_v1, dot_v2), 1.0f * 3.0f + 2.0f * 4.0f); // 3 + 8 = 11
    test_end();

    // gm_vec2_spline_point (Requires a full spline implementation, skip for now or mock)
    // test_start("gm_vec2_spline_point");
    // // ...
    // test_end();

    test_suite_end();
}

void test_vec3_functions() {
    test_suite_start("GM_Vec3 Functions");

    test_start("gm_vec3_create");
    GM_Vec3 v = gm_vec3_create(1.0f, 2.0f, 3.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    ASSERT_NEAR(v.z, 3.0f);
    test_end();

    test_start("gm_vec3_add");
    GM_Vec3 v1 = GM_Vec3Lit(1, 2, 3);
    GM_Vec3 v2 = GM_Vec3Lit(4, 5, 6);
    GM_Vec3 sum = gm_vec3_add(v1, v2);
    ASSERT_VEC3_NEAR(sum, GM_Vec3Lit(5, 7, 9));
    test_end();

    test_start("gm_vec3_sub");
    GM_Vec3 diff = gm_vec3_sub(v2, v1);
    ASSERT_VEC3_NEAR(diff, GM_Vec3Lit(3, 3, 3));
    test_end();

    test_start("gm_vec3_scale");
    GM_Vec3 scaled = gm_vec3_scale(v1, 2.0f);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3Lit(2, 4, 6));
    test_end();

    test_start("gm_vec3_scale_non_uniform");
    GM_Vec3 s_nu = GM_Vec3Lit(2.0f, 3.0f, 0.5f);
    scaled = gm_vec3_scale_non_uniform(v1, s_nu);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3Lit(2.0f, 6.0f, 1.5f));
    test_end();

    test_start("gm_vec3_scale_xyz");
    scaled = gm_vec3_scale_xyz(v1, 2.0f, 3.0f, 0.5f);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3Lit(2.0f, 6.0f, 1.5f));
    test_end();

    test_start("gm_vec3_lerp");
    GM_Vec3 a = GM_Vec3Lit(0, 0, 0);
    GM_Vec3 b = GM_Vec3Lit(10, 10, 10);
    GM_Vec3 lerp_res = gm_vec3_lerp(a, b, 0.5f);
    ASSERT_VEC3_NEAR(lerp_res, GM_Vec3Lit(5, 5, 5));
    test_end();

    test_start("gm_vec3_projection");
    GM_Vec3 proj_a = GM_Vec3Lit(1, 1, 1);
    GM_Vec3 proj_b = GM_Vec3Lit(5, 0, 0);
    GM_Vec3 proj_res = gm_vec3_projection(proj_a, proj_b);
    ASSERT_VEC3_NEAR(proj_res, GM_Vec3Lit(1, 0, 0));
    test_end();

    test_start("gm_vec3_normalize");
    GM_Vec3 norm_v = GM_Vec3Lit(3, 4, 0); // Magnitude 5
    GM_Vec3 normalized = gm_vec3_normalize(norm_v);
    ASSERT_VEC3_NEAR(normalized, GM_Vec3Lit(0.6f, 0.8f, 0.0f));
    GM_Vec3 zero_v = GM_Vec3Lit(0,0,0);
    normalized = gm_vec3_normalize(zero_v);
    ASSERT_VEC3_NEAR(normalized, GM_Vec3Lit(0,0,0)); // Should return zero vector for zero input
    test_end();

    test_start("gm_vec3_distance");
    GM_Vec3 d1 = GM_Vec3Lit(0, 0, 0);
    GM_Vec3 d2 = GM_Vec3Lit(3, 0, 4);
    ASSERT_NEAR(gm_vec3_distance(d1, d2), 5.0f);
    test_end();

    test_start("gm_vec3_magnitude and gm_vec3_magnitude_squared");
    GM_Vec3 mag_v = GM_Vec3Lit(2, -2, 1);
    ASSERT_NEAR(gm_vec3_magnitude(mag_v), 3.0f);
    ASSERT_NEAR(gm_vec3_magnitude_squared(mag_v), 9.0f);
    test_end();

    test_start("gm_vec3_dot");
    GM_Vec3 dot_v1 = GM_Vec3Lit(1, 2, 3);
    GM_Vec3 dot_v2 = GM_Vec3Lit(4, -5, 6);
    ASSERT_NEAR(gm_vec3_dot(dot_v1, dot_v2), 1.0f * 4.0f + 2.0f * -5.0f + 3.0f * 6.0f); // 4 - 10 + 18 = 12
    test_end();

    test_start("gm_vec3_cross");
    GM_Vec3 cross_v1 = GM_Vec3Lit(1, 0, 0); // X-axis
    GM_Vec3 cross_v2 = GM_Vec3Lit(0, 1, 0); // Y-axis
    GM_Vec3 cross_res = gm_vec3_cross(cross_v1, cross_v2); // X x Y = Z
    ASSERT_VEC3_NEAR(cross_res, GM_Vec3Lit(0, 0, 1));

    cross_v1 = GM_Vec3Lit(1, 2, 3);
    cross_v2 = GM_Vec3Lit(4, 5, 6);
    cross_res = gm_vec3_cross(cross_v1, cross_v2);
    // (2*6 - 3*5) = 12 - 15 = -3
    // (3*4 - 1*6) = 12 - 6 = 6
    // (1*5 - 2*4) = 5 - 8 = -3
    ASSERT_VEC3_NEAR(cross_res, GM_Vec3Lit(-3, 6, -3));
    test_end();

    test_suite_end();
}

void test_vec4_functions() {
    test_suite_start("GM_Vec4 Functions");

    test_start("gm_vec4_create");
    GM_Vec4 v = gm_vec4_create(1.0f, 2.0f, 3.0f, 4.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    ASSERT_NEAR(v.z, 3.0f);
    ASSERT_NEAR(v.w, 4.0f);
    test_end();

    test_start("gm_vec4_add");
    GM_Vec4 v1 = GM_Vec4Lit(1, 2, 3, 4);
    GM_Vec4 v2 = GM_Vec4Lit(5, 6, 7, 8);
    GM_Vec4 sum = gm_vec4_add(v1, v2);
    ASSERT_VEC4_NEAR(sum, GM_Vec4Lit(6, 8, 10, 12));
    test_end();

    test_start("gm_vec4_sub");
    GM_Vec4 diff = gm_vec4_sub(v2, v1);
    ASSERT_VEC4_NEAR(diff, GM_Vec4Lit(4, 4, 4, 4));
    test_end();

    test_start("gm_vec4_scale");
    GM_Vec4 scaled = gm_vec4_scale(v1, 2.0f);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4Lit(2, 4, 6, 8));
    test_end();

    test_start("gm_vec4_scale_non_uniform");
    GM_Vec4 s_nu = GM_Vec4Lit(2.0f, 3.0f, 0.5f, 1.0f);
    scaled = gm_vec4_scale_non_uniform(v1, s_nu);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4Lit(2.0f, 6.0f, 1.5f, 4.0f));
    test_end();

    test_start("gm_vec4_scale_xyz"); // (Note: Should be _xyzw for Vec4, but testing as declared)
    scaled = gm_vec4_scale_xyz(v1, 2.0f, 3.0f, 0.5f, 1.0f);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4Lit(2.0f, 6.0f, 1.5f, 4.0f));
    test_end();

    test_start("gm_vec4_lerp");
    GM_Vec4 a = GM_Vec4Lit(0, 0, 0, 0);
    GM_Vec4 b = GM_Vec4Lit(10, 10, 10, 10);
    GM_Vec4 lerp_res = gm_vec4_lerp(a, b, 0.5f);
    ASSERT_VEC4_NEAR(lerp_res, GM_Vec4Lit(5, 5, 5, 5));
    test_end();

    test_start("gm_vec4_normalize");
    GM_Vec4 norm_v = GM_Vec4Lit(1, 2, 2, 0); // Magnitude 3
    GM_Vec4 normalized = gm_vec4_normalize(norm_v);
    ASSERT_VEC4_NEAR(normalized, GM_Vec4Lit(1.0f/3.0f, 2.0f/3.0f, 2.0f/3.0f, 0.0f));
    test_end();

    test_start("gm_vec4_distance");
    GM_Vec4 d1 = GM_Vec4Lit(0, 0, 0, 0);
    GM_Vec4 d2 = GM_Vec4Lit(1, 2, 2, 0);
    ASSERT_NEAR(gm_vec4_distance(d1, d2), 3.0f);
    test_end();

    test_start("gm_vec4_magnitude and gm_vec4_magnitude_squared");
    GM_Vec4 mag_v = GM_Vec4Lit(1, -2, 2, 0);
    ASSERT_NEAR(gm_vec4_magnitude(mag_v), 3.0f);
    ASSERT_NEAR(gm_vec4_magnitude_squared(mag_v), 9.0f);
    test_end();

    test_start("gm_vec4_dot");
    GM_Vec4 dot_v1 = GM_Vec4Lit(1, 2, 3, 4);
    GM_Vec4 dot_v2 = GM_Vec4Lit(5, -6, 7, -8);
    ASSERT_NEAR(gm_vec4_dot(dot_v1, dot_v2), 1.0f*5.0f + 2.0f*-6.0f + 3.0f*7.0f + 4.0f*-8.0f); // 5 - 12 + 21 - 32 = -18
    test_end();

    test_suite_end();
}

// MATRIX tests
// NOTE: These tests assume your GM_Matrix4 operations are implemented for a ROW-MAJOR system
// as you stated: "i'm using a row major system".
// This means the matrix data[i] should be interpreted as m[row * 4 + col]
// and multiplications will be done with row vectors (v' = v * M).

void test_matrix4_identity() {
    test_suite_start("GM_Matrix4 Identity");

    test_start("gm_mat4_identity produces correct identity matrix");
    GM_Matrix4 identity = gm_mat4_identity();
    GM_Matrix4 expected = {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    }; // Row-major representation of identity
    ASSERT_MATRIX4_NEAR(identity, expected);
    test_end();

    test_suite_end();
}

void test_matrix4_multiply() {
    test_suite_start("GM_Matrix4 Multiply");

    test_start("gm_mat4_mult with identity");
    GM_Matrix4 A = {
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 10, 11, 12,
        13, 14, 15, 16
    }; // Arbitrary row-major matrix
    GM_Matrix4 identity = gm_mat4_identity();
    GM_Matrix4 result = gm_mat4_mult(A, identity); // A * I = A
    ASSERT_MATRIX4_NEAR(result, A);
    result = gm_mat4_mult(identity, A); // I * A = A
    ASSERT_MATRIX4_NEAR(result, A);
    test_end();

    test_start("gm_mat4_mult with simple matrices (row-major)");
    // Example matrices (simple, for easy manual calculation)
    // A = [[1, 2], [3, 4]] (extend to 4x4 with identity)
    // B = [[5, 6], [7, 8]] (extend to 4x4 with identity)
    // A * B = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]]
    //       = [[5+14, 6+16], [15+28, 18+32]]
    //       = [[19, 22], [43, 50]]

    GM_Matrix4 matA = gm_mat4_identity();
    matA.data[0] = 1; matA.data[1] = 2; // Row 0
    matA.data[4] = 3; matA.data[5] = 4; // Row 1

    GM_Matrix4 matB = gm_mat4_identity();
    matB.data[0] = 5; matB.data[1] = 6; // Row 0
    matB.data[4] = 7; matB.data[5] = 8; // Row 1

    GM_Matrix4 matC = gm_mat4_mult(matA, matB);

    GM_Matrix4 expectedC = gm_mat4_identity();
    expectedC.data[0] = 19; expectedC.data[1] = 22;
    expectedC.data[4] = 43; expectedC.data[5] = 50;

    ASSERT_MATRIX4_NEAR(matC, expectedC);
    test_end();

    test_suite_end();
}

void test_matrix4_transpose() {
    test_suite_start("GM_Matrix4 Transpose");

    test_start("gm_mat4_transpose produces correct transpose");
    GM_Matrix4 original = {
        1,  2,  3,  4,  // Row 0
        5,  6,  7,  8,  // Row 1
        9, 10, 11, 12, // Row 2
        13, 14, 15, 16  // Row 3
    }; // Row-major representation

    GM_Matrix4 transposed = gm_mat4_transpose(original);

    // Expected transpose (will be column-major if original was row-major,
    // but stored back into a row-major `GM_Matrix4` struct)
    GM_Matrix4 expected = {
        1,  5,  9, 13, // First column of original becomes first row
        2,  6, 10, 14, // Second column of original becomes second row
        3,  7, 11, 15, // Third column of original becomes third row
        4,  8, 12, 16  // Fourth column of original becomes fourth row
    };
    ASSERT_MATRIX4_NEAR(transposed, expected);
    test_end();

    test_start("Double transpose returns original");
    GM_Matrix4 double_transposed = gm_mat4_transpose(transposed);
    ASSERT_MATRIX4_NEAR(double_transposed, original);
    test_end();

    test_suite_end();
}

void test_matrix4_inverse() {
    test_suite_start("GM_Matrix4 Inverse");

    test_start("gm_mat4_inverse - Identity matrix");
    GM_Matrix4 identity = gm_mat4_identity();
    bool success;
    GM_Matrix4 inverse_identity = gm_mat4_inverse(identity, &success);
    ASSERT_TRUE(success);
    ASSERT_MATRIX4_NEAR(inverse_identity, identity);
    test_end();

    test_start("gm_mat4_inverse - Simple invertible matrix");
    // Example from https://math.stackexchange.com/questions/112613/how-to-find-the-inverse-of-a-4x4-matrix
    // This example is in row-major
    GM_Matrix4 m = {
        1.0f, 0.0f, 2.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    GM_Matrix4 expected_inv = {
        1.0f, 0.0f, -2.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    GM_Matrix4 inv_m = gm_mat4_inverse(m, &success);
    ASSERT_TRUE(success);
    ASSERT_MATRIX4_NEAR(inv_m, expected_inv);
    test_end();

    test_start("gm_mat4_inverse - Test M * inv(M) = Identity");
    // A slightly more complex, but still simple invertible matrix
    GM_Matrix4 test_mat = {
        2.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 2.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    }; // Determinant = 2*2 - 1*1 = 3
    GM_Matrix4 test_mat_inv = gm_mat4_inverse(test_mat, &success);
    ASSERT_TRUE(success);
    GM_Matrix4 product = gm_mat4_mult(test_mat, test_mat_inv);
    ASSERT_MATRIX4_NEAR(product, gm_mat4_identity());
    test_end();

    test_start("gm_mat4_inverse - Singular matrix (no inverse)");
    GM_Matrix4 singular_matrix = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, // All zeros row makes it singular
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    GM_Matrix4 inv_singular = gm_mat4_inverse(singular_matrix, &success);
    ASSERT_FALSE(success); // Should report failure
    // When singular, the function should return identity or some default error matrix.
    // The specific return value when `success` is false might vary, but for this test,
    // we'll assume it returns identity as a safe fallback.
    ASSERT_MATRIX4_NEAR(inv_singular, gm_mat4_identity());
    test_end();

    test_suite_end();
}

void test_matrix4_translation_scale_rotation() {
    test_suite_start("GM_Matrix4 Translation, Scale, Rotation");

    test_start("gm_mat4_translate_xyz");
    GM_Matrix4 mat_t = gm_mat4_identity();
    mat_t = gm_mat4_translate_xyz(mat_t, 5.0f, -2.0f, 3.0f);
    GM_Matrix4 expected_t = gm_mat4_identity();
    expected_t.v[0].w = 5.0f; // m03
    expected_t.v[1].w = -2.0f; // m13
    expected_t.v[2].w = 3.0f; // m23
    ASSERT_MATRIX4_NEAR(mat_t, expected_t);
    test_end();

    test_start("gm_mat4_scale_xyz");
    GM_Matrix4 mat_s = gm_mat4_identity();
    mat_s = gm_mat4_scale_xyz(mat_s, 2.0f, 3.0f, 0.5f);
    GM_Matrix4 expected_s = gm_mat4_identity();
    expected_s.data[0] = 2.0f; // m00
    expected_s.data[5] = 3.0f; // m11
    expected_s.data[10] = 0.5f; // m22
    ASSERT_MATRIX4_NEAR(mat_s, expected_s);
    test_end();

    test_start("gm_mat4_rotate_xyz - z-axis 90 deg (Row-Major)");
    GM_Matrix4 mat_r = gm_mat4_identity();
    mat_r = gm_mat4_rotate_xyz(mat_r, 90.0f, 0.0f, 0.0f, 1.0f);
    GM_Matrix4 expected_r = gm_mat4_identity();
    expected_r.data[0] = 0.0f;  expected_r.data[1] = -1.0f; // cos(90), -sin(90)
    expected_r.data[4] = 1.0f; expected_r.data[5] = 0.0f;   // sin(90), cos(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_start("gm_mat4_rotate_xyz - X-axis 90 deg (Row-Major)");
    mat_r = gm_mat4_identity();
    mat_r = gm_mat4_rotate_xyz(mat_r, 90.0f, 1.0f, 0.0f, 0.0f);
    
    // Row-major rotation matrix structure for X-rotation:
    // [ 1   0       0     0 ]
    // [ 0 cos(a)  sin(a)  0 ]
    // [ 0 -sin(a) cos(a)  0 ]
    // [ 0   0       0     1 ]
    expected_r = gm_mat4_identity();
    expected_r.data[5] = 0.0f;  expected_r.data[6] = -1.0f; // cos(90), -sin(90)
    expected_r.data[9] = 1.0f; expected_r.data[10] = 0.0f;  // sin(90), cos(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_start("gm_mat4_rotate_xyz - Y-axis 90 deg (Row-Major)");
    mat_r = gm_mat4_identity();
    mat_r = gm_mat4_rotate_xyz(mat_r, 90.0f, 0.0f, 1.0f, 0.0f);
    
    // Row-major rotation matrix structure for Y-rotation:
    // [ cos(a) 0 -sin(a) 0 ]
    // [   0    1   0     0 ]
    // [ sin(a) 0 cos(a)  0 ]
    // [   0    0   0     1 ]
    expected_r = gm_mat4_identity();
    expected_r.data[0] = 0.0f;  expected_r.data[2] = 1.0f; // cos(90),  sin(90)
    expected_r.data[8] = -1.0f; expected_r.data[10] = 0.0f;  // -sin(90), cos(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_suite_end();
}

void test_matrix4_transform() {
    test_suite_start("GM_Matrix4 Transform");

    test_start("gm_mat4_transform - Simple Scale, Rotate, Translate");
    GM_Vec3 s = GM_Vec3Lit(2.0f, 2.0f, 2.0f);
    float angle = 90.0f; // 90 degree Z rotation
    GM_Vec3 axis = GM_Vec3Lit(0.0f, 0.0f, 1.0f);
    GM_Vec3 t = GM_Vec3Lit(10.0f, 0.0f, 0.0f);

    GM_Matrix4 transform_mat = gm_mat4_transform(s, angle, axis, t);
    GM_Matrix4 expected_transform = {
        0.0f, -2.0f, 0.0f, 10.0f,
        2.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 2.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    ASSERT_MATRIX4_NEAR(transform_mat, expected_transform);
    test_end();

    test_start("gm_mat4_transform_inverse - M * inv(M) = Identity");
    GM_Vec3 s_inv = GM_Vec3Lit(1.0f, 1.0f, 1.0f);
    float angle_inv = 45.0f; // Arbitrary angle
    GM_Vec3 axis_inv = GM_Vec3Lit(0.0f, 1.0f, 0.0f); // Y-axis rotation
    GM_Vec3 t_inv = GM_Vec3Lit(1.0f, 2.0f, 3.0f);

    GM_Matrix4 M = gm_mat4_transform(s_inv, angle_inv, axis_inv, t_inv);
    GM_Matrix4 M_inv = gm_mat4_inverse_transform(s_inv, angle_inv, axis_inv, t_inv);
    
    GM_Matrix4 product = gm_mat4_mult(M, M_inv);
    ASSERT_MATRIX4_NEAR(product, gm_mat4_identity());
    test_end();

    test_suite_end();
}

void test_shapes3d_creation() {
    test_suite_start("GM_Shapes 3D Functions");

    test_start("gm_rectangle3d_create");
    GM_Rectangle3D r3d = gm_rectangle3d_create(1, 2, 3, 4, 5, 6);
    ASSERT_VEC3_NEAR(r3d.center, GM_Vec3Lit(1, 2, 3));
    ASSERT_NEAR(r3d.width, 4.0f);
    ASSERT_NEAR(r3d.height, 5.0f);
    ASSERT_NEAR(r3d.length, 6.0f);
    test_end();

    test_start("gm_rectangle_reference3d_create");
    GM_Vec3 center_ref = GM_Vec3Lit(7, 8, 9);
    GM_RectangleReference3D r_ref3d = gm_rectangle_reference3d_create(&center_ref, 10, 11, 12);
    ASSERT_VEC3_NEAR((*r_ref3d.center), GM_Vec3Lit(7, 8, 9));
    ASSERT_NEAR(r_ref3d.width, 10.0f);
    ASSERT_NEAR(r_ref3d.height, 11.0f);
    ASSERT_NEAR(r_ref3d.length, 12.0f);
    test_end();

    test_start("gm_circle3d_create");
    GM_Circle3D c3d = gm_circle3d_create(1.0f, 2.0f, 3.0f, 5.0f);
    ASSERT_VEC3_NEAR(c3d.center, GM_Vec3Lit(1.0f, 2.0f, 3.0f));
    ASSERT_NEAR(c3d.radius, 5.0f);
    test_end();

    test_start("gm_circle_reference3d_create");
    GM_Vec3 circle_center_ref = GM_Vec3Lit(11.0f, 22.0f, 33.0f);
    GM_CircleReference3D c_ref3d = gm_circle_reference3d_create(&circle_center_ref, 15.0f);
    ASSERT_VEC3_NEAR((*c_ref3d.center), GM_Vec3Lit(11.0f, 22.0f, 33.0f));
    ASSERT_NEAR(c_ref3d.radius, 15.0f);
    test_end();

    test_suite_end();
}

// INTERSECTION tests
void test_intersection2d_line_line() {
    test_suite_start("GM_Intersection 2D Line-Line");

    test_start("gm_intersection2d_line_line - Intersecting");
    GM_Vec2 A = GM_Vec2Lit(0, 0);
    GM_Vec2 B = GM_Vec2Lit(2, 2);
    GM_Vec2 C = GM_Vec2Lit(0, 2);
    GM_Vec2 D = GM_Vec2Lit(2, 0);
    GM_Vec2 intersection_point;
    ASSERT_TRUE(gm_intersection2d_line_line(A, B, C, D, &intersection_point));
    ASSERT_VEC2_NEAR(intersection_point, GM_Vec2Lit(1, 1));
    test_end();

    test_start("gm_intersection2d_line_line - Parallel");
    A = GM_Vec2Lit(0, 0);
    B = GM_Vec2Lit(1, 0);
    C = GM_Vec2Lit(0, 1);
    D = GM_Vec2Lit(1, 1);
    ASSERT_FALSE(gm_intersection2d_line_line(A, B, C, D, NULL));
    test_end();

    test_start("gm_intersection2d_line_line - Collinear (no overlap)");
    A = GM_Vec2Lit(0, 0);
    B = GM_Vec2Lit(1, 0);
    C = GM_Vec2Lit(2, 0);
    D = GM_Vec2Lit(3, 0);
    ASSERT_FALSE(gm_intersection2d_line_line(A, B, C, D, NULL));
    test_end();

    test_start("gm_intersection2d_line_line - Collinear (overlap)");
    A = GM_Vec2Lit(0, 0);
    B = GM_Vec2Lit(2, 0);
    C = GM_Vec2Lit(1, 0);
    D = GM_Vec2Lit(3, 0);
    // This case might return true or false depending on exact implementation,
    // often false for segment-segment intersection to avoid infinite solutions.
    // Given that your function uses a determinant-based approach (den),
    // it's likely to return false for collinear segments due to `den` being near zero.
    ASSERT_FALSE(gm_intersection2d_line_line(A, B, C, D, NULL));
    test_end();

    test_suite_end();
}

void test_intersection2d_line_aabb() {
    test_suite_start("GM_Intersection 2D Line-AABB");

    test_start("gm_intersection2d_line_aabb - Intersecting with specific example");
    GM_Vec2 p0 = GM_Vec2Lit(0, 0);
    GM_Vec2 p1 = GM_Vec2Lit(10, 10);
    GM_Vec2 aabb_position = GM_Vec2Lit(1.5f, 3.6f);
    GM_RectangleReference2D aabb = gm_rectangle_reference2d_create(&aabb_position, 3, 3);
    GM_Vec2 inPoint, outPoint;

    // Expected values from previous logic:
    // AABB: x_min = 0, x_max = 3; y_min = 2.1, y_max = 5.1
    // Line y=x
    // Entry at y=2.1 => (2.1, 2.1)
    // Exit at x=3 => (3.0, 3.0)
    bool intersected = gm_intersection2d_line_aabb(p0, p1, aabb, &inPoint, &outPoint);
    ASSERT_TRUE(intersected);
    ASSERT_VEC2_NEAR(inPoint, GM_Vec2Lit(2.1f, 2.1f));
    ASSERT_VEC2_NEAR(outPoint, GM_Vec2Lit(3.0f, 3.0f));
    test_end();

    test_start("gm_intersection2d_line_aabb - No intersection");
    p0 = GM_Vec2Lit(0, 0);
    p1 = GM_Vec2Lit(1, 1);
    aabb_position = GM_Vec2Lit(10, 10);
    aabb = gm_rectangle_reference2d_create(&aabb_position, 2, 2);
    GM_Vec2 dummy_in, dummy_out;
    ASSERT_FALSE(gm_intersection2d_line_aabb(p0, p1, aabb, &dummy_in, &dummy_out));
    test_end();

    test_start("gm_intersection2d_line_aabb - Line inside AABB entirely");
    p0 = GM_Vec2Lit(1, 1);
    p1 = GM_Vec2Lit(2, 2);
    aabb_position = GM_Vec2Lit(0, 0);
    aabb = gm_rectangle_reference2d_create(&aabb_position, 10, 10); // AABB from -5 to 5
    // Should intersect, entry point is p0, exit is p1
    intersected = gm_intersection2d_line_aabb(p0, p1, aabb, &inPoint, &outPoint);
    ASSERT_TRUE(intersected);
    ASSERT_VEC2_NEAR(inPoint, p0);
    ASSERT_VEC2_NEAR(outPoint, p1);
    test_end();

    test_suite_end();
}

void test_intersection3d_line_aabb() {
    test_suite_start("GM_Intersection 3D Line-AABB");

    test_start("gm_intersection3d_line_aabb - Intersecting");
    GM_Vec3 p0_3d = GM_Vec3Lit(0, 0, 0);
    GM_Vec3 p1_3d = GM_Vec3Lit(5, 5, 5);
    GM_Vec3 aabb_pos_3d = GM_Vec3Lit(1.0f, 1.0f, 1.0f);
    GM_RectangleReference3D aabb_3d = gm_rectangle_reference3d_create(&aabb_pos_3d, 1.0f, 1.0f, 1.0f); // AABB from (0.5,0.5,0.5) to (1.5,1.5,1.5)
    GM_Vec3 in_3d, out_3d;

    ASSERT_TRUE(gm_intersection3d_line_aabb(p0_3d, p1_3d, aabb_3d, &in_3d, &out_3d));
    ASSERT_VEC3_NEAR(in_3d, GM_Vec3Lit(0.5f, 0.5f, 0.5f));
    ASSERT_VEC3_NEAR(out_3d, GM_Vec3Lit(1.5f, 1.5f, 1.5f));
    test_end();

    test_start("gm_intersection3d_line_aabb - No Intersection");
    p0_3d = GM_Vec3Lit(0, 0, 0);
    p1_3d = GM_Vec3Lit(0.1f, 0.1f, 0.1f);
    aabb_pos_3d = GM_Vec3Lit(10, 10, 10);
    aabb_3d = gm_rectangle_reference3d_create(&aabb_pos_3d, 1, 1, 1);
    GM_Vec3 dummy_in_3d, dummy_out_3d;
    ASSERT_FALSE(gm_intersection3d_line_aabb(p0_3d, p1_3d, aabb_3d, &dummy_in_3d, &dummy_out_3d));
    test_end();
    
    test_start("gm_intersection3d_line_aabb - Line inside AABB entirely");
    p0_3d = GM_Vec3Lit(1, 1, 1);
    p1_3d = GM_Vec3Lit(2, 2, 2);
    aabb_pos_3d = GM_Vec3Lit(0, 0, 0);
    aabb_3d = gm_rectangle_reference3d_create(&aabb_pos_3d, 10, 10, 10); // AABB from -5 to 5
    // Should intersect, entry point is p0, exit is p1
    bool intersected = gm_intersection3d_line_aabb(p0_3d, p1_3d, aabb_3d, &in_3d, &out_3d);
    ASSERT_TRUE(intersected);
    ASSERT_VEC3_NEAR(in_3d, p0_3d);
    ASSERT_VEC3_NEAR(out_3d, p1_3d);
    test_end();

    test_suite_end();
}

// QUATERNION tests
void test_quaternion_functions() {
    test_suite_start("GM_Quaternion Functions");

    test_start("gm_quat_create (Identity)");
    GM_Quaternion q_identity = gm_quat_create(0.0f, GM_Vec3Lit(1.0f, 0.0f, 0.0f)); // 0-degree rotation
    ASSERT_QUAT_NEAR(q_identity, GM_QuaternionLit(1.0f, 0.0f, 0.0f, 0.0f));
    test_end();

    test_start("gm_quat_create (X-axis 180 deg)");
    GM_Quaternion q_x180 = gm_quat_create(180.0f, GM_Vec3Lit(1.0f, 0.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_x180, GM_QuaternionLit(0.0f, 1.0f, 0.0f, 0.0f));
    test_end();

    test_start("gm_quat_create (Z-axis 90 deg)");
    GM_Quaternion q_z90 = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 0.0f, 1.0f));
    // cos(45) = sqrt(2)/2 = 0.7071
    // sin(45) = sqrt(2)/2 = 0.7071
    ASSERT_QUAT_NEAR(q_z90, GM_QuaternionLit(0.7071f, 0.0f, 0.0f, 0.7071f));
    test_end();

    test_start("gm_quat_mult");
    GM_Quaternion q_x90 = gm_quat_create(90.0f, GM_Vec3Lit(1.0f, 0.0f, 0.0f)); // w=0.7071, x=0.7071, y=0, z=0
    GM_Quaternion q_y90 = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 1.0f, 0.0f)); // w=0.7071, x=0, y=0.7071, z=0

    GM_Quaternion q_result = gm_quat_mult(q_y90, q_x90);
    ASSERT_QUAT_NEAR(q_result, GM_QuaternionLit(0.5f, 0.5f, 0.5f, -0.5f));
    test_end();

    test_start("gm_quat_to_mat4 (Identity)");
    GM_Matrix4 mat_from_identity_quat = gm_quat_to_mat4(gm_quat_create(0.0f, GM_Vec3Lit(1.0f, 0.0f, 0.0f)));
    ASSERT_MATRIX4_NEAR(mat_from_identity_quat, gm_mat4_identity());
    test_end();

    test_start("gm_quat_to_mat4 (Z-axis 90 deg)");
    GM_Matrix4 mat_from_z90_quat = gm_quat_to_mat4(gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 0.0f, 1.0f)));
    GM_Matrix4 expected_mat_z90 = gm_mat4_identity();
    expected_mat_z90.data[0] = 0.0f;  expected_mat_z90.data[1] = -1.0f;
    expected_mat_z90.data[4] = 1.0f; expected_mat_z90.data[5] = 0.0f;
    ASSERT_MATRIX4_NEAR(mat_from_z90_quat, expected_mat_z90);
    test_end();

    test_start("gm_quat_slerp");
    GM_Quaternion q_from = GM_QuaternionLit(1.0f, 0.0f, 0.0f, 0.0f); // Identity
    GM_Quaternion q_to = gm_quat_create(180.0f, GM_Vec3Lit(0.0f, 1.0f, 0.0f)); // 180 deg around Y (0, 0, 1, 0)
    GM_Quaternion slerp_half = gm_quat_slerp(q_from, q_to, 0.5f);
    // Halfway between identity and Y-180 should be Y-90
    GM_Quaternion expected_y90 = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 1.0f, 0.0f));
    ASSERT_QUAT_NEAR(slerp_half, expected_y90);
    test_end();

    test_start("gm_quat_transform_vec3");
    GM_Vec3 vec_to_rotate = GM_Vec3Lit(1.0f, 0.0f, 0.0f); // X-axis
    GM_Quaternion rotate_around_z = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 0.0f, 1.0f));
    GM_Vec3 rotated_vec = gm_quat_vector_mult(rotate_around_z, vec_to_rotate);
    ASSERT_VEC3_NEAR(rotated_vec, GM_Vec3Lit(0.0f, 1.0f, 0.0f)); // X rotated 90 deg around Z becomes Y
    test_end();

    test_suite_end();
}

// UTILITY tests
void test_utility_functions() {
    test_suite_start("GM_Utility Functions");

    test_start("gm_utility_clamp");
    ASSERT_NEAR(CLAMP(5.0f, 0.0f, 10.0f), 5.0f);
    ASSERT_NEAR(CLAMP(-1.0f, 0.0f, 10.0f), 0.0f);
    ASSERT_NEAR(CLAMP(11.0f, 0.0f, 10.0f), 10.0f);
    test_end();

    test_start("degrees to radians");
    ASSERT_NEAR(DEGREES_TO_RAD(0.0f), 0.0f);
    ASSERT_NEAR(DEGREES_TO_RAD(180.0f), PI);
    ASSERT_NEAR(DEGREES_TO_RAD(90.0f), PI / 2.0f);
    test_end();

    test_start("gm_utility_radians_to_degrees");
    ASSERT_NEAR(RAD_TO_DEGREES(0.0f), 0.0f);
    ASSERT_NEAR(RAD_TO_DEGREES(PI), 180.0f);
    ASSERT_NEAR(RAD_TO_DEGREES(PI / 2.0f), 90.0f);
    test_end();

    test_start("gm_utility_lerp_float");
    ASSERT_NEAR(gm_lerp(0.0f, 10.0f, 0.5f), 5.0f);
    ASSERT_NEAR(gm_lerp(0.0f, 10.0f, 0.0f), 0.0f);
    ASSERT_NEAR(gm_lerp(0.0f, 10.0f, 1.0f), 10.0f);
    ASSERT_NEAR(gm_lerp(10.0f, 0.0f, 0.5f), 5.0f);
    test_end();

    test_start("gm_utility_is_near_zero");
    ASSERT_TRUE(NEAR_ZERO(0.0f));
    ASSERT_TRUE(NEAR_ZERO(EPSILON / 2.0f));
    ASSERT_FALSE(NEAR_ZERO(EPSILON * 2.0f));
    test_end();

    test_suite_end();
}

// EULER tests
void test_euler_functions() {
    test_suite_start("GM_Euler Functions");

    test_start("gm_euler_to_quat (Identity)");
    GM_Vec3 euler_zero = GM_Vec3Lit(0.0f, 0.0f, 0.0f); // Roll, Pitch, Yaw
    GM_Quaternion q_from_euler_zero = gm_quat_from_euler(euler_zero);
    ASSERT_QUAT_NEAR(q_from_euler_zero, GM_QuaternionLit(1.0f, 0.0f, 0.0f, 0.0f));
    test_end();

    test_start("gm_euler_to_quat (90 deg Yaw)");
    GM_Vec3 euler_yaw90 = GM_Vec3Lit(0.0f, 0.0f, 90.0f); // Yaw (Z-axis rotation)
    GM_Quaternion q_from_euler_yaw90 = gm_quat_from_euler(euler_yaw90);
    GM_Quaternion expected_z90 = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 0.0f, 1.0f));
    ASSERT_QUAT_NEAR(q_from_euler_yaw90, expected_z90);
    test_end();

    test_start("gm_euler_to_quat (90 deg Pitch)");
    GM_Vec3 euler_pitch90 = GM_Vec3Lit(0.0f, 90.0f, 0.0f); // Pitch (Y-axis rotation)
    GM_Quaternion q_from_euler_pitch90 = gm_quat_from_euler(euler_pitch90);
    GM_Quaternion expected_y90 = gm_quat_create(90.0f, GM_Vec3Lit(0.0f, 1.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_from_euler_pitch90, expected_y90);
    test_end();

    test_start("gm_euler_to_quat (90 deg Roll)");
    GM_Vec3 euler_roll90 = GM_Vec3Lit(90.0f, 0.0f, 0.0f); // Roll (X-axis rotation)
    GM_Quaternion q_from_euler_roll90 = gm_quat_from_euler(euler_roll90);
    GM_Quaternion expected_x90 = gm_quat_create(90.0f, GM_Vec3Lit(1.0f, 0.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_from_euler_roll90, expected_x90);
    test_end();

    test_suite_end();
}