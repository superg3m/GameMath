#include <gm.hpp>

// --- Global Test Counters ---
static int g_tests_run = 0;
static int g_tests_passed = 0;
static const char* g_current_test_name = NULL;
static const char* g_current_suite_name = NULL;

// --- Helper Functions for Assertions ---
bool float_equals_near(float a, float b) {
    return abs(a - b) <= EPSILON;
}

/*
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
*/

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
    do { if (actual != expected) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f), got (%.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, actual.x, actual.y); \
        return; \
    } } while (0)

#define ASSERT_VEC3_NEAR(actual, expected) \
    do { if (actual != expected) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f, %.4f), got (%.4f, %.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, expected.z, actual.x, actual.y, actual.z); \
        return; \
    } } while (0)

#define ASSERT_VEC4_NEAR(actual, expected) \
    do { if (actual != expected) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Expected (%.4f, %.4f, %.4f, %.4f), got (%.4f, %.4f, %.4f, %.4f)\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__, expected.x, expected.y, expected.z, expected.w, actual.x, actual.y, actual.z, actual.w); \
        return; \
    } } while (0)

#define ASSERT_MATRIX4_NEAR(actual, expected) \
    do { if (actual != expected) { \
        printf("    \033[31mFAIL\033[0m: %s.%s at %s:%d - Matrix mismatch.\n", g_current_suite_name, g_current_test_name, __FILE__, __LINE__); \
        return; \
    } } while (0)

#define ASSERT_QUAT_NEAR(actual, expected) \
    do { if (actual != expected) { \
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

// --- Main Test Runner ---
int main() {
    printf("Starting all GM Math Unit Tests...\n");

    // Call all test suites
    test_vec2_functions();
    test_vec3_functions();
    test_vec4_functions();
    test_matrix4_multiply();
    test_matrix4_transpose();
    test_matrix4_inverse();
    test_matrix4_translation_scale_rotation();
    test_matrix4_transform(); // and transform_inverse
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

// VECTOR tests
void test_vec2_functions() {
    test_suite_start("GM_Vec2 Functions");

    test_start("gm_vec2_create");
    GM_Vec2 v = GM_Vec2(1.0f, 2.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    test_end();

    test_start("gm_vec2_add");
    GM_Vec2 v1 = GM_Vec2(1, 2);
    GM_Vec2 v2 = GM_Vec2(3, 4);
    GM_Vec2 sum = v1 + v2;
    ASSERT_VEC2_NEAR(sum, GM_Vec2(4, 6));
    test_end();

    test_start("gm_vec2_sub");
    GM_Vec2 diff = v2 - v1;
    ASSERT_VEC2_NEAR(diff, GM_Vec2(2, 2));
    test_end();

    test_start("gm_vec2_scale");
    GM_Vec2 scaled = v1.scale(5.0f);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2(5, 10));
    test_end();

    test_start("gm_vec2_scale_non_uniform");
    GM_Vec2 s_nu = GM_Vec2(2.0f, 3.0f);
    scaled = v1.scale(s_nu.x, s_nu.y);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2(2.0f, 6.0f));
    test_end();

    test_start("gm_vec2_scale_xyz"); // (Note: Should be _xy for Vec2, but testing as declared)
    scaled = v1.scale(2.0f, 3.0f);
    ASSERT_VEC2_NEAR(scaled, GM_Vec2(2.0f, 6.0f));
    test_end();

    test_start("gm_vec2_lerp");
    GM_Vec2 a = GM_Vec2(0, 0);
    GM_Vec2 b = GM_Vec2(10, 10);
    GM_Vec2 lerp_res = GM_Vec2::lerp(a, b, 0.5f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2(5, 5));
    lerp_res = GM_Vec2::lerp(a, b, 0.0f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2(0, 0));
    lerp_res = GM_Vec2::lerp(a, b, 1.0f);
    ASSERT_VEC2_NEAR(lerp_res, GM_Vec2(10, 10));
    test_end();

    test_start("gm_vec2_normalize");
    GM_Vec2 norm_v = GM_Vec2(3, 4); // Magnitude 5
    GM_Vec2 normalized = norm_v.normalize();
    ASSERT_VEC2_NEAR(normalized, GM_Vec2(0.6f, 0.8f));
    test_end();

    test_start("gm_vec2_distanfce");
    GM_Vec2 d1 = GM_Vec2(0, 0);
    GM_Vec2 d2 = GM_Vec2(3, 4);
    ASSERT_NEAR(GM_Vec2::distance(d1, d2), 5.0f);
    test_end();

    test_start("gm_vec2_magnitude and gm_vec2_magnitude_squared");
    GM_Vec2 mag_v = GM_Vec2(3, -4);
    ASSERT_NEAR(mag_v.magnitude(), 5.0f);
    ASSERT_NEAR(mag_v.magnitudeSquared(), 25.0f);
    test_end();

    test_start("gm_vec2_dot");
    GM_Vec2 dot_v1 = GM_Vec2(2, 2);
    GM_Vec2 dot_v2 = GM_Vec2(3, 4);
    ASSERT_NEAR(GM_Vec2::dot(dot_v1, dot_v2), 2.0f * 3.0f + 2.0f * 4.0f); // 6 + 8 = 14
    test_end();

    test_suite_end();
}

void test_vec3_functions() {
    test_suite_start("GM_Vec3 Functions");

    test_start("gm_vec3_create");
    GM_Vec3 v = GM_Vec3(1.0f, 2.0f, 3.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    ASSERT_NEAR(v.z, 3.0f);
    test_end();

    test_start("gm_vec3_add");
    GM_Vec3 v1 = GM_Vec3(1, 2, 3);
    GM_Vec3 v2 = GM_Vec3(4, 5, 6);
    GM_Vec3 sum = v1 + v2;
    ASSERT_VEC3_NEAR(sum, GM_Vec3(5, 7, 9));
    test_end();

    test_start("gm_vec3_sub");
    GM_Vec3 diff = v2 - v1;
    ASSERT_VEC3_NEAR(diff, GM_Vec3(3, 3, 3));
    test_end();

    test_start("gm_vec3_scale");
    GM_Vec3 scaled = v1.scale(2.0f);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3(2, 4, 6));
    test_end();

    test_start("gm_vec3_scale_non_uniform");
    GM_Vec3 s_nu = GM_Vec3(2.0f, 3.0f, 0.5f);
    scaled = v1.scale(s_nu);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3(2.0f, 6.0f, 1.5f));
    test_end();

    test_start("gm_vec3_scale_xyz");
    scaled = v1.scale(2.0f, 3.0f, 0.5f);
    ASSERT_VEC3_NEAR(scaled, GM_Vec3(2.0f, 6.0f, 1.5f));
    test_end();

    test_start("gm_vec3_lerp");
    GM_Vec3 a = GM_Vec3(0, 0, 0);
    GM_Vec3 b = GM_Vec3(10, 10, 10);
    GM_Vec3 lerp_res = GM_Vec3::lerp(a, b, 0.5f);
    ASSERT_VEC3_NEAR(lerp_res, GM_Vec3(5, 5, 5));
    test_end();

    test_start("gm_vec3_normalize");
    GM_Vec3 norm_v = GM_Vec3(3, 4, 0); // Magnitude 5
    GM_Vec3 normalized = norm_v.normalize();
    ASSERT_VEC3_NEAR(normalized, GM_Vec3(0.6f, 0.8f, 0.0f));
    GM_Vec3 zero_v = GM_Vec3(0,0,0);
    normalized = zero_v.normalize();
    ASSERT_VEC3_NEAR(normalized, GM_Vec3(0, 0, 0)); // Should return zero vector for zero input
    test_end();

    test_start("gm_vec3_distanfce");
    GM_Vec3 d1 = GM_Vec3(0, 0, 0);
    GM_Vec3 d2 = GM_Vec3(3, 0, 4);
    ASSERT_NEAR(GM_Vec3::distance(d1, d2), 5.0f);
    test_end();

    test_start("gm_vec3_magnitude and gm_vec3_magnitude_squared");
    GM_Vec3 mag_v = GM_Vec3(2, -2, 1);
    ASSERT_NEAR(mag_v.magnitude(), 3.0f);
    ASSERT_NEAR(mag_v.magnitudeSquared(), 9.0f);
    test_end();

    test_start("gm_vec3_dot");
    GM_Vec3 dot_v1 = GM_Vec3(1, 2, 3);
    GM_Vec3 dot_v2 = GM_Vec3(4, -5, 6);
    ASSERT_NEAR(GM_Vec3::dot(dot_v1, dot_v2), 1.0f * 4.0f + 2.0f * -5.0f + 3.0f * 6.0f); // 4 - 10 + 18 = 12
    test_end();

    test_start("gm_vec3_cross");
    GM_Vec3 cross_v1 = GM_Vec3(1, 0, 0); // X-axis
    GM_Vec3 cross_v2 = GM_Vec3(0, 1, 0); // Y-axis
    GM_Vec3 cross_res = GM_Vec3::cross(cross_v1, cross_v2); // X x Y = Z
    ASSERT_VEC3_NEAR(cross_res, GM_Vec3(0, 0, 1));

    cross_v1 = GM_Vec3(1, 2, 3);
    cross_v2 = GM_Vec3(4, 5, 6);
    cross_res = GM_Vec3::cross(cross_v1, cross_v2);
    // (2*6 - 3*5) = 12 - 15 = -3
    // (3*4 - 1*6) = 12 - 6 = 6
    // (1*5 - 2*4) = 5 - 8 = -3
    ASSERT_VEC3_NEAR(cross_res, GM_Vec3(-3, 6, -3));
    test_end();

    test_suite_end();
}

void test_vec4_functions() {
    test_suite_start("GM_Vec4 Functions");

    test_start("GM_Vec4");
    GM_Vec4 v = GM_Vec4(1.0f, 2.0f, 3.0f, 4.0f);
    ASSERT_NEAR(v.x, 1.0f);
    ASSERT_NEAR(v.y, 2.0f);
    ASSERT_NEAR(v.z, 3.0f);
    ASSERT_NEAR(v.w, 4.0f);
    test_end();

    test_start("gm_vec4_add");
    GM_Vec4 v1 = GM_Vec4(1, 2, 3, 4);
    GM_Vec4 v2 = GM_Vec4(5, 6, 7, 8);
    GM_Vec4 sum = v1 + v2;
    ASSERT_VEC4_NEAR(sum, GM_Vec4(6, 8, 10, 12));
    test_end();

    test_start("gm_vec4_sub");
    GM_Vec4 diff = v2 - v1;
    ASSERT_VEC4_NEAR(diff, GM_Vec4(4, 4, 4, 4));
    test_end();

    test_start("gm_vec4_scale");
    GM_Vec4 scaled = v1.scale(2.0f);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4(2, 4, 6, 8));
    test_end();

    test_start("gm_vec4_scale_non_uniform");
    GM_Vec4 s_nu = GM_Vec4(2.0f, 3.0f, 0.5f, 1.0f);
    scaled = v1.scale(s_nu);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4(2.0f, 6.0f, 1.5f, 4.0f));
    test_end();

    test_start("gm_vec4_scale_xyz"); // (Note: Should be _xyzw for Vec4, but testing as declared)
    scaled = v1.scale(2.0f, 3.0f, 0.5f, 1.0f);
    ASSERT_VEC4_NEAR(scaled, GM_Vec4(2.0f, 6.0f, 1.5f, 4.0f));
    test_end();

    test_start("gm_vec4_lerp");
    GM_Vec4 a = GM_Vec4(0, 0, 0, 0);
    GM_Vec4 b = GM_Vec4(10, 10, 10, 10);
    GM_Vec4 lerp_res = GM_Vec4::lerp(a, b, 0.5f);
    ASSERT_VEC4_NEAR(lerp_res, GM_Vec4(5, 5, 5, 5));
    test_end();

    test_start("gm_vec4_normalize");
    GM_Vec4 norm_v = GM_Vec4(1, 2, 2, 0); // Magnitude 3
    GM_Vec4 normalized = norm_v.normalize();
    ASSERT_VEC4_NEAR(normalized, GM_Vec4(1.0f/3.0f, 2.0f/3.0f, 2.0f/3.0f, 0.0f));
    test_end();

    test_start("gm_vec4_distanfce");
    GM_Vec4 d1 = GM_Vec4(0, 0, 0, 0);
    GM_Vec4 d2 = GM_Vec4(1, 2, 2, 0);
    ASSERT_NEAR(GM_Vec4::distance(d1, d2), 3.0f);
    test_end();

    test_start("gm_vec4_magnitude and gm_vec4_magnitude_squared");
    GM_Vec4 mag_v = GM_Vec4(1, -2, 2, 0);
    ASSERT_NEAR(mag_v.magnitude(), 3.0f);
    ASSERT_NEAR(mag_v.magnitudeSquared(), 9.0f);
    test_end();

    test_start("gm_vec4_dot");
    GM_Vec4 dot_v1 = GM_Vec4(1, 2, 3, 4);
    GM_Vec4 dot_v2 = GM_Vec4(5, -6, 7, -8);
    ASSERT_NEAR(GM_Vec4::dot(dot_v1, dot_v2), 1.0f*5.0f + 2.0f*-6.0f + 3.0f*7.0f + 4.0f*-8.0f); // 5 - 12 + 21 - 32 = -18
    test_end();

    test_suite_end();
}

// MATRIX tests
// NOTE: These tests assume your GM_Matrix4 operations are implemented for a ROW-MAJOR system
// as you stated: "i'm usinfg a row major system".
// This means the matrix data[i] should be interpreted as m[row * 4 + col]
// and multiplications will be done with row vectors (v' = v * M).

void test_matrix4_multiply() {
    test_suite_start("GM_Matrix4 Multiply");

    test_start("gm_mat4_mult with identity");
    GM_Matrix4 A = {
        1,  2,  3,  4,
        5,  6,  7,  8,
        9,  10, 11, 12,
        13, 14, 15, 16
    }; // Arbitrary row-major matrix

    GM_Matrix4 identity = GM_Matrix4::identity();
    GM_Matrix4 result = A * identity; // A * I = A
    ASSERT_MATRIX4_NEAR(result, A);
    result = identity * A; // I * A = A
    ASSERT_MATRIX4_NEAR(result, A);
    test_end();

    test_start("gm_mat4_mult with simple matrices (row-major)");
    // Example matrices (simple, for easy manual calculation)
    // A = [[1, 2], [3, 4]] (extend to 4x4 with identity)
    // B = [[5, 6], [7, 8]] (extend to 4x4 with identity)
    // A * B = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]]
    //       = [[5+14, 6+16], [15+28, 18+32]]
    //       = [[19, 22], [43, 50]]

    GM_Matrix4 matA = GM_Matrix4::identity();
    matA.v[0].x = 1; matA.v[0].y = 2; // Row 0
    matA.v[1].x = 3; matA.v[1].y = 4; // Row 1

    GM_Matrix4 matB = GM_Matrix4::identity();
    matB.v[0].x = 5; matB.v[0].y = 6; // Row 0
    matB.v[1].x = 7; matB.v[1].y = 8; // Row 1

    GM_Matrix4 matC = matA * matB;

    GM_Matrix4 expectedC = GM_Matrix4::identity();
    expectedC.v[0].x = 19; expectedC.v[0].y = 22;
    expectedC.v[1].x = 43; expectedC.v[1].y = 50;

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
        9,  10, 11, 12, // Row 2
        13, 14, 15, 16  // Row 3
    }; // Row-major representation

    GM_Matrix4 transposed = original.transpose();

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
    GM_Matrix4 float_transposed = transposed.transpose();
    ASSERT_MATRIX4_NEAR(float_transposed, original);
    test_end();

    test_suite_end();
}

void test_matrix4_inverse() {
    test_suite_start("GM_Matrix4 Inverse");

    test_start("GM_Matrix4::inverse - Identity matrix");
    GM_Matrix4 identity = GM_Matrix4::identity();
    bool success;
    GM_Matrix4 inverse_identity = identity.inverse(success);
    ASSERT_TRUE(success);
    ASSERT_MATRIX4_NEAR(inverse_identity, identity);
    test_end();

    test_start("GM_Matrix4::inverse - Simple invertible matrix");
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
    GM_Matrix4 inv_m = m.inverse(success);
    ASSERT_TRUE(success);
    ASSERT_MATRIX4_NEAR(inv_m, expected_inv);
    test_end();

    test_start("GM_Matrix4::inverse - Test M * inv(M) = Identity");
    // A slightly more complex, but still simple invertible matrix
    GM_Matrix4 test_mat = {
        GM_Vec4(2.0f, 1.0f, 0.0f, 0.0f),
        GM_Vec4(1.0f, 2.0f, 0.0f, 0.0f),
        GM_Vec4(0.0f, 0.0f, 1.0f, 0.0f),
        GM_Vec4(0.0f, 0.0f, 0.0f, 1.0f)
    }; // Determinant = 2*2 - 1*1 = 3
    GM_Matrix4 test_mat_inv = test_mat.inverse(success);
    ASSERT_TRUE(success);
    GM_Matrix4 product = test_mat * test_mat_inv;
    ASSERT_MATRIX4_NEAR(product, GM_Matrix4::identity());
    test_end();

    test_start("GM_Matrix4::inverse - Singular matrix (no inverse)");
    GM_Matrix4 sinfgular_matrix = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 0.0f, // All zeros row makes it sinfgular
        0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    GM_Matrix4 inv_sinfgular = sinfgular_matrix.inverse(success);
    ASSERT_FALSE(success); // Should report failure
    // When sinfgular, the function should return identity or some default error matrix.
    // The specific return value when `success` is false might vary, but for this test,
    // we'll assume it returns identity as a safe fallback.
    ASSERT_MATRIX4_NEAR(inv_sinfgular, GM_Matrix4::identity());
    test_end();

    test_suite_end();
}

void test_matrix4_translation_scale_rotation() {
    test_suite_start("GM_Matrix4 Translation, Scale, Rotation");

    test_start("GM_Matrix4::translate");
    GM_Matrix4 mat_t = GM_Matrix4::identity();
    mat_t = GM_Matrix4::translate(mat_t, 5.0f, -2.0f, 3.0f);
    GM_Matrix4 expected_t = GM_Matrix4::identity();
    expected_t.v[0].w = 5.0f; // m03
    expected_t.v[1].w = -2.0f; // m13
    expected_t.v[2].w = 3.0f; // m23
    ASSERT_MATRIX4_NEAR(mat_t, expected_t);
    test_end();

    test_start("GM_Matrix4::scale");
    GM_Matrix4 mat_s = GM_Matrix4::identity();
    mat_s = GM_Matrix4::scale(mat_s, 2.0f, 3.0f, 0.5f);
    GM_Matrix4 expected_s = GM_Matrix4::identity();
    expected_s.v[0].x = 2.0f; // m00
    expected_s.v[1].y = 3.0f; // m11
    expected_s.v[2].z = 0.5f; // m22
    ASSERT_MATRIX4_NEAR(mat_s, expected_s);
    test_end();

    test_start("GM_Matrix4::rotate - z-axis 90 deg (Row-Major)");
    GM_Matrix4 mat_r = GM_Matrix4::identity();
    mat_r = GM_Matrix4::rotate(mat_r, 90.0f, 0.0f, 0.0f, 1.0f);
    GM_Matrix4 expected_r = GM_Matrix4::identity();
    expected_r.v[0].x = 0.0f;  expected_r.v[0].y = -1.0f; // cosf(90), -sinf(90)
    expected_r.v[1].x = 1.0f; expected_r.v[1].y = 0.0f;   // sinf(90), cosf(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_start("GM_Matrix4::rotate - X-axis 90 deg (Row-Major)");
    mat_r = GM_Matrix4::identity();
    mat_r = GM_Matrix4::rotate(mat_r, 90.0f, 1.0f, 0.0f, 0.0f);
    
    // Row-major rotation matrix structure for X-rotation:
    // [ 1   0       0     0 ]
    // [ 0 cosf(a)  sinf(a)  0 ]
    // [ 0 -sinf(a) cosf(a)  0 ]
    // [ 0   0       0     1 ]
    expected_r = GM_Matrix4::identity();
    expected_r.v[1].y = 0.0f; expected_r.v[1].z = -1.0f; // cosf(90), -sinf(90)
    expected_r.v[2].y = 1.0f; expected_r.v[2].z = 0.0f;  // sinf(90), cosf(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_start("GM_Matrix4::rotate - Y-axis 90 deg (Row-Major)");
    mat_r = GM_Matrix4::identity();
    mat_r = GM_Matrix4::rotate(mat_r, 90.0f, 0.0f, 1.0f, 0.0f);
    
    // Row-major rotation matrix structure for Y-rotation:
    // [ cosf(a) 0 -sinf(a) 0 ]
    // [   0    1   0     0 ]
    // [ sinf(a) 0 cosf(a)  0 ]
    // [   0    0   0     1 ]
    expected_r = GM_Matrix4::identity();
    expected_r.v[0].x = 0.0f;  expected_r.v[0].z = 1.0f; // cosf(90),  sinf(90)
    expected_r.v[2].x = -1.0f; expected_r.v[2].z = 0.0f;  // -sinf(90), cosf(90)
    
    ASSERT_MATRIX4_NEAR(mat_r, expected_r);
    test_end();

    test_suite_end();
}

void test_matrix4_transform() {
    test_suite_start("GM_Matrix4 Transform");

    test_start("GM_Matrix4::transform - Simple Scale, Rotate, Translate");
    GM_Vec3 s = GM_Vec3(2.0f, 2.0f, 2.0f);
    float angle = 90.0f; // 90 degree Z rotation
    GM_Vec3 axis = GM_Vec3(0.0f, 0.0f, 1.0f);
    GM_Vec3 t = GM_Vec3(10.0f, 0.0f, 0.0f);

    GM_Matrix4 transform_mat = GM_Matrix4::transform(s, angle, axis, t);
    GM_Matrix4 expected_transform = {
        0.0f, -2.0f, 0.0f, 10.0f,
        2.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 2.0f, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    ASSERT_MATRIX4_NEAR(transform_mat, expected_transform);
    test_end();

    test_start("gm_mat4_transform_inverse - M * inv(M) = Identity");
    GM_Vec3 s_inv = GM_Vec3(1.0f, 1.0f, 1.0f);
    float angle_inv = 45.0f; // Arbitrary angle
    GM_Vec3 axis_inv = GM_Vec3(0.0f, 1.0f, 0.0f); // Y-axis rotation
    GM_Vec3 t_inv = GM_Vec3(1.0f, 2.0f, 3.0f);

    GM_Matrix4 M = GM_Matrix4::transform(s_inv, angle_inv, axis_inv, t_inv);
    GM_Matrix4 M_inv = GM_Matrix4::inverse_transform(s_inv, angle_inv, axis_inv, t_inv);
    
    GM_Matrix4 product = M * M_inv;
    ASSERT_MATRIX4_NEAR(product, GM_Matrix4::identity());
    test_end();

    test_suite_end();
}



// QUATERNION tests
void test_quaternion_functions() {
    test_suite_start("GM_Quaternion Functions");

    test_start("gm_quat_create (Identity)");
    GM_Quaternion q_identity = GM_Quaternion(0, 1, 0, 0); // 0-degree rotation
    ASSERT_QUAT_NEAR(q_identity, GM_Quaternion::identity());
    test_end();

    test_start("gm_quat_create (X-axis 180 deg)");
    GM_Quaternion q_x180 = GM_Quaternion(180.0f, GM_Vec3(1.0f, 0.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_x180, GM_Quaternion::literal(0.0f, 1.0f, 0.0f, 0.0f));
    test_end();

    test_start("gm_quat_create (Z-axis 90 deg)");
    GM_Quaternion q_z90 = GM_Quaternion(90.0f, GM_Vec3(0.0f, 0.0f, 1.0f));
    // cosf(45) = sqrtf(2)/2 = 0.7071
    // sinf(45) = sqrtf(2)/2 = 0.7071
    ASSERT_QUAT_NEAR(q_z90, GM_Quaternion::literal(0.7071f, 0.0f, 0.0f, 0.7071f));
    test_end();

    test_start("gm_quat_mult");
    GM_Quaternion q_x90 = GM_Quaternion(90.0f, GM_Vec3(1.0f, 0.0f, 0.0f)); // w=0.7071, x=0.7071, y=0, z=0
    GM_Quaternion q_y90 = GM_Quaternion(90.0f, GM_Vec3(0.0f, 1.0f, 0.0f)); // w=0.7071, x=0, y=0.7071, z=0

    GM_Quaternion q_result = q_y90 * q_x90;
    ASSERT_QUAT_NEAR(q_result, GM_Quaternion(120.0f, GM_Vec3(1.0f, 1.0f, -1.0f)));
    test_end();

    test_start("gm_quat_to_mat4 (Identity)");
    GM_Matrix4 mat_from_identity_quat = GM_Quaternion::literal(1, 0, 0, 0).toMatrix4();
    ASSERT_MATRIX4_NEAR(mat_from_identity_quat, GM_Matrix4::identity());
    test_end();

    test_start("gm_quat_to_mat4 (Z-axis 90 deg)");
    GM_Matrix4 mat_from_z90_quat = GM_Quaternion(90.0f, GM_Vec3(0.0f, 0.0f, 1.0f)).toMatrix4();
    GM_Matrix4 expected_mat_z90 = GM_Matrix4::identity();
    expected_mat_z90.v[0].x = 0.0f;  expected_mat_z90.v[0].y = -1.0f;
    expected_mat_z90.v[1].x = 1.0f; expected_mat_z90.v[1].y = 0.0f;
    ASSERT_MATRIX4_NEAR(mat_from_z90_quat, expected_mat_z90);
    test_end();

    test_start("gm_quat_slerp");
    GM_Quaternion q_from = GM_Quaternion::identity(); // Identity
    GM_Quaternion q_to = GM_Quaternion(180.0f, GM_Vec3(0.0f, 1.0f, 0.0f)); // 180 deg around Y (0, 0, 1, 0)
    GM_Quaternion slerp_half = GM_Quaternion::slerp(q_from, q_to, 0.5f);
    // Halfway between identity and Y-180 should be Y-90
    GM_Quaternion expected_y90 = GM_Quaternion(90.0f, GM_Vec3(0.0f, 1.0f, 0.0f));
    ASSERT_QUAT_NEAR(slerp_half, expected_y90);
    test_end();

    test_start("gm_quat_transform_vec3");
    GM_Vec3 vec_to_rotate = GM_Vec3(1.0f, 0.0f, 0.0f); // X-axis
    GM_Quaternion rotate_around_z = GM_Quaternion(90.0f, GM_Vec3(0.0f, 0.0f, 1.0f));
    GM_Vec3 rotated_vec = rotate_around_z * vec_to_rotate;
    ASSERT_VEC3_NEAR(rotated_vec, GM_Vec3(0.0f, 1.0f, 0.0f));
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
    GM_Vec3 euler_zero = GM_Vec3(0.0f, 0.0f, 0.0f); // Roll, Pitch, Yaw
    GM_Quaternion q_from_euler_zero = GM_Quaternion::fromEuler(euler_zero);
    ASSERT_QUAT_NEAR(q_from_euler_zero, GM_Quaternion::literal(1.0f, 0.0f, 0.0f, 0.0f));
    test_end();

    test_start("gm_euler_to_quat (90 deg Yaw)");
    GM_Vec3 euler_yaw90 = GM_Vec3(0.0f, 0.0f, 90.0f); // Yaw (Z-axis rotation)
    GM_Quaternion q_from_euler_yaw90 = GM_Quaternion::fromEuler(euler_yaw90);
    GM_Quaternion expected_z90 = GM_Quaternion(90.0f, GM_Vec3(0.0f, 0.0f, 1.0f));
    ASSERT_QUAT_NEAR(q_from_euler_yaw90, expected_z90);
    test_end();

    test_start("gm_euler_to_quat (90 deg Pitch)");
    GM_Vec3 euler_pitch90 = GM_Vec3(0.0f, 90.0f, 0.0f); // Pitch (Y-axis rotation)
    GM_Quaternion q_from_euler_pitch90 = GM_Quaternion::fromEuler(euler_pitch90);
    GM_Quaternion expected_y90 = GM_Quaternion(90.0f, GM_Vec3(0.0f, 1.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_from_euler_pitch90, expected_y90);
    test_end();

    test_start("gm_euler_to_quat (90 deg Roll)");
    GM_Vec3 euler_roll90 = GM_Vec3(90.0f, 0.0f, 0.0f); // Roll (X-axis rotation)
    GM_Quaternion q_from_euler_roll90 = GM_Quaternion::fromEuler(euler_roll90);
    GM_Quaternion expected_x90 = GM_Quaternion(90.0f, GM_Vec3(1.0f, 0.0f, 0.0f));
    ASSERT_QUAT_NEAR(q_from_euler_roll90, expected_x90);
    test_end();

    test_suite_end();
}