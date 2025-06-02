#include <gm.h>

void test_intersection() {
    GM_Vec2 p0 = {0, 0};
    GM_Vec2 p1 = {10, 10};
    GM_Rectangle2D aabb = gm_rectangle2d_create(2, 2, 3, 3);
    GM_Vec2 inPoint, outPoint;

    if (gm_intersection2d_line_aabb(p0, p1, aabb, &inPoint, &outPoint)) {
        printf("Intersection occurred!\n");
        printf("Entry point: (%f, %f)\n", inPoint.x, inPoint.y);
        printf("Exit point:  (%f, %f)\n", outPoint.x, outPoint.y);
    } else {
        printf("No intersection.\n");
    }
}

int main() {
    GM_Matrix4 A = {
        .data = {
            1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12,
            13, 14, 15, 16
        }
    };

    GM_Matrix4 B = {
        .data = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            1, 1, 1, 1
        }
    };

    GM_Matrix4 expected = {
        .data = {
            5, 6, 7, 4,
            13, 14, 15, 8,
            21, 22, 23, 12,
            29, 30, 31, 16
        }
    };

    GM_Matrix4 result = gm_mat4_mult(A, B);

    if (gm_mat4_equal(&expected, &result)) {
        printf("Test Passed\n");
    } else {
        printf("Test Failed\n");
    }

    result = gm_mat4_transpose(B);
    GM_Matrix4 expected2 = {
        .data = {
            1, 0, 0, 1,
            0, 1, 0, 1,
            0, 0, 1, 1,
            0, 0, 0, 1
        }
    };



    if (gm_mat4_equal(&expected2, &result)) {
        printf("Test Passed\n");
    } else {
        printf("Test Failed\n");
    }

    test_intersection();

    return 0;
}