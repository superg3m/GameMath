#include <gm.h>

int mat_equal(GM_Matrix4* A, GM_Matrix4* B) {
    for (int i = 0; i < 16; i++) {
        if (fabsf(A->data[i] - B->data[i]) > 1e-5f) {
            return 0;
        }
    }
    return 1;
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

    if (mat_equal(&expected, &result)) {
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



    if (mat_equal(&expected2, &result)) {
        printf("Test Passed\n");
    } else {
        printf("Test Failed\n");
    }

    return 0;
}