#include "./gm.hpp"

GM_Vec2::GM_Vec2(float fill) {
    this->x = fill;
    this->y = fill;
}

GM_Vec2::GM_Vec2(float x, float y) {
    this->x = x;
    this->y = y;
}

float GM_Vec2::magnitude() {
    return sqrtf((this->x * this->x) + (this->y * this->y));
}

float GM_Vec2::magnitudeSquared() {
    return (this->x * this->x) + (this->y * this->y);
}

GM_Vec2 GM_Vec2::normalize() {
    GM_Vec2 ret(0, 0);
    const float magnitude = this->magnitude();
    if (magnitude == 0) {
        return GM_Vec2(0,0);
    }

    ret.x = this->x / magnitude;
    ret.y = this->y / magnitude;

    return ret;
}

GM_Vec2 GM_Vec2::scale(float scale) {
    return GM_Vec2(this->x * scale, this->y * scale);
}

GM_Vec2 GM_Vec2::scale(GM_Vec2 s) {
    return GM_Vec2(this->x * s.x, this->y * s.y);
}

GM_Vec2 GM_Vec2::scale(float scale_x, float scale_y) {
    return GM_Vec2(this->x * scale_x, this->y * scale_y);
}

float GM_Vec2::dot(GM_Vec2 a, GM_Vec2 b) {
    return (a.x * b.x) + (a.y * b.y);
}

float GM_Vec2::distanfce(GM_Vec2 a, GM_Vec2 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y));
}

float GM_Vec2::distanfceSquared(GM_Vec2 a, GM_Vec2 b) {
    return SQUARED(b.x - a.x) + SQUARED(b.y - a.y);
}

GM_Vec2 GM_Vec2::lerp(GM_Vec2 a, GM_Vec2 b, float t) {
    GM_Vec2 ab = (b - a);
    return a + (ab.scale(t));
}

GM_Vec2 GM_Vec2::operator+(const GM_Vec2 &right) {
    return GM_Vec2(this->x + right.x, this->y + right.y);
}
GM_Vec2& GM_Vec2::operator+=(const GM_Vec2 &right) {
    this->x += right.x;
    this->y += right.y;

    return *this;
}

GM_Vec2 GM_Vec2::operator-(const GM_Vec2 &right) {
    return GM_Vec2(this->x - right.x, this->y - right.y);
}
GM_Vec2& GM_Vec2::operator-=(const GM_Vec2 &right) {
    this->x -= right.x;
    this->y -= right.y;

    return *this;
}

GM_Vec2 GM_Vec2::operator*(const GM_Vec2 &right) {
    return GM_Vec2(this->x * right.x, this->y * right.y);
}

GM_Vec2& GM_Vec2::operator*=(const GM_Vec2 &right) {
    this->x *= right.x;
    this->y *= right.y;

    return *this;
}

GM_Vec2 GM_Vec2::operator/(const GM_Vec2 &right) {
    return GM_Vec2(this->x / right.x, this->y / right.y);
}
GM_Vec2& GM_Vec2::operator/=(const GM_Vec2 &right) {
    this->x /= right.x;
    this->y /= right.y;

    return *this;
}

bool GM_Vec2::operator==(const GM_Vec2 &right) {
    return NEAR_ZERO(this->x - right.x) && NEAR_ZERO(this->y - right.y);
}

bool GM_Vec2::operator!=(const GM_Vec2 &right) {
    return !(*this == right);
}

GM_Vec3::GM_Vec3(float fill) {
    this->x = fill;
    this->y = fill;
    this->z = fill;
}

GM_Vec3::GM_Vec3(float x, float y, float z) {
    this->x = x;
    this->y = y;
    this->z = z;
}

float GM_Vec3::magnitude() {
    return sqrtf(SQUARED(this->x) + SQUARED(this->y) + SQUARED(this->z));
}

float GM_Vec3::magnitudeSquared() {
    return SQUARED(this->x) + SQUARED(this->y) + SQUARED(this->z);
}

GM_Vec3 GM_Vec3::normalize() {
    GM_Vec3 ret(0, 0, 0);
    const float magnitude = this->magnitude();
    if (magnitude == 0) {
        return GM_Vec3(0, 0, 0);
    }

    ret.x = this->x / magnitude;
    ret.y = this->y / magnitude;
    ret.z = this->z / magnitude;

    return ret;
}

GM_Vec3 GM_Vec3::scale(float scale) {
    return GM_Vec3(this->x * scale, this->y * scale, this->z * scale);
}

GM_Vec3 GM_Vec3::scale(GM_Vec3 s) {
    return GM_Vec3(this->x * s.x, this->y * s.y, this->z * s.z);
}

GM_Vec3 GM_Vec3::scale(float scale_x, float scale_y, float scale_z) {
    return GM_Vec3(this->x * scale_x, this->y * scale_y, this->z * scale_z);
}


float GM_Vec3::distanfce(GM_Vec3 a, GM_Vec3 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z));
}

float GM_Vec3::distanfceSquared(GM_Vec3 a, GM_Vec3 b) {
    return SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z);
}

float GM_Vec3::dot(GM_Vec3 a, GM_Vec3 b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

GM_Vec3 GM_Vec3::lerp(GM_Vec3 a, GM_Vec3 b, float t) {
    GM_Vec3 ab = (b - a);
    return a + (ab.scale(t));
}

GM_Vec3 GM_Vec3::cross(GM_Vec3 A, GM_Vec3 B) {
    GM_Vec3 ret(0, 0, 0);
    ret.x = A.y * B.z - A.z * B.y;
    ret.y = A.z * B.x - A.x * B.z;
    ret.z = A.x * B.y - A.y * B.x;

    return ret;
}

GM_Vec3 GM_Vec3::operator+(const GM_Vec3 &right) {
    return GM_Vec3(this->x + right.x, this->y + right.y, this->z + right.z);
}
GM_Vec3& GM_Vec3::operator+=(const GM_Vec3 &right) {
    this->x += right.x;
    this->y += right.y;
    this->z += right.z;

    return *this;
}

GM_Vec3 GM_Vec3::operator-(const GM_Vec3 &right) {
    return GM_Vec3(this->x - right.x, this->y - right.y, this->z - right.z);
}
GM_Vec3& GM_Vec3::operator-=(const GM_Vec3 &right) {
    this->x -= right.x;
    this->y -= right.y;
    this->z -= right.z;

    return *this;
}


GM_Vec3 GM_Vec3::operator*(const GM_Vec3 &right) {
    return GM_Vec3(this->x * right.x, this->y * right.y, this->z * right.z);
}
GM_Vec3& GM_Vec3::operator*=(const GM_Vec3 &right) {
    this->x *= right.x;
    this->y *= right.y;
    this->z *= right.z;

    return *this;
}


GM_Vec3 GM_Vec3::operator/(const GM_Vec3 &right) {
    return GM_Vec3(this->x / right.x, this->y / right.y, this->z / right.z);
}
GM_Vec3& GM_Vec3::operator/=(const GM_Vec3 &right) {
    this->x /= right.x;
    this->y /= right.y;
    this->z /= right.z;

    return *this;
}

bool GM_Vec3::operator==(const GM_Vec3 &right) {
    return NEAR_ZERO(this->x - right.x) && NEAR_ZERO(this->y - right.y) && NEAR_ZERO(this->z - right.z);
}
bool GM_Vec3::operator!=(const GM_Vec3 &right) {
    return !(*this == right);
}

GM_Vec4::GM_Vec4(float fill) {
    this->x = fill;
    this->y = fill;
    this->z = fill;
    this->w = fill;
}

GM_Vec4::GM_Vec4(float x, float y, float z, float w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}

float GM_Vec4::magnitude() {
    return sqrtf(SQUARED(this->x) + SQUARED(this->y) + SQUARED(this->z) + SQUARED(this->w));
}

float GM_Vec4::magnitudeSquared() {
    return SQUARED(this->x) + SQUARED(this->y) + SQUARED(this->z) + SQUARED(this->w);
}

GM_Vec4 GM_Vec4::normalize() {
    GM_Vec4 ret(0, 0, 0, 0);
    const float magnitude = this->magnitude();
    if (magnitude == 0) {
        return GM_Vec4(0, 0, 0, 0);
    }

    ret.x = this->x / magnitude;
    ret.y = this->y / magnitude;
    ret.z = this->z / magnitude;
    ret.w = this->w / magnitude;

    return ret;
}

GM_Vec4 GM_Vec4::scale(float scale) {
    return GM_Vec4(this->x * scale, this->y * scale, this->z * scale, this->w * scale);
}

GM_Vec4 GM_Vec4::scale(GM_Vec4 s) {
    return GM_Vec4(this->x * s.x, this->y * s.y, this->z * s.z, this->w * s.w);
}

GM_Vec4 GM_Vec4::scale(float scale_x, float scale_y, float scale_z, float scale_w) {
    return GM_Vec4(this->x * scale_x, this->y * scale_y, this->z * scale_z, this->w * scale_w);
}

float GM_Vec4::distanfce(GM_Vec4 a, GM_Vec4 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z) + SQUARED(b.w - a.w));
}

float GM_Vec4::distanfceSquared(GM_Vec4 a, GM_Vec4 b) {
    return SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z) + SQUARED(b.w - a.w);
}

float GM_Vec4::dot(GM_Vec4 a, GM_Vec4 b) {
    return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

GM_Vec4 GM_Vec4::lerp(GM_Vec4 a, GM_Vec4 b, float t) {
    GM_Vec4 ab = (b - a);
    return a + (ab.scale(t));
}


GM_Vec4 GM_Vec4::operator+(const GM_Vec4 &right) {
    return GM_Vec4(this->x + right.x, this->y + right.y, this->z + right.z, this->w + right.w);
}
GM_Vec4& GM_Vec4::operator+=(const GM_Vec4 &right) {
    this->x += right.x;
    this->y += right.y;
    this->z += right.z;
    this->w += right.w;

    return *this;
}

GM_Vec4 GM_Vec4::operator-(const GM_Vec4 &right) {
    return GM_Vec4(this->x - right.x, this->y - right.y, this->z - right.z, this->w - right.w);
}
GM_Vec4& GM_Vec4::operator-=(const GM_Vec4 &right) {
    this->x -= right.x;
    this->y -= right.y;
    this->z -= right.z;
    this->w -= right.w;

    return *this;
}


GM_Vec4 GM_Vec4::operator*(const GM_Vec4 &right) {
    return GM_Vec4(this->x * right.x, this->y * right.y, this->z * right.z, this->w * right.w);
}
GM_Vec4& GM_Vec4::operator*=(const GM_Vec4 &right) {
    this->x *= right.x;
    this->y *= right.y;
    this->z *= right.z;
    this->w *= right.w;

    return *this;
}


GM_Vec4 GM_Vec4::operator/(const GM_Vec4 &right) {
    return GM_Vec4(this->x / right.x, this->y / right.y, this->z / right.z, this->w / right.w);
}
GM_Vec4& GM_Vec4::operator/=(const GM_Vec4 &right) {
    this->x /= right.x;
    this->y /= right.y;
    this->z /= right.z;
    this->w /= right.w;

    return *this;
}

bool GM_Vec4::operator==(const GM_Vec4 &right) {
    return NEAR_ZERO(this->x - right.x) && NEAR_ZERO(this->y - right.y) && NEAR_ZERO(this->z - right.z) && NEAR_ZERO(this->w - right.w);
}
bool GM_Vec4::operator!=(const GM_Vec4 &right) {
    return !(*this == right);
}

GM_Matrix4::GM_Matrix4() {
    v[0] = GM_Vec4(0, 0, 0, 0);
    v[1] = GM_Vec4(0, 0, 0, 0);
    v[2] = GM_Vec4(0, 0, 0, 0);
    v[3] = GM_Vec4(0, 0, 0, 0);
}

GM_Matrix4::GM_Matrix4(GM_Vec4 r0, GM_Vec4 r1, GM_Vec4 r2, GM_Vec4 r3) {
    this->v = {
        r0, 
        r1, 
        r2, 
        r3 
    };
}

GM_Matrix4::GM_Matrix4(float m00, float m01, float m02, float m03,
        float m10, float m11, float m12, float m13,
        float m20, float m21, float m22, float m23,
        float m30, float m31, float m32, float m33) {
    v[0] = GM_Vec4(m00, m01, m02, m03);
    v[1] = GM_Vec4(m10, m11, m12, m13);
    v[2] = GM_Vec4(m20, m21, m22, m23);
    v[3] = GM_Vec4(m30, m31, m32, m33);
}

GM_Matrix4 GM_Matrix4::transpose() {
    GM_Matrix4 ret;

    ret.v[0].x = this->v[0].x;
    ret.v[0].y = this->v[1].x;
    ret.v[0].z = this->v[2].x;
    ret.v[0].w = this->v[3].x;

    ret.v[1].x = this->v[0].y;
    ret.v[1].y = this->v[1].y;
    ret.v[1].z = this->v[2].y;
    ret.v[1].w = this->v[3].y;

    ret.v[2].x = this->v[0].z;
    ret.v[2].y = this->v[1].z;
    ret.v[2].z = this->v[2].z;
    ret.v[2].w = this->v[3].z;

    ret.v[3].x = this->v[0].w;
    ret.v[3].y = this->v[1].w;
    ret.v[3].z = this->v[2].w;
    ret.v[3].w = this->v[3].w;      

    return ret;
}

GM_Matrix4 GM_Matrix4::identity() {
    GM_Matrix4 ret;
    ret.v = {
        GM_Vec4(1, 0, 0, 0),
        GM_Vec4(0, 1, 0, 0),
        GM_Vec4(0, 0, 1, 0),
        GM_Vec4(0, 0, 0, 1)
    };

    return ret;
}

GM_Matrix4 GM_Matrix4::scale(GM_Matrix4 mat, float scale) {
    return GM_Matrix4::scale(mat, GM_Vec3(scale, scale, scale));
}

GM_Matrix4 GM_Matrix4::scale(GM_Matrix4 mat, GM_Vec3 s) {
    GM_Matrix4 scale_matrix;
    scale_matrix.v = {
        GM_Vec4(s.x,  0.0f, 0.0f, 0.0f),
        GM_Vec4(0.0f, s.y,  0.0f, 0.0f),
        GM_Vec4(0.0f, 0.0f, s.z,  0.0f),
        GM_Vec4(0.0f, 0.0f, 0.0f, 1.0f) 
    };

    return scale_matrix * mat;
}

GM_Matrix4 GM_Matrix4::scale(GM_Matrix4 mat, float scale_x, float scale_y, float scale_z) {
    return GM_Matrix4::scale(mat, GM_Vec3(scale_x, scale_y, scale_z));
}


GM_Matrix4 GM_Matrix4::rotate(GM_Matrix4 mat, float theta, GM_Vec3 axis) {
    float rad = DEGREES_TO_RAD(theta);
    float c = cosf(rad);
    float s = sinf(rad);
    float t = 1.0f - c;

    axis = axis.normalize();
    float x = axis.x;
    float y = axis.y;
    float z = axis.z;

    GM_Matrix4 rot;
    rot.v = {
        GM_Vec4(t * x * x + c,     t * x * y - z * s, t * x * z + y * s, 0.0f),
        GM_Vec4(t * x * y + z * s, t * y * y + c,     t * y * z - x * s, 0.0f),
        GM_Vec4(t * x * z - y * s, t * y * z + x * s, t * z * z + c,     0.0f),
        GM_Vec4(0.0f,              0.0f,              0.0f,              1.0f)
    };

    return rot * mat;
}

GM_Matrix4 GM_Matrix4::rotate(GM_Matrix4 mat, float theta, float rot_x, float rot_y, float rot_z) {
    return GM_Matrix4::rotate(mat, theta, GM_Vec3(rot_x, rot_y, rot_z));
}

GM_Matrix4 GM_Matrix4::translate(GM_Matrix4 mat, GM_Vec3 t) {
    GM_Matrix4 translate_matrix;
    translate_matrix.v = {
        GM_Vec4(1.0f, 0.0f, 0.0f, t.x),
        GM_Vec4(0.0f, 1.0f, 0.0f, t.y),
        GM_Vec4(0.0f, 0.0f, 1.0f, t.z),
        GM_Vec4(0.0f, 0.0f, 0.0f, 1.0f)
    };

    return translate_matrix * mat;
}

GM_Matrix4 GM_Matrix4::translate(GM_Matrix4 mat, float x, float y, float z) {
    return GM_Matrix4::translate(mat, GM_Vec3(x, y, z));
}

GM_Matrix4 GM_Matrix4::transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t) {
    GM_Matrix4 scale_matrix = GM_Matrix4::scale(GM_Matrix4::identity(), s);
    GM_Matrix4 rotation_matrix = GM_Matrix4::rotate(GM_Matrix4::identity(), theta, axis);
    GM_Matrix4 translation_matrix = GM_Matrix4::translate(GM_Matrix4::identity(), t);

    return translation_matrix * rotation_matrix * scale_matrix;
}

GM_Matrix4 GM_Matrix4::inverse_transform(GM_Vec3 s, float theta, GM_Vec3 axis, GM_Vec3 t) {
    GM_Matrix4 inverse_scale_matrix = GM_Matrix4::scale(GM_Matrix4::identity(), s.scale(1 / s.x, 1 / s.y, 1 / s.z));
    GM_Matrix4 inverse_rotation_matrix = GM_Matrix4::rotate(GM_Matrix4::identity(), theta, axis).transpose();
    GM_Matrix4 inverse_translation_matrix = GM_Matrix4::translate(GM_Matrix4::identity(), t.scale(-1));

    return inverse_scale_matrix * inverse_rotation_matrix * inverse_translation_matrix;
}

GM_Matrix4 GM_Matrix4::perspective(float fov_degrees, float aspect, float near_plane, float far_plane) {
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

    GM_Matrix4 ret;
    ret.v = {
        GM_Vec4(A,  0,  0,  0),
        GM_Vec4(0,  B,  0,  0),
        GM_Vec4(0,  0,  C,  D),
        GM_Vec4(0,  0, -1,  0)
    };

    return ret;
}

// Found at: https://en.wikipedia.org/wiki/Orthographic_projection
GM_Matrix4 GM_Matrix4::orthographic(float left, float right, float bottom, float top, float near_plane, float far_plane) {
    const float A = 2.0f / (right - left);
    const float B = 2.0f / (top - bottom);
    const float C = -2.0f / (far_plane - near_plane);
    const float D = -(right + left) / (right - left);
    const float E = -(top + bottom) / (top - bottom);
    const float F = -(far_plane + near_plane) / (far_plane - near_plane);

    GM_Matrix4 ret;
    ret.v = {
        GM_Vec4(A,  0,  0,  D),
        GM_Vec4(0,  B,  0,  E),
        GM_Vec4(0,  0,  C,  F),
        GM_Vec4(0,  0, -1,  0)
    };

    return ret;
}

// Found at: https://www.khronos.org/opengl/wiki/GluLookAt_code
GM_Matrix4 GM_Matrix4::lookat(GM_Vec3 position, GM_Vec3 target, GM_Vec3 world_up) {
    GM_Vec3 forward = (position - target).normalize();
    GM_Vec3 right   = GM_Vec3::cross(world_up, forward).normalize();
    GM_Vec3 up      = GM_Vec3::cross(forward, right).normalize();

    GM_Matrix4 rotation;
    rotation.v = {
        GM_Vec4(right.x,   right.y,   right.z,   0),
        GM_Vec4(up.x,      up.y,      up.z,      0),
        GM_Vec4(forward.x, forward.y, forward.z, 0),
        GM_Vec4(0,         0,         0,         1)
    };
    
    GM_Matrix4 translation = GM_Matrix4::translate(GM_Matrix4::identity(), -position.x, -position.y, -position.z);

    return rotation * translation;
}

#define M_ELEM(mat, row, col) (mat.data[((row) * 4) + (col)])

float gm_mat3_determinant_helper(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
    return a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g);
}

GM_Matrix4 GM_Matrix4::inverse(GM_Matrix4 mat, bool* success) {
    if (success) {
        *success = false;
    }

    float m00 = mat.v[0].x, m01 = mat.v[0].y, m02 = mat.v[0].z, m03 = mat.v[0].w;
    float m10 = mat.v[1].x, m11 = mat.v[1].y, m12 = mat.v[1].z, m13 = mat.v[1].w;
    float m20 = mat.v[2].x, m21 = mat.v[2].y, m22 = mat.v[2].z, m23 = mat.v[2].w;
    float m30 = mat.v[3].x, m31 = mat.v[3].y, m32 = mat.v[3].z, m33 = mat.v[3].w;

    float c00 = gm_mat3_determinant_helper(m11, m12, m13, m21, m22, m23, m31, m32, m33);
    float c01 = gm_mat3_determinant_helper(m10, m12, m13, m20, m22, m23, m30, m32, m33);
    float c02 = gm_mat3_determinant_helper(m10, m11, m13, m20, m21, m23, m30, m31, m33);
    float c03 = gm_mat3_determinant_helper(m10, m11, m12, m20, m21, m22, m30, m31, m32);

    float det = m00 * c00 - m01 * c01 + m02 * c02 - m03 * c03;
    if (NEAR_ZERO(det)) {
        return GM_Matrix4::identity();
    }

    float invDet = 1.0f / det;
    GM_Matrix4 inv;

    // Row 0
    inv.v[0].x = invDet * c00;
    inv.v[0].y = invDet * (-gm_mat3_determinant_helper(m01, m02, m03, m21, m22, m23, m31, m32, m33));
    inv.v[0].z = invDet * gm_mat3_determinant_helper(m01, m02, m03, m11, m12, m13, m31, m32, m33);
    inv.v[0].w = invDet * (-gm_mat3_determinant_helper(m01, m02, m03, m11, m12, m13, m21, m22, m23));

    // Row 1
    inv.v[1].x = invDet * (-c01);
    inv.v[1].y = invDet * gm_mat3_determinant_helper(m00, m02, m03, m20, m22, m23, m30, m32, m33);
    inv.v[1].z = invDet * (-gm_mat3_determinant_helper(m00, m02, m03, m10, m12, m13, m30, m32, m33));
    inv.v[1].w = invDet * gm_mat3_determinant_helper(m00, m02, m03, m10, m12, m13, m20, m22, m23);

    // Row 2
    inv.v[2].x = invDet * c02;
    inv.v[2].y = invDet * (-gm_mat3_determinant_helper(m00, m01, m03, m20, m21, m23, m30, m31, m33));
    inv.v[2].z = invDet * gm_mat3_determinant_helper(m00, m01, m03, m10, m11, m13, m30, m31, m33);
    inv.v[2].w = invDet * (-gm_mat3_determinant_helper(m00, m01, m03, m10, m11, m13, m20, m21, m23));

    // Row 3
    inv.v[3].x = invDet * (-c03);
    inv.v[3].y = invDet * gm_mat3_determinant_helper(m00, m01, m02, m20, m21, m22, m30, m31, m32);
    inv.v[3].z = invDet * (-gm_mat3_determinant_helper(m00, m01, m02, m10, m11, m12, m30, m31, m32));
    inv.v[3].w = invDet * gm_mat3_determinant_helper(m00, m01, m02, m10, m11, m12, m20, m21, m22);

    if (success) {
        *success = true;
    }

    return inv;
}

GM_Matrix4 GM_Matrix4::operator*(const GM_Matrix4 &right) {
    GM_Matrix4 C;

    for (int i = 0; i < 4; i++) {
        C.v[i].x += this->v[i].x * right.v[0].x;
        C.v[i].x += this->v[i].y * right.v[1].x;
        C.v[i].x += this->v[i].z * right.v[2].x;
        C.v[i].x += this->v[i].w * right.v[3].x;

        C.v[i].y += this->v[i].x * right.v[0].y;
        C.v[i].y += this->v[i].y * right.v[1].y;
        C.v[i].y += this->v[i].z * right.v[2].y;
        C.v[i].y += this->v[i].w * right.v[3].y;

        C.v[i].z += this->v[i].x * right.v[0].z;
        C.v[i].z += this->v[i].y * right.v[1].z;
        C.v[i].z += this->v[i].z * right.v[2].z;
        C.v[i].z += this->v[i].w * right.v[3].z;
        
        C.v[i].w += this->v[i].x * right.v[0].w;
        C.v[i].w += this->v[i].y * right.v[1].w;
        C.v[i].w += this->v[i].z * right.v[2].w;
        C.v[i].w += this->v[i].w * right.v[3].w;
    }
        
    return C;
}
GM_Matrix4& GM_Matrix4::operator*=(const GM_Matrix4 &right) {
    *this = *this * right;
    return *this;
}

bool GM_Matrix4::operator==(const GM_Matrix4 &right) {
    return (this->v[0] == right.v[0]) && (this->v[1] == right.v[1]) && (this->v[2] == right.v[2]) && (this->v[3] == right.v[3]);
}

bool GM_Matrix4::operator!=(const GM_Matrix4 &right) {
    return !(*this == right);
}