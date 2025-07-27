#include "./gm.hpp"

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

GM_Vec2::GM_Vec2() {
    this->x = 0.0f;
    this->y = 0.0f;
}

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

GM_Vec2 GM_Vec2::scale(float scale) const {
    return GM_Vec2(this->x * scale, this->y * scale);
}

GM_Vec2 GM_Vec2::scale(GM_Vec2 s) const {
    return GM_Vec2(this->x * s.x, this->y * s.y);
}

GM_Vec2 GM_Vec2::scale(float scale_x, float scale_y) const {
    return GM_Vec2(this->x * scale_x, this->y * scale_y);
}

float GM_Vec2::dot(GM_Vec2 a, GM_Vec2 b) {
    return (a.x * b.x) + (a.y * b.y);
}

float GM_Vec2::distance(GM_Vec2 a, GM_Vec2 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y));
}

float GM_Vec2::distanceSquared(GM_Vec2 a, GM_Vec2 b) {
    return SQUARED(b.x - a.x) + SQUARED(b.y - a.y);
}

GM_Vec2 GM_Vec2::lerp(GM_Vec2 a, GM_Vec2 b, float t) {
    GM_Vec2 ab = (b - a);
    return a + (ab.scale(t));
}

GM_Vec2 GM_Vec2::euler(float yaw, float pitch) {
    GM_Vec2 ret = GM_Vec2(0, 0);
    ret.x = cosf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));
    ret.y = sinf(DEGREES_TO_RAD(pitch));

    return ret;
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

GM_Vec3::GM_Vec3() {
    this->x = 0.0f;
    this->y = 0.0f;
    this->z = 0.0f;
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

GM_Vec3::GM_Vec3(GM_Vec2 v, float z) {
    this->x = v.x;
    this->y = v.y;
    this->z = z;
}

GM_Vec3::GM_Vec3(GM_Vec4 v) {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
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

GM_Vec3 GM_Vec3::scale(float scale) const {
    return GM_Vec3(this->x * scale, this->y * scale, this->z * scale);
}

GM_Vec3 GM_Vec3::scale(GM_Vec3 s) const {
    return GM_Vec3(this->x * s.x, this->y * s.y, this->z * s.z);
}

GM_Vec3 GM_Vec3::scale(float scale_x, float scale_y, float scale_z) const {
    return GM_Vec3(this->x * scale_x, this->y * scale_y, this->z * scale_z);
}

float GM_Vec3::distance(GM_Vec3 a, GM_Vec3 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z));
}

float GM_Vec3::distanceSquared(GM_Vec3 a, GM_Vec3 b) {
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

GM_Vec3 GM_Vec3::euler(float yaw, float pitch) {
    GM_Vec3 ret = GM_Vec3(0, 0, 0);
    ret.x = cosf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));
    ret.y = sinf(DEGREES_TO_RAD(pitch));
    ret.z = sinf(DEGREES_TO_RAD(yaw)) * cosf(DEGREES_TO_RAD(pitch));

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

GM_Vec4::GM_Vec4() {
    this->x = 0.0f;
    this->y = 0.0f;
    this->z = 0.0f;
    this->w = 0.0f;
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

GM_Vec4::GM_Vec4(GM_Vec3 v, float w) {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
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

GM_Vec4 GM_Vec4::scale(float scale) const {
    return GM_Vec4(this->x * scale, this->y * scale, this->z * scale, this->w * scale);
}

GM_Vec4 GM_Vec4::scale(GM_Vec4 s) const {
    return GM_Vec4(this->x * s.x, this->y * s.y, this->z * s.z, this->w * s.w);
}

GM_Vec4 GM_Vec4::scale(float scale_x, float scale_y, float scale_z, float scale_w) const {
    return GM_Vec4(this->x * scale_x, this->y * scale_y, this->z * scale_z, this->w * scale_w);
}

float GM_Vec4::distance(GM_Vec4 a, GM_Vec4 b) {
    return sqrtf(SQUARED(b.x - a.x) + SQUARED(b.y - a.y) + SQUARED(b.z - a.z) + SQUARED(b.w - a.w));
}

float GM_Vec4::distanceSquared(GM_Vec4 a, GM_Vec4 b) {
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

GM_AABB GM_AABB::fromCenterExtents(GM_Vec3 center, GM_Vec3 extents) {
    return GM_AABB(center - extents, center + extents);
}

bool GM_AABB::intersection(GM_AABB aabb, GM_Vec3 p0, GM_Vec3 p1) {
    float tmin = -10000;
    float tmax = 10000;
 
    // X coordinate
    if (fabs(p1.x) > EPSILON) {
        float t1 = (aabb.min.x - p0.x) / p1.x;
        float t2 = (aabb.max.x - p0.x) / p1.x;
 
        tmin = MAX(tmin, MIN(t1, t2));
        tmax = MIN(tmax, MAX(t1, t2));
    }
 
    // Y coordinate
    if (fabs(p1.y) > EPSILON) {
        float t1 = (aabb.min.y - p0.y) / p1.y;
        float t2 = (aabb.max.y - p0.y) / p1.y;
 
        tmin = MAX(tmin, MIN(t1, t2));
        tmax = MIN(tmax, MAX(t1, t2));
    }
 
    // Z coordinate
    if (fabs(p1.z) > EPSILON) {
        float t1 = (aabb.min.z - p0.z) / p1.z;
        float t2 = (aabb.max.z - p0.z) / p1.z;
 
        tmin = MAX(tmin, MIN(t1, t2));
        tmax = MIN(tmax, MAX(t1, t2));
    }
 
 
    if (tmax > tmin && tmax > 0.0) {
        return true;
    } else {
        return false;
    }
}

GM_AABB::GM_AABB() {
    this->min = GM_Vec3(0, 0, 0);
    this->max = GM_Vec3(0, 0, 0);
}

GM_AABB::GM_AABB(GM_Vec3 min, GM_Vec3 max) {
    this->min = min;
    this->max = max;
}

GM_AABB::GM_AABB(float min_x, float min_y, float min_z, float max_x, float max_y, float max_z) {
    this->min.x = min_x;
    this->min.y = min_y;
    this->min.z = min_z;

    this->max.x = max_x;
    this->max.y = max_y;
    this->max.z = max_z;
}

GM_Vec3 GM_AABB::getCenter() {
    GM_Vec3 extents = this->getExtents();
    return GM_Vec3(min.x + extents.x, min.y + extents.y, min.z + extents.z);
}

GM_Vec3 GM_AABB::getExtents() {
    return (max - min).scale(0.5f);
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

GM_Matrix4 GM_Matrix4::fromColumnMajor(const float mat[16]) {
    GM_Matrix4 ret = {
        GM_Vec4{mat[0], mat[4], mat[8], mat[12]},
        GM_Vec4{mat[1], mat[5], mat[9], mat[13]},
        GM_Vec4{mat[2], mat[6], mat[10], mat[14]},
        GM_Vec4{mat[3], mat[7], mat[11], mat[15]},
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

GM_Matrix4 GM_Matrix4::rotate(GM_Matrix4 mat, GM_Quaternion quat) {
    float theta;
    GM_Vec3 axis;
    quat.toAngleAxis(theta, axis);
    return GM_Matrix4::rotate(mat, theta, axis);
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

void GM_Matrix4::decompose(GM_Matrix4 mat, GM_Vec3* out_position, GM_Quaternion* out_orientation, GM_Vec3* out_scale) {
    GM_Vec3 translation = GM_Vec3(mat.v[0].w, mat.v[1].w, mat.v[2].w);

    GM_Vec3 scale = GM_Vec3(0);
    {
        GM_Vec3 column1 = GM_Vec3(mat.v[0].x, mat.v[1].x, mat.v[2].x);
        GM_Vec3 column2 = GM_Vec3(mat.v[0].y, mat.v[1].y, mat.v[2].y);
        GM_Vec3 column3 = GM_Vec3(mat.v[0].z, mat.v[1].z, mat.v[2].z);
        
        float scale_x = column1.magnitude();
        float scale_y = column2.magnitude();
        float scale_z = column3.magnitude();

        scale = GM_Vec3(scale_x, scale_y, scale_z);
    }

    GM_Quaternion orientation = GM_Quaternion::identity();
    {
        GM_Vec3 column1 = GM_Vec3(mat.v[0].x, mat.v[1].x, mat.v[2].x);
        GM_Vec3 column2 = GM_Vec3(mat.v[0].y, mat.v[1].y, mat.v[2].y);
        GM_Vec3 column3 = GM_Vec3(mat.v[0].z, mat.v[1].z, mat.v[2].z);

        // Normalize columns to get pure rotation basis
        column1 = column1.scale(scale.x);
        column2 = column2.scale(scale.y);
        column3 = column3.scale(scale.z);

        GM_Matrix4 rotation_matrix = GM_Matrix4(
            GM_Vec4{column1.x, column2.x, column3.x, 0},
            GM_Vec4{column1.y, column2.y, column3.y, 0},
            GM_Vec4{column1.z, column2.z, column3.z, 0},
            GM_Vec4{0,         0,         0,         0}
        );

        orientation = GM_Quaternion::fromRotationMatrix(rotation_matrix);
    }

    if (out_position) {
        *out_position = translation;
    }

    if (out_scale) {
        *out_scale = scale;
    }

    if (out_orientation) {
        *out_orientation = orientation;
    }
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

    GM_Matrix4 ret = {
        A,  0,  0,  D,
        0,  B,  0,  E,
        0,  0,  C,  F,
        0,  0,  0,  1
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

float gm_mat3_determinant_helper(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
    return a * (e * i - f * h) -
            b * (d * i - f * g) +
            c * (d * h - e * g);
}

GM_Matrix4 GM_Matrix4::inverse(bool &success) {
    success = false;

    float m00 = this->v[0].x, m01 = this->v[0].y, m02 = this->v[0].z, m03 = this->v[0].w;
    float m10 = this->v[1].x, m11 = this->v[1].y, m12 = this->v[1].z, m13 = this->v[1].w;
    float m20 = this->v[2].x, m21 = this->v[2].y, m22 = this->v[2].z, m23 = this->v[2].w;
    float m30 = this->v[3].x, m31 = this->v[3].y, m32 = this->v[3].z, m33 = this->v[3].w;

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

    success = true;

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

GM_Vec4 GM_Matrix4::operator*(const GM_Vec4 &right) {
    GM_Vec4 ret;
    ret.x += this->v[0].x * right.x;
    ret.x += this->v[0].y * right.y;
    ret.x += this->v[0].z * right.z;
    ret.x += this->v[0].w * right.w;

    ret.y += this->v[1].x * right.x;
    ret.y += this->v[1].y * right.y;
    ret.y += this->v[1].z * right.z;
    ret.y += this->v[1].w * right.w;

    ret.z += this->v[2].x * right.x;
    ret.z += this->v[2].y * right.y;
    ret.z += this->v[2].z * right.z;
    ret.z += this->v[2].w * right.w;
    
    ret.w += this->v[3].x * right.x;
    ret.w += this->v[3].y * right.y;
    ret.w += this->v[3].z * right.z;
    ret.w += this->v[3].w * right.w;
        
    return ret;
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

GM_Quaternion::GM_Quaternion() {
    this->w = 1;
    this->v = GM_Vec3(0, 0, 0);
}

GM_Quaternion::GM_Quaternion(float theta, GM_Vec3 axis) {
    float radians = DEGREES_TO_RAD(theta);
    this->w = cosf(radians / 2.0f);
    if (NEAR_ZERO(this->w)) {
        this->w = 0.0f;
    }

    axis = axis.normalize();
    this->v = axis.scale(sinf(radians / 2.0f));
}

GM_Quaternion::GM_Quaternion(float theta, float x, float y, float z) {
    *this = GM_Quaternion(theta, GM_Vec3(x, y, z));
}

GM_Quaternion GM_Quaternion::inverse() {
    GM_Quaternion ret(1, 0, 0, 0);

    float magnitude_squared = SQUARED(this->w) + GM_Vec3::dot(this->v, this->v);
    if (magnitude_squared == 0.0f) { 
        return GM_Quaternion::identity();
    }

    ret.w = this->w / magnitude_squared;
    ret.v = this->v.scale(-1.0f / magnitude_squared);

    return ret;
}

GM_Quaternion GM_Quaternion::scale(float scale) {
    GM_Quaternion ret;

    ret.w   = this->w   * scale;
    ret.v.x = this->v.x * scale;
    ret.v.y = this->v.y * scale;
    ret.v.z = this->v.z * scale;

    return ret;
}

GM_Quaternion GM_Quaternion::normalize() {
    GM_Vec4 temp = GM_Vec4(this->w, this->v.x, this->v.y, this->v.z).normalize();
    
    GM_Quaternion ret;
    ret.w = temp.x;
    ret.v.x = temp.y;
    ret.v.y = temp.z;
    ret.v.z = temp.w;

    return ret;
}

GM_Quaternion GM_Quaternion::identity() {
    return {1, 0, 0, 0};
}

GM_Quaternion GM_Quaternion::literal(float w, GM_Vec3 axis) {
    GM_Quaternion ret;
    ret.w = w;
    ret.v = axis;

    return ret;
}

GM_Quaternion GM_Quaternion::literal(float w, float x, float y, float z) {
    GM_Quaternion ret;
    ret.w = w;
    ret.v = GM_Vec3(x, y, z);

    return ret;
}

GM_Quaternion GM_Quaternion::fromEuler(GM_Vec3 euler_angles_degrees) {
    float roll_rad_half = DEGREES_TO_RAD(euler_angles_degrees.x) * 0.5f;
    float pitch_rad_half = DEGREES_TO_RAD(euler_angles_degrees.y) * 0.5f;
    float yaw_rad_half = DEGREES_TO_RAD(euler_angles_degrees.z) * 0.5f;

    float cx = cosf(roll_rad_half);
    float sx = sinf(roll_rad_half);
    float cy = cosf(pitch_rad_half);
    float sy = sinf(pitch_rad_half);
    float cz = cosf(yaw_rad_half);
    float sz = sinf(yaw_rad_half);

    GM_Quaternion q = GM_Quaternion::identity();
    q.w = cx * cy * cz + sx * sy * sz;
    q.v.x = sx * cy * cz - cx * sy * sz;
    q.v.y = cx * sy * cz + sx * cy * sz;
    q.v.z = cx * cy * sz - sx * sy * cz;

    return q;
}

GM_Quaternion GM_Quaternion::fromAngleAxis(float angle, GM_Vec3 axis) {
    float half_angle = DEGREES_TO_RAD(angle) * 0.5f;
    float sinf_half = sinf(half_angle);
    float cosf_half = cosf(half_angle);

    GM_Quaternion q = GM_Quaternion::identity();
    axis = axis.normalize();

    q.w     = cosf_half;
    q.v.x   = axis.x * sinf_half;
    q.v.y   = axis.y * sinf_half;
    q.v.z   = axis.z * sinf_half;

    return q;
}

GM_Quaternion GM_Quaternion::fromRotationMatrix(const float m[16]) {
    GM_Quaternion q;

    float m00 = m[0],  m01 = m[1],  m02 = m[2];
    float m10 = m[4],  m11 = m[5],  m12 = m[6];
    float m20 = m[8],  m21 = m[9],  m22 = m[10];
    
    float t;
    if (m22 < 0) {
        if (m00 > m11) {
            t = 1 + m00 - m11 - m22;
            q.v.x = t;
            q.v.y = m01 + m10;
            q.v.z = m02 + m20;
            q.w = m21 - m12;
        } else {
            t = 1 - m00 + m11 - m22;
            q.v.x = m01 + m10;
            q.v.y = t;
            q.v.z = m12 + m21;
            q.w = m02 - m20;
        }
    } else {
        if (m00 < -m11) {
            t = 1 - m00 - m11 + m22;
            q.v.x = m02 + m20;
            q.v.y = m12 + m21;
            q.v.z = t;
            q.w = m10 - m01;
        } else {
            t = 1 + m00 + m11 + m22;
            q.v.x = m21 - m12;
            q.v.y = m02 - m20;
            q.v.z = m10 - m01;
            q.w = t;
        }
    }
    
    float s = 0.5f / sqrtf(t);
    q.v.x *= s;
    q.v.y *= s;
    q.v.z *= s;
    q.w   *= s;

    return q;
}

GM_Quaternion GM_Quaternion::fromRotationMatrix(GM_Matrix4 mat) {
    return GM_Quaternion::fromRotationMatrix(&mat.v[0].x);
}

GM_Matrix4 GM_Quaternion::toMatrix4() {
    GM_Matrix4 result = GM_Matrix4::identity();

    float x2 = this->v.x * this->v.x;
    float y2 = this->v.y * this->v.y;
    float z2 = this->v.z * this->v.z;

    float xy = this->v.x * this->v.y;
    float xz = this->v.x * this->v.z;
    float yz = this->v.y * this->v.z;
    float xw = this->v.x * this->w;
    float yw = this->v.y * this->w;
    float zw = this->v.z * this->w;

    result.v[0].x = 1.0f - 2.0f * (y2 + z2);  // m00
    result.v[0].y = 2.0f * (xy - zw);         // m01
    result.v[0].z = 2.0f * (xz + yw);         // m02
    result.v[0].w = 0.0f;                     // m03

    result.v[1].x = 2.0f * (xy + zw);         // m10
    result.v[1].y = 1.0f - 2.0f * (x2 + z2);  // m11
    result.v[1].z = 2.0f * (yz - xw);         // m12
    result.v[1].w = 0.0f;                     // m13

    result.v[2].x  = 2.0f * (xz - yw);        // m20
    result.v[2].y  = 2.0f * (yz + xw);        // m21
    result.v[2].z  = 1.0f - 2.0f * (x2 + y2); // m22
    result.v[2].w  = 0.0f;                    // m23

    result.v[3].x = 0.0f;                     // m30
    result.v[3].y = 0.0f;                     // m31
    result.v[3].z = 0.0f;                     // m32
    result.v[3].w = 1.0f;                     // m33

    return result;
}

void GM_Quaternion::toAngleAxis(float &theta, GM_Vec3 &vec) {
    GM_Quaternion quat = this->normalize();
    float sinf_half_theta = quat.v.magnitude();

    if (sinf_half_theta < EPSILON) {
        vec = GM_Vec3(1, 0, 0);
    } else {
        vec = quat.v.scale(1.0f / sinf_half_theta);
    }

    // Clamp w to [-1, 1] to avoid NaNs due to precision issues
    float w = CLAMP(quat.w, -1.0f, 1.0f);
    theta = 2.0f * acosf(w);
    theta = RAD_TO_DEGREES(theta);
}

GM_Quaternion GM_Quaternion::slerp(GM_Quaternion q, GM_Quaternion r, float t) {
    q = q.normalize();
    r = r.normalize();
    float dot = GM_Quaternion::dot(q, r);

    if (dot < 0.0f) {
        r = r.scale(-1.0f);
        dot = -dot;
    }

    if (dot > 0.9995f) {
        GM_Quaternion lerp = q + (r - q).scale(t);
        return lerp.normalize();
    }

    float theta_0 = acosf(dot);
    float theta = theta_0 * t;

    GM_Quaternion q3 = r - q.scale(dot);
    q3 = q3.normalize();

    GM_Quaternion term1 = q.scale(cosf(theta));
    GM_Quaternion term2 = q3.scale(sinf(theta));
    return term1 + term2;
}

float GM_Quaternion::dot(GM_Quaternion a, GM_Quaternion b) {
    float dot = a.w   * b.w   +
                a.v.x * b.v.x +
                a.v.y * b.v.y +
                a.v.z * b.v.z;

    return dot;
}

GM_Quaternion GM_Quaternion::operator+(const GM_Quaternion &right) {
    GM_Quaternion ret = GM_Quaternion::identity();

    ret.w = this->w + right.w;
    ret.v = this->v + right.v;

    return ret;
}
GM_Quaternion& GM_Quaternion::operator+=(const GM_Quaternion &right) {
    *this = *this + right;
    return *this;
}

GM_Quaternion GM_Quaternion::operator-(const GM_Quaternion &right) {
    GM_Quaternion ret = GM_Quaternion::identity();

    ret.w = this->w - right.w;
    ret.v = this->v - right.v;

    return ret;
}
GM_Quaternion& GM_Quaternion::operator-=(const GM_Quaternion &right) {
    *this = *this - right;
    return *this;
}

GM_Quaternion GM_Quaternion::operator*(const GM_Quaternion &right) {
    GM_Quaternion ret = GM_Quaternion::identity();
    ret.w = (this->w * right.w) - GM_Vec3::dot(this->v, right.v);
    ret.v = (this->v.scale(right.w) + right.v.scale(this->w)) + GM_Vec3::cross(this->v, right.v);
    
    return ret;
}
GM_Quaternion& GM_Quaternion::operator*=(const GM_Quaternion &right) {
    *this = *this * right;
    return *this;
}

GM_Vec3 GM_Quaternion::operator*(const GM_Vec3 &right) {
    GM_Quaternion q = *this;
    GM_Quaternion p = GM_Quaternion::literal(0, right);

    return ((q * p) * q.inverse()).v;
}

bool GM_Quaternion::operator==(const GM_Quaternion &right) {
    return NEAR_ZERO(this->w - right.w) && (this->v == right.v);
}
bool GM_Quaternion::operator!=(const GM_Quaternion &right) {
    return !(*this == right);
}