// Vector
// Matrix

// row major
typedef struct GM_Matix4 {
    float data[16]; 
} Matrix4;

typedef struct GM_Vec2 {
    float x, y;
} GM_Vec2;

typedef struct GM_Vec3 {
    float x, y, z;
} GM_Vec3;

typedef struct GM_Vec4 {
    float x, y, z, w;
} GM_Vec4;


// Creation
Matrix4 mat4_identity();
Matrix4 mat4_translation(GM_Vec3 t);
Matrix4 mat4_scale(GM_Vec3 s);

// Convert to radians in the implementation
Matrix4 mat4_rotation_x(float degrees);
Matrix4 mat4_rotation_y(float degress);
Matrix4 mat4_rotation_z(float degress);


Matrix4 mat4_perspective(float fov_rad, float aspect, float near, float far);
Matrix4 mat4_orthographic(float left, float right, float bottom, float top, float near, float far);
Matrix4 mat4_look_at(GM_Vec3 eye, GM_Vec3 center, GM_Vec3 up);

// Math
Matrix4 gm_mat4_mult(Matrix4 a, Matrix4 b);
GM_Vec4 mat4_mult_vec4(Matrix4 m, GM_Vec4 v);
Matrix4 mat4_inverse(Matrix4 m);
Matrix4 mat4_transpose(Matrix4 m);

// float  gm_dot_product()
