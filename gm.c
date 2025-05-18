#define GM_IMPL
#include "gm.h"

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


float gm_v3_dot_product(GM_Vec3 A, GM_Vec3 B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

float gm_v4_dot_product(GM_Vec4 A, GM_Vec4 B) {
    return (A.x * B.x) + (A.y * B.y) + (A.z * B.z) + (A.w * B.w);
}

float gm_v4_magnitude(GM_Vec4 A) {
    return sqrt((A.x*A.x) + (A.x*A.x) + (A.z*A.z));
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

GM_Matix4 mat4_identity() {
    GM_Matix4 ret = {
        .data = {
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_mult(GM_Matix4 A, GM_Matix4 B) {
    GM_Matix4 ret = {0};
    
    for (int i = 0; i < 16; i++) {
        const row = i / 4;
        const column = i % 4;

        const GM_Vec4 a_row = A.v[i / 4];
        const GM_Vec4 b_column = {
            B.data[column + (0 * 4)], 
            B.data[column + (1 * 4)], 
            B.data[column + (2 * 4)], 
            B.data[column + (3 * 4)]
        };

        ret.data[i] = gm_v4_dot_product(a_row, b_column);
    }

    return ret;
}

GM_Matix4 mat4_translation(GM_Vec3 t) {
    GM_Matix4 ret = {
        .data = {
            1, 0, 0, t.x,
            0, 1, 0, t.y,
            0, 0, 1, t.z,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 mat4_scale(GM_Vec3 s) {
    GM_Matix4 ret = {
        .data = {
            s.x, 0, 0, 0,
            0, s.y, 0, 0,
            0, 0, s.z, 0,
            0, 0, 0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_x(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
        .data = {
            1, 0,  0, 0,
            0, c, -s, 0,
            0, s,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_y(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
        .data = {
            c, 0, -s, 0,
            0, 1,  0, 0,
            s, 0,  c, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

GM_Matix4 gm_mat4_rotation_z(float degrees) {
    float radians = DEGREES_TO_RAD(degrees);
    float c = cosf(radians);
    float s = sinf(radians);
    
    GM_Matix4 ret = {
        .data = {
            c, -s, 0, 0,
            s, c,  0, 0,
            0, 0,  1, 0,
            0, 0,  0, 1
        }
    };
    
    return ret;
}

// Z is negative going away from the viewer here so Right-handed coordinate system OpenGL
GM_Matix4 gm_mat4_perspective(float fov_degrees, float aspect, float near_plane, float far_plane) {
    float fov_radians = DEGREES_TO_RAD(fov_degrees);
    float p = 1 / (tan(fov_radians) / 2.0f);

    const range = near_plane - far_plane;
    const A = (-far_plane - near_plane) / range; 
    const B = (2 * far_plane * near_plane) / range; 

    GM_Matix4 ret = {
        .data = {
            p / aspect, 0, 0, 0,
            0, p, 0, 0,
            0, 0, A, B,
            0, 0, 1, 0
        }
    };
    
    return ret;
}

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

u32 gm_rgba_to_u32(GM_RGBA color) {
	u32 alpha = ((int)(color.rgba.a) << 24);
	u32 red = ((int)(color.rgba.r) << 16);
	u32 green = ((int)(color.rgba.g) << 8);
	u32 blue = ((int)(color.rgba.b) << 0);
					
	u32 rgba = alpha|red|green|blue;

	return rgba;
}

GM_RGBA gm_rgba_from_u32(u32 color) {
	GM_RGBA ret = {0};

	ret.rgba.b = ((color >> 0) & 0xFF); 
	ret.rgba.g = ((color >> 8) & 0xFF); 
	ret.rgba.r = ((color >> 16) & 0xFF); 
	ret.rgba.a = ((color >> 24) & 0xFF); 
				
	return ret;
}

GM_RGBA gm_rgba_multiply(GM_RGBA color, float value) {
	GM_RGBA ret = {0};
	ret.rgba.r = color.rgba.r * value;
	ret.rgba.g = color.rgba.g * value;
	ret.rgba.b = color.rgba.b * value;
	ret.rgba.a = color.rgba.a * value;

	return ret;
}

u32 gm_rgba_u32_multiply(u32 color, float value) {
	u8 b = (u8)(((color >> 0) & 0xFF)  * value);
	u8 g = (u8)(((color >> 8) & 0xFF)  * value);
	u8 r = (u8)(((color >> 16) & 0xFF) * value);
	u8 a = (u8)(((color >> 24) & 0xFF) * value);

	color = r|g|b|a;

	return color;
}

GM_RGBA gm_rgba_alpha_blend(GM_RGBA front_color, GM_RGBA back_color) {
	GM_RGBA ret = {0};

	float normalized_back_alpha = (float)back_color.a / 255.0f;

	ret.rgba.a = back_color.rgba.a;
	ret.rgba.r = (back_color.rgba.r * normalized_back_alpha) + ((u32)front_color.rgba.r * (1 - normalized_back_alpha));
	ret.rgba.g = (back_color.rgba.g * normalized_back_alpha) + ((u32)front_color.rgba.g * (1 - normalized_back_alpha));
	ret.rgba.b = (back_color.rgba.b * normalized_back_alpha) + ((u32)front_color.rgba.b * (1 - normalized_back_alpha));

	return ret;
}

GM_RGBA gm_rgba_u32_alpha_blend(u32 front_color_u32, u32 back_color_u32) {
	GM_RGBA front_color = {0};
	front_color.rgba.a = (u8)(((u32)front_color_u32 >> 24) & 0xFF);
	front_color.rgba.r = (u8)(((u32)front_color_u32 >> 16) & 0xFF);
	front_color.rgba.g = (u8)(((u32)front_color_u32 >> 8) & 0xFF);
	front_color.rgba.b = (u8)(((u32)front_color_u32 >> 0) & 0xFF);

	GM_RGBA back_color = {0};
	back_color.rgba.a = (u8)(((u32)back_color_u32 >> 24) & 0xFF);
	back_color.rgba.r = (u8)(((u32)back_color_u32 >> 16) & 0xFF);
	back_color.rgba.g = (u8)(((u32)back_color_u32 >> 8) & 0xFF);
	back_color.rgba.b = (u8)(((u32)back_color_u32 >> 0) & 0xFF);
				
	return gm_rgba_alpha_blend(front_color, back_color);
}


typedef struct GM_Bitmap {
	u32 width;
	u32 height;
	u8* memory;
} GM_Bitmap;


#pragma pack(push, 1)
typedef struct BmpHeader { 	// total 40 bytes
	char signature[2]; 			// 2 bytes
	u32 file_size; 				// 4 bytes
	u32 reserved;  				// 4 bytes
	u32 data_offset; 			// 4 bytes

	u32 size;					// 4 bytes
	u32 width;					// 4 bytes
	u32 height;					// 4 bytes
	u16 planes;					// 2 bytes
	u16 bits_per_pixel; 		// 2 bytes
	u32 compression;			// 4 bytes
	u32 compressed_image_size;	// 4 bytes
	u32 x_pixels_per_meter; 	// 4 bytes // horizontal resolution: Pixels/meter
	u32 y_pixels_per_meter; 	// 4 bytes // vertical resolution: Pixels/meter
	u32 colors_used;            // 4 bytes // Number of actually used colors. For a 8-bit / pixel bitmap this will be 100h or 256.
	u32 important_colors;		// 4 bytes // Number of important colors
} BmpHeader;
#pragma pack(pop)

GM_API GM_Bitmap ckit_parser_load_bmp(u8* bmp_file_data, size_t file_size);

GM_Bitmap ckit_parser_load_bmp(u8* bmp_file_data, size_t file_size) {
	GM_Bitmap ret = {0};

	BmpHeader bmp_header;
	ckit_memory_copy(bmp_file_data, &bmp_header, sizeof(bmp_header), file_size);
	ret.bytes_per_pixel = (u8)(bmp_header.bits_per_pixel / 8);
	u8* bmp_data = bmp_file_data + bmp_header.data_offset;

	// size_t bmp_size = file_size - sizeof(BmpHeader);
	ret.width = bmp_header.width;
	ret.height = bmp_header.height;
	ret.memory = bmp_data;
	
	return ret;
}