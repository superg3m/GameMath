#define GM_IMPL
#include "gm.h"






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

CKIT_Vector2 ckit_vector2_spline_point(CKIT_Vector2* spline_points, u32 spline_points_count, double t) {
    if (spline_points_count == 1) {
        return spline_points[0];
    }

    CKIT_Vector2* points_vector = NULLPTR;
    for (u32 i = 0; i < spline_points_count - 1; i++) { // get lerp points
        CKIT_Vector2 lerp_point = ckit_vector2_lerp(spline_points[i + 1], spline_points[i], t);
        ckit_vector_push(points_vector, lerp_point);
    }

    CKIT_Vector2 ret = ckit_vector2_spline_point(points_vector, ckit_vector_count(points_vector), t); // feed back points
    ckit_vector_free(points_vector);
    return ret;
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