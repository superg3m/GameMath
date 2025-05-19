#define GM_IMPL
#include "gm.h"

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