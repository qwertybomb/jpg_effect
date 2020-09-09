/* standard headers */
#include <cstdint>
#include <vector>
#include <array>
#include <cmath>
#include <numbers>
#include <numeric>
#include <algorithm>

/* png++ headers */
#pragma clang  diagnostic push
/* prevent warnings */
#pragma clang diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#pragma clang diagnostic ignored "-Wunused-parameter"
    #include "png++/png.hpp"
#pragma clang  diagnostic pop

namespace 
{
	constexpr std::uint32_t block_size = 8;
	
	auto image_to_grid(std::array<std::uint8_t, 3> const* image, std::uint32_t const width, std::uint32_t const height)
	{
		/* create an array of blocks */
		std::vector<std::array<std::array<std::uint8_t, 3>, block_size *  block_size>> result_grid 
			{ (width / block_size) * (height / block_size) };
	
		/* copy each block of the image to the array of blocks */	
		for(std::uint32_t i = 0; i < width / block_size; ++i)
		{
			for(std::uint32_t j = 0; j < height / block_size; ++j)
			{
				for(std::uint32_t col_idx = i * block_size; col_idx < width && col_idx < (i + 1) * block_size; ++col_idx)
				{
					for(std::uint32_t row_idx = j * block_size; row_idx < height && row_idx < (j + 1) * block_size; ++row_idx)
					{
						result_grid[j * (width / block_size) + i]
						[(row_idx % block_size) * block_size + (col_idx % block_size)] = image[row_idx * width + col_idx];
					}	
				}
			}
		}
	
		return result_grid;
	}

	auto grid_to_image(std::vector<std::array<std::array<std::uint8_t, 3>, block_size *  block_size>> const &grid
	,std::array<std::uint8_t, 3> * image, std::uint32_t width, std::uint32_t height)
	{
		/* copy each block of the image to the array of blocks */	
		for(std::uint32_t i = 0; i < width / block_size; ++i)
		{
			for(std::uint32_t j = 0; j < height / block_size; ++j)
			{
				for(std::uint32_t col_idx = i * block_size; col_idx < width && col_idx < (i + 1) * block_size; ++col_idx)
				{
					for(std::uint32_t row_idx = j * block_size; row_idx < height && row_idx < (j + 1) * block_size; ++row_idx)
					{
						image[row_idx * width + col_idx] = grid[j * (width / block_size) + i]
						[(row_idx % block_size) * block_size + (col_idx % block_size)];
					}	
				}
			}
		}
	}
	
	/* discret cosine transform */
	auto compute_dct(std::array<std::array<std::uint8_t, 3>, block_size * block_size> const &in)
	{
		std::array<std::array<double, 3>, block_size * block_size> result;
	
		for (std::uint32_t i = 0; i < block_size; i++)
		{
			for (std::uint32_t j = 0; j < block_size; j++) 
	   		{	
	   			std::array<double, 3> s = {};
	      		for (std::uint32_t  u = 0; u < block_size; u++)
	      		{
	        		for (std::uint32_t v = 0; v < block_size; v++)
	        		{
	          			s[0] += ((double)in[v * block_size + u][0] - 127.0) * cos((2.0 * (double)u + 1.0) * (double)i * std::numbers::pi / ((double)block_size * 2.0)) * cos((2.0 * (double)v + 1.0) * (double)j * std::numbers::pi / ((double)block_size * 2.0)) * ((i == 0) ? 1.0 / sqrt(2) : 1.0) * ((j == 0) ? 1.0 / sqrt(2) : 1.0);
	          			s[1] += ((double)in[v * block_size + u][1] - 127.0) * cos((2.0 * (double)u + 1.0) * (double)i * std::numbers::pi / ((double)block_size * 2.0)) * cos((2.0 * (double)v + 1.0) * (double)j * std::numbers::pi / ((double)block_size * 2.0)) * ((i == 0) ? 1.0 / sqrt(2) : 1.0) * ((j == 0) ? 1.0 / sqrt(2) : 1.0);
	          			s[2] += ((double)in[v * block_size + u][2] - 127.0) * cos((2.0 * (double)u + 1.0) * (double)i * std::numbers::pi / ((double)block_size * 2.0)) * cos((2.0 * (double)v + 1.0) * (double)j * std::numbers::pi / ((double)block_size * 2.0)) * ((i == 0) ? 1.0 / sqrt(2) : 1.0) * ((j == 0) ? 1.0 / sqrt(2) : 1.0);
					}
				}
	      		s[0] /= sqrt(2.0 * (double)block_size);
	      		s[1] /= sqrt(2.0 * (double)block_size);
	      		s[2] /= sqrt(2.0 * (double)block_size);
				result[j * block_size + i] = s;
	       }
	   }
	   
	   return result;
	}
	
	/* inverse discret cosine transform */
	auto compute_idct(std::array<std::array<double, 3>, block_size * block_size> const &in)
	{
		std::array<std::array<std::uint8_t, 3>, block_size * block_size> result;
	
		for (std::uint32_t i = 0; i < block_size; i++)
		{
			for (std::uint32_t j = 0; j < block_size; j++) 
	   		{	
	   			std::array<double, 3> s = {};
	      		for (std::uint32_t  u = 0; u < block_size; u++)
	      		{
	        		for (std::uint32_t v = 0; v < block_size; v++)
	        		{
	          			s[0] += in[v * block_size + u][0] * cos((2.0 * (double)i + 1.0) * (double)u * std::numbers::pi / (double(block_size) * 2.0)) * cos((2.0 * (double)j + 1.0) * (double)v * std::numbers::pi / (double(block_size) * 2.0)) * ((u == 0) ? 1.0 / sqrtf(2) : 1.0) * ((v == 0) ? 1.0 / sqrtf(2) : 1.0);
						s[1] += in[v * block_size + u][1] * cos((2.0 * (double)i + 1.0) * (double)u * std::numbers::pi / (double(block_size) * 2.0)) * cos((2.0 * (double)j + 1.0) * (double)v * std::numbers::pi / (double(block_size) * 2.0)) * ((u == 0) ? 1.0 / sqrtf(2) : 1.0) * ((v == 0) ? 1.0 / sqrtf(2) : 1.0);
						s[2] += in[v * block_size + u][2] * cos((2.0 * (double)i + 1.0) * (double)u * std::numbers::pi / (double(block_size) * 2.0)) * cos((2.0 * (double)j + 1.0) * (double)v * std::numbers::pi / (double(block_size) * 2.0)) * ((u == 0) ? 1.0 / sqrtf(2) : 1.0) * ((v == 0) ? 1.0 / sqrtf(2) : 1.0);
					}
				}
	      		s[0] /= sqrt(2.0 * (double)block_size);
	      		s[1] /= sqrt(2.0 * (double)block_size);
	      		s[2] /= sqrt(2.0 * (double)block_size);
				result[j * block_size + i][0] = std::clamp(std::abs(s[0] + 127.0), 0.0, 255.0);
				result[j * block_size + i][1] = std::clamp(std::abs(s[1] + 127.0), 0.0, 255.0);
				result[j * block_size + i][2] = std::clamp(std::abs(s[2] + 127.0), 0.0, 255.0);					
	       }
	   }
	   
	   return result;
	}

	constexpr std::array<double, block_size * block_size> compute_quant_block(std::uint32_t quality_level)
	{
		std::array<double, block_size * block_size> result { 16, 11, 10, 16, 24,  40,  51,  61,  12, 12, 14, 19, 26,  58,  60,  55, 14, 13, 16, 24, 40,  57,  69,  56,  14, 17, 22, 29, 51,  87,  80,  62, 18, 22, 37, 56, 68,  109, 103, 77,  24, 35, 55, 64, 81,  104, 113, 92, 49, 64, 78, 87, 103, 121, 120, 101, 72, 92, 95, 98, 112, 100, 103, 99};
		double s = quality_level < 50 ? 5000 / quality_level : 200 - 2 * (quality_level);

		for(auto& value : result)
		{
			value = (std::int64_t)((s * value + 50.0) / 100.0);
			if (value == 0) value = 1;
		}

		return result;
    }

	void quantize(std::array<double, block_size * block_size> const &quant_block, std::array<std::array<double, 3>, block_size * block_size>  &block )
	{
		for(std::uint32_t i = 0; i < block_size * block_size; ++i)
		{
			block[i][0] = std::round(block[i][0] / quant_block[i]) * quant_block[i];
			block[i][1] = std::round(block[i][1] / quant_block[i]) * quant_block[i];
			block[i][2] = std::round(block[i][2] / quant_block[i]) * quant_block[i];
		}
	}
	
	void dct_compress(std::array<std::uint8_t, 3> *image, std::array<double, block_size * block_size> const &quant_block, std::uint32_t width, std::uint32_t height)
	{
		auto grid = image_to_grid(image, width, height);
		for(auto& block : grid)
		{
			std::array<std::array<double, 3>, block_size * block_size> dct_block = compute_dct(block);

			/* quantize */
			quantize(quant_block, dct_block);

			block = compute_idct(dct_block);
		}
		grid_to_image(grid, image, width, height);
	}

	auto readpng(char const *filepath)
	{
		/* load the png at filepath */
		png::image<png::rgb_pixel> input_image { filepath };

		/* create an array to store the image */
		std::vector<std::array<std::uint8_t, 3>> image_result { input_image.get_width() * input_image.get_height() }; 

		for(std::uint32_t i = 0; i < input_image.get_width(); ++i)
		{
			for(std::uint32_t j = 0; j < input_image.get_height(); ++j)
			{
				/* copy the pixels from input_image to image_result */
				image_result[j * input_image.get_width() + i][0] = input_image.get_pixel(i, j).red;
				image_result[j * input_image.get_width() + i][1] = input_image.get_pixel(i, j).green;
				image_result[j * input_image.get_width() + i][2] = input_image.get_pixel(i, j).blue;
			}
		}
		
		struct result_t
		{
			std::vector<std::array<std::uint8_t, 3>> image;	
			std::uint32_t width;
			std::uint32_t height;
		} result{ std::move(image_result), input_image.get_width(), input_image.get_height() };
		return result;
	}

	
	void writepng(char const *filepath, std::array<std::uint8_t, 3> const *image, std::uint32_t width, std::uint32_t height)
	{
		/* load the png at filepath */
		png::image<png::rgb_pixel> output_image { width, height };
		
		for(std::uint32_t i = 0; i < width; ++i)
		{
			for(std::uint32_t j = 0; j < height; ++j)
			{
				/* copy the pixels from input_image to image_result */
				output_image.set_pixel(i, j,  png::rgb_pixel{image[j * width + i][0], image[j * width + i][1], image[j * width + i][2]});
			}
		}

		output_image.write(filepath);
	}
}	

int main(int argc, char *argv[])
{
	std::uint32_t quality_level = 50;
	auto quant_block = compute_quant_block(quality_level);
	if (argc % 2 == 1)
	{
		for(int i = 1; i < argc; i += 2)
		{
			if(strcmp(argv[i], "-q") == 0)
			{
				quality_level = 0;
				for(; argv[i + 1][0] != '\0'; argv[i + 1]++)
				{
					if(argv[i + 1][0] < '0' || argv[i + 1][0] > '9')
					{
						goto invalid_arg;
					}
					else 
					{
						quality_level *= 10;
						quality_level += argv[i + 1][0] - '0';
					}
				}
				quant_block = compute_quant_block(quality_level);
			}
			else
			{
				auto [image, width, height] = readpng(argv[i]);
				dct_compress(image.data(), quant_block, width, height);
				writepng(argv[i + 1], image.data(), width, height);
			}
		}
	}
	else
	{
		invalid_arg:
		std::cerr << "Error: invalid arguments\n";
		std::cerr << "Usage: " << argv[0] << " [flags] [input_png] [output_png] ...\n";
		std::cerr << "Flags: [-q n] sets the quailty level to n; the default level is set to 50\n";
		return -1;
	}
}
