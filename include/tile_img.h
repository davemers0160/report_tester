// This file was created by Victoria Lockridge, NSWC Crane
// July 18, 2023
#pragma once

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui.hpp"
#include <opencv2/opencv.hpp>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "cell_information.h"

// Purpose: To know where each tile image falls in overall image to relate image to DEF cells
// Need this information to find DEF cells to overlay and process to new locations.

class Tile_img
{
public:
	// Attributes
	std::string image_file{ "" };
	std::string txt_file{ "" };
	int32_t height = 0; // Height of image in pixels
	int32_t width = 0; // Width of image in pixels
	double pixel_size_x = 1.0; // from PixelSize nm/pixel in.txt (resolution)
	double pixel_size_y = 1.0;
	int32_t initLocColX = 0; // X initial location in pixels in stitched image
	int32_t initLocRowY = 0; // Y initial location in pixels in stitched image
	int32_t imgLocColX = 0; // X final location in pixels in stitched image after image defined from(0, 0) at upper left corner
	int32_t imgLocRowY = 0; // Y final location in pixels in stitched image after image defined from(0, 0) at upper left corner
	std::vector<Cell_info> tile_cells; // standard cells/gates associated with this tile
	bool aligned = false; // After tilecells are transformed by operator, "aligned" becomes true
	bool binarized = false;
	uint8_t struct_elem_size = 5;
	uint16_t binary_threshold = 164;
	std::vector<cv::Point> a_matrix;
	std::vector<cv::Point> b_matrix;
	std::vector<Cell_info> corner_pts;
	std::vector<Cell_info> orig_corner_pts;
	std::vector<double> gnd_pixels_y;
	std::vector<double> pwr_pixels_y;
	bool scored = false;
	std::vector<cv::Scalar> used_colors;

	std::string img_filename = "";

	//Methods/Functions
	Tile_img(std::string img_file, std::string txt_file);
	Tile_img(const cv::Mat& current_tile_image, std::string img_file, int32_t core_box_dx, int32_t core_box_dy);
	Tile_img(int32_t imgLocColX_, int32_t imgLocRowY_, int32_t h_, int32_t w_) : imgLocColX(imgLocColX_), imgLocRowY(imgLocRowY_), height(h_), width(w_) {};
	Tile_img(int32_t imgLocColX_, int32_t imgLocRowY_, int32_t h_, int32_t w_, bool al_, bool bin_) : imgLocColX(imgLocColX_), imgLocRowY(imgLocRowY_), height(h_), width(w_), aligned(al_), binarized(bin_) {};

	~Tile_img();
	void read_txtdata(std::string txt_file);
	void find_pixel_size(const cv::Mat& rotated_img, std::string file_name, int32_t core_box_dx, int32_t core_box_dy);
	void find_ImgLoc(Tile_img& tile0);
	void get_tile_cells(const std::vector<Cell_info>& def_cells);
	void plot_DEFpts(const std::vector<Cell_info>& def_cells);
	cv::Scalar get_color(std::string cellname);
};
