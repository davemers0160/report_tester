#define _CRT_SECURE_NO_WARNINGS

//#include <ryml_all.hpp>

#if defined(_WIN32) | defined(__WIN32__) | defined(__WIN32) | defined(_WIN64) | defined(__WIN64)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <cstdint>
#include <cstdlib>
#include <string>
#include <iostream>
#include <atomic>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <fstream>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// project specific includes
#include "cell_information.h"
//#include "utilities.h"
#include "utils.h"
//#include "display.h"
#include "tile_img.h"
#include "matched_template.h"

#include "html_tags.h"


//#include "ocv_threshold_functions.h"
//#include "get_directory_listing.h"
//#include "file_ops.h"

//#include "utilities.h"
//#include "image_tile.h"
//#include "get_tiled_images.h"
////#include "parse_cell_input.h"
//#include "cell_helper.h"
//#include "num2string.h"


//#include "temp_functions.h"





//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
    uint32_t idx, jdx;

    char key = 0;

    //uint32_t img_w, img_h, img_c;
    double threshold = 30.0;
    uint64_t x = 0;
    uint64_t y = 0;
    uint32_t tile_x, tile_y;
    uint32_t tile_w, tile_h;
    uint32_t tile_index = 0;

    uint64_t cell_w, cell_h;
    int32_t cells_max_x, cells_max_y;
    uint64_t overlap_x = 0, overlap_y = 0;
    uint32_t img_w = 9000;
    uint32_t img_h = 3000;


    std::vector<int32_t> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(8);


    //if (argc < 2)
    //{
    //    std::cout << "Enter the location of the image..." << std::endl;
    //    std::cin.ignore();
    //}


    //----------------------------------------------------------------------------
    // random image generation piece
    std::string cv_window = "Image";
    cv::namedWindow(cv_window, cv::WINDOW_NORMAL | cv::WINDOW_KEEPRATIO);
    int bp = 0;
    
    cv::RNG rng;
    cv::Mat current_tile_image = cv::Mat(img_h, img_w, CV_8UC3);
    rng.fill(current_tile_image, cv::RNG::UNIFORM, 0, 255);

    //----------------------------------------------------------------------------
    // random matched template generation piece
    uint32_t mt1_w = 100;
    uint32_t mt1_h = 30;
    uint32_t mt2_w = 60;
    uint32_t mt2_h = 30;

    cv::Mat mt1 = cv::Mat(mt1_h, mt1_w, CV_8UC1);
    rng.fill(mt1, cv::RNG::UNIFORM, 0, 255);

    cv::Mat mt2 = cv::Mat(mt2_h, mt2_w, CV_8UC1);
    rng.fill(mt2, cv::RNG::UNIFORM, 0, 255);

    std::vector<matched_template> matched_templates;

    // put the template 1 image into the matched template container
    matched_templates.push_back(matched_template("template 1"));
    matched_templates[0].imgs.push_back(mt1);
    matched_templates[0].imgs.push_back(mt1);
    matched_templates[0].imgs.push_back(mt1);
    matched_templates[0].imgs.push_back(mt1);

    // put the template 2 image into the matched template container
    matched_templates.push_back(matched_template("template 2"));
    matched_templates[1].imgs.push_back(mt2);
    matched_templates[1].imgs.push_back(mt2);
    matched_templates[1].imgs.push_back(mt2);
    matched_templates[1].imgs.push_back(mt2);


    //----------------------------------------------------------------------------
    // random cell location generation piece
    std::vector<Pt_info> ct1_pts;
    std::vector<Pt_info> ct2_pts;

    uint32_t num_cells1 = 1000;
    uint32_t num_cells2 = 600;

     


    //----------------------------------------------------------------------------
    std::vector<Tile_img> tiles;
    bool whole = false;


    //----------------------------------------------------------------------------
    std::vector<Cell_info> cell_list;
    cell_list.push_back(Cell_info("cell_type_1", mt1_h, mt1_w, { 0,1,2,3 }, ct1_pts));
    cell_list.push_back(Cell_info("cell_type_2", mt2_h, mt2_w, { 0,1,2,3 }, ct2_pts));


    //----------------------------------------------------------------------------
    // detects
    std::vector<std::vector<Pt_info>> detects;



    //----------------------------------------------------------------------------
    // conversion info
    Pixel_to_DEF_conversion def_conversion_info;



    //----------------------------------------------------------------------------
    // report piece
    // directory to save report
    std::string report_dir = "../report";

    // Creates data directory if it doesn't exist
    std::string html_img_dir = report_dir + "/html_images";
    std::string templ_img_dir = report_dir + "/template_images";
    std::string h_html_dir = "html_images";
    std::string h_templ_dir = "template_images";
    std::filesystem::create_directory(html_img_dir);
    std::filesystem::create_directory(templ_img_dir);


    std::ofstream out_file(report_dir + "/report.html");
    write_toplines(out_file);
    write_template_imgs(matched_templates, cell_list, templ_img_dir);

    // Write heading
    Html_top_heading top_heading;
    out_file << top_heading;

    //******************* IMAGE PART
    //Html_img img_whole(cell_list, detects, tiles, whole, current_tile_image, html_img_dir, templ_img_dir, def_conversion_info, matched_templates);
    Html_img img_whole(cell_list, detects, tiles, whole, current_tile_image, html_img_dir, h_html_dir, h_templ_dir, def_conversion_info, matched_templates);
    out_file << img_whole;;
    //**********************

    write_ending(out_file);

    out_file.close();


    std::cout << "Press Enter to close..." << std::endl;

    std::cin.ignore();

    cv::destroyAllWindows();

    return 0;
}
