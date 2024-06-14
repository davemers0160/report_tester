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

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

// project specific includes
#include "cell_information.h"
#include "utilities.h"
#include "utils.h"
//#include "display.h"
#include "html_tags.h"
#include "tile_img.h"
#include "matched_template.h"


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
    uint32_t = img_w = 9000;
    uint32_t = img_h = 3000;


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




    //----------------------------------------------------------------------------
    // random cell location generation piece



    //----------------------------------------------------------------------------



    //----------------------------------------------------------------------------
    // report piece
    // directory to save report
    std::string report_dir = "../ report";


    std::ofstream out_file(report_dir + "/report.html");
    write_toplines(out_file);
    write_template_imgs(matched_templates, cell_list, templ_img_dir);

    // Write heading
    Html_top_heading top_heading;
    out_file << top_heading;



    std::cout << "Press Enter to close..." << std::endl;

    std::cin.ignore();

    cv::destroyAllWindows();

    return 0;
}
