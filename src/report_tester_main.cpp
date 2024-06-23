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
#include "test_functions.h"

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
    uint32_t tile_index = 0;

    uint64_t cell_w, cell_h;
    int32_t cells_max_x, cells_max_y;
    uint64_t overlap_x = 0, overlap_y = 0;
    uint32_t img_w = 9000;
    uint32_t img_h = 3000;


    std::vector<int32_t> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(8);

    cv::Mat tmp_img;

    //if (argc < 2)
    //{
    //    std::cout << "Enter the location of the image..." << std::endl;
    //    std::cin.ignore();
    //}


    //----------------------------------------------------------------------------
    // random image generation piece
    std::string cv_window = "Image";
    int bp = 0;
    
    cv::RNG rng(1234567890);
    cv::Mat current_tile_image = cv::Mat(img_h, img_w, CV_8UC1);
    rng.fill(current_tile_image, cv::RNG::UNIFORM, 0, 255);
    cv::cvtColor(current_tile_image, current_tile_image, cv::COLOR_GRAY2BGR);

    //----------------------------------------------------------------------------
    // setup cell parameters
    uint32_t min_cell_w = 30;
    uint32_t min_cell_h = 30;

    uint32_t cell1_w = 4* min_cell_w;
    uint32_t cell1_h = min_cell_h;
    uint32_t cell2_w = 2* min_cell_w;
    uint32_t cell2_h = min_cell_h;

    double pixel_size_x = 20.0; // from PixelSize nm/pixel in.txt (resolution)
    double pixel_size_y = 20.0;

    //----------------------------------------------------------------------------
    // random matched template generation piece


    cv::Mat mt1 = cv::Mat(cell1_h, cell1_w, CV_8UC1);
    rng.fill(mt1, cv::RNG::UNIFORM, 0, 2);
    //cv::cvtColor(mt1, mt1, cv::COLOR_GRAY2BGR);

    cv::Mat mt2 = cv::Mat(cell2_h, cell2_w, CV_8UC1);
    rng.fill(mt2, cv::RNG::UNIFORM, 0, 2);
    //cv::cvtColor(mt2, mt2, cv::COLOR_GRAY2BGR);

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

    uint32_t buffer_w = 60;
    uint32_t buffer_h = 30;

    uint32_t num_rows = floor((img_h - (2 * buffer_h)) / (double)min_cell_h);
    uint32_t num_cols = floor((img_w - (2 * buffer_w)) / (double)min_cell_w);

    uint32_t col_idx, row_idx;
    // rows
    for (idx = 0; idx < num_rows; ++idx)
    {
        row_idx = (idx * min_cell_h) + buffer_h;
        col_idx = buffer_w;

        while (col_idx < (img_w - buffer_w))
        {
            // randomly select one of the cells to put in into the list
            uint32_t cell_t = rng.uniform(0, 3);
            int32_t status = rng.uniform(0, 3);

            switch (cell_t)
            {
            case 0:
                ct1_pts.push_back(Pt_info(0, cv::Point(col_idx, row_idx), status));
                ct1_pts[ct1_pts.size()-1].on_tile = 0;
                ct1_pts[ct1_pts.size() - 1].defpt = cv::Point(col_idx* pixel_size_x, row_idx* pixel_size_y);
                col_idx += cell1_w;
                break;
            case 1:
                ct2_pts.push_back(Pt_info(0, cv::Point(col_idx, row_idx), status));
                ct2_pts[ct2_pts.size() - 1].on_tile = 0;
                ct2_pts[ct2_pts.size() - 1].defpt = cv::Point(col_idx * pixel_size_x, row_idx * pixel_size_y);
                col_idx += cell2_w;
                break;
            default:
                col_idx += min_cell_w*rng.uniform(0, 5);
                break;
            }
        }
    }

    // draw the locations of each of the cells for verification
    tmp_img = current_tile_image.clone();

    for (idx = 0; idx < ct1_pts.size(); ++idx)
    {
        cv::rectangle(tmp_img, cv::Rect(ct1_pts[idx].pt, cv::Size(cell1_w, cell1_h)), cv::Scalar(255, 0, 0), 2, 8);
    }

    for (idx = 0; idx < ct2_pts.size(); ++idx)
    {
        cv::rectangle(tmp_img, cv::Rect(ct2_pts[idx].pt, cv::Size(cell2_w, cell2_h)), cv::Scalar(0, 255, 0), 2, 8);
    }

    // view the results of the random selection
    cv::namedWindow(cv_window, cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
    cv::imshow(cv_window, tmp_img);
    cv::waitKey(0);

    //----------------------------------------------------------------------------
    // create the cell lists
    std::vector<Cell_info> cell_list;
    cell_list.push_back(Cell_info("cell_type_1", cell1_h, cell1_w, { 0,1,2,3 }, ct1_pts));
    cell_list.push_back(Cell_info("cell_type_2", cell2_h, cell2_w, { 0,1,2,3 }, ct2_pts));

    cell_list[0].selected = true;
    cell_list[0].height_nm = cell1_h * pixel_size_y;
    cell_list[0].width_nm = cell1_w * pixel_size_x;

    cell_list[1].selected = true;
    cell_list[1].height_nm = cell2_h * pixel_size_y;
    cell_list[1].width_nm = cell2_w * pixel_size_x;

    //----------------------------------------------------------------------------
    std::vector<Tile_img> tiles;
    bool whole = true;

    tiles.push_back(Tile_img(0, 0, img_h, img_w));
    tiles[0].aligned = true;
    tiles[0].binarized = true;
    tiles[0].pixel_size_x = pixel_size_x;
    tiles[0].pixel_size_y = pixel_size_y;

    //----------------------------------------------------------------------------
    // detects
    std::vector<std::vector<Pt_info>> detects(2);

    for (idx = 0; idx < ct1_pts.size(); ++idx)
    {
        if (ct1_pts[idx].matched == 0)
            detects[0].push_back(Pt_info(0, ct1_pts[idx].pt, -1));
    }

    for (idx = 0; idx < ct2_pts.size(); ++idx)
    {
        if (ct2_pts[idx].matched == 0)
            detects[1].push_back(Pt_info(0, ct2_pts[idx].pt, -1));
    }


    //----------------------------------------------------------------------------
    // conversion info
    Pixel_to_DEF_conversion def_conversion_info;
    cv::Mat tmp_tr = (cv::Mat_<double>(3, 3) << 1, 0, img_w>>1, 0, 1, img_h>>1, 0, 0, 1);
    def_conversion_info.transf_matrices.push_back(tmp_tr);  // cv::Mat::eye(3,3,CV_64FC1)


    //----------------------------------------------------------------------------
    // report piece
    // directory to save report
    std::string report_dir = "../report";

    // Creates data directory if it doesn't exist
    std::string html_img_dir = report_dir + "/html_images";
    std::string templ_img_dir = report_dir + "/template_images";
    std::string h_html_dir = "html_images";
    std::string h_templ_dir = "template_images";

    std::filesystem::create_directory(report_dir);
    std::filesystem::create_directory(html_img_dir);
    std::filesystem::create_directory(templ_img_dir);


    std::ofstream out_file(report_dir + "/report.html");
    write_toplines(out_file);
    write_template_imgs(matched_templates, cell_list, templ_img_dir);

    // Write heading
    Html_top_heading top_heading;
    out_file << top_heading;

    //----------------------------------------------------------------------------
    // start here!!!!!!
    std::vector<cv::Rect> report_grid;
    int32_t tile_w = 1800;
    int32_t tile_h = 1800;
    int32_t max_cell_width = 7 * min_cell_w;     // <-- this needs to be calculated based on the largest cell in cell_list
    int32_t max_cell_height = min_cell_h;     // <-- this needs to be calculated based on the largest cell in cell_list

    // only need to call this once and feed to the Html_img report piece
    report_grid_generation(img_w, img_h, tile_w, tile_h, max_cell_width, max_cell_height, report_grid);

    // temp code just to view the tile layout
    tmp_img = current_tile_image.clone();

    for (idx = 0; idx < report_grid.size(); ++idx)
    {
        cv::rectangle(tmp_img, report_grid[idx], cv::Scalar(10*rng.uniform(12,26), 10 * rng.uniform(12, 26), 10 * rng.uniform(12, 26)), 4, 8);
    }

    // view the results of the random selection
    cv::namedWindow(cv_window, cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
    cv::imshow(cv_window, tmp_img);
    cv::waitKey(0);

    // here's how to crop an image using the report grid rect
    cv::Mat cropped_img = current_tile_image(report_grid[0]).clone();

    cv::namedWindow("cropped_image", cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
    cv::imshow("cropped_image", cropped_img);
    cv::waitKey(0);

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
