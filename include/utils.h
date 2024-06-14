// This file was created by Victoria Lockridge, NSWC Crane
// October 2, 2023
#pragma once

//#include <ryml_all.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <map>
#include "cell_information.h"
#include "tile_img.h"

class Core_box
{
public:
    double ll_x = 0.0;
    double ur_x = 0.0;
    double ll_y = 0.0;
    double ur_y = 0.0;
    double area = 0.0; // nanometers
};

class Pixel_to_DEF_conversion
{
public:
    // Attributes
    double xcenter = 0; // nanometers
    double ycenter = 0; // nanometers
    cv::Point user_transl = cv::Point(0, 0); // pixels
    std::vector<cv::Mat> transf_matrices;
    // Note:  need pixel_size_x and pixel_size_y from tiles
};

//************************ PARSERS
//////////////////////////////////////////////////////
// Parses .yml file
//void parse_cell_input(std::string param_filename, uint32_t& tile_x, uint32_t& tile_y, uint32_t& tile_w, uint32_t& tile_h, std::vector<Cell_info>& cell_list)
//{
//  uint32_t h, w;
//  std::string name;
//  std::vector<uint16_t> orientation;
//  cell_list.clear();
//
//  try
//  {
//    std::ifstream tmp_stream(param_filename);
//    std::stringstream buffer;
//    buffer << tmp_stream.rdbuf();
//    std::string contents = buffer.str();
//    tmp_stream.close();
//
//    ryml::Tree config = ryml::parse_in_arena(ryml::to_csubstr(contents));
//
//    //config["scale"] >> scale;
//
//    ryml::NodeRef tiles = config["tiles"];
//    tiles["x"] >> tile_x;
//    tiles["y"] >> tile_y;
//    tiles["width"] >> tile_w;
//    tiles["height"] >> tile_h;
//
//    ryml::NodeRef cell_node = config["cell_list"];
//
//    for (ryml::ConstNodeRef n : cell_node.children())
//    {
//      n["name"] >> name;
//      n["height"] >> h;
//      n["width"] >> w;
//      //n["orientations"] >> orientation;
//      //Cell_info cell_tmp(name, std::floor(h / scale + 0.5), std::floor(w / scale + 0.5), orientation);
//      //Cell_info cell_tmp(name, h, w, orientation);
//      Cell_info cell_tmp(name, h, w);
//
//      cell_list.push_back(cell_tmp);
//    }
//  }
//  catch (std::exception& e)
//  {
//    std::cout << "Error: " << e.what() << std::endl;
//  }
//}   // end of parse_cell_input

//////////////////////////////////////////////////////////
// Finds number for symmetry
int symnum(std::string symmetry)
{
  std::map<std::string, int> m{{"N", 0}, {"FN", 1}, {"S", 2},{"FS", 3}};
  return m[symmetry];
}

/////////////////////////////////////////////////////////
// Finds letters for symmetry
std::string symname(int symnum)
{
  std::map<int, std::string> m{{0, "N"}, {1, "FN"}, {2, "S"},{3, "FS"}};
  return m[symnum];
}

//////////////////////////////////////////////////////////
// Gets number of lines in file
std::int32_t get_num_lines(std::string file_name)
{
  std::ifstream infile_test(file_name);
  std::string aline;
  std::int32_t ipr = 0;
  while (std::getline(infile_test, aline))
  {
    ipr++;
  }
  infile_test.close();
  return ipr - 1;
}

////////////////////////////////////////////////////////////
// Parse line of file
void parse_line(std::string line, std::vector<Cell_info>& cells, Core_box& corebox)
{
  // Alternate ways to parse
  // Get substring of each cell name in cell list and search def file over and over 
  // OR pull out cell name in def file and compare to selected list (USED THIS METHOD)
  
  // Get cellname
  auto const regex = std::regex("scc9gena_[a-z0-9_]+ ");
  auto searchResults = std::smatch{};
  bool const txtcontainsregex = std::regex_search(line, searchResults, regex);
  if (txtcontainsregex)
  {
    // Get matched string position in the line
    int64_t ridx = searchResults.prefix().length();
    // Pull substring of interest
    std::string rstring = line.substr(ridx, line.size());
    // To get the end of the cell name, find the first + after the "scc9" part
    std::size_t found = rstring.find_first_of("+");
    std::string cellnm = rstring.substr(0, found - 1);

    // Compare the found cell name in def file to the cell names in cell_list from .yml file
    for (int idx = 0; idx < cells.size(); ++idx)
    {
      // If found a cell name in cells, then parse out x,y,orientation information
      if (cellnm == cells[idx].name)
      {
        // Get x and y, find "PLACED ( [0-9] [0-9] ) [FNS]"
        auto const regex2 = std::regex(" \\( [0-9]+ [0-9]+ \\) [FNS]+");
        auto searchResults2 = std::smatch{};
        bool const txtcontainsregex2 = std::regex_search(rstring, searchResults2, regex2);

        if (txtcontainsregex2)
        {
          // Get matched string position in the line
          int64_t ridx2 = searchResults2.prefix().length();
          // Pull substring of interest
          std::string rstring2 = rstring.substr(ridx2, rstring.size());
          std::istringstream iss(rstring2);
          std::string word;
          std::vector<std::string> parts;
          while (std::getline(iss, word, ' '))
          {
            parts.push_back(word);
          }
          // Store cell information in the vector cells
          Pt_info tmp;
          tmp.pt.x = stoi(parts[2]);
          tmp.pt.y = stoi(parts[3]);
          tmp.defpt = tmp.pt;
          // Symmetry N, FN, S, FS to 0,1,2,or 3
          tmp.orientation = symnum(parts[5]);
          cells[idx].pt_locations.push_back(tmp);
        }

        break; //if cellnm found in suggested cells and processed, then break out of the for loop to next line of file
      }
    }// End of for loop that checks cell names in cells for a match
  }

  int index;
  // Check line for substrings related to FE_CORE_BOX_
  std::vector<std::string> cb_search{"FE_CORE_BOX_LL_X REAL ", "FE_CORE_BOX_UR_X REAL ", "FE_CORE_BOX_LL_Y REAL ", "FE_CORE_BOX_UR_Y REAL "};
  for (int idx = 0; idx < cb_search.size(); ++idx)
  {
    std::size_t foundcb = line.find(cb_search[idx]);
    if (foundcb != std::string::npos)
    {
      index = (int)(foundcb + cb_search[idx].size());
      std::istringstream line_iss(line.substr(index, line.length() - 1));
      std::string lword;
      std::vector<std::string> parts;
      while (std::getline(line_iss, lword, ' '))
      {
        parts.push_back(lword);
      }
      switch (idx)
      {
      case(0):
        corebox.ll_x = stod(parts[0]);
        break;
      case(1):
        corebox.ur_x = stod(parts[0]);
        break;
      case(2):
        corebox.ll_y = stod(parts[0]);
        break;
      case(3):
        corebox.ur_y = stod(parts[0]);
        break;
      }
      break; //after processing the found cb, break for loop (don't continue checking that line for cb_search)
    }
  }
}


//************************ ALIGNMENT
//////////////////////////////////////////////////////////
// Flip one cell across horizontal or vertical line
cv::Point flip_cell(int horiz_or_vert, cv::Point pt, double centerV)
{
  cv::Point refpt;
  if (horiz_or_vert == 0) //horizontal
  {
    refpt.x = pt.x;
    refpt.y = 2 * centerV - pt.y;
    //refpt.y = (refpr.y <= centerV) ? refpt.y + 2 * (centerV - refpt.y) : refpt.y - 2 * (refpt.y - centerV);
  }
  else if (horiz_or_vert == 1) //vertical
  {
    refpt.x = 2 * centerV - pt.x;
    refpt.y = pt.y;
  }
  return refpt;
}

//////////////////////////////////////////////////////////
// Flip one cell across horizontal or vertical line
cv::Point flip_cell(int horiz_or_vert, cv::Point2d pt, double centerV)
{
  cv::Point refpt;
  if (horiz_or_vert == 0) //horizontal
  {
    refpt.x = pt.x;
    refpt.y = 2.0 * centerV - pt.y;
    //refpt.y = (refpr.y <= centerV) ? refpt.y + 2 * (centerV - refpt.y) : refpt.y - 2 * (refpt.y - centerV);
  }
  else if (horiz_or_vert == 1) //vertical
  {
    refpt.x = 2.0 * centerV - pt.x;
    refpt.y = pt.y;
  }
  return refpt;
}

double flip_cell(int horiz_or_vert, double ypt, double centerV)
{
  double rypt;
  
  rypt = 2.0 * centerV - ypt;
  
  return rypt;
}

//////////////////////////////////////////////////////////
// Mirrors cells across horizontal and vertical center lines
// Converts cell width, height, and point locations to pixels from nanometers
void mirror_cell_locations(std::vector<Cell_info>& cell_list, double pixel_size_x, double pixel_size_y, double xcenter, double ycenter)
{
  for (Cell_info& cell : cell_list)
  {
    // Convert cell width and height from nanometers to pixels
    //cell.width = (uint32_t)floor((double)cell.width / tile0.pixel_size_x + 0.5);
    //cell.height = (uint32_t)floor((double)cell.height / tile0.pixel_size_y + 0.5);
    cell.width_nm = cell.width;
    cell.width = (uint32_t)ceil((double)cell.width / pixel_size_x);
    cell.height_nm = cell.height;
    cell.height = (uint32_t)ceil((double)cell.height / pixel_size_y);

    for (Pt_info& ptloc : cell.pt_locations)
    {
      // Move pt_locations(anchor) of all cells in cell_list to upper left corner of cell from upper right corner of cell
      ptloc.pt = flip_cell(0, ptloc.pt, ycenter); //Flip around horizontal center (nanometers)
      ptloc.pt = flip_cell(1, ptloc.pt, xcenter); //Flip around vertical center (nanometers)
      ptloc.pt.x = (int32_t)floor((double)ptloc.pt.x / pixel_size_x);  // Change nanometers to pixels
      ptloc.pt.y = (int32_t)floor((double)ptloc.pt.y / pixel_size_y); // Change nanometers to pixels
      // Move pt_locations(anchor) of all cells in cell_list to upper left corner of cell from upper right corner of cell
      ptloc.pt.x = ptloc.pt.x - (int32_t)cell.width;
    }
  }
}

////////////////////////////////////////////////////////////
//// Mirrors power/ground points across horizontal and vertical center lines
//// Converts point locations to pixels from nanometers
//void mirror_pwgnd_locations(std::vector<cv::Point2d>& points, double pixel_size_x, double pixel_size_y, double xcenter, double ycenter)
//{
//  for (cv::Point2d& pt : points)
//  {
//    pt = flip_cell(0, pt, ycenter); //Flip around horizontal center (nanometers)
//    pt = flip_cell(1, pt, xcenter); //Flip around vertical center (nanometers)
//    //pt.x = (int32_t)floor((double)pt.x / pixel_size_x);  // Change nanometers to pixels
//    //pt.y = (int32_t)floor((double)pt.y / pixel_size_y); // Change nanometers to pixels
//    pt.x = (double)pt.x / pixel_size_x;  // Change nanometers to pixels
//    pt.y = (double)pt.y / pixel_size_y; // Change nanometers to pixels
//  }
//}


void translate_cells(std::vector<Cell_info>& cell_list, cv::Point transl)
{
  for (Cell_info& cell : cell_list)
  {
    for (Pt_info& ptloc : cell.pt_locations)
    {
      ptloc.pt.x = ptloc.pt.x + transl.x;
      ptloc.pt.y = ptloc.pt.y + transl.y;
    }
  }
}

//void translate_pwgnd(std::vector<cv::Point2d>& points, cv::Point transl)
//{
//  for (cv::Point2d& pt : points)
//  {
//    pt.x = pt.x + transl.x;
//    pt.y = pt.y + transl.y;
//  }
//}



///////////////////////////////////////////////////////
// Calculates transformation matrix
void calculate_transformation_matrix(std::vector<cv::Point> a_matrix, std::vector<cv::Point> b_matrix, cv::Mat &M)
{
  // Find number of corners with points
  int count = 0;
  for (int i = 0; i < a_matrix.size(); i++)
  {
    if (a_matrix[i].x != 0)
      count++;
  }
  cv::Mat A = cv::Mat::ones(3, count, CV_64FC1);
  cv::Mat B = cv::Mat::ones(3, count, CV_64FC1);

  int j = 0;
  for (int i = 0; i < a_matrix.size(); i++)
  {
    if (a_matrix[i].x != 0)
    {
      A.at<double>(0, j) = a_matrix[i].x;
      A.at<double>(1, j) = a_matrix[i].y;
      B.at<double>(0, j) = b_matrix[i].x;
      B.at<double>(1, j) = b_matrix[i].y;
      j++;
    }

  }

  cv::Mat A_inv(A.cols, A.rows, CV_64FC1, cv::Scalar::all(1));
  cv::invert(A, A_inv, cv::DecompTypes::DECOMP_SVD);

  M = B * A_inv;

  // Debug test
  //cv::Mat Bnew(B.rows, B.cols, CV_64FC1);
  //Bnew = M * A;
}


//////////////////////////////////////////////////////
// Apply transformation matrix to all tile cells
//void apply_transformation(cv::Mat M, Tile_img& tile)
void apply_transformation_matrix(std::vector<Cell_info> &tile_cells, cv::Mat M)
{
  // Transformation, need 3xlength matrix (row 1: x values, row 2: y values, row 3: ones)
  // Create matrix for each cell in corner_pts and multiply
  // Then return corner_pts
  // Apply transformation to original corner points 
  //    and create new_corner_pts for display and to update b_matrix

  // Apply the latest transformation matrix M to all tile cells and save
  for (int i = 0; i < tile_cells.size(); i++)
  {
    for (int j = 0; j < tile_cells[i].pt_locations.size(); j++)
    {
      cv::Mat mtile = cv::Mat::ones(3, 1, CV_64FC1);
      mtile.at<double>(0, 0) = tile_cells[i].pt_locations[j].pt.x; 
      mtile.at<double>(1, 0) = tile_cells[i].pt_locations[j].pt.y;

      cv::Mat new_mtile = M * mtile;

      tile_cells[i].pt_locations[j].pt.x = new_mtile.at<double>(0, 0);  
      tile_cells[i].pt_locations[j].pt.y = new_mtile.at<double>(1, 0);
    }
  }
}   // end of apply_transformation_matrix

//void apply_transformation_matrix_pwgnd(std::vector<cv::Point2d> points, cv::Mat M)
//{
//  //Converted whole, not by tiles !!!!
//  for (int j = 0; j < points.size(); j++)
//  {
//    cv::Mat mtile = cv::Mat::ones(3, 1, CV_64FC1);
//    mtile.at<double>(0, 0) = points[j].x;
//    mtile.at<double>(1, 0) = points[j].y;
//
//    cv::Mat new_mtile = M * mtile;
//
//    points[j].x = new_mtile.at<double>(0, 0);
//    points[j].y = new_mtile.at<double>(1, 0);
//  }
//}

//////////////////////////////////////////////////////
// Apply translation to all tile cells when only one corner was used for alignment
void apply_transformation_1c(Tile_img& tile)
{
  // Do only translation if one corner
  // Find translation between a and b matrix
  std::int32_t transx = 0;
  std::int32_t transy = 0;

  for (int i = 0; i < tile.b_matrix.size(); i++)
  {
    if (tile.b_matrix[i].x > 0)
    {
      transx = tile.b_matrix[i].x - tile.a_matrix[i].x;
      transy = tile.b_matrix[i].y - tile.a_matrix[i].y;
    }
  }

  // Apply translation to all tile cells and save
  for (int i = 0; i < tile.tile_cells.size(); i++)
  {
    for (int j = 0; j < tile.tile_cells[i].pt_locations.size(); j++)
    {
      tile.tile_cells[i].pt_locations[j].pt.x += transx;
      tile.tile_cells[i].pt_locations[j].pt.y += transy;
    }
  }
}   // end of apply_transformation_matrix_1c

//////////////////////////////////////////////////////
// Apply translation to all tile cells when only two corners were used for alignment
// First finds average translation
void apply_transformation_2c(Tile_img& tile)
{
  // Do only translation if only 2 corners
  // Find average x and y translation between a and b matrix
  std::int32_t transx = 0;
  std::int32_t transy = 0;

  for (int i = 0; i < tile.b_matrix.size(); i++)
  {
    if (tile.b_matrix[i].x > 0)
    {
      transx = transx + (tile.b_matrix[i].x - tile.a_matrix[i].x);
      transy = transy + (tile.b_matrix[i].y - tile.a_matrix[i].y);
    }
  }

  transx = (std::int32_t)floor(transx / 2);
  transy = (std::int32_t)floor(transy / 2);

  // Apply the latest translation to all tile cells and save
  for (int i = 0; i < tile.tile_cells.size(); i++)
  {
    for (int j = 0; j < tile.tile_cells[i].pt_locations.size(); j++)
    {
      tile.tile_cells[i].pt_locations[j].pt.x += transx;
      tile.tile_cells[i].pt_locations[j].pt.y += transy;
    }
  }
}   // end of apply_transformation_matrix_1c


///////////////////////////////////////////////////////
// After applying transformation matrix to tile cells, this method updates cell in cell_list (global)
void update_cell_list(Tile_img tile, std::vector<Cell_info>& cell_list)
{
  for (Cell_info c : tile.tile_cells)
  {
    for (Pt_info ptloc : c.pt_locations)
    {
      // glob_indx[0] is cell index, glob_indx[1] is pt_location index in cell
      //if (!cell_list[ptloc.glob_indx[0]].pt_locations[ptloc.glob_indx[1]].glcell_aligned)
      //{
        cell_list[ptloc.glob_indx[0]].pt_locations[ptloc.glob_indx[1]].pt = ptloc.pt;
        //cell_list[ptloc.glob_indx[0]].pt_locations[ptloc.glob_indx[1]].glcell_aligned = true;
        //ptloc.glcell_aligned = true;
      //}
    }
  }
}


//************************ HELPFUL TOOLS
///////////////////////////////////////////////////////
// 
bool is_still_on_tile(int32_t tilex, int32_t tiley, int32_t tile_width, int32_t tile_height, Pt_info ptloc)
{
  bool is_on_tile = false;
  if ((ptloc.pt.x > tilex) & (ptloc.pt.x < tilex + tile_width))
  {
    if ((ptloc.pt.y > tiley) & (ptloc.pt.y < tiley + tile_height))
      is_on_tile = true;
  }
    
  return is_on_tile;
}


///////////////////////////////////////////////////////
// Assign tile indices to point locations on cell_list
void assign_tile_idx_to_global_cells(std::vector<Cell_info>& tile_cells, std::vector<Cell_info>& cell_list, std::int32_t idxtile, int32_t tilex, int32_t tiley, int32_t tile_width, int32_t tile_height)
{
  // Check if cells are still on tile after alignment
  std::int32_t cidx, pidx;

  for (int i = 0; i < tile_cells.size(); i++)
  {
    for (int j = 0; j < tile_cells[i].pt_locations.size(); j++)
    {
      if (is_still_on_tile(tilex, tiley, tile_width, tile_height, tile_cells[i].pt_locations[j]))
      {
        cidx = tile_cells[i].pt_locations[j].glob_indx[0];
        pidx = tile_cells[i].pt_locations[j].glob_indx[1];
        cell_list[cidx].pt_locations[pidx].on_tile = idxtile;
      }
      else
      {
        tile_cells[i].pt_locations.erase(tile_cells[i].pt_locations.begin() + j);
        j = j - 1;
      }
    }
  }
}


///////////////////////////////////////////////////////
// Find widest cell in cell_list
int32_t get_widest_cell_width(const std::vector<Cell_info>& cells)
{
  int32_t cell_w = 0;
  for (Cell_info c : cells)
    cell_w = (cell_w > (int32_t)c.width ? cell_w : (int32_t)c.width);

  return cell_w;
}

///////////////////////////////////////////////////////
// Find tallest cell in cell_list
int32_t get_tallest_cell_height(const std::vector<Cell_info>& cells)
{
  int32_t cell_h = 0;
  for (Cell_info c : cells)
    cell_h = (cell_h > (int32_t)c.height ? cell_h : (int32_t)c.height);

  return cell_h;
}


//************************ FOR REPORT
///////////////////////////////////////////////////////
// For report, get all cell names that have been processed
std::vector<std::string> get_processed_cell_names(const std::vector<Cell_info> &cells)
{
  std::vector<std::string> names;
  for (int i = 0; i < cells.size(); i++)
  {
    if (cells[i].selected == true)
      names.push_back(cells[i].name);
  }
  return names;
};


///////////////////////////////////////////////////////
// For report, get all image file names that have been processed
std::vector<std::string> get_processed_image_files(const std::vector<Tile_img> &tiles, bool whole_img)
{
  std::vector<std::string> files;
  for (int i = 0; i < tiles.size(); i++)
  {
    if (tiles[i].scored)
    {
      std::int32_t last = (std::int32_t)tiles[i].image_file.length();
      std::size_t found = tiles[i].image_file.find_last_of("/\\");
      std::string tmp = tiles[i].image_file.substr(found+1, last);
      // Add image global pixel range to string
      if (!whole_img)
        tmp = tmp + " - pixel location on image is (" + std::to_string(tiles[i].imgLocColX) + ":" + std::to_string(tiles[i].imgLocColX + tiles[i].width) + ", " + std::to_string(tiles[i].imgLocRowY) + ":" + std::to_string(tiles[i].imgLocRowY + tiles[i].height) + ")";
      files.push_back(tmp);
    } 
  }
  return files;
};


///////////////////////////////////////////////////////
// Remove cells that are cutoff by right or lower edge of image
void remove_partial_edge_cells(std::vector<Cell_info> &tile_cells, int32_t tile_width, int32_t tile_height, int32_t tile_imgX, int32_t tile_imgY)
{
  // For each of the 10 types of cells in tile_cells, check the pt_locations to see if they are cutoff by the right or lower edge of the tile image
  for (int i = 0; i < tile_cells.size(); i++)
  {
    // Check right edge
    for (int j = 0; j < tile_cells[i].pt_locations.size(); j++)
      if (tile_cells[i].pt_locations[j].pt.x > (tile_imgX + (tile_width - tile_cells[i].width)))
      {
        tile_cells[i].pt_locations.erase(tile_cells[i].pt_locations.begin() + j);
        j = j - 1;
      } 

    // Of those left, check lower edge
    for (int j = 0; j < tile_cells[i].pt_locations.size(); j++)
      if (tile_cells[i].pt_locations[j].pt.y > (tile_imgY + (tile_height - tile_cells[i].height)))
      {
        tile_cells[i].pt_locations.erase(tile_cells[i].pt_locations.begin() + j);
        j = j - 1;
      }
  }
}


////////////////////////////////////////////////////////
// Translate pixel coordinates to DEF file coordinates in nanometers - for report of errant cells
// Steps:
// Remove affine transformation performed during alignment???
// Remove translation of overall image points by user 
// Remove scaling
// Remove mirroring
cv::Point convert_pixel_to_DEF_pt(Pixel_to_DEF_conversion conv_info, cv::Point pixel_pt, int32_t cell_width, int32_t tile_idx, double pixel_size_x, double pixel_size_y, bool errant)
{
  // Reverse process to get pixels from DEF points in nanometers:
  // Reverse affine transformation (in pixels)
  // Add back user translation (in pixels)
  // Add back cell width to x dimension (in pixels)
  // Convert pixels to nanometers using scale
  // Flip points across center of corebox

  cv::Point new_def_pt;

  // Transformation matrix
  // M * A_orig = pixel_pt_matrix_B;
  // A_orig = M_inverse * pixel_pt_B;
  cv::Mat M = conv_info.transf_matrices[tile_idx].clone();
  cv::Mat M_inv(M.cols, M.rows, CV_64FC1, cv::Scalar::all(1));
  cv::invert(M, M_inv, cv::DecompTypes::DECOMP_SVD);
  cv::Mat TestM = M_inv * M;

  cv::Mat B = cv::Mat::ones(3, 1, CV_64FC1);
  B.at<double>(0, 0) = pixel_pt.x;
  B.at<double>(1, 0) = pixel_pt.y;
  cv::Mat A;
  A = M_inv * B;
  new_def_pt.x = A.at<double>(0, 0);
  new_def_pt.y = A.at<double>(1, 0);

  // Remove user translation to point conv_info.user_transl
  new_def_pt = new_def_pt - conv_info.user_transl;

  // Move point x to right by width of cell
  new_def_pt.x = new_def_pt.x + cell_width;
  
  // Convert from pixels to nanometers
  // double pixel_size_x nm/pixel (resolution)
  new_def_pt.x = new_def_pt.x * pixel_size_x;
  new_def_pt.y = new_def_pt.y * pixel_size_y;

  //Flip cells across corebox centers 
  new_def_pt = flip_cell(1, new_def_pt, conv_info.xcenter); //vertical flip
  new_def_pt = flip_cell(0, new_def_pt, conv_info.ycenter); //horizontal filp

  return new_def_pt;
}


////////////////////////////////////////////////////////
// DEBUG
void debug_calculate_DEFcoords_cells(std::vector<Cell_info> cell_list, Pixel_to_DEF_conversion def_conversion_info, std::vector<Tile_img> tiles )
{
    // DEBUG - write out all cells on tiles that have been aligned
    // Test for aligned tiles, then test if cells in cell_list are on this tile
    //  Then write cell, coords, and calculate def coords - write to csv file

  std::ofstream ofile("C:/Data/Trust_ME/Whole Image/cell_data/location_all_cells_01292024.csv");
  for (int i = 0; i < tiles.size(); i++)
  {
    if (tiles[i].aligned)
    {
      for (Cell_info cell : cell_list)
      {
        for (Pt_info ptloc : cell.pt_locations)
        {
          if (ptloc.on_tile == i)
          {
            // convert_pixel_to_DEF_pt(Pixel_to_DEF_conversion conv_info, cv::Point pixel_pt, int32_t cell_width, int32_t tile_idx, double pixel_size_x, double pixel_size_y, bool errant)
            cv::Point test_pt;
            test_pt = convert_pixel_to_DEF_pt(def_conversion_info, ptloc.pt, cell.width, i, tiles[i].pixel_size_x, tiles[i].pixel_size_y, false);
            ofile << cell.name << "," << ptloc.pt.x << "," << ptloc.pt.y << "," << ptloc.defpt.x << "," << ptloc.defpt.y << "," <<  test_pt.x << "," << test_pt.y << std::endl;
          }
        }
      }
    }
  }
}


//************************ POWER/GROUND WORK
////////////////////////////////////////////////////////
// Parsing def file for power and ground locations

std::string word_search(bool power)
{
  std::map<int, std::string> m{{0, "vgnd"}, {1, "vpwr"}};
  return m[power];
}


void get_pwgnd_locations(bool power, std::string def_file, std::vector<cv::Point2d> &pin_locations)
{
  std::ifstream infile(def_file);
  std::string line;
  bool quit = false;

  std::string search_str = word_search(power);
  
  while (std::getline(infile, line))
  {
    if (line.find(search_str, 0) != std::string::npos)
      while (std::getline(infile, line)) // get the lines after it
      {
        if (line.find(";", 0) != std::string::npos)
        {
          quit = true;
          break;
        }

        else
        {
          // Break line up by spaces
          std::istringstream iss(line);
          std::string word;
          std::vector<std::string> parts;
          while (std::getline(iss, word, ' '))
          {
            parts.push_back(word);
          }

          // Find the words to save in Power_data format
          cv::Point2d tmp_pt;
          for (int i = 0; i < parts.size(); i++)
          {
            //tmp_pd.pwer = saved_str(power);
            if (parts[i] == "(")
            {
              int new_idx = i + 1;
              int new_idx2 = i + 2;
              tmp_pt.x = stod(parts[new_idx]);
              tmp_pt.y = stod(parts[new_idx2]);
              pin_locations.push_back(tmp_pt);
              break; // break from for loop that is evaluating one line of parts
            }
          }
        }
      }
    if (quit)
      break;
  }
  infile.close();
}

////////////////////////////////////////////////////////
// Pull power/ground info out of core.final.def
void parse_def_power_ground(std::vector<cv::Point2d> &power_pins, std::vector<cv::Point2d> &ground_pins, std::string def_file)
{
  get_pwgnd_locations(0, def_file, ground_pins);
  get_pwgnd_locations(1, def_file, power_pins);
}

//////////////////////////////////////////////////////////////////////////////
// Find pixel values of unique power/ground y values using cell_list def y locations - average value of pixel values in cell_list
void find_pixel_vals_pwrgnd(const std::vector<Cell_info>& cell_list, std::vector<double> pins_defy, std::vector<double>& pixels_y, int img_width)
{
  pixels_y.clear();
  double avg;
  int img_start = 500;
  for (int i = 0; i < pins_defy.size(); i++)
  {
    double pwr_sum = 0;
    int32_t count_pwr = 0;
    std::vector<double> pwr_pts;
    for (Cell_info cell : cell_list)
    {
      for (Pt_info& ptloc : cell.pt_locations)
      {
        if (pins_defy[i] == ptloc.defpt.y) // find cells with the same def y location
          if ((ptloc.pt.x > img_start) &(ptloc.pt.x < img_width))  
          {
            pwr_sum = pwr_sum + ptloc.pt.y;
            pwr_pts.push_back(ptloc.pt.y);
            count_pwr++;
            if (count_pwr > 50)
              break;
          }
      }
      if (count_pwr > 50)
        break;
    }
    if (count_pwr > 0)
      avg = pwr_sum / (double)pwr_pts.size();
    else
      avg = 0;
    pixels_y.push_back(avg);
  }

}


////////////////////////////////////////////////////////
// Find the unique y values of the power and ground pins
void unique_y_pins(std::vector<cv::Point2d>& pins, std::vector<double>& ypoints)
{
  for (cv::Point2d pin : pins)
  {
    ypoints.push_back(pin.y);
  }
  // Sort ypoints
  std::sort(ypoints.begin(), ypoints.end());
  // Remove consective, adjacent duplicates
  auto last = std::unique(ypoints.begin(), ypoints.end());
  ypoints.erase(last, ypoints.end());
}

void closer_to_power_ground(cv::Point cell_location, std::vector<double>& power_pins_y, std::vector<double>& ground_pins_y, bool& power)//, double& dist)
{
  double max_pwr_pin = *std::max_element(power_pins_y.begin(), power_pins_y.end());
  double min_pwr_pin = *std::min_element(power_pins_y.begin(), power_pins_y.end());
  double max_gnd_pin = *std::max_element(ground_pins_y.begin(), ground_pins_y.end());
  double min_gnd_pin = *std::min_element(ground_pins_y.begin(), ground_pins_y.end());

  double min_pwr_dist = max_pwr_pin - min_pwr_pin;
  double min_gnd_dist = max_gnd_pin - min_gnd_pin;

  for (double pw_y : power_pins_y)
  {
    if (abs(cell_location.y - pw_y) < min_pwr_dist)
      min_pwr_dist = abs(cell_location.y - pw_y);
  }

  for (double gnd_y : ground_pins_y)
  {
    if (abs(cell_location.y - gnd_y) < min_gnd_dist)
      min_gnd_dist = abs(cell_location.y - gnd_y);
  }

  power = (min_pwr_dist < min_gnd_dist) ? 1:0;
    
}

void debug_cells_close_to_power(std::vector<Cell_info> cell_list)
{
  std::string filename;
  filename = "C:/Data/Trust_ME/mydata/cpp_parsed/parse_pwgnd.csv";
  std::ofstream ofile(filename);

  for (Cell_info cell : cell_list)
    for (Pt_info ptloc : cell.pt_locations)
      //ofile << cell.name << "," << ptloc.pt.x << "," << ptloc.pt.y << "," << ptloc.orientation << "," << ptloc.power << "," <<  ptloc.pwr_distance << std::endl;
      ofile << cell.name << "," << ptloc.pt.x << "," << ptloc.pt.y << "," << ptloc.orientation << "," << ptloc.power << ","  << std::endl;
}

////////////////////////////////////////////////////////
// DEBUG
void debug_cells_pwgnd_yvalues(std::vector<cv::Point2d> points, std::vector<cv::Point2d> points_def, bool power, Pixel_to_DEF_conversion conv_info, double pixel_size_x, double pixel_size_y)
{
  // DEBUG - write out all cells on tiles that have been aligned
  // Test for aligned tiles, then test if cells in cell_list are on this tile
  //  Then write cell, coords, and calculate def coords - write to csv file

  //std::ofstream ofile("C:/Data/Trust_ME/mydata/cpp_parsed/whole_pwrgnd_pixels.csv");
  std::string pwgnd;
  std::string filename;
  if (power)
  {
    pwgnd = "power";
    filename = "C:/Data/Trust_ME/mydata/cpp_parsed/def_conv/whole_pwrg_pixels_01092024.csv";
  }
    
  else
  {
    pwgnd = "ground";
    filename = "C:/Data/Trust_ME/mydata/cpp_parsed/def_conv/whole_gnd_pixels_01092024.csv";
  }
    
  std::ofstream ofile(filename);
  for (int i = 0; i < points.size(); i++)
  {
    cv::Point tmp = cv::Point((int32_t)points[i].x, (int32_t)points[i].y);
    cv::Point defpt = convert_pixel_to_DEF_pt(conv_info, tmp, 0, 0, pixel_size_x, pixel_size_y, false);
    ofile << pwgnd << "," << points[i].x << "," << points[i].y << "," << points_def[i].x << "," << points_def[i].y << "," << defpt.x << "," << defpt.y << std::endl;
  }
}


//************************ SELECTED CELLS
///////////////////////////////////////////////
// 
//void set_selected_switch_on_cells(QModelIndexList selected_cells, std::vector<Cell_info>& cell_list)
//{
//  for (int indx = 0; indx < cell_list.size(); indx++)
//  {
//    // Clear out previous data - reinitialize "selected"
//    cell_list[indx].selected = 0;
//
//    if (selected_cells.size() == 0)
//      cell_list[indx].selected = 1;
//    else if (selected_cells.size() > 0)
//    {
//      for (int i = 0; i < selected_cells.size(); i++)
//      {
//        if (indx == (int)selected_cells[i].row())
//        {
//          cell_list[indx].selected = 1;
//          break;
//        }
//      }
//    }
//  }
//}


//************************ DISTINGUISHING SIMILAR CELLS WHEN CHECKING ERRANT CELLS
///////////////////////////////////////////////
// ASSUMPTION:  All North orientations are close to ground, FN: close to ground, S: close to power, FS: close to power
bool cell_type_member(uint16_t orientation, bool power)
{
  if (power)
  {
    std::map<uint16_t, bool> m{{0, false}, {1, false}, {2, true},{3, true}};
    return m[orientation];
  }
  else
  {
    std::map<uint16_t, bool> m{{0, true}, {1, true}, {2, false},{3, false}};
    return m[orientation];
  }
}

//////////////////////////////////////////////////////////////////////////////
// Check errant detects as last step in score_detects() function
// For each cell type in detects, for each detection in celltype, if matched == -1 (errant),
//    for the orientation, is the power/ground closeness correct
//             yes:  errant cell of that cell type
//             no:  not that cell type
//                  Is there another celltype with that same width?
//                        no:  unknown errant cell
//                        yes:  check the cell types of same width for cells located at same location as errant cell
//                              Same location?
//                                Was this cell matched?
//                                    no:  it's an errant cell
//                                   yes:  not an errant cell (it's a cell of cell type with same width)
void check_errant_detects_for_similar_cell(std::vector<Cell_info> cell_list, std::vector<std::vector<Pt_info>>& detects, std::vector<double> pwr_pixels_y, std::vector<double> gnd_pixels_y)
{
  int count_errant_members = 0;
  int count_located_but_missed = 0;
  int break_2loops = 0;
  bool power = 0;
  for (int cnum = 0; cnum < detects.size(); cnum++)
  {
    count_errant_members = 0;
    count_located_but_missed = 0;
    for (int i = 0; i < detects[cnum].size(); i++)
    {
      if (detects[cnum][i].matched == -1) // errant detection
      {
        // Check power/ground
        closer_to_power_ground(detects[cnum][i].pt, pwr_pixels_y, gnd_pixels_y, power);
        detects[cnum][i].power = power;
        bool member;
        member = cell_type_member(detects[cnum][i].orientation, power);
        if (member)
        {
          count_errant_members++;
          // Errant cell of this cell type 
          break; // Go to next detect
        }
          
        else // not a member
        {
          break_2loops = 0;
          for (int k = 0; k < cell_list.size(); k++)
          {
            
            if ((cnum != k) & (cell_list[k].width_nm == cell_list[cnum].width_nm))
            {
              for (int j = 0; j < cell_list[k].pt_locations.size(); j++)
              {
                
                // Check if there is a cell located at the same location as detects[cnum][i] with distance of 10 pixels in x and in y
                // TODO: remove hard coded 10 in favor of a dynamic value that tracks with image/cell size
                if ((detects[cnum][i].pt.x > cell_list[k].pt_locations[j].pt.x - 10) & (detects[cnum][i].pt.x < cell_list[k].pt_locations[j].pt.x + 10))
                {
                  if ((detects[cnum][i].pt.y > cell_list[k].pt_locations[j].pt.y - 10) & (detects[cnum][i].pt.y < cell_list[k].pt_locations[j].pt.y + 10))
                  {
                    if (cell_list[k].pt_locations[j].matched == 1)  
                    {
                      // remove this detect
                      detects[cnum].erase(detects[cnum].begin() + i);
                      i = i - 1; // fix vector location
                      // break out of 2 for loops
                      break_2loops = 1;
                      break;
                    }
                    else if (cell_list[k].pt_locations[j].matched == 2)// If it was missed?????
                    {
                      // still an errant cell - do what?
                      // Found cell that should be there but it wasn't matched with this cell type so leave it errant
                      count_located_but_missed++;
                      break_2loops = 1;
                      break;
                    }
                  }
                }
              }
            }
            if (break_2loops == 1)
              break;
          }
        }
      }
    }
  }
}


void find_max_cell_dims(std::vector<Cell_info> cell_list, uint32_t& max_width, uint32_t& max_height)
{
  for (int i = 0; i < cell_list.size(); i++)
  {
    if (cell_list[i].width > max_width)
      max_width = cell_list[i].width;
    
    if (cell_list[i].height > max_height)
      max_height = cell_list[i].height;
  }
}


void find_cellrow_yvalues(std::vector<Cell_info> cell_list, std::vector<int>& row_yvalues, std::vector<std::vector<int>>& x_values)
{
  std::vector<cv::Point> all_points;
  // Get unique y values in cell_list in pixels
  for (int i = 0; i < cell_list.size(); i++)
  {
    for (int j = 0; j < cell_list[i].pt_locations.size(); j++)
    {
      row_yvalues.push_back(cell_list[i].pt_locations[j].pt.y);
      all_points.push_back(cell_list[i].pt_locations[j].pt);
    }
      

    // Sort row_yvalues
    std::sort(row_yvalues.begin(), row_yvalues.end());
    // Remove consective, adjacent duplicates
    auto last = std::unique(row_yvalues.begin(), row_yvalues.end());
    row_yvalues.erase(last, row_yvalues.end());
  }

  // Get range of x values with that y pixel value
  //std::vector<std::vector<int>> x_values;
  for (int yvalue : row_yvalues)
  {
    std::vector<int> xvals;
    for (int i = 0; i < all_points.size(); i++)
    {
      if (yvalue == all_points[i].y)
        xvals.push_back(all_points[i].x);
    }
    x_values.push_back(xvals);
  }
  

}

////////////////////////////////////////////////////////////////////////////////
// // Function finds closest cell row - used to see if cell is in a vertical cell row - not misaligned
//ASSUMPTION:  CELLS ARE VERTICALLY ALIGNED IN GRIDDED ROWS(in y - direction)
int find_closest_row_diff(cv::Point detect_pt, std::vector<int>& row_yvalues, std::vector<std::vector<int>>& x_values, uint32_t max_cell_height)
{
  int diff_y = 2*max_cell_height;

  for (int j = 0; j < x_values.size(); j++)
  {
    // Get min and max x for this row j - due to alignment issues, need to find y values based on x values
    std::vector<int>::iterator resultn = std::min_element(x_values[j].begin(), x_values[j].end());
    int minx = *resultn;
    std::vector<int>::iterator resultx = std::max_element(x_values[j].begin(), x_values[j].end());
    int maxx = *resultx;
    if ((detect_pt.x <= maxx) &(detect_pt.x >= minx))
    {
      if (abs(row_yvalues[j] - detect_pt.y) < diff_y)
      {
        diff_y = abs(row_yvalues[j] - detect_pt.y);
      }
    }
  }
  /*for (int i = 0; i < row_yvalues.size(); i++)
  {
    if (abs(row_yvalues[i] - detect_pt.y) < diff_y)
    {
      diff_y = abs(row_yvalues[i] - detect_pt.y);
    }
  }*/

  return diff_y;
}

struct Check_cell
{
  cv::Point pt;
  int8_t matched;
  uint32_t cell_width;
};

/////////////////////////////////////////////////////////
// Test if errant cell is actually 2 or more cells
// Only check cells in cell_list in same row (use y-value)
void check_cell_overlap(std::vector<Cell_info> cell_list, cv::Point detect_pt, uint32_t detect_width, int cnum, uint32_t max_cell_width, bool& remove)
{
  std::vector<Check_cell> test_cells;

  // Find all cells 9called test_cells) that are in the same row and with x value range AND were matched
  for (int cidx = 0; cidx < cell_list.size(); cidx++)
  {
    for (int pidx = 0; pidx < cell_list[cidx].pt_locations.size(); pidx++)
    {
      //if ((cnum != cidx) & (cell_list[cidx].pt_locations[pidx].matched == 1)) // Doesn't work for fill_16 cells
      if (cell_list[cidx].pt_locations[pidx].matched == 1)
      {
        int test_val = 10; // in pixels
        if ((cell_list[cidx].pt_locations[pidx].pt.y >= detect_pt.y - test_val) & (cell_list[cidx].pt_locations[pidx].pt.y <= detect_pt.y + test_val))
          if ((cell_list[cidx].pt_locations[pidx].pt.x >= detect_pt.x - max_cell_width) & (cell_list[cidx].pt_locations[pidx].pt.x <= detect_pt.x + max_cell_width))
          {
            Check_cell tmp;
            tmp.cell_width = cell_list[cidx].width;
            tmp.pt = cell_list[cidx].pt_locations[pidx].pt;
            tmp.matched = cell_list[cidx].pt_locations[pidx].matched;
            // Save cells to be checked for overlaps
            test_cells.push_back(tmp);
          }
      }
    }
  }

  // Check cells for overlaps on a matched cell (matched == 1)
  // Detect cell overlaps the cell if they are in same y row and if any part of cell (x to x + cell_width) falls between detect x and detect x + detect's cell_width
  for (int i = 0; i < test_cells.size(); i++)
  {
    int min_x = test_cells[i].pt.x;
    int max_x = test_cells[i].pt.x + test_cells[i].cell_width;
    if (((detect_pt.x >= min_x) & (detect_pt.x <= max_x)) |
      ((detect_pt.x + detect_width >= min_x) & (detect_pt.x + detect_width <= max_x)))
    {
        remove = true;
        break;
    }
  }
}


/////////////////////////////////////////////////////////
// 1.  Checks if errant cell is actually on a cell row - if not, not an errant cell
// 2.  Test if errant cell is actually 2 or more cells - if errant cell overlaps 
//     defined cells and these cells were matched (matched == 1), then not errant
//     Only check for overlapping cells in cell_list in same row (use y-value)
// NOTE:  If defined cell was not matched (matched == 2), then errant cell is still errant.
void check_errant_detects_against_matched_cells(std::vector<Cell_info> cell_list, std::vector<std::vector<Pt_info>>& detects)
{
  // Find max cell width of all cells in cell_list
  uint32_t max_cell_width = 1;
  uint32_t max_cell_height = 1;
  find_max_cell_dims(cell_list, max_cell_width, max_cell_height);

  std::vector<int> row_yvalues;
  std::vector<std::vector<int>> x_values;
  find_cellrow_yvalues(cell_list, row_yvalues, x_values);

  for (int cnum = 0; cnum < detects.size(); cnum++)
  {
    for (int i = 0; i < detects[cnum].size(); i++)
    {
      if (detects[cnum][i].matched == -1) // errant cell detection
      {
        // ASSUMPTION:  CELLS ARE VERTICALLY ALIGNED IN GRIDDED ROWS (in y-direction)
        // Test if detect's y-value is near cell rows in cell_list 
        // Near - continue testing
        // Far - not errant, remove from detect list

        // Find closest row
        int diff_y;
        diff_y = find_closest_row_diff(detects[cnum][i].pt, row_yvalues, x_values, max_cell_height);
        int test_val = (int)(0.10 * max_cell_height);  //% of cell height
        if (diff_y > test_val) // Cell is not vertically aligned with cell rows
        {
          // Not errant, delete from detects list
          detects[cnum].erase(detects[cnum].begin() + i);
          i = i - 1; // fix vector location
        }
        else
        {
          // Test if errant cell part of another cell
          // Only check cells in cell_list in same row (use y-value)
          bool remove = false;
          check_cell_overlap(cell_list, detects[cnum][i].pt, cell_list[cnum].width, cnum, max_cell_width, remove);
          if (remove)
          {
            // Not errant, delete from detects list
            detects[cnum].erase(detects[cnum].begin() + i);
            i = i - 1; // fix vector location
          }
        }
      }
    }
  }
}


/*
//////////////////////////////////////////
// Get cell template of one cell type
void find_one_matched_template(std::vector<Tile_img>& tiles, matched_template &mtemplate, int celltype_idx, cv::Mat src = cv::Mat(), uint32_t num_samples = 7)
{
  // This work is mostly Dave's work in utilities.h, find_cell_templates()

  uint32_t num_tiles = tiles.size();
  uint32_t num_cells = tiles[0].tile_cells.size();
  std::vector<size_t> count(num_tiles);
  //std::vector<uint32_t> cell_type_index;

  cv::Mat img, full_image;
  cv::Mat binary_img;
  cv::Mat mt;

  cv::Rect r;

  std::string win_name = "Cell Type";
  std::string cell_name;

  int tile_max_idx;

  // Step 1: go through the tiles and cells and find the tile that has the max number of cells for a given cell type
  if (num_tiles > 1)
  {
    // cycle through each tile
    for (int jdx = 0; jdx < num_tiles; ++jdx)
    {
      // get the number of cells for each tile
      if (tiles[jdx].aligned & tiles[jdx].binarized)
        count[jdx] = tiles[jdx].tile_cells[celltype_idx].pt_locations.size();
      else
        count[jdx] = 0;
    }

    // get the iterator to the max value and store the index of the tile that has the most cells
    auto result = std::max_element(count.begin(), count.end());
    tile_max_idx = (int)(std::distance(count.begin(), result)); 
    //cell_type_index.push_back(std::distance(count.begin(), result));
    //}
  }
  else
    tile_max_idx = 0;
     
  // Step 2: get the matched template for cell type if selected
  cv::namedWindow(win_name, cv::WINDOW_GUI_EXPANDED);
  cv::setWindowProperty(win_name, cv::WND_PROP_ASPECT_RATIO, cv::WINDOW_KEEPRATIO);

  // convert the image to a gray scale and 32-bit float
  if (!src.empty())
  {
    if (src.channels() > 1)
    {
      cv::cvtColor(src, full_image, cv::COLOR_BGR2GRAY);
      full_image.convertTo(full_image, CV_32FC1);
    }
    else
    {
      src.convertTo(full_image, CV_32FC1);
    }
  }

  // read in image or use the dimensions from tiles to get the right image
  if (src.empty())
  {
    img = cv::imread(tiles[tile_max_idx].image_file, cv::IMREAD_GRAYSCALE | cv::IMREAD_ANYDEPTH);
    img.convertTo(img, CV_32FC1);
  }
  else
  {
    r = cv::Rect(tiles[tile_max_idx].imgLocColX, tiles[tile_max_idx].imgLocRowY, tiles[tile_max_idx].width, tiles[tile_max_idx].height);
    img = full_image(r);
  }

  // binarize image
  binarize_image(img, binary_img, tiles[tile_max_idx].binary_threshold, tiles[tile_max_idx].struct_elem_size);

  // get the matching template
  get_matching_template_(binary_img, tiles[tile_max_idx].tile_cells[celltype_idx], tiles[tile_max_idx], mtemplate, num_samples);

  cell_name = "Cell: " + tiles[tile_max_idx].tile_cells[tile_max_idx].name;

  cv::setWindowTitle(win_name, cell_name);
  //if (!matched_templates[idx].img.empty())
  //{
  cv::Mat tmp_mosiac = mtemplate.get_mosaic();
  cv::imshow(win_name, tmp_mosiac);
  cv::waitKey(500);

  cv::destroyWindow(win_name);
}
*/


void setup_cell_list(std::set<std::string> uniq_cellnames, std::vector<Cell_info>& cell_list, std::vector<std::string> nonf_cells)
{
  cell_list.clear();
  for (std::string cellname: uniq_cellnames)
  {
    if (!nonf_cells.empty() & (std::find(nonf_cells.begin(), nonf_cells.end(), cellname) != nonf_cells.end()))
      continue;
    else
    {
      Cell_info tmp_cell(cellname);
      cell_list.push_back(tmp_cell);
    }
  }

}

void add_cell_width_height(std::vector<Cell_info>& cell_list, const std::vector<uint32_t>& cell_width, const std::vector<uint32_t>& cell_height)
{
  for (int i = 0; i < cell_list.size(); i++)
  {
    cell_list[i].width = cell_width[i];
    cell_list[i].height = cell_height[i];
  }
}














