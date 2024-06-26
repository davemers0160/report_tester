// This file was created by Victoria Lockridge, NSWC Crane
// September 5, 2023
#include "tile_img.h"
#include <random>

//////////////////////////////////////////////////////////
void onMouse2(int evt, int x, int y, int flags, void* param)
{
  if (evt == cv::EVENT_LBUTTONDOWN)
  {
    std::vector<cv::Point>* ptPtr = (std::vector<cv::Point>*)param;
    ptPtr->push_back(cv::Point(x, y));
  }
}

void Tile_img::find_pixel_size(const cv::Mat &img, std::string file_name, int32_t core_box_dx, int32_t core_box_dy)
{
  cv::Mat climg = img.clone();
  // Size of image to display to see UL and LR corners
  const int32_t simg_width = 3000;
  const int32_t simg_height = 2000;
  
  // Display upper left corner and prompt user to select point
  std::string windowname("Select upper left point of cell area. Then hit ESC ");
  cv::namedWindow(windowname, cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
  cv::Rect roi_ul = cv::Rect(0,0, simg_width, simg_height);
  cv::Mat ul_img = climg(roi_ul);
  cv::Mat orig_ul_img = ul_img.clone();
  cv::imshow(windowname, ul_img);

  //Prompt user for UL boundary point of image core box
  std::vector<cv::Point> points_ul;
  cv::setMouseCallback(windowname, onMouse2, (void*)&points_ul);
  while (true)
  {
    for (cv::Point& point : points_ul)
    {
      ul_img = orig_ul_img.clone();
      // Draw points on image
      cv::circle(ul_img, point, 5, cv::Scalar{0, 200, 0}, -1);
      cv::imshow(windowname, ul_img);
    }
    
    if (cv::waitKey(1) == 27)break;
    
  }
  
  cv::destroyAllWindows();
  cv::Point point1 = cv::Point(0,0);
  // Define point to find scale or pixel size
  if (!points_ul.empty())
    point1 = points_ul.back();
  
  // Display lower right corner and prompt user to select point
  windowname = "Select lower right point of cell area. Then hit ESC ";
  cv::namedWindow(windowname, cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
  height = climg.rows;
  width = climg.cols;
  int32_t start_x = width - simg_width - 1;
  int32_t start_y = height - simg_height - 1;
  cv::Rect roi_lr = cv::Rect(start_x, start_y, simg_width, simg_height);
  cv::Mat lr_img = climg(roi_lr);
  cv::Mat orig_lr_img = lr_img.clone();
  cv::imshow(windowname, lr_img);
  
  //Prompt user for LR boundary point of image core box
  std::vector<cv::Point> points_lr;
  cv::setMouseCallback(windowname, onMouse2, (void*)&points_lr);
  while (true)
  {
    for (cv::Point& point : points_lr)
    {
      lr_img = orig_lr_img.clone();
      // Draw points on image
      cv::circle(lr_img, point, 5, cv::Scalar{0, 200, 0}, -1);
      cv::imshow(windowname, lr_img);
    }
    
    if (cv::waitKey(1) == 27)break;
  }
  cv::destroyAllWindows();

  if (!points_lr.empty() && !points_ul.empty())
  {
    cv::Point point2 = points_lr.back();
    point2.x = start_x + point2.x;
    point2.y = start_y + point2.y;

    // Calculate scale
    // core_box_dx and core_box_dy are found from the core box data in the def file (in micrometers)
    double userdx = (double)point2.x - (double)point1.x;
    double userdy = (double)point2.y - (double)point1.y;
    pixel_size_x = (double)core_box_dx * 1000.0 / (double)userdx;
    pixel_size_y = (double)core_box_dy * 1000.0 / (double)userdy;
  }
  else
  {
    // Needs improvement - indicate to user that points were not selected!
    pixel_size_x = 1;
    pixel_size_y = 1;
  }
  
}

/////////////////////////////////
//Constructor for whole image
Tile_img::Tile_img(const cv::Mat &current_tile_image, std::string img_file, int32_t core_box_dx, int32_t core_box_dy) : image_file(img_file)
{
  //imgLocColX and imgLocRowY already equal (0,0)
  // Need pixel_size
  find_pixel_size(current_tile_image, img_file, core_box_dx, core_box_dy); // width, height, and pixel sizes are defined
}
//Constructor for individual images 
Tile_img::Tile_img(std::string img_file, std::string txt_file) : image_file(img_file), txt_file(txt_file)
{
  read_txtdata(txt_file);
}


/////////////////////////////////
// Divides the string at dlm and selects the section to return based on loc
std::string split_string(std::string line, const char* dlm, int32_t loc)
{
  std::stringstream ssline = std::stringstream(line);
  std::string substring;
  std::vector<std::string> substrlist;

  while (std::getline(ssline, substring, *dlm))
  {
    substrlist.push_back(substring);
  }
  return substrlist[loc];
}

/////////////////////////////////
// Parses the tile .txt files to get PixelSize (resolutions), StagePositionX, StagePositionY, DataSize(rowsxcols)
void Tile_img::read_txtdata(std::string txt_file)
{
  double stagePositionX = 0;
  double stagePositionY = 0;
  std::ifstream infile(txt_file);
  std::string line;

  // Search .txt file for "PixelSize", "StagePositionX", "StagePositionY", "DataSize"
  while (std::getline(infile, line))
  {
    if (line.find("PixelSize") != std::string::npos)
    {
      pixel_size_x = std::stod(split_string(line, "=", 1));
      pixel_size_y = pixel_size_x;
    }

    if (line.find("StagePositionX") != std::string::npos)
    {
      stagePositionX = abs(std::stod(split_string(line, "=", 1)));
    }

    if (line.find("StagePositionY") != std::string::npos)
    {
      stagePositionY = abs(std::stod(split_string(line, "=", 1)));
    }

    if (line.find("DataSize") != std::string::npos)
    {
      std::string datasz = split_string(line, "=", 1);
      height = std::stoi(split_string(datasz, "x", 1));
      width = std::stoi(split_string(datasz, "x", 0));
    }
  }
  initLocColX = (int32_t)floor(stagePositionX / pixel_size_x); // Drops all decimals to result in pixels
  initLocRowY = (int32_t)floor(stagePositionY / pixel_size_y); // Drops all decimals to result in pixels
}

/////////////////////////////////
// Returns the x and y location of this tile image in the entire stitched image (in pixels)
void Tile_img::find_ImgLoc(Tile_img& tile0)
{
  this->imgLocColX = abs(this->initLocColX - tile0.initLocColX);
  this->imgLocRowY = abs(this->initLocRowY - tile0.initLocRowY);
}


/////////////////////////////////
// Assigns colors to different standard cells for plotting
cv::Scalar Tile_img::get_color(std::string cellname)
{
  // Note that rand is seeded to 1 by default
  cv::Scalar color;
  int32_t a, b, c;
  int32_t m = 27;

  // Check if color was already used
  // Check color generated and make sure it's not close to red or yellow
  // Don't want dark colors
  // Don't want near white
  
  bool test = false;
  bool used_colors_test = false;
  int32_t count_loop = 0;

  // Some cells will be kept as specific colors
  if (cellname == "scc9gena_dfrtp_1")
  {
    color = cv::Scalar(0, 250, 180); //Green Yellow
    used_colors.push_back(color);
  }
  else if (cellname == "scc9gena_xorlp2_m")
  {
    color = cv::Scalar(255, 200, 150); //Blue
    used_colors.push_back(color);
  }
  else if (cellname == "scc9gena_xnorlp2_m")
  {
    color = cv::Scalar(0, 150, 255); //Orange
    used_colors.push_back(color);
  }
  else
  {
    test = true;
    count_loop = 0;
    while (test)
    {
      a = (rand() % 10) * m;
      b = (rand() % 10) * m;
      c = (rand() % 10) * m;
      if ((a < 70) && (b < 70) && (c > 150)) // close to red color
        continue;
      if ((a < 50) && (b > 150) && (c > 150)) // close to yellow color
        continue;
      if ((a < 75) && (b < 75) && (c < 75)) // color too dark
        continue;
      if ((a > 225) && (b > 225) && (c > 225)) // color too white
        continue;
      // Test if near a color in used_colors
      used_colors_test = false;
      for (int i = 0; i < used_colors.size(); i++)
      {
        if ((a < used_colors[i][0] + 20) && (a > used_colors[i][0] - 20) &&
            (b < used_colors[i][1] + 20) && (b > used_colors[i][1] - 20) &&
            (c < used_colors[i][2] + 20) && (c > used_colors[i][2] - 20))
        {
          used_colors_test = true;
          break;
        }
      }
      if (used_colors_test == true)
      {
        if (count_loop > 100) // If while loop has iterated 100 times, quit
          test = false; // push the current color onto used_colors
        else
        {
          count_loop++;
          continue;
        }
      }
      else
        test = false;
    }
    color = cv::Scalar(a, b, c);
    used_colors.push_back(color);
  }
     
  return color;
}

/////////////////////////////////
// Returns the cells associated with a specific ("this") tile image
void Tile_img::get_tile_cells(const std::vector<Cell_info> &def_cells)
{
  // Find the DEF cells associated with the given tile.
  // Clean out anything in tile.tilecells
  tile_cells.clear();

  for (int i = 0; i < def_cells.size(); i++ )
  {
    Cell_info tmp;
    tmp.name = def_cells[i].name;
    tmp.height = def_cells[i].height;
    tmp.width = def_cells[i].width;
    tmp.orientations = def_cells[i].orientations;
    tmp.color = get_color(tmp.name);
    for (int j = 0; j < def_cells[i].pt_locations.size(); j++)
    {
      if ((def_cells[i].pt_locations[j].pt.x <= imgLocColX + this->width) && (def_cells[i].pt_locations[j].pt.x >= imgLocColX))
      {
        if ((def_cells[i].pt_locations[j].pt.y <= imgLocRowY + this->height) && (def_cells[i].pt_locations[j].pt.y >= imgLocRowY))
        {

          Pt_info pttmp;
          pttmp.pt = def_cells[i].pt_locations[j].pt;
          pttmp.orientation = def_cells[i].pt_locations[j].orientation;
          pttmp.matched = def_cells[i].pt_locations[j].matched;
          pttmp.glob_indx = def_cells[i].pt_locations[j].glob_indx;
          tmp.pt_locations.push_back(pttmp);
        }
      }
    }
    tile_cells.push_back(tmp);  
  }
}



/////////////////////////////////
// Plot the standard cells by color onto the tile image
void Tile_img::plot_DEFpts(const std::vector<Cell_info> &def_cells)
{
  get_tile_cells(def_cells);

  cv::Mat mat;
  mat = cv::imread(image_file);
  cv::namedWindow("Display Points", cv::WINDOW_GUI_EXPANDED | cv::WINDOW_KEEPRATIO);
  cv::resizeWindow("Display Points", 1000, 700);

  for (Cell_info c : tile_cells)
  {
    for (Pt_info pl : c.pt_locations)
    {
      cv::Point tmpt;
      tmpt.x = pl.pt.x - imgLocColX;
      tmpt.y = pl.pt.y - imgLocRowY;
      cv::circle(mat, tmpt, 5, c.color, -1);
    }
  }

  cv::imshow("Display Points", mat);
  cv::waitKey(1);
}


////////////////////////////////////
// Destructor
Tile_img::~Tile_img()
{
}


