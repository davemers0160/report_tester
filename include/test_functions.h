#ifndef _TEST_FUNCTIONS_H_
#define	_TEST_FUNCTIONS_H_

#define _CRT_SECURE_NO_WARNINGS 1

#include <limits>


//-----------------------------------------------------------------------------
void report_grid_generation(int64_t img_w, int64_t img_h, int32_t tile_w, int32_t tile_h, int32_t cell_w, int32_t cell_h, std::vector<Tile_img>& grid)
{
    int32_t idx, jdx;

    int32_t t_w, t_h;
    int32_t row = 0, col = 0;

    // clear out the vector
    grid.clear();

    // run through the image height and width in tile_h and tile_w increments 
    row = 0;
    while (row < img_h)
    {
        col = 0;
        while (col < img_w)
        {

            t_w = (col + tile_w >= img_w) ? img_w - col : tile_w;
            t_h = (row + tile_h >= img_h) ? img_h - row : tile_h;

            // push the rect back into the vector
            grid.push_back(Tile_img(col, row, t_h, t_w, true, true));

            col += (tile_w - cell_w);
        }

        row += (tile_h - cell_h);
    }

}   // end of report_grid_generation


//-----------------------------------------------------------------------------
void calculate_cell_maximums(std::vector<Cell_info>& cell_list, int32_t& max_cell_width, int32_t& max_cell_height)
{
    uint32_t idx;

    max_cell_width = 0;
    max_cell_height = 0;
    
    // go through each cell and find the max size
    for (Cell_info& c : cell_list)
    {
        max_cell_width = std::max(max_cell_width, (int32_t)c.width);
        max_cell_height = std::max(max_cell_height, (int32_t)c.height);
    }

}   // end of calculate_cell_maximums



#endif // _TEST_FUNCTIONS_H_
