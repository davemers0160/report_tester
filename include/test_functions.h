#ifndef _TEST_FUNCTIONS_H_
#define	_TEST_FUNCTIONS_H_

#define _CRT_SECURE_NO_WARNINGS 1


//-----------------------------------------------------------------------------
void report_grid_generation(int64_t img_w, int64_t img_h, int32_t tile_w, int32_t tile_h, int32_t cell_w, int32_t cell_h, std::vector<cv::Rect>& grid)
{
    int32_t idx, jdx;

    int32_t t_w, t_h;
    int32_t row = 0, col = 0;

    //int32_t xl, yl, xr, yr;

    cv::Rect r;

    // clear out the vector that 
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

            r = cv::Rect(cv::Point(col, row), cv::Size(t_w, t_h));

            // push the rect back into the vector
            grid.push_back(r);

            col += (tile_w - cell_w);
        }

        row += (tile_h - cell_h);
    }

}   // end of report_grid_generation


#endif // _TEST_FUNCTIONS_H_
