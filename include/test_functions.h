#ifndef _TEST_FUNCTIONS_H_
#define	_TEST_FUNCTIONS_H_

#define _CRT_SECURE_NO_WARNINGS 1


//-----------------------------------------------------------------------------
void report_grid_generation(int64_t img_w, int64_t img_h, int32_t tile_w, int32_t tile_h, int32_t cell_w, int32_t cell_h, std::vector<cv::Rect>& grid)
{
    int32_t idx, jdx;

    int32_t xl, yl, xr, yr;

    cv::Rect r;

    // clear out the vector that 
    grid.clear();

    // run through the image height and width in tile_h and tile_w increments 
    for (idx = -cell_h; idx < img_h; idx += (tile_h))
    {
        for (jdx = -cell_w; jdx < img_w; jdx += (tile_w))
        {
            // create a +/- cell_w and cell_h buffer for the image
            xl = std::max<int64_t>(0, jdx);
            xr = std::min<int64_t>(img_w, xl+tile_w);

            yl = std::max<int64_t>(0, idx);
            yr = std::min<int64_t>(img_h, yl+tile_h);

            // create the rect
            r = cv::Rect(cv::Point(xl, yl), cv::Point(xr, yr));

            // push the rect back into the vector
            grid.push_back(r);
        }
    }

}   // end of report_grid_generation


#endif // _TEST_FUNCTIONS_H_
