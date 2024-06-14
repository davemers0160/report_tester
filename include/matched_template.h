#ifndef _MATCHTED_TEMPLATE_CLASS_H_
#define _MATCHTED_TEMPLATE_CLASS_H_

#include <cstdint>
#include <cstdlib>
#include <string>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>

//-----------------------------------------------------------------------------
class matched_template
{
public:
	std::string name;
    std::vector<cv::Mat> imgs;
    uint16_t orientation;
    
    //-----------------------------------------------------------------------------
	matched_template() = default;

    matched_template(std::string n_) : name(n_), orientation(0)
    {
        imgs.clear();
    }

    matched_template(uint16_t o_) : name(""), orientation(o_) 
    {
        imgs.clear();
    }

    matched_template(std::string n_, uint16_t o_) : name(n_), orientation(o_) 
    {
        imgs.clear();
    }

    matched_template(std::string n_, cv::Mat img_, uint16_t o_) : name(n_), orientation(o_)
    {
        imgs.clear();
        imgs.push_back(img_);
    }

    //-----------------------------------------------------------------------------
    cv::Mat get_mosaic(void)
    {
        uint32_t idx;
        cv::Mat mosaic = cv::Mat::ones(2*imgs[0].rows + 4, 2*imgs[0].cols + 4, imgs[0].type());
        mosaic = 255 * mosaic;

        for (idx = 0; idx < imgs.size(); ++idx)
        {
            switch (idx)
            {
            case 0:
                imgs[idx].copyTo(mosaic(cv::Rect(0, 0, imgs[idx].cols, imgs[idx].rows)));
                //mosiac(cv::Rect(0, 0, imgs[idx].cols, imgs[idx].rows)) = imgs[idx].clone();
                break;

            case 1:
                imgs[idx].copyTo(mosaic(cv::Rect(imgs[idx].cols + 3, 0, imgs[idx].cols, imgs[idx].rows)));
                break;

            case 2:
                imgs[idx].copyTo(mosaic(cv::Rect(0, imgs[idx].rows + 3, imgs[idx].cols, imgs[idx].rows)));
                break;

            case 3:
                imgs[idx].copyTo(mosaic(cv::Rect(imgs[idx].cols + 3, imgs[idx].rows + 3, imgs[idx].cols, imgs[idx].rows)));
                break;
            }
        }

        return mosaic;

    }   // end of get_mosaic


private:

};	// end of matched_template

#endif	// _MATCHTED_TEMPLATE_CLASS_H_
