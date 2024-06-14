#pragma once

#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>


//////////////////////////////////////////////////////////////////////////////////////////
// Class Pt_info for pt_locations in Cell_information
class Pt_info
{
public:
	uint16_t orientation = 99;
	cv::Point pt;
	cv::Point defpt;
	std::vector<std::int32_t> glob_indx; //0: cell index, 1: pt_location index
	int8_t matched = 0;		// -1: errant cell with no match, 0: not searched yet, 1: found, 2: not found
	double distance;
  std::int32_t on_tile = -1;
	bool power;
	//double pwr_distance;

	// Constructors
	Pt_info() = default;
	Pt_info(uint16_t o_, cv::Point pt_, int8_t m_) : orientation(o_), pt(pt_), matched(m_) 
	{
		distance = -1.0;		// std::numeric_limits<double>::max();
	}
	Pt_info(const Pt_info& p, double d_) : orientation(p.orientation), pt(p.pt), matched(p.matched), distance(d_) { }

private:

};


//////////////////////////////////////////////////////////////////////////////////////////
// Class Cell_information to store all standard cells listed in DEF file
class Cell_info
{
public:
	// Attributes
	std::string name;
	uint32_t height = 0;
	uint32_t width = 0;
	uint32_t height_nm = 0;
	uint32_t width_nm = 0;
	std::vector<uint16_t> orientations; // Possible orientations for this cell name
	std::vector<Pt_info> pt_locations;  // Cell locations listed in DEF file for this cell name
	cv::Scalar color;
	bool selected = 0;

	// Constructors
	Cell_info() = default;

	Cell_info(std::string n_, uint32_t h_, uint32_t w_) : name(n_), height(h_), width(w_), orientations({0,1,2,3})
	{		
		pt_locations.clear();
		color = cv::Scalar(255, 0, 0);
	}

	Cell_info(std::string n_, uint32_t h_, uint32_t w_, std::vector<uint16_t> o_) : name(n_), height(h_), width(w_), orientations(o_)
	{
		pt_locations.clear();
		color = cv::Scalar(255, 0, 0);
	}

	Cell_info(std::string n_, uint32_t h_, uint32_t w_, std::vector<uint16_t> o_, std::vector<Pt_info> pt_) : name(n_), height(h_), width(w_), orientations(o_), pt_locations(pt_)
	{
		color = cv::Scalar(255, 0, 0);
	}

	Cell_info(const Cell_info& ci) : name(ci.name), height(ci.height), width(ci.width), orientations(ci.orientations), pt_locations(ci.pt_locations), color(ci.color), height_nm(ci.height_nm), width_nm(ci.width_nm)
	{
		//pt_locations.clear();
	}

	Cell_info(std::string n_) : name(n_)
	{
		pt_locations.clear();
		color = cv::Scalar(255, 0, 0);
	}

private:

};
