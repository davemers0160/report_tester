// This file was created by Victoria Lockridge, NSWC Crane
// October 3, 2023
#pragma once
#define _CRT_SECURE_NO_WARNINGS 1
#define IMAGE_LIMIT 25

#include <vector>
#include <string>
#include <time.h>
#include "cell_information.h"
#include "tile_img.h"
#include "utils.h"
#include "matched_template.h"


////////////////////////////////////////////////
////////////////////////////////////////////////
/// Unordered Lists (Images Processed, Cells Processed)
class Html_ul
{
private:
	std::vector<std::string> data_list;
	std::uint32_t is_cells; 

public:
	Html_ul(const std::vector<std::string> &data_list, bool is_cells);
	~Html_ul();
	void write_heading(std::ostream& out, bool is_cells);
	void write_rows(std::ostream& out);
	friend std::ostream& operator << (std::ostream& out, Html_ul ul_tmp);
};

//******************************************************
Html_ul::Html_ul(const std::vector<std::string> &data_list, bool is_cells)
{
	this->data_list = data_list;
	this->is_cells = is_cells;
}
Html_ul::~Html_ul()
{
};

/////////////////////////////////////////////////////////////////////////////////////////
//
void Html_ul::write_heading(std::ostream& out, bool is_cells)
{
	std::string ul_heading;
	if (is_cells)
		ul_heading = "<h2>Cells Processed </h2>";
	else
		ul_heading ="<h2>Images Processed </h2>";
	out << ul_heading << std::endl;
};

/////////////////////////////////////////////////////////////////////////////////////////
//
void Html_ul::write_rows(std::ostream& out)
{
	for (int i = 0; i < this->data_list.size(); i++)
	{
		out << "<li>" << data_list[i] << "</li>" << std::endl;
	}

}






//******************************************************
// Overloaded stream output operator << for ul         *
//******************************************************
	std::ostream& operator << (std::ostream& out, Html_ul ul_tmp)
	{
		out << "<div class=\"section_list\">" << std::endl;
		ul_tmp.write_heading(out, ul_tmp.is_cells);
		// Write beginning tag for ul
		out << "<ul>" << std::endl;
		ul_tmp.write_rows(out);
		// Write end tag for ul
		out << "</ul>" << std::endl;
		// Write end tag for div
		out << "</div>" << std::endl;
		return out;
	}
	

////////////////////////////////////////////////
////////////////////////////////////////////////
// Write heading and date for heading to .html
class Html_top_heading
{
private:
	std::string top_heading = "Analysis Report";
	std::string date_heading;
	std::string get_date_string();
public:
	Html_top_heading();
	void write_heading(std::ostream& out);
	friend std::ostream& operator << (std::ostream& out, Html_top_heading heading_tmp);

};

//******************************************************
Html_top_heading::Html_top_heading()
{
	this->date_heading = get_date_string();
}

std::string Html_top_heading::get_date_string()
{
	std::array<char, 16> buffer;
	buffer.fill(0);
	time_t rawtime;
	time(&rawtime);
	const auto timeinfo = localtime(&rawtime);
	strftime(buffer.data(), sizeof(buffer), "%d-%m-%Y", timeinfo);
	std::string month[] = {"", "January", "February", "March", "April", "May", "June",
													"July", "August", "September", "October", "November", "December"};
	std::string str(buffer.data());
	std::string mo = str.substr(3,2);
	std::string day = str.substr(0, 2);
	std::string year = str.substr(6,4);
	std::string mword = month[stoi(mo)];
	std::string date_words = day + " " + mword + " " + year;

	return date_words;
}

//******************************************************
// Overloaded stream output operator << for heading          *
//******************************************************
std::ostream& operator << (std::ostream& out, Html_top_heading heading_tmp)
{
	out << "<div class=\"top_heading\">" << std::endl;
	// Write top_heading
	out << "<h1>" << heading_tmp.top_heading << "</h1>" << std::endl;;
	// Write date_heading
	out << "<h1>" << heading_tmp.date_heading << "</h1>" << std::endl;;
	// Write end tag for div
	out << "</div>" << std::endl;
	return out;
}



////////////////////////////////////////////////
////////////////////////////////////////////////
// Create section for 2 tables of parameters
class Params_table
{
public:
	std::string section_heading;
	std::vector<std::string> table_headings;
	std::vector<std::string> labels;
	std::vector<std::string> values;
	std::vector<std::string> values2;
	std::vector<std::string> values3;

	Params_table();
};

Params_table::Params_table() {};

class Html_params
{
public:
	std::string heading = "Parameter Settings";
	Params_table search_params; 
	Params_table tile_binary_params;

	Html_params();
	void get_params(std::vector<Tile_img> tiles, double match_thresh, double overlap_threshold, double distance_error);
	void setup_search_params(double overlap_threshold, double distance_error);
	void setup_tile_params(std::vector<Tile_img> tiles, double match_thresh);
	friend std::ostream& operator << (std::ostream& out, Html_params params_section);
};

Html_params::Html_params() {};

void Html_params::get_params(std::vector<Tile_img> tiles, double match_thresh, double overlap_threshold, double distance_error)
{
	setup_search_params(overlap_threshold, distance_error); 
	setup_tile_params( tiles, match_thresh);
}

void Html_params::setup_search_params(double overlap_threshold, double distance_error)
{
	search_params.section_heading = "Search Parameter Settings";
	search_params.table_headings = {"Parameter", "Setting"};
	search_params.labels = {"Overlap Threshold", "Distance Error"};
	search_params.values = {std::to_string(overlap_threshold) , std::to_string(distance_error)};

}

void Html_params::setup_tile_params(std::vector<Tile_img> tiles, double match_thresh)
{
	std::string tmp_val, tmp_val2, tmp_val3, tmp_label;
	
	tile_binary_params.section_heading = "Tile Binary Parameter Settings";
	tile_binary_params.table_headings = {"Tile Number", "Structural Element Size","Binary Threshold", "Match Threshold"};

	for (int i = 0; i < tiles.size(); i++)
	{
		tmp_val = std::to_string(tiles[i].struct_elem_size);
		tile_binary_params.values.push_back(tmp_val);
		tmp_val2 = std::to_string(tiles[i].binary_threshold);
		tile_binary_params.values2.push_back(tmp_val2);
		tmp_val3 = std::to_string(match_thresh);
		tile_binary_params.values3.push_back(tmp_val3);
		tmp_label = std::to_string(i);
		tile_binary_params.labels.push_back(tmp_label);

	}
}


//******************************************************
// Overloaded stream output operator << for Html_params          *
//******************************************************
std::ostream& operator << (std::ostream& out, Html_params params)
{
	// Place both global and tile parameter settings in one section
	// 2 tables required
	// Pass in class with global params and tile params
	out << "<div class=\"section_params_table\">" << std::endl;
	// Write  heading
	out << "<h2>" << params.heading << "</h2>" << std::endl;
	// Write search table heading
	out << "<h3>" << params.search_params.section_heading << "</h3>" << std::endl;
	// Write search table
	out << "<table class=\"table2\">" << std::endl;
	out << "<thead>" << std::endl;
	out << "<tr>" << std::endl;
	out << "<th>" << params.search_params.table_headings[0] << "</th>" << std::endl;
	for (int i = 1; i < params.search_params.table_headings.size(); i++)
	{
		out << "<th scope=\"col\">" << params.search_params.table_headings[i] << "</th>" << std::endl;
	}
	out << "</tr>" << std::endl;
	out << "</thead>" << std::endl;

	out << "<tbody>" << std::endl;
	for (int i = 0; i < params.search_params.labels.size(); i++)
	{
		out << "<tr>" << std::endl;
		out << "<th scope = \"row\">" << params.search_params.labels[i] << "</th>" << std::endl;
		out << "<td>" << params.search_params.values[i] << "</td>" << std::endl;
		out << "</tr>" << std::endl;
	}

	out << "</tbody>" << std::endl;
	out << "</table>" << std::endl;

	// Write  binary parameters table
	// Write binary parameters table heading
	out << "<h3>" << params.tile_binary_params.section_heading << "</h3>" << std::endl;
	// Write search table
	out << "<table class=\"table3\">" << std::endl;
	out << "<thead>" << std::endl;
	out << "<tr>" << std::endl;
	out << "<th>" << params.tile_binary_params.table_headings[0] << "</th>" << std::endl;
	for (int i = 1; i < params.tile_binary_params.table_headings.size(); i++)
	{
		out << "<th scope=\"col\">" << params.tile_binary_params.table_headings[i] << "</th>" << std::endl;
	}
	out << "</tr>" << std::endl;
	out << "</thead>" << std::endl;

	out << "<tbody>" << std::endl;
	for (int i = 0; i < params.tile_binary_params.labels.size(); i++)
	{
		out << "<tr>" << std::endl;
		out << "<th scope = \"row\">" << params.tile_binary_params.labels[i] << "</th>" << std::endl;
		out << "<td>" << params.tile_binary_params.values[i] << "</td>" << std::endl;
		out << "<td>" << params.tile_binary_params.values2[i] << "</td>" << std::endl;
		out << "<td>" << params.tile_binary_params.values3[i] << "</td>" << std::endl;
		out << "</tr>" << std::endl;
	}

	out << "</tbody>" << std::endl;
	out << "</table>" << std::endl;

	out << "</div>" << std::endl;
	return out;
}



////////////////////////////////////////////////
////////////////////////////////////////////////
// Create a vector of table rows, each row tallies one cell type 
class Table_row
{
public:
	std::string cell;
	std::int32_t design_num;
	std::int32_t found_num;
	std::int32_t missed_num;
	std::int32_t errant_num;
	std::int32_t on_edge_num;
};

////////////////////////////////////////////////
////////////////////////////////////////////////
/// Table for results
class Html_table
{
private:
	std::vector<Table_row> table_rows;
	std::string heading = "High Level Analysis Results";
	void get_table_rows(const std::vector<Cell_info> &cells, const std::vector<std::vector<Pt_info>>& detects, const std::vector<Tile_img>& tiles);
	std::vector<std::string> table_headings = {"Cell", "Design", "Found", "Missed", "Errant"};
	void calc_percent_area_found(const std::vector<Cell_info>& cells, const std::vector<std::vector<Pt_info>>& detects, const std::vector<Tile_img>& tiles, double corebox_area);
	std::int32_t get_total_cells(const std::vector<Tile_img>& tiles, int32_t num_cells);
	double proc_time;

public:
	Html_table(const std::vector<Cell_info> &cells, const std::vector<std::vector<Pt_info>> &detects, const std::vector<Tile_img>& tiles, double corebox_area, double corr_time);
	double percent_area_found = 0;  // percentage of area of cells found to area of cells selected to find
	std::int32_t total_design = 0, total_matched = 0, total_missed = 0, total_errant = 0;
	friend std::ostream& operator << (std::ostream& out, Html_table table_tmp);
};

//******************************************************
Html_table::Html_table(const std::vector<Cell_info> &cells, const std::vector<std::vector<Pt_info>> &detects, const std::vector<Tile_img>& tiles, double corebox_area, double corr_time)
{
	proc_time = corr_time;
	get_table_rows(cells, detects, tiles); // , missed, errant);
	// Report calculated area used
	// Calculate used cell area now - before converting nanometers to pixels (mirror_cell_locations)
	// Only report tiles looked at
	// Pass in tile_cells and matched info
	
	// Pass in matched and errant cells
	calc_percent_area_found(cells, detects, tiles, corebox_area);
	std::int32_t num_cells = (int32_t)cells.size();
	total_design = get_total_cells(tiles, num_cells);

};

int32_t Html_table::get_total_cells(const std::vector<Tile_img>& tiles, std::int32_t num_cells)
{
	std::int32_t total_cells = 0;
	for (int i = 0; i < tiles.size(); i++)
		for (int j = 0; j < num_cells; j++)
			if (tiles[i].tile_cells[j].selected)
				total_cells = total_cells + (int32_t)tiles[i].tile_cells[j].pt_locations.size();
	return total_cells;
}


void Html_table::get_table_rows(const std::vector<Cell_info> &cells, const std::vector<std::vector<Pt_info>> &detects, const std::vector<Tile_img>& tiles)
{
	for (int i = 0; i < cells.size(); i++)
	{
		if (cells[i].selected == true)
		{
			Table_row tmp_row;
			tmp_row.cell = cells[i].name;
			// Find number of found, missed, errant cells based on cell type
			int32_t count_matched = 0;
			int32_t count_missed = 0;
			int32_t count_errant = 0;
			int32_t count_design = 0;
			int32_t count_on_edge = 0;

			for (int jdx = 0; jdx < tiles.size(); ++jdx)
			{
				for (int kdx = 0; kdx < cells[i].pt_locations.size(); ++kdx)
				{
					if (cells[i].selected)
					{
						if ((cells[i].pt_locations[kdx].on_tile == jdx) & (cells[i].pt_locations[kdx].matched == 1))
						{
							count_matched++;
							total_matched++;
						}
						else if ((cells[i].pt_locations[kdx].on_tile == jdx) & (cells[i].pt_locations[kdx].matched == 2))
						{
							count_missed++;
							total_missed++;
							if ((cells[i].pt_locations[kdx].pt.x > (tiles[jdx].width - cells[i].width)) | (cells[i].pt_locations[kdx].pt.y > (tiles[jdx].height - cells[i].height)))
								count_on_edge++;
						}
					}
				}
			}

			// Find errants for this cell in detects
			// Assumption:  detects vector includes all detections from all tiles that have been searched during this one instance of processing
			for (int k = 0; k < detects[i].size(); k++)
			{
				if (detects[i][k].matched == -1)
				{
					count_errant++;
					total_errant++;
				}
			}

			tmp_row.missed_num = count_missed;
			tmp_row.found_num = count_matched;
			tmp_row.errant_num = count_errant;
			tmp_row.on_edge_num = count_on_edge;

			// Get design cell count from tiles
			for (int idx = 0; idx < tiles.size(); idx++)
			{
				if (tiles[idx].aligned)
					count_design = count_design + tiles[idx].tile_cells[i].pt_locations.size();
			}
			tmp_row.design_num = count_design;
			table_rows.push_back(tmp_row);
		}
	}
}

///////////////////////////////////////////////////////
// Calculate the percentage of (area of cells found) to (area of cells searched for)
// This calculation used pixels - could use nanometers
void Html_table::calc_percent_area_found(const std::vector<Cell_info>& cells, const std::vector<std::vector<Pt_info>>& detects, const std::vector<Tile_img>& tiles, double corebox_area)
{
	double cell_area = 0;
	double design_area = 0;
	double found_area = 0;

	for (int i = 0; i < cells.size(); i++)
	{
		if (cells[i].selected == true)
		{
			cell_area = (double)cells[i].width * (double)cells[i].height;
			for (int jdx = 0; jdx < tiles.size(); ++jdx)
			{
				for (int kdx = 0; kdx < cells[i].pt_locations.size(); ++kdx)
				{
					if ((cells[i].pt_locations[kdx].on_tile == jdx) & (cells[i].pt_locations[kdx].matched == 1))
					{
						design_area += cell_area;
						found_area += cell_area;
					}
					else if ((cells[i].pt_locations[kdx].on_tile == jdx) & (cells[i].pt_locations[kdx].matched == 2))
						design_area += cell_area;
				}
			}

			// Negelected errant cells for this calculation
		}
	}

	percent_area_found = (found_area / design_area) * 100;

}


//******************************************************
// Overloaded stream output operator << for table          *
//******************************************************
std::ostream& operator << (std::ostream& out, Html_table table_tmp)
{
	out << "<div class=\"section_table\">" << std::endl;
	// Write table heading
	out << "<h2 id=\"Table1\">" << table_tmp.heading << "</h2>" << std::endl;
	// Write table
	out << "<table class=\"table1\">" << std::endl;
	out << "<thead>" << std::endl;
	out << "<tr>" << std::endl;
	out << "<th>"<< table_tmp.table_headings[0] << "</th>" << std::endl;
	for (int i = 1; i < table_tmp.table_headings.size(); i++)
	{
		out << "<th scope=\"col\">" << table_tmp.table_headings[i] << "</th>" << std::endl;
	}
	out << "</tr>" << std::endl;
	out << "</thead>" << std::endl;

	out << "<tbody>" << std::endl;
	for (int i = 0; i < table_tmp.table_rows.size(); i++)
	{
		out << "<tr>" << std::endl;
		out << "<th scope = \"row\"><a href=\"#" << table_tmp.table_rows[i].cell << "\">" << table_tmp.table_rows[i].cell << "</a>" << "</th>" << std::endl;
		out << "<td>" << table_tmp.table_rows[i].design_num << "</td>" << std::endl;
		out << "<td scope = \"col\">" << table_tmp.table_rows[i].found_num << "</td>" << std::endl;
		out << "<td scope = \"col\">" << table_tmp.table_rows[i].missed_num << "</td>" << std::endl;
		out << "<td scope = \"col\">" << table_tmp.table_rows[i].errant_num << "</td>" << std::endl;
		out << "</tr>" << std::endl;
	}
	out << "</tbody>" << std::endl;

	
	// Add percentage of design cells
	out << "<tfoot>" << std::endl;
	out << "<tr>" << std::endl;
	out << "<th>Percentage of Design</th>" << std::endl;
	out << "<td scope=\"col\">" << std::to_string((double)table_tmp.total_design * 100.0/ (double)table_tmp.total_design).substr(0, 5) << "</td>" << std::endl;
	out << "<td scope=\"col\">" << std::to_string((double)table_tmp.total_matched * 100.0 / (double)table_tmp.total_design).substr(0, 5) << "</td>" << std::endl;
	out << "<td scope=\"col\">" << std::to_string((double)table_tmp.total_missed * 100.0 / (double)table_tmp.total_design).substr(0, 4) << "</td>" << std::endl;
	out << "<td scope=\"col\">" << std::to_string((double)table_tmp.total_errant * 100.0 / (double)table_tmp.total_design).substr(0, 4) << "</td>" << std::endl;
	out << "</tr>" << std::endl;
	out << "</tfoot>" << std::endl;

	// Add totals in table tfoot
	out << "<tfoot>" << std::endl;
	out << "<tr>" << std::endl;
	out << "<th>Totals</th>" << std::endl;
	out << "<td scope=\"col\">" << table_tmp.total_design << "</td>" << std::endl;
	out << "<td scope=\"col\">" << table_tmp.total_matched << "</td>" << std::endl;
	out << "<td scope=\"col\">" << table_tmp.total_missed << "</td>" << std::endl;
	out << "<td scope=\"col\">" << table_tmp.total_errant << "</td>" << std::endl;
	out << "</tr>" << std::endl;
	out << "</tfoot>" << std::endl;

	out << "</table>" << std::endl;

	// Write percentage of found area
	std::string s = std::to_string(table_tmp.percent_area_found).substr(0,4);
	std::string str_darea = "Percent of area of found cells to area of design cells : " + s + "%";
	out << "<h4>" << str_darea << "</h4>" << std::endl;

	std::string proc_string1 = "Time to process finding cells: " + std::to_string(table_tmp.proc_time).substr(0, 6) + " seconds.";
	std::string proc_string2 = "Time to process finding cells: " + std::to_string(table_tmp.proc_time/3600.0).substr(0, 5) + " hours.";
	out << "<h4>" << proc_string1 << "</h4>" << std::endl;
	out << "<h4>" << proc_string2 << "</h4>" << std::endl;

	out << "</div>" << std::endl;
	return out;
}


////////////////////////////////////////////////
////////////////////////////////////////////////
// Images
// Html_img is made up of Img_sections.  Each Img_section (based on cell) is made up of 2 vectors (missed, errant) of Img_subs.

// Contains all information needed to printout a missed or errant cell (text and image)
class Img_sub 
{
public:
	cv::Point coords; //x,y
	cv::Point def_coords; //x,y
	std::string img_file_name;
	std::string himg_filename;
	std::int32_t sub_img_width;
	std::int32_t sub_img_height;
	std::int32_t html_img_width;
	std::int32_t html_img_height = 150;

	Img_sub(Pt_info ptloc, std::string cell_name, std::string type, std::string html_img_dir, std::string h_img_dir); // type is "m" or "e"
	void write_image(int ptloc_x, int ptloc_y, int cell_width, int cell_height, cv::Mat img, std::string sub_img_filename, std::int32_t tilex, std::int32_t tiley, std::string type, std::int32_t large_cell_width);
};

// Constructor:  added image zoom information to create image
// Used missed and errant to setup separate sections
// Note: Type:  m for missed, e for errant
Img_sub::Img_sub(Pt_info ptloc, std::string cell_name, std::string type, std::string html_img_dir, std::string h_img_dir)
{
	coords = ptloc.pt;
	def_coords = ptloc.defpt;
	img_file_name = html_img_dir + "/" + cell_name + "_" + std::to_string(coords.x) + "_" + std::to_string(coords.y) + "_" + type + ".png";
	himg_filename = h_img_dir + "/" + cell_name + "_" + std::to_string(coords.x) + "_" + std::to_string(coords.y) + "_" + type + ".png";
}


void Img_sub::write_image(int ptloc_x, int ptloc_y, int cell_width, int cell_height, cv::Mat img, std::string sub_img_filename, std::int32_t tilex, std::int32_t tiley, std::string type, std::int32_t large_cell_width)
{
	cv::Mat imgW = img.clone();
	// Draw rectangle on image (blue missed, red errant)
	cv::Point display_pt = cv::Point(ptloc_x - tilex, ptloc_y - tiley);
	cv::Rect roi(display_pt.x, display_pt.y, (int)cell_width, (int)cell_height);
	if (type == "m")
		cv::rectangle(imgW, roi, cv::Scalar{200, 0, 0}, 4);
	else if (type == "e")
		cv::rectangle(imgW, roi, cv::Scalar{0, 0, 200}, 4);

	// Define rectangle for subsection image
	sub_img_width = 4*large_cell_width; // 4 * cell_width;
	sub_img_height = 4*cell_height;
	cv::Size sz = imgW.size();
	int img_width = sz.width;
	int img_height = sz.height;
	// Find center of cell and subtract to get start of image
	cv::Point cell_center; 
	cell_center.x = ptloc_x - tilex + (int)floor(cell_width / 2);
	cell_center.y = ptloc_y - tiley + (int)floor(cell_height / 2);
	int sx = cell_center.x - (int)floor(sub_img_width / 2);
	int sy = cell_center.y - (int)floor(sub_img_height / 2);
	// Test for left edge/top of img
	sx = std::max(1, sx);
	sy = std::max(1, sy);
	
	// Test for right edge/bottom of img
	if (img_width - sx < sub_img_width)
		sx = img_width - sub_img_width;
	if (img_height - sy < sub_img_height)
		sy = img_height - sub_img_height;

	cv::Rect rect = cv::Rect(sx, sy, sub_img_width, sub_img_height);  
	cv::Mat sub_image = imgW(rect);

	// Write/Save images of missed and errant cells
	std::vector<int> cmp(2);
	cmp[0] = cv::IMWRITE_PNG_COMPRESSION;
	cmp[1] = 4; // compression factor - use values of 0 to 9 - default is 1

	// Reduce the number of rows and cols
	cv::Mat resized_img;
	double reduction = 0.50;
	cv::resize(sub_image, resized_img, cv::Size(sub_image.cols * reduction, sub_image.rows * reduction), 0, 0, cv::INTER_LINEAR);

	cv::imwrite(sub_img_filename, resized_img, cmp);

	// Find width for html
	html_img_width = (int32_t)floor((double)sub_img_width * ((double)html_img_height/ (double)sub_img_height));
}


/////////////////////////////////////////////////////////////////////////////////////////
// Write all matched_templates template images for each orientation to dir template_images
void write_template_imgs(std::vector<matched_template> matched_templates, std::vector<Cell_info> cell_list, std::string templ_dir)
{
	for (int i = 0; i < matched_templates.size(); i++)
	{
		for (int j = 0; j < matched_templates[i].imgs.size(); j++)
		{
			std::string cname = cell_list[i].name;
			std::string orient = symname(j);
			std::string tfilename = templ_dir + "/" + "template_" + cname + "_" + orient + ".png";
			cv::Mat timg;
			matched_templates[i].imgs[j].convertTo(timg, CV_8UC1, 255.0);
			std::vector<int> cmp(2);
			cmp[0] = cv::IMWRITE_PNG_COMPRESSION;
			cmp[1] = 4; // compression factor - use values of 0 to 9 - default is 1
			cv::imwrite(tfilename, timg), cmp;
		}
	}
}

// Contains all information to write out missed and errant cells attributed to one cell type
class Img_section
{
public:
	std::string cell_name;
	cv::Mat cell_templ_img;
	std::string templ_filename;
	std::string cell_orientation;
	std::string section_heading;
	std::string section_subheading_width;
	std::string section_subheading_height;
	std::vector<Img_sub> missed;
	std::vector<Img_sub> errant;

	Img_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, std::int32_t large_cell_width, bool whole_img, const cv::Mat& current_tile_image, std::string html_id, std::string h_id, std::string tmp_ld, Pixel_to_DEF_conversion convert_info, cv::Mat cell_template_img, int orientation)
		: html_img_dir(html_id), h_img_dir(h_id), templ_dir(tmp_ld)
	{
		setup_section(cell, cell_detect, tiles, large_cell_width, whole_img, current_tile_image, convert_info, cell_template_img, orientation);
	}
	
	void setup_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, std::int32_t large_cell_width, bool whole_img, const cv::Mat& current_tile_image, Pixel_to_DEF_conversion convert_info, cv::Mat cell_template_img, int orientation);

	void setup_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, bool whole_img, const cv::Mat& current_tile_image, cv::Mat cell_template_img, int orientation);

private:
	std::string html_img_dir;
	std::string h_img_dir;
	std::string templ_dir;

	int64_t max_tile_height = 1000;
	int64_t max_tile_width = 1000;

	double tile_height_percent = 0.1;
	double tile_width_percent = 0.1;

	//-----------------------------------------------------------------------------
	void create_template_header(Cell_info cell, cv::Mat cell_template_img, int orientation)
	{
		cell_template_img.convertTo(cell_templ_img, CV_8UC1, 255.0);
		cell_name = cell.name;
		cell_orientation = symname((int)orientation);
		templ_filename = templ_dir + "/" + "template_" + cell.name;

		// One cell type gets passed from cell list and from detects list
		section_heading = " id=\"" + cell.name + "\">Cell: " + cell.name;
		section_subheading_width = std::to_string((int)cell.width_nm);
		section_subheading_height = std::to_string((int)cell.height_nm);

	}	// end of create_template_header

	//-----------------------------------------------------------------------------
	void create_report_tiles(int64_t img_w, int64_t img_h, int64_t tile_w, int64_t tile_h, int64_t cell_w, int64_t cell_h, std::vector<cv::Rect>& vr)
	{
		int64_t idx, jdx;

		int64_t xl, yl, xr, yr;

		cv::Rect r;

		vr.clear();

		for (idx = tile_h; idx < img_h; idx += tile_h)
		{

			for (jdx = tile_w; jdx < img_w; jdx += tile_w)
			{
				xl = std::max<int64_t>(0, jdx - cell_w - tile_w);
				yl = std::max<int64_t>(0, idx - cell_h - tile_h);

				xr = std::min<int64_t>(img_w, jdx + (2 * cell_w));
				yr = std::min<int64_t>(img_h, idx + (2 * cell_h));

				r = cv::Rect(cv::Point(xl, yl), cv::Point(xr, yr));

				vr.push_back(r);

			}

		}

	}
};

//Img_section::Img_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, std::int32_t large_cell_width, bool whole_img, const cv::Mat& current_tile_image, std::string html_img_dir, std::string h_img_dir, std::string templ_dir, Pixel_to_DEF_conversion convert_info, cv::Mat cell_template_img, int orientation)
//{
//
//
//}


void Img_section::setup_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, std::int32_t large_cell_width, bool whole_img, const cv::Mat& current_tile_image, Pixel_to_DEF_conversion convert_info, cv::Mat cell_template_img, int orientation)
{
	// Creates subsections inside each cell section, one for missed cells and one for errant cells
	// Steps:  searches for tiles that are aligned
	//         searches cell_list for cells on that tile that is aligned
	//         if there are points for that cell type that are marked "missed", create image and add to vector missed
	//         if there are points for that cell type that are marked "errant", create image and add to vector errant

	// Added DEBUG test to stop program if too many images are created - limit
	bool limit_images = true;
	int32_t num_images_created = 0;
	int32_t num_images_created_e = 0;

	// create the cell header info display
	create_template_header(cell, cell_template_img, orientation);

	// loop through the tiles
	for (int i = 0; i < tiles.size(); i++)
	{
		if (tiles[i].aligned)
		{
			// Read in image for tiles[i]
			cv::Mat img;
			if (whole_img)
				img = current_tile_image.clone();
			else
				img = cv::imread(tiles[i].image_file);
			
			for (Pt_info ptloc : cell.pt_locations)
			{
				// Make sure point is still on tile - otherwise can't display it
				if (ptloc.on_tile == i)
				{
					if ((ptloc.pt.y > tiles[i].imgLocRowY) & (ptloc.pt.y < tiles[i].imgLocRowY + tiles[i].height))
					{
						if ((ptloc.pt.x > tiles[i].imgLocColX) & (ptloc.pt.x < tiles[i].imgLocColX + tiles[i].width))
						{
							if (ptloc.matched == 2) // misssed cell
							{
								// Make sub_img here and push to vector "missed"
								Img_sub tmp_sub(ptloc, cell.name, "m", html_img_dir, h_img_dir);
								// Create image on html
								// LIMIT, only write image if it's not passed limit
								if (limit_images)
								{
									if (num_images_created < IMAGE_LIMIT)
										tmp_sub.write_image(ptloc.pt.x, ptloc.pt.y, cell.width, cell.height, img, tmp_sub.img_file_name, tiles[i].imgLocColX, tiles[i].imgLocRowY, "m", large_cell_width);
								}
								else
									tmp_sub.write_image(ptloc.pt.x, ptloc.pt.y, cell.width, cell.height, img, tmp_sub.img_file_name, tiles[i].imgLocColX, tiles[i].imgLocRowY, "m", large_cell_width);
								num_images_created++;
								missed.push_back(tmp_sub);
							}
							
						}
					}
				}
			}

			// ISSUES FOR ERRANT CELLS - WHAT TILE DO THEY FALL ON? 
			// The cell type passed in is the same as the detect cell type passed in
			
			for (int d = 0; d < cell_detect.size(); d++)
			{
				if (cell_detect[d].matched == -1)
				{
					if ((cell_detect[d].pt.x > tiles[i].imgLocColX) & (cell_detect[d].pt.x < tiles[i].imgLocColX + tiles[i].width))
					{
						if ((cell_detect[d].pt.y > tiles[i].imgLocRowY) & (cell_detect[d].pt.y < tiles[i].imgLocRowY + tiles[i].height))
						{
							Img_sub tmp_sub2(cell_detect[d], cell.name, "e", html_img_dir, h_img_dir);
							if (limit_images)
							{
								if (num_images_created_e < IMAGE_LIMIT)
									tmp_sub2.write_image(cell_detect[d].pt.x, cell_detect[d].pt.y, cell.width, cell.height, img, tmp_sub2.img_file_name, tiles[i].imgLocColX, tiles[i].imgLocRowY, "e", large_cell_width);
							}
							else
								tmp_sub2.write_image(cell_detect[d].pt.x, cell_detect[d].pt.y, cell.width, cell.height, img, tmp_sub2.img_file_name, tiles[i].imgLocColX, tiles[i].imgLocRowY, "e", large_cell_width);
							
							num_images_created_e++;
							
							// Convert coords of tmp_sub2 to DEF point
							cv::Point pixel_pt = cv::Point(tmp_sub2.coords.x, tmp_sub2.coords.y);
							
							// DON"T KNOW THE CORRECT CELL WIDTH SINCE IT"S ERRANT SO WON"T GET EXACT DEF LOCATION
							tmp_sub2.def_coords = convert_pixel_to_DEF_pt(convert_info, pixel_pt, cell.width, i, tiles[i].pixel_size_x, tiles[i].pixel_size_y, true);
							
							// Debug
							//cv::Point test_pt1;
							//test_pt1 = convert_pixel_to_DEF_pt(convert_info, cv::Point(2461, 2945), 99, i, tiles[i].pixel_size_x, tiles[i].pixel_size_y, true);
							
							errant.push_back(tmp_sub2);
						}
					}
				}
			}
		}
	}
}

void Img_section::setup_section(Cell_info cell, std::vector<Pt_info> cell_detect, const std::vector<Tile_img>& tiles, bool whole_img, const cv::Mat& current_tile_image, cv::Mat cell_template_img, int orientation)
{
	// Creates subsections inside each cell section, one for missed cells and one for errant cells
	// Steps:  searches for tiles that are aligned
	//         searches cell_list for cells on that tile that is aligned
	//         if there are points for that cell type that are marked "missed", create image and add to vector missed
	//         if there are points for that cell type that are marked "errant", create image and add to vector errant

	uint32_t idx, jdx;

	int32_t num_images_created = 0;
	int32_t num_images_created_e = 0;

	int64_t report_tile_w;
	int64_t report_tile_h;

	cv::Mat img;
	std::vector<cv::Mat> report_tiles;
	std::vector<cv::Rect> report_rects;

	// create the cell header info display
	create_template_header(cell, cell_template_img, orientation);

	// loop through the tiles
	for (idx = 0; idx < tiles.size(); ++idx)
	{
		// Read in image for tiles[idx]
		if (whole_img)
			img = current_tile_image.clone();
		else
			img = cv::imread(tiles[idx].image_file);

		// generate a set up cv::Rect for the report tiles from the image tile
	    report_tile_w = std::min((int64_t)(tile_width_percent*img.cols), max_tile_width);
	    report_tile_h = std::min((int64_t)(tile_height_percent*img.rows), max_tile_height);		
		create_report_tiles(img.cols, img.rows, report_tile_w, report_tile_h, cell.width, cell.height, report_rects);

		// go through each report rect and get the list of missed and errant cells for the image
		for (jdx = 0; jdx < report_rects.size(); ++jdx)
		{
			// missed


			// errants

		}


	}	// end for loop

}	// end of setup_section





//-----------------------------------------------------------------------------
//  Contains all information to write out html document
class Html_img
{
public:
	std::vector<Img_section> sections;
	std::string heading = "Individual Cell Analysis Results";

	Html_img(const std::vector<Cell_info> &cell_list, const std::vector<std::vector<Pt_info>>& detects, const std::vector<Tile_img> &tiles, bool whole_img, const cv::Mat& current_tile_image, std::string html_img_dir, std::string h_img_dir, std::string templ_dir, Pixel_to_DEF_conversion convert_info, std::vector<matched_template> matched_templates);
	friend std::ostream& operator << (std::ostream& out, Html_img img_tmp);
};

Html_img::Html_img(const std::vector<Cell_info> &cell_list, const std::vector<std::vector<Pt_info>> &detects, const std::vector<Tile_img> &tiles, bool whole_img, const cv::Mat& current_tile_image, std::string html_img_dir, std::string h_img_dir, std::string templ_dir, Pixel_to_DEF_conversion convert_info, std::vector<matched_template> matched_templates)
{
	// Make sections based on cell type
	for (int cnum = 0; cnum < cell_list.size(); cnum ++)
	{
		if (cell_list[cnum].selected == true) 
		{
			if (matched_templates[cnum].imgs.size() > 0)
			{
				// Creates section of images (missed and errant) by cell type (scc9gena_dfrtp_1, sccgena_mux2_2, etc)
				Img_section tmp_sect(cell_list[cnum], detects[cnum], tiles, cell_list[0].width, whole_img, current_tile_image, html_img_dir, h_img_dir, templ_dir, convert_info, matched_templates[cnum].get_mosaic(), matched_templates[cnum].orientation);
				sections.push_back(tmp_sect);
			}
			//else
				// ERROR
		}
	}
}

//******************************************************
// Write template images in a table
//******************************************************
void write_template_table(std::ostream& out, std::string templ_filename)
{
	
	out << "<table style = \" width: fit-content; border-collapse: collapse; border: 3px solid rgb(140 140 140):\">" << std::endl;
	out << "<tr style = \"height: 50px;\">" << std::endl;
	std::vector<std::string> snames{"N", "FN"};
	for (int i = 0; i < snames.size(); i++)
	{

		out << "<td style=\"text - align: center;  padding: 0px; border: 3px solid rgb(140, 140, 140);\">" << std::endl;
		out << "<img style = \"display: block;\" width = \"100%\" height = \"100%\" src = \"" + templ_filename + "_" + snames[i] + ".png\" alt = \"template image\" />" << std::endl;

		out << "</td>" << std::endl;
	}

	// Add back commented rows if you want 2 rows with 2 images each
	//out << "</tr>" << std::endl;
	//out << "<tr style = \"height: 50px;\">" << std::endl;

	std::vector<std::string> snames2{"S", "FS"};
	for (int i = 0; i < snames2.size(); i++)
	{
		//out << "<table style = \" width: fit-content; border-collapse: collapse; border: 2px solid rgb(140 140 140):\">";
		//out << "<tr style = \"height: 100px;\">";
		out << "<td style=\"text - align: center;  padding: 0px; border: 3px solid rgb(140, 140, 140);\">" << std::endl;
		out << "<img style = \"display: block;\" width = \"100%\" height = \"100%\" src = \"" + templ_filename + "_" + snames2[i] + ".png\" alt = \"test image\" />" << std::endl;

		out << "</td>" << std::endl;
	}

	out << "</tr>" << std::endl;
	out << "</table>" << std::endl;

}

//******************************************************
// Overloaded stream output operator << for image section          *
//******************************************************
std::ostream& operator << (std::ostream& out, Html_img img_tmp)
{
	//int limit_output = 25;
	out << "<div class=\"section_img_html\">" << std::endl;
	// Write image heading
	out << "<h2>" << img_tmp.heading << "</h2>" << std::endl;

	for (int i = 0; i < img_tmp.sections.size(); i++)
	{

			// Save image of template 
			//cv::imwrite(img_tmp.sections[i].templ_img_filename, img_tmp.sections[i].cell_templ_img);

			out << "<div class=\"img_section\">" << std::endl;
			out << "<h2" << img_tmp.sections[i].section_heading << "</h2>" << std::endl;
			out << "<h4>" << "Cell width: " << img_tmp.sections[i].section_subheading_width << "</h4>" << std::endl;
			out << "<h4>" << "Cell height: " << img_tmp.sections[i].section_subheading_height << "</h4>" << std::endl;
			// Add NEW template images
			write_template_table(out, img_tmp.sections[i].templ_filename);
			out << "<h4>" << "Binary Cell Templates:  Orientations N, FN, S, FS" << "</h4>" << std::endl;
			if (img_tmp.sections[i].missed.empty() & img_tmp.sections[i].errant.empty())
				out << "<h4>" << "Note:  No missed or errant cells." << "</h4>" << std::endl;
			// end template image

			// Magnifying factor for images - fix this based on reduction above
			int magnify = 1;  // Change if the missed and errant images on html are too small.
			if (!img_tmp.sections[i].missed.empty())
			{
				std::vector<Img_sub> missed = img_tmp.sections[i].missed;
				out << "<div class=\"img_sub\">" << std::endl;
				out << "<h3>Missed Cells</h3>" << std::endl;

				int limit_missed = IMAGE_LIMIT;
				limit_missed = (limit_missed < missed.size()) ? limit_missed : missed.size();

				for (int j = 0; j < limit_missed; j++)
				{
					std::string glob_heading = "Global Pixel Coordinates(x,y): " + std::to_string(missed[j].coords.x) + ", " + std::to_string(missed[j].coords.y);
					std::string def_heading = "DEF Coordinates(x,y): " + std::to_string(missed[j].def_coords.x) + ", " + std::to_string(missed[j].def_coords.y);
					out << "<h4>" << glob_heading << "</h4>" << std::endl;
					out << "<h4>" << def_heading << "</h4>" << std::endl;
					out << "<img" << std::endl;
					out << "class = \"fit-picture\"" << std::endl;
					out << "src = \"" << missed[j].himg_filename << "\"" << std::endl;
					out << "alt = \"Missed standard cell\"" << std::endl;
					// Defines image height/width on html
					out << "width = \"" << std::to_string(magnify * missed[j].html_img_width) << "\"" << std::endl;
					out << "height = \"" << std::to_string(magnify * missed[j].html_img_height) << "\"" << std::endl;
					out << "/>" << std::endl;
					out << "<form action=\"#Table1\"> <input type=\"submit\" value=\"Back to Table\" /> </form>" << std::endl;
				}
				out << "</div>" << std::endl;
			}

			
			if (!img_tmp.sections[i].errant.empty())
			{
				std::vector<Img_sub> errant = img_tmp.sections[i].errant;
				out << "<div class=\"img_sub_e\">" << std::endl;
				out << "<h3>Errant Cells</h3>" << std::endl;

				int limit_errant = IMAGE_LIMIT;
				limit_errant = (limit_errant < errant.size()) ? limit_errant : errant.size();

				for (int j = 0; j < limit_errant; j++)
				{
					std::string glob_heading = "Global Pixel Coordinates(x,y): " + std::to_string(errant[j].coords.x) + ", " + std::to_string(errant[j].coords.y);
					std::string def_heading = "Approximate DEF Coordinates(x,y): " + std::to_string(errant[j].def_coords.x) + ", " + std::to_string(errant[j].def_coords.y);
					out << "<h4>" << glob_heading << "</h4>" << std::endl;
					out << "<h4>" << def_heading << "</h4>" << std::endl;
					out << "<img" << std::endl;
					out << "class = \"fit-picture\"" << std::endl;
					out << "src = \"" << errant[j].himg_filename << "\"" << std::endl;
					out << "alt = \"Errant standard cell\"" << std::endl;
					out << "width = \"" << std::to_string(magnify * errant[j].html_img_width) << "\"" << std::endl;
					out << "height = \"" << std::to_string(magnify * errant[j].html_img_height) << "\"" << std::endl;
					out << "/>" << std::endl;
					out << "<form action=\"#Table1\"> <input type=\"submit\" value=\"Back to Table\" /> </form>" << std::endl;
				}
				out << "</div>" << std::endl;
			}
			out << "</div>" << std::endl;
		
	}
	out << "</div>" << std::endl;

	return out;
}

void write_toplines(std::ostream& out)
{
	out << "<html>" << std::endl;
	out << "<head>" << std::endl;
	out << "<style>" << std::endl;

	out << "div {" << std::endl;
  out << "border: 5px solid #555;" << std::endl;
	out << "background-color: #eee;" << std::endl;
	out << "padding: 0.5rem;" << std::endl;
	out << "display: flex;" << std::endl;
	out << "flex-direction: column;" << std::endl;
	out << "margin: 20px;" << std::endl;
	out << "}" << std::endl;

	// Headings
	out << "h1 {" << std::endl;
	out << "margin: 0.2rem 0.2rem 0 0;" << std::endl;
	out << "color: #071094cb;" << std::endl;
	out << "font-size: 1.5rem;" << std::endl;
	out << "font-family:Arial, Helvetica, sans-serif;" << std::endl;
	out << "text-align: center;" << std::endl;
	out << "}" << std::endl;

	out << "h2 {" << std::endl;
	out << "margin: 0.2rem 0.2rem 0 0;" << std::endl;
	out << "font-size: 1.0rem;" << std::endl;
	out << "font-family:Arial, Helvetica, sans-serif;" << std::endl;
	out << "text-align: left;" << std::endl;
	out << "padding: 0.5rem;" << std::endl;
	out << "}" << std::endl;

	out << "h3 {" << std::endl;
	out << "margin: 0.2rem 0.2rem 0 0;" << std::endl;
	out << "font-size: 0.9rem;" << std::endl;
	out << "font-family:Arial, Helvetica, sans-serif;" << std::endl;
	out << "text-align: left;" << std::endl;
	out << "font-weight: 900;" << std::endl;
	out << "}" << std::endl;

	out << "h4 {" << std::endl;
	out << "margin: 0.2rem 0.2rem 0 0;" << std::endl;
	out << "font-size: 0.80rem;" << std::endl;
	out << "font-family:'Courier New', Courier, monospace;" << std::endl;
	out << "text-align: left;" << std::endl;
	out << "}" << std::endl;

	// List items for unordered lists
	out << "li{" << std::endl;
	out << "list-style-type: disc;" << std::endl;
	out << "font-family:Arial, Helvetica, sans-serif;" << std::endl;
	out << "font-size: 0.9rem;" << std::endl;
	out << "}" << std::endl;

	// Top heading with dates
	out << ".top_heading {" << std::endl;
	out << "background-color: #fff;" << std::endl;
	out << "border: none;" << std::endl;
	out << "}" << std::endl;

	// Parameter Table sections
	out << ".section_params_table > h3 {" << std::endl;
	out << "margin: 0.4rem 0.2rem 0.2rem 0;" << std::endl;
	out << "color: #385985" << std::endl;
	out << "}" << std::endl;

	// Tables
	out << "thead {" << std::endl;
	out << "background-color: #3f5ba6;" << std::endl;
	out << "color: #fff;" << std::endl;
	out << "}" << std::endl;
	out << "tfoot {" << std::endl;
	out << "background-color: #3f5ba6;" << std::endl;
	out << "color: #fff;" << std::endl;
	out << "}" << std::endl;
	out << "tbody {" << std::endl;
	out << "background-color: #c9dbf7;" << std::endl;
	out << "}" << std::endl;

	out << "table {" << std::endl;
	out << "border-collapse: collapse;" << std::endl;
	out << "border: 2px solid rgb(200, 200, 200);" << std::endl;
	out << "letter-spacing: 1px;" << std::endl;
	out << "font-family: sans-serif;" << std::endl;
	out << "font-size: 0.8rem;" << std::endl;
	out << "width: fit-content;" << std::endl;
	out << "}" << std::endl;

	out << "td {" << std::endl;
	out << "border: 1px solid rgb(190, 190, 190);" << std::endl;
	out << "padding: 5px 10px;" << std::endl;
	out << "text-align: center;" << std::endl;
	out << "}" << std::endl;

	out << "th {" << std::endl;
	out << "border: 1px solid rgb(190, 190, 190);" << std::endl;
	out << "padding: 5px 25px;" << std::endl;
	out << "text-align: left;" << std::endl;
	out << "}" << std::endl;

	out << "tr:nth-child(even) {" << std::endl;
	out << "background-color: #dce7fa;" << std::endl;
	out << "}" << std::endl;

	/* Image section */
	out << ".img_section {" << std::endl;
	out << "border: 5px solid #555;" << std::endl;
	out << "background-color: #d5f1b575;" << std::endl;
	out << "}" << std::endl;

	out << ".img_section > h2 {" << std::endl;
	out << "color: #075314;" << std::endl;
	out << "}" << std::endl;

	out << ".img_sub {" << std::endl;
	out << "background-color: #f5f7f4e3;" << std::endl;
	out << "}" << std::endl;
	out << ".img_sub > h3 {" << std::endl;
	out << "color: #1205c2;" << std::endl;
	out << "}" << std::endl;

	out << ".img_sub {" << std::endl;
	out << "background-color: #f5f7f4e3;" << std::endl;
	out << "}" << std::endl;
	out << ".img_sub_e > h3 {" << std::endl;
	out << "color: #8f0808;" << std::endl;
	out << "}" << std::endl;

	out << ".fit-picture{" << std::endl;
	out << "margin: 0.2rem;" << std::endl;
	out << "}" << std::endl;

  // Last lines of style area
	out << "</style>" << std::endl;
	out << "</head>" << std::endl;
	out << "<body>" << std::endl;
}

void write_ending(std::ostream& out)
{
	out << "</body>" << std::endl;
	out << "</html>" << std::endl;
}



