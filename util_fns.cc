#include <algorithm>
#include <cmath>
#include <sstream>
#include "util_fns.h"

namespace utils {

  std::vector< std::vector<int> > read_data(std::ifstream& input_file, const char delimiter) {
    /*

      Saves each row of data as a column in output_data
    
      The data is expected to be of the form:
      x1 y1 z1 ...
      x2 y2 z2 ...
      ...
      thus each "column" of output_data contains all
      input data of one variable

    */
    // determine number of columns by reading first line
    // count rows, not columns
    std::string line;
    int nrows = 0;
    while(std::getline(input_file, line)) {
      nrows++;
    }
    // move back to beginning of file
    // need to clear EOF flag first
    input_file.clear();
    input_file.seekg(0);
    std::string val;
    std::vector< std::vector<int> > output_data(nrows);
    int j = 0;
    while(std::getline(input_file, line)) {
      std::stringstream ss(line);
      while(std::getline(ss, val, delimiter)) {
	output_data[j].push_back(std::atoi(val.c_str()));
      }
      j++;
    }
    return output_data;
  }

} 
