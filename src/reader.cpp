// reader.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "reader.h"

/*
cv::Mat readRGB(std::string imagePath) {
  cv::Mat image = cv::imread(imagePath);

  // Check for failure
  if (image.empty()) {
    std::string errorMsg = "Could not open or find " + imagePath +
                           "\nPress any key to continue.";
    std::cout << errorMsg << std::endl;
    std::cin.get();
  }

  return image;
};


cv::Mat readDepth(std::string depthPath) {
  int flags = CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_ANYCOLOR;
  cv::Mat depth = cv::imread(depthPath, flags);
  depth.convertTo(depth, CV_32F);

  // Check for failure
  if (depth.empty()) {
    std::string errorMsg = "Could not open or find " + depthPath +
                           "\nPress any key to continue.";
    std::cout << errorMsg << std::endl;
    std::cin.get();
  }

  return depth;
};
*/

std::vector<float> readObj(const char *filename) {

  std::string line;
  std::ifstream myfile;
  const char * ptr;
  std::vector<float> vec = {};
  float vx, vy, vz;

  myfile.open(filename);

  if(!myfile.is_open()) {
    perror("ERROR OPENING .OBJ FILE.");
    exit(EXIT_FAILURE);
  }

  while(getline(myfile, line)) {
    ptr = line.c_str();
    if (*ptr == 'v' && *(ptr+1) == ' ') {
      std::istringstream iss(++ptr);
      iss >> vx >> vy >> vz;
      vec.insert(vec.end(), vx);
      vec.insert(vec.end(), vy);
      vec.insert(vec.end(), vz);
    }
  }
  return vec;
};
