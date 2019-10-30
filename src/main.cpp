// main.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <iostream>
#include <vector>

#include "gpa.h"
#include "reader.h"

int main()
{
  int nClouds = 3;
  Eigen::MatrixXf* clouds;

  float s01;
  float s02;
  Eigen::Matrix3f R01;
  Eigen::Matrix3f R02;
  Eigen::Vector3f t01;
  Eigen::Vector3f t02;

  R01 = rotMat(10,0,0);
  s01 = 1;
  t01 << 2, 0, 0;
  srtTransformation transf01 = {s01, R01, t01};

  R02 = rotMat(0,-12,0);
  s02 = 1;
  t02 << 0, 1, 0;
  srtTransformation transf02 = {s02, R02, t02};

  clouds = new Eigen::MatrixXf [nClouds];
  std::vector<float> cloud0_vec;
  cloud0_vec = readObj("../data/male_head.obj");

  clouds[0] = vectorToEigen(cloud0_vec);
  clouds[1] = srtWarp(clouds[0], transf01);
  clouds[2] = srtWarp(clouds[1], transf02);

  gpaipc(clouds, 3, 50, true);

  delete [] clouds;
  return 0;
}
