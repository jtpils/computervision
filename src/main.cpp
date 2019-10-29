// main.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <ANN/ANN.h>
#include <map>

#include "gpa.h"
#include "reader.h"


int main()
{
  int dim = 3;
  int nClouds = 3;
  Eigen::MatrixXf* clouds;

  float s01;
  float s02;
  Eigen::Matrix3f R01;
  Eigen::Matrix3f R02;
  Eigen::Vector3f t01;
  Eigen::Vector3f t02;

  R01 = rotmat(10,0,0);
  s01 = 1;
  t01 << 2, 0, 0;
  struct srtTransformation transf01 = {s01, R01, t01};

  R02 = rotmat(0,-12,0);
  s02 = 1;
  t02 << 0, 1, 0;
  struct srtTransformation transf02 = {s02, R02, t02};


  clouds = new Eigen::MatrixXf [nClouds];
  std::vector<float> cloud0_vec;
  cloud0_vec = readObj("../data/male_head.obj");

  clouds[0] = vector2Eigen(cloud0_vec);
  clouds[1] = srtWarp(clouds[0], transf01);
  clouds[2] = srtWarp(clouds[1], transf02);


  /*
   * Generalized Procrustes Analysis starts here
   */
  int MAX_ITER = 50;
  for (int iiter = 0; iiter < MAX_ITER; ++iiter) {
    std::cout << "Epoch: "  << iiter << " \n";
    Eigen::MatrixXf* indexMatrix;
    indexMatrix = nearestNeighbors(clouds, nClouds);

    std::map<std::vector<int>, int> groupMap;
    findGroups(indexMatrix, nClouds, groupMap);
    delete [] indexMatrix;
    Centroid centroid(clouds, nClouds, groupMap);

    Eigen::MatrixXf current;
    Eigen::MatrixXf target;
    float mse;
    struct srtTransformation srt;
    for (int i = 0; i < nClouds; ++i) {
      target = centroid.centroidByCloudNum(i);
      current = centroid.cloudByCloudNum(i);
      mse = eigenMSE(target, current);
      std::cout << " Loss (MSE) (cloud" << i << "): "  << mse;
      std::cout << " (using " << target.rows() << " mutual neighbors)" << '\n';
      srt = procrustes(current, target);
      clouds[i] = srtWarp(clouds[i], srt);
    }
  }
  /*
  std::cout << "C0:" << std::endl;
  std::cout << c0.rows() << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "C1:" << std::endl;
  std::cout << c1.rows() << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "C2:" << std::endl;
  std::cout << c2.rows() << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "PC0:" << std::endl;
  std::cout << pc0.rows() << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "PC1:" << std::endl;
  std::cout << pc1.rows() << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << "PC2:" << std::endl;
  std::cout << pc2.rows() << std::endl;Vector3d
  int count;
  for(std::map<std::vector<int>, int>::iterator iter = groupMap.begin(); iter != groupMap.end(); ++iter) {
    std::vector<int> vec = iter->first;
    count = iter->second;
    for (std::vector<int>::const_iterator j = vec.begin(); j != vec.end(); ++j) {
      std::cout << *j << ' ';
    }
    std::cout << count;
    std::cout << '\n';
  }
  */
  delete [] clouds;
  return 0;
}
