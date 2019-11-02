#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <chrono>

#include "reader.h"

#define DIM 3

Eigen::MatrixXf vectorToEigen(std::vector<float> vec) {
  // Converts std::vector<float> to Eigen::MatriXf assuming a row-major
  // storage order, i.e., [x0, y0, z0, x1, ..., zn] gets converted to
  // [[x0, y0, z0], [x1, y1, z1], ...] as an Eigen matrix.
  int j = 0;
  int i = 0;
  int nPts;

  Eigen::MatrixXf mat;
  nPts = ((int) vec.size() / DIM);
  mat = Eigen::MatrixXf::Zero(nPts, DIM);

  for(int k = 0; k < vec.size(); ++k) {
    j = k % DIM;
    i = ((int) k / DIM);
    mat(i, j) = vec[k];
  }
  return mat;
}


std::vector<float> eigenToVector(Eigen::MatrixXf mat) {
  // Converts std::vector<float> to Eigen::MatriXf assuming a row-major
  Eigen::MatrixXf auxMat(mat.rows(), mat.cols());
  auxMat = mat.transpose();
  std::vector<float> vec(auxMat.data(), auxMat.data() + auxMat.size());
  return vec;
}





int main() {
  std::vector<float> cloud0_vec;
  cloud0_vec = readObj("male_head.obj");

  Eigen::Matrix3f mat;
  mat << 1, 1, 1, 0, 1, 1, 0, 1, 1;

  Eigen::Matrix3f mat2;
  mat2 = mat.transpose();
  mat2(0,0) = 0;

  std::cout << mat << std::endl;
  std::cout << mat2 << "\n\n";
  std::cout << mat << std::endl;

  return 0;
}
