// gpa.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#ifndef gpa_h
#define gpa_h

#include <iostream>
#include <map>

#include <Eigen/Dense>
#include <ANN/ANN.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

typedef std::map<std::vector<int>, int> groupmap;


//------------------------------------------------------------------------------
typedef struct srtTransformation {
  float s;                          // Scale
  Eigen::Matrix3f R;                // Rotation
  Eigen::Vector3f t;                // Translation
}srtTransformation;
//-------------------------------------------------------------------------------


//------------------------------------------------------------------------------
float radiansToDegrees(float d);
//------------------------------------------------------------------------------
float eigenMSE(Eigen::MatrixXf target, Eigen::MatrixXf current);
//------------------------------------------------------------------------------
Eigen::Matrix3f rotMat(float alpha_deg, float beta_deg, float gamma_deg);
//------------------------------------------------------------------------------
Eigen::MatrixXf* nearestNeighbors(Eigen::MatrixXf* clouds, int nClouds);
//------------------------------------------------------------------------------
Eigen::MatrixXf vectorToEigen(std::vector<float> vec);
//------------------------------------------------------------------------------
std::vector<float> eigenToVector(Eigen::MatrixXf mat);
//------------------------------------------------------------------------------
Eigen::MatrixXf srtWarp(Eigen::MatrixXf A, srtTransformation transf);
//------------------------------------------------------------------------------
ANNpointArray EigenMatrixToANNpointArray(Eigen::MatrixXf cloud);
//------------------------------------------------------------------------------
void findGroups(Eigen::MatrixXf* corr, int nClouds, groupmap& m);
//------------------------------------------------------------------------------
srtTransformation procrustes(Eigen::MatrixXf A, Eigen::MatrixXf B);
//------------------------------------------------------------------------------
void centerClouds(Eigen::MatrixXf* clouds, int nClouds);
//------------------------------------------------------------------------------
void gpaipc(Eigen::MatrixXf* clouds, int nClouds, int maxIter, glm::mat4* model, bool verbose);
//------------------------------------------------------------------------------


class Centroid {
  public:
    //--------------------------------------------------------------------------
    Centroid(Eigen::MatrixXf* clouds, int nClouds, groupmap m);
    //--------------------------------------------------------------------------
    Eigen::MatrixXf centroidByCloudNum(int cloudNumber);
    //--------------------------------------------------------------------------
    Eigen::MatrixXf cloudByCloudNum(int cloudNumber);
    //--------------------------------------------------------------------------
    bool isGroupValid(std::vector<int> group, int count);
    //--------------------------------------------------------------------------

  private:
    std::vector<Eigen::MatrixXf>* centralClouds;
    std::vector<Eigen::MatrixXf>* filteredClouds;
};

#endif
