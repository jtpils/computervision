// gpa.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <vector>
#include <map>
#include <math.h>
#include <Eigen/Dense>

#include <ANN/ANN.h>
#include "gpa.h"


float eigenMSE(Eigen::MatrixXf target, Eigen::MatrixXf current) {
  // Computes the mean-square error between two point matrices
  Eigen::MatrixXf delta;
  Eigen::MatrixXf ddt;
  float mse;

  delta = target - current;
  ddt = delta * delta.transpose();
  mse = ddt.trace() / delta.rows();
  return mse;
}


std::vector<float> eigenToVector(Eigen::MatrixXf mat) {
  // Converts Eigen::Matrix to std::vector<float> assuming a row-major
  // storage order. The data is fully copied.
  Eigen::MatrixXf auxMat(mat.rows(), mat.cols());
  auxMat = mat.transpose();
  std::vector<float> vec(auxMat.data(), auxMat.data() + auxMat.size());
  return vec;
}


Eigen::MatrixXf vectorToEigen(std::vector<float> vec) {
  // Converts std::vector<float> to Eigen::MatriXf assuming a row-major
  // storage order, i.e., [x0, y0, z0, x1, ..., zn] gets converted to
  // [[x0, y0, z0], [x1, y1, z1], ...] as an Eigen matrix. Data is fully
  // copied (deep copy).
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


float degreesToRadians(float d) {
  // Converts an input in degrees to radians
  return (d / 180.0) * ((float) M_PI);
};


Eigen::Matrix3f rotMat(float alpha_deg, float beta_deg, float gamma_deg) {
  // Creates a 3x3 rotation matrix according to the provided rotation angles
  // about the z-axis, y-axis and x-axis (according to the order of the arguments)
  Eigen::Matrix3f Rx;
  Eigen::Matrix3f Ry;
  Eigen::Matrix3f Rz;
  float alpha_rad = degreesToRadians(alpha_deg);
  float beta_rad = degreesToRadians(beta_deg);
  float gamma_rad = degreesToRadians(gamma_deg);

  Rx << 1, 0, 0,
        0, cos(gamma_rad), sin(gamma_rad),
        0, sin(gamma_rad), cos(gamma_rad);

  Ry << cos(beta_rad), 0, sin(beta_rad),
        0, 1, 0,
        -sin(beta_rad), 0, cos(beta_rad);

  Rz << cos(alpha_rad), -sin(alpha_rad), 0,
        sin(alpha_rad), cos(alpha_rad), 0,
        0, 0, 1;

  Eigen::Matrix3f R = Rz*Ry*Rx;

  return R;
};


Eigen::MatrixXf srtWarp(Eigen::MatrixXf A, srtTransformation transf) {
  // Transforms a nPts x 3 point matrix A according to a SRT transformation
  Eigen::MatrixXf srtA;
  Eigen::MatrixXf transfA;
  srtA = transf.s * transf.R * A.transpose();
  srtA.colwise() += transf.t;
  transfA = srtA.transpose();
  return transfA;
};


srtTransformation procrustes(Eigen::MatrixXf A, Eigen::MatrixXf B) {
  // Solves the least squares Procrustes problem of finding the SRT
  // srtTransformation that minimizes the error E = sAR + tj - B
  srtTransformation transf = {};

  Eigen::Vector3f centroid_A = A.colwise().mean();
  Eigen::Vector3f centroid_B = B.colwise().mean();

  A.rowwise() -= centroid_A.transpose();
  B.rowwise() -= centroid_B.transpose();

  Eigen::MatrixXf C = A.transpose() * B;

  Eigen::JacobiSVD<Eigen::MatrixXf> svd(C, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3f U = svd.matrixU();
  Eigen::Matrix3f V = svd.matrixV();

  transf.R = V * U.transpose();

  if (transf.R.determinant() < 0) {
    Eigen::Matrix3f reflex;
    reflex << 1, 0, 0,
              0, 1, 0,
              0, 0, transf.R.determinant();
    transf.R = U * reflex * V.transpose();
  };

  transf.t = centroid_B - transf.R * centroid_A;
  transf.s = 1;

  return transf;
};


ANNpointArray EigenMatrixToANNpointArray(Eigen::MatrixXf cloud) {
  // Converts an EigenMatrix to an ANNpointArray to be used when
  // computing nearest neighbors.
  int nPts = cloud.rows();
  int dim = cloud.cols();
  ANNpointArray dataPts;
  dataPts = annAllocPts(nPts, dim);
  for(int i = 0; i < nPts; i++){
    ANNpoint p = annAllocPt(dim);
    p[0] = ((float) cloud(i,0));
    p[1] = ((float) cloud(i,1));
    p[2] = ((float) cloud(i,2));
    dataPts[i] = p;
  };
  return dataPts;
};


Eigen::MatrixXf* nearestNeighbors(Eigen::MatrixXf* clouds, int nClouds) {
  // Creates nClouds correspondence matrices with the indices of the
  // nearest neighbors of each cloud, in all the other clouds.
  int nPts = clouds[0].rows();
  int dim = clouds[0].cols();
  Eigen::MatrixXf* corr;
  corr = new Eigen::MatrixXf [nClouds];           // correspondence matrix

  // Initialize the correspondence matrices
  for(int i = 0; i < nClouds; i++) {
    corr[i] = Eigen::MatrixXf::Zero(nPts, nClouds);
  }

  // Allocate stuff before the search loop
  ANNidxArray nnIdx = new ANNidx[1];					    // near neighbor index
  ANNdistArray dists = new ANNdist[1];            // near neighbor distance
  ANNpoint queryPt = annAllocPt(dim);             // point we lookin for
  ANNpointArray dataPts = annAllocPts(nPts, dim); // search space

  for (int ki = 0; ki < nClouds; ki++) {
    ANNkd_tree* kdTree;
    dataPts = EigenMatrixToANNpointArray(clouds[ki]);
    kdTree = new ANNkd_tree(dataPts, nPts, dim);  // k-d tree init

    for (int kj = 0; kj < nClouds; kj++) {
      for (int i = 0; i < nPts; i++) {
          queryPt[0] = clouds[kj](i,0);
          queryPt[1] = clouds[kj](i,1);
          queryPt[2] = clouds[kj](i,2);
          kdTree->annkSearch(queryPt, 1, nnIdx, dists, 0);
          corr[kj].operator()(i,ki) = ((int) nnIdx[0]);
      }
    }
    delete kdTree;                                // clean things up
    annClose();
  }
  delete [] nnIdx;			                          // clean things up
  delete [] dists;                                // clean things up
  return corr;
};


void findGroups(Eigen::MatrixXf* corr, int nClouds, groupmap& m) {
  // Finds groups of points which are mutually nearest neighbors. A map m holds
  // the information about the groups. The keys are the std::vectors describing
  // the group, and the values are the number of occurences of each group.
  int nPts = corr[0].rows();
  int ind_in_i;                        // cloud-i index of the nearest point
  int ind_in_j;                        // cloud-j index of the nearest point
  int nMutuals;                        // no. of mutually nearest nbrs in a group

  for (int i = 0; i < nClouds; i++) {
    for (int p = 0; p < nPts; p++) {

      std::vector<int> ind = {};       // mutually nearest neighbors indices
      nMutuals = 0;

      for (int j = 0; j < nClouds; j++) {
        ind_in_j = corr[i](p,j);
        ind_in_i = corr[j](ind_in_j, i);

        if (p == ind_in_i) {           // mutually nearest neighbor condition
          nMutuals++;
          ind.insert(ind.end(), ind_in_j);
        }
        else {
          ind.insert(ind.end(), -1);   // if no mutually nearest point, add -1
        };
      };

      if (nMutuals > 1)
        ++m[ind];                      // add to the map if there's at least 2

    };
  };
};


Centroid::Centroid(Eigen::MatrixXf* clouds, int nClouds, groupmap m) {
  int count;
  int cloudNumber;
  int groupNumber;
  std::vector<int> validClouds;
  std::vector<int> validIndexes;
  std::vector<int> group;
  Eigen::MatrixXf groupCentroid(1,DIM);

  centralClouds = new std::vector<Eigen::MatrixXf> [nClouds];
  filteredClouds = new std::vector<Eigen::MatrixXf> [nClouds];
  groupNumber = 0;

  // iterate over all the groups in the map
  for(std::map<std::vector<int>, int>::iterator iter = m.begin(); iter != m.end(); ++iter) {
    group = iter->first;           // get the group itself
    count = iter->second;          // get the number of occurences of this group
    cloudNumber = 0;
    validClouds = {};
    groupCentroid = Eigen::MatrixXf::Zero(1,DIM);

    // check if group occurences equals number of clouds in the group
    if (isGroupValid(group, count)) {
      // iterate over all the indices in the group
      for (std::vector<int>::const_iterator j = group.begin(); j != group.end(); ++j) {
        if (*j != -1) {
          validClouds.insert(validClouds.end(), cloudNumber);
          groupCentroid.row(0) += clouds[cloudNumber].row(*j);
          filteredClouds[cloudNumber].insert(filteredClouds[cloudNumber].end(), clouds[cloudNumber].row(*j));
        }
        ++cloudNumber;
      };

      groupCentroid = groupCentroid / validClouds.size();

      for (std::vector<int>::const_iterator k = validClouds.begin(); k != validClouds.end(); ++k)
        centralClouds[*k].insert(centralClouds[*k].end(), groupCentroid);
    }
    ++groupNumber;
  };
};


bool Centroid::isGroupValid(std::vector<int> group, int count) {
  int numCloudsInGroup;
  // Get number of clouds in the group (number of entries different from -1)
  numCloudsInGroup = group.size() - std::count(group.begin(), group.end(), -1);

  if (numCloudsInGroup == count) {
    return true;
  }
  else {
    return false;
  }
};


Eigen::MatrixXf Centroid::centroidByCloudNum(int cloudNumber) {
  // Get a matrix with the centroid points calculated with cloud cloudsNumber
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(centralClouds[cloudNumber].size(), DIM);
  float* pt;

  for (int i = 0; i < centralClouds[cloudNumber].size(); i++) {
    pt = centralClouds[cloudNumber][i].data();
    Eigen::MatrixXf mat_row = Eigen::Map<Eigen::Matrix<float, 1, 3>>(pt);
    mat.row(i) += mat_row;
  }
  return mat;
};


Eigen::MatrixXf Centroid::cloudByCloudNum(int cloudNumber) {
  // Get a matrix with the points from cloud cloudNumber which were used to
  // compute the centroid
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(filteredClouds[cloudNumber].size(), DIM);
  float* pt;

  for (int i = 0; i < filteredClouds[cloudNumber].size(); i++) {
    pt = filteredClouds[cloudNumber][i].data();
    Eigen::MatrixXf mat_row = Eigen::Map<Eigen::Matrix<float,1,DIM> >(pt);
    mat.row(i) += mat_row;
  }
  return mat;
};


void gpaipc(Eigen::MatrixXf* clouds, int nClouds, int maxIter, bool verbose) {
  float meanSquaredError = 0;
  srtTransformation srt;
  Eigen::MatrixXf* indexMatrix;

  for (int iiter = 0; iiter < maxIter; ++iiter) {
    if (verbose == true)
      std::cout << "Epoch: "  << iiter << " \n";

    indexMatrix = nearestNeighbors(clouds, nClouds);
    groupmap groups;
    findGroups(indexMatrix, nClouds, groups);
    delete [] indexMatrix;
    Centroid centroid(clouds, nClouds, groups);

    Eigen::MatrixXf current;
    Eigen::MatrixXf target;

    // run procrustes for all the clouds with the centroid
    for (int i = 0; i < nClouds; ++i) {
      target = centroid.centroidByCloudNum(i);
      current = centroid.cloudByCloudNum(i);
      meanSquaredError = eigenMSE(target, current);

      if (verbose == true) {
        std::cout << " Loss (MSE) (cloud" << i << "): "  << meanSquaredError;
        std::cout << " (using " << target.rows() << " mutual neighbors)" << '\n';
      }

      srt = procrustes(current, target);
      clouds[i] = srtWarp(clouds[i], srt);
    }
  }
};
