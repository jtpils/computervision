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
  Eigen::MatrixXf delta;
  Eigen::MatrixXf ddt;
  float mse;

  delta = target - current;
  ddt = delta * delta.transpose();
  mse = ddt.trace() / delta.rows();
  return mse;
}


Eigen::MatrixXf vector2Eigen(std::vector<float> vec) {
  int j = 0;
  int i = 0;
  int nPts;

  Eigen::MatrixXf mat;
  nPts = ((int) vec.size() / 3);
  mat = Eigen::MatrixXf::Zero(nPts, 3);

  for(int k = 0; k < vec.size(); ++k) {
    j = k % 3;
    i = ((int) k/3);
    mat(i, j) = vec[k];
  }
  return mat;
}


float r2d(float d) {
  return (d / 180.0) * ((float) M_PI);
};


Eigen::Matrix3f rotmat(float alpha_deg, float beta_deg, float gamma_deg) {
  Eigen::Matrix3f Rx;
  Eigen::Matrix3f Ry;
  Eigen::Matrix3f Rz;
  float alpha_rad = r2d(alpha_deg);
  float beta_rad = r2d(beta_deg);
  float gamma_rad = r2d(gamma_deg);

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
  Eigen::MatrixXf srtA;
  Eigen::MatrixXf transfA;
  srtA = transf.s * transf.R * A.transpose();
  srtA.colwise() += transf.t;
  transfA = srtA.transpose();
  return transfA;
};


srtTransformation procrustes(Eigen::MatrixXf A, Eigen::MatrixXf B) {
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


ANNpointArray matrix2ANN(Eigen::MatrixXf cloud) {
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
  int nPts = clouds[0].rows();
  int dim = clouds[0].cols();
  Eigen::MatrixXf* corr;
  corr = new Eigen::MatrixXf [nClouds];

  for(int i = 0; i < nClouds; i++) {
    corr[i] = Eigen::MatrixXf::Zero(nPts, nClouds);
  }

  for (int ki = 0; ki < nClouds; ki++){
    ANNkd_tree* kdTree;
    ANNpointArray dataPts = annAllocPts(nPts, dim);
    //std::cout<< "Converting point cloud..." << std::endl;
    dataPts = matrix2ANN(clouds[ki]);
    //std::cout<< "NN-search (two point-clouds)." << std::endl;
    kdTree = new ANNkd_tree(dataPts, nPts, dim);
    ANNidxArray nnIdx;					// near neighbor indices
    ANNdistArray dists;
    ANNpoint queryPt = annAllocPt(dim);
    for (int kj = 0; kj < nClouds; kj++){
      for (int i = 0; i < nPts; i++){
          queryPt[0] = clouds[kj](i,0);
          queryPt[1] = clouds[kj](i,1);
          queryPt[2] = clouds[kj](i,2);
          nnIdx = new ANNidx[1];
          dists = new ANNdist[1];
          kdTree->annkSearch(queryPt, 1, nnIdx, dists, 0);
          corr[kj].operator()(i,ki) = ((int) nnIdx[0]);
          delete [] nnIdx;			// clean things up
          delete [] dists;
        }
    }
    delete kdTree;
    annClose();
  }
  return corr;
};



void findGroups(Eigen::MatrixXf* corr, int nClouds, std::map<std::vector<int>, int>& m) {
  //std::cout << "Finding MUTUAL neighbor pairs..." << '\n';
  int nPts = corr[0].rows();
  int ind_in_i;
  int ind_in_j;

  for (int i = 0; i < nClouds; i++) {
    for (int p = 0; p < nPts; p++) {
      std::vector<int> ind = {};
      int nMutuals = 0;
      for (int j = 0; j < nClouds; j++) {
        ind_in_j = corr[i](p,j);
        ind_in_i = corr[j](ind_in_j, i);

        if (p == ind_in_i) {
          nMutuals++;
          ind.insert(ind.end(), ind_in_j);
        }
        else {
          ind.insert(ind.end(), -1);
        };
      };
      if (nMutuals > 1) {
        ++m[ind];
      };
    };
  };
};


Centroid::Centroid(Eigen::MatrixXf* clouds, int nClouds, std::map<std::vector<int>, int> m) {
  //std::cout << "Computing centroid..." << '\n';
  int count;
  int cloudNumber;
  int groupNumber;
  std::vector<int> validClouds;
  std::vector<int> validIndexes;
  Eigen::MatrixXf groupCentroid(1,3);

  centralClouds = new std::vector<Eigen::MatrixXf> [nClouds];
  filteredClouds = new std::vector<Eigen::MatrixXf> [nClouds];
  groupNumber = 0;

  // iterate over all the groups in the map
  for(std::map<std::vector<int>, int>::iterator iter = m.begin(); iter != m.end(); ++iter) {
    std::vector<int> group = iter->first; // get the group itself
    count = iter->second; // get the number of occurences of this group
    cloudNumber = 0;
    validClouds = {};
    groupCentroid = Eigen::MatrixXf::Zero(1,3);

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
  numCloudsInGroup = group.size() - std::count(group.begin(), group.end(), -1);

  if (numCloudsInGroup == count) {
    return true;
  }
  else {
    return false;
  }
};


Eigen::MatrixXf Centroid::centroidByCloudNum(int cloudNumber) {
  //std::cout << "Getting centroid generated by cloud No. " << cloudNumber << ".\n";
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(centralClouds[cloudNumber].size(), 3);
  float* pt;

  for (int i = 0; i < centralClouds[cloudNumber].size(); i++) {
    pt = centralClouds[cloudNumber][i].data();
    Eigen::MatrixXf mat_row = Eigen::Map<Eigen::Matrix<float, 1, 3>>(pt);
    mat.row(i) += mat_row;
  }
  return mat;
};


Eigen::MatrixXf Centroid::cloudByCloudNum(int cloudNumber) {
  //std::cout << "Retrieving the points of cloud" << cloudNumber << " w/ mutual NNs." << "\n";
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(filteredClouds[cloudNumber].size(), 3);
  float* pt;

  for (int i = 0; i < filteredClouds[cloudNumber].size(); i++) {
    pt = filteredClouds[cloudNumber][i].data();
    Eigen::MatrixXf mat_row = Eigen::Map<Eigen::Matrix<float, 1, 3>>(pt);
    mat.row(i) += mat_row;
  }
  return mat;
};
