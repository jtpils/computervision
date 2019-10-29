#ifndef gpa_h
#define gpa_h


typedef struct srtTransformation {
  float s;
  Eigen::Matrix3f R;
  Eigen::Vector3f t;
}srtTransformation;


float r2d(float d);

float eigenMSE(Eigen::MatrixXf target, Eigen::MatrixXf current);

Eigen::Matrix3f rotmat(float alpha_deg, float beta_deg, float gamma_deg);

Eigen::MatrixXf* nearestNeighbors(Eigen::MatrixXf* clouds, int nClouds);

Eigen::MatrixXf vector2Eigen(std::vector<float> vec);

Eigen::MatrixXf srtWarp(Eigen::MatrixXf A, srtTransformation transf);

ANNpointArray matrix2ANN(Eigen::MatrixXf cloud);

void findGroups(Eigen::MatrixXf* corr, int nClouds, std::map<std::vector<int>, int>& m);

srtTransformation procrustes(Eigen::MatrixXf A, Eigen::MatrixXf B);



class Centroid {

  public:
    Centroid(Eigen::MatrixXf* clouds, int nClouds, std::map<std::vector<int>, int> m);
    Eigen::MatrixXf centroidByCloudNum(int cloudNumber);
    Eigen::MatrixXf cloudByCloudNum(int cloudNumber);
    bool isGroupValid(std::vector<int> group, int count);

  private:
    std::vector<Eigen::MatrixXf>* centralClouds;
    std::vector<Eigen::MatrixXf>* filteredClouds;

};


#endif
