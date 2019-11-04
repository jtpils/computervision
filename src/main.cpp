#include <vector>
#include <string>

#include "reader.h"
#include "viz.h"
#include "gpa.h"
#include "common.h"


std::mutex g_display_mutex; // Declared in common.h, needs a definition here

int main()
{
    for (int i = 0; i < NCLOUDS; ++i) {
      genModelMat(i, glm::mat4(1.0f));
    }

    Eigen::MatrixXf* clouds;
    clouds = new Eigen::MatrixXf [NCLOUDS];
    std::vector<float> cloud_vec;
    std::string root = "../data/";
    std::string path;
    std::string view_no;

    for(int i = 0; i < NCLOUDS; ++i) {
      view_no = std::to_string(i);
      path = root + "bunnyview" + view_no + ".obj";
      cloud_vec = readObj(path.c_str());
      clouds[i] = vectorToEigen(cloud_vec);
      clouds[i] /= 20;
    }


    /*
    cv::Mat c = readDepth("../data/heads/seq-01/frame-000000.depth.png");
    cv::Mat K = generateCalibrationMat(585, 585, 320, 240);
    cv::Mat cloud;
    cv::rgbd::depthTo3d(c, K, cloud);
    */

    int nVertices = clouds[0].rows();

    float *vertices[NCLOUDS];
    std::vector<float> vecs[NCLOUDS];

    for (int i = 0; i < NCLOUDS; ++i) {
      vecs[i] = eigenToVector(clouds[i]);
      vertices[i] = vecs[i].data();
    }

    std::future<void> drawResult(std::async(drawThread, vertices, nVertices, NCLOUDS));
    std::cin.get();
    std::future<void> updateModelsMatrixResult(std::async(updateModelsMatrix, clouds, NCLOUDS));

    drawResult.get();
    updateModelsMatrixResult.get();

    return 0;
};
