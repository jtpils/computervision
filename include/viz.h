// viz.h
// Date: 2019-29-10
// Created by: Gabriel Moreira

#ifndef VIZ_H
#define VIZ_H

#include <future>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <Eigen/Dense>

//#include <opencv2/opencv.hpp>

#include "shader.h"
#include "common.h"


//------------------------------------------------------------------------------
// VERTEX AND ATTRIBUTE BUFFERS
static unsigned int VBO;
static unsigned int VAO;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// SCREEN GLOBALS
const unsigned int SCR_WIDTH  = 800;
const unsigned int SCR_HEIGHT = 500;
const float ASPECT_RATIO      = (float) SCR_WIDTH / (float) SCR_HEIGHT;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// CAMERA'S COORDINATE SYSTEM AXIS INITIALIZATION
typedef struct Camera {
  glm::vec3 pos;
  glm::vec3 front;
  glm::vec3 up;
  float fov;
  float nearField;
  float farField;
}Camera;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
// MOUSE INPUT CONTROL
typedef struct Mouse {
  bool firstMouse;
  bool leftButtonPressed;
  float lastX;
  float lastY;
}Mouse;
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void genViewMat(Camera camera);
//------------------------------------------------------------------------------
void genModelMat(int modelNum, glm::mat4 newModel);
//------------------------------------------------------------------------------
void genProjectionMat(Camera camera);
//------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
//------------------------------------------------------------------------------
void leftButtonPressedRotation();
//------------------------------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
//------------------------------------------------------------------------------
void processInput(GLFWwindow *window);
//------------------------------------------------------------------------------
GLFWwindow* windowInit();
//------------------------------------------------------------------------------
void windowLoop(GLFWwindow* window, float** vertices, int nVertices, int nClouds, cloudShader shader);
//------------------------------------------------------------------------------
void drawThread(float **vertices, int nVertices, int nClouds);
//------------------------------------------------------------------------------
void clearScreen();
//------------------------------------------------------------------------------
void updateModelsMatrix(Eigen::MatrixXf* clouds, int nClouds);
//------------------------------------------------------------------------------
//cv::Mat generateCalibrationMat(float fx, float fy, float x0, float y0);
//------------------------------------------------------------------------------
Camera cameraInit();
//------------------------------------------------------------------------------
Mouse mouseInit();
//------------------------------------------------------------------------------



#endif
