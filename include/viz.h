// viz.h
// Date: 2019-29-10
// Created by: Gabriel Moreira

#ifndef VIZ_H
#define VIZ_H

#define NCLOUDS 10

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"
#include <mutex>
#include <future>

std::mutex g_display_mutex;

//------------------------------------------------------------------------------
// VERTEX AND ATTRIBUTE BUFFERS
unsigned int VBO;
unsigned int VAO;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// SCREEN
const unsigned int SCR_WIDTH  = 800;
const unsigned int SCR_HEIGHT = 500;
float ASPECT_RATIO = (float) SCR_WIDTH / (float) SCR_HEIGHT;
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

Camera camera = {glm::vec3(0.0f, 0.0f, 35.0f),
                 glm::vec3(0.0f, 0.0f, -1.0f),
                 glm::vec3(0.0f, 1.0f,  0.0f),
                 45.0f,
                 0.1f,
                 100.0f};
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// MOUSE CONTROL GLOBAL VARIABLES
bool firstMouse        = true;
bool leftButtonPressed = false;
float lastX            =  SCR_WIDTH / 2.0;
float lastY            =  SCR_HEIGHT / 2.0;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// COORDINATE SYSTEM TRANSFORMATION MATRICES
//------------------------------------------------------------------------------
glm::mat4 model[NCLOUDS]  = {glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f),
                             glm::mat4(1.0f)};

glm::mat4 view            =  glm::mat4(1.0f);
glm::mat4 projection      =  glm::mat4(1.0f);
//------------------------------------------------------------------------------

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
void windowLoop(GLFWwindow* window, float** vertices, int nVertices, int nClouds, int buffersize, cloudShader shader);
//------------------------------------------------------------------------------
void drawThread(float **vertices, int nVertices, int nClouds, int buffersize);
//------------------------------------------------------------------------------

#endif
