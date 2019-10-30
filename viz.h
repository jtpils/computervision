// viz.h
// Date: 2019-29-10
// Created by: Gabriel Moreira

#ifndef VIZ_H
#define VIZ_H

#define NCLOUDS 5

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"
#include <mutex>
#include <future>

std::mutex g_display_mutex;

const unsigned int SCR_WIDTH  = 1800;
const unsigned int SCR_HEIGHT = 1400;

typedef struct Camera {
  glm::vec3 pos;
  glm::vec3 front;
  glm::vec3 up;
}Camera;

Camera camera = {glm::vec3(0.0f, 0.0f, 35.0f),
                 glm::vec3(0.0f, 0.0f, -1.0f),
                 glm::vec3(0.0f, 1.0f,  0.0f)};

bool firstMouse = true;
float yaw       = -90.0f;
float pitch     =  0.0f;
float lastX     =  800.0f / 2.0;
float lastY     =  600.0 / 2.0;
float fov       =  45.0f;

glm::mat4 model[NCLOUDS]  = {glm::mat4(1.0f), glm::mat4(1.0f), glm::mat4(1.0f), glm::mat4(1.0f), glm::mat4(1.0f)};
glm::mat4 view       =  glm::mat4(1.0f);
glm::mat4 projection =  glm::mat4(1.0f);


void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void processInput(GLFWwindow *window);
GLFWwindow* windowInit();
void windowLoop(GLFWwindow* window, float** vertices, int nVertices, int nClouds, int buffersize, cloudShader shader) ;
void drawThread(float **vertices, int nVertices, int nClouds, int buffersize);

#endif
