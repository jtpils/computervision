// viz.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <iostream>
#include <vector>

#include "reader.h"
#include "gpa.h"
#include "viz.h"
#include <future>
#include <mutex>

unsigned int VBO, VAO;

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    float sensitivity = 0.05;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw   += xoffset;
    pitch += yoffset;

    if(pitch > 89.0f)
        pitch = 89.0f;
    if(pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 front;
    front.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    front.y = sin(glm::radians(pitch));
    front.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    camera.front = glm::normalize(front);
}


void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    float cameraSpeed = 0.25f; // adjust accordingly
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
      camera.pos += cameraSpeed * camera.front;
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
      camera.pos -= cameraSpeed * camera.front;
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
      camera.pos -= glm::normalize(glm::cross(camera.front, camera.up)) * cameraSpeed;
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
      camera.pos += glm::normalize(glm::cross(camera.front, camera.up)) * cameraSpeed;
}


GLFWwindow* windowInit() {
  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Nabla", NULL, NULL);
  if (window == NULL) {
      std::cout << "Failed to create GLFW window" << std::endl;
      glfwTerminate();
  }
  glfwMakeContextCurrent(window);
  glfwSetCursorPosCallback(window, mouse_callback);
  glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
      std::cout << "Failed to initialize GLAD" << std::endl;
  }

  glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
  return window;
}


void windowLoop(GLFWwindow* window, float** vertices, int nVertices, int nClouds, int buffersize, cloudShader shader) {
  while(!glfwWindowShouldClose(window)) {
    // Check for user input
    processInput(window);

    // Rendering commands here
    glClearColor(0.0f, 0.f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    shader.use();

    view = glm::lookAt(camera.pos, camera.pos + camera.front, camera.up);
    projection = glm::perspective(glm::radians(55.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);

    // Pass matrices to the shader program
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);

    for (int i = 0; i < nClouds; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, VBO);
      glBufferData(GL_ARRAY_BUFFER, buffersize, vertices[i], GL_DYNAMIC_DRAW);
      glBindVertexArray(VAO);

      g_display_mutex.lock();
      shader.setMat4("model", model[i]);
      glDrawArrays(GL_POINTS, 0, nClouds*nVertices);
      g_display_mutex.unlock();
    }
    // Callbacks and swap the rendering buffers
    glfwPollEvents();
    glfwSwapBuffers(window);
  }
}


void drawThread(float **vertices, int nVertices, int nClouds, int buffersize) {
  GLFWwindow* window = windowInit();
  cloudShader shader("vertexShader.vs", "fragmentShader.fs");

  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);

  // bind the Vertex Array Object first, then bind and set vertex buffer(s),
  // and then configure vertex attributes(s).
  glBindVertexArray(VAO);
  glBindBuffer(GL_ARRAY_BUFFER, VBO);

  glBufferData(GL_ARRAY_BUFFER, buffersize, vertices[0], GL_DYNAMIC_DRAW);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  shader.use();
  windowLoop(window, vertices, nVertices, nClouds, buffersize, shader);

  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
  glfwTerminate();
}


void updateModelsMatrix(Eigen::MatrixXf* clouds, int nClouds) {
  float meanSquaredError = 0;
  srtTransformation srt;
  Eigen::MatrixXf* indexMatrix;
  int maxIter = 300;
  for (int iiter = 0; iiter < maxIter; ++iiter) {
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
      srt = procrustes(current, target);
      clouds[i] = srtWarp(clouds[i], srt);
      float meanSquaredError = eigenMSE(target, current);

      std::cout << " Loss (MSE) (cloud" << i << "): "  << meanSquaredError;
      std::cout << " (using " << target.rows() << " mutual neighbors)" << '\n';

      float newModel[16] = {srt.R(0,0), srt.R(1,0), srt.R(2,0), 0,
                            srt.R(0,1), srt.R(1,1), srt.R(2,1), 0,
                            srt.R(0,2), srt.R(1,2), srt.R(2,2), 0,
                            srt.t(0), srt.t(1), srt.t(2), 1};

      g_display_mutex.lock();
      model[i] *= glm::make_mat4(newModel);
      g_display_mutex.unlock();
    }
  }
}



int main()
{
    int nClouds = 3;
    Eigen::MatrixXf* clouds;

    float s01;
    float s02;
    Eigen::Matrix3f R01;
    Eigen::Matrix3f R02;
    Eigen::Vector3f t01;
    Eigen::Vector3f t02;

    R01 = rotMat(10,0,0);
    s01 = 1;
    t01 << 2, 0, 0;
    srtTransformation transf01 = {s01, R01, t01};

    R02 = rotMat(0,-12,0);
    s02 = 1;
    t02 << 0, 1, 0;
    srtTransformation transf02 = {s02, R02, t02};

    clouds = new Eigen::MatrixXf [nClouds];
    std::vector<float> cloud0_vec;
    cloud0_vec = readObj("male_head.obj");

    clouds[0] = vectorToEigen(cloud0_vec);
    clouds[1] = srtWarp(clouds[0], transf01);
    clouds[2] = srtWarp(clouds[0], transf02);
    clouds[0]; ///= 20;
    clouds[1]; ///= 20;
    clouds[2];// /= 20;
    int nVertices = clouds[0].rows();

    float *vertices[3];
    std::vector<float> vecs[3];

    for (int i = 0; i < nClouds; ++i) {
      vecs[i] = eigenToVector(clouds[i]);
      vertices[i] = vecs[i].data();
      std::cout << eigenToVector(clouds[i]).data() << std::endl;
    }



    std::future<void> drawResult(std::async(drawThread, vertices, nVertices, nClouds, sizeof(float)*nVertices*3));
    std::future<void> updateModelsMatrixResult(std::async(updateModelsMatrix, clouds, nClouds));

    drawResult.get();
    //updateModelsMatrixResult.get();

    return 0;

}
