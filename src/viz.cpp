// viz.cpp
// Date: 2019-29-10
// Created by: Gabriel Moreira

#include <iostream>
#include <vector>
#include <math.h>

#include "reader.h"
#include "gpa.h"
#include "viz.h"


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}


void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
        leftButtonPressed = true;
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
        leftButtonPressed = false;
}


void leftButtonPressedRotation(double xpos, double ypos) {
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

  float xSensitivity = 0.01;
  float ySensitivity = 0.25;

  xoffset *= xSensitivity;
  yoffset *= ySensitivity;

  if (xoffset != 0 || yoffset != 0) {
    float norm = glm::length(camera.pos);
    glm::vec3 delta_pos = glm::normalize(camera.up)*yoffset + glm::cross(camera.up, camera.pos)*xoffset;
    glm::vec3 normal_vec = glm::cross(camera.pos, camera.up);
    camera.pos -= delta_pos;
    camera.pos = glm::normalize(camera.pos)*norm;
    camera.up = glm::normalize(glm::cross(normal_vec, camera.pos));
  };
};


void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
  if (leftButtonPressed == false)
    firstMouse = true;

  if (leftButtonPressed == true) {
    leftButtonPressedRotation(xpos, ypos);
  };
};


void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    float cameraSpeed = 0.25f; // adjust accordingly

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
      camera.pos += cameraSpeed * camera.front;

    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
      camera.pos -= cameraSpeed * camera.front;

    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS)
      camera.pos += camera.up;

    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS)
      camera.pos -= camera.up;

    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS)
      camera.pos -= glm::cross(camera.front, camera.up);

    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS)
      camera.pos += glm::cross(camera.front, camera.up);
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
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
      std::cout << "Failed to initialize GLAD" << std::endl;
  }

  glViewport(0, 0, SCR_WIDTH, SCR_HEIGHT);
  return window;
}


void windowLoop(GLFWwindow* window, float** vertices, int nVertices, int nClouds, int buffersize, cloudShader shader) {
  glm::vec4 cloudColor [nClouds];
  cloudColor[0] = glm::vec4(1.0, 1.0, 1.0, 1.0);
  cloudColor[1] = glm::vec4(1.0, 0.0, 0.0, 1.0);
  cloudColor[2] = glm::vec4(0.0, 1.0, 0.0, 1.0);
  cloudColor[3] = glm::vec4(0.0, 0.0, 1.0, 1.0);
  cloudColor[4] = glm::vec4(1.0, 1.0, 0.0, 1.0);
  cloudColor[5] = glm::vec4(0.0, 1.0, 1.0, 1.0);
  cloudColor[6] = glm::vec4(1.0, 0.0, 1.0, 1.0);
  cloudColor[7] = glm::vec4(0.2, 0.5, 0.5, 1.0);
  cloudColor[8] = glm::vec4(0.5, 0.5, 2.0, 1.0);
  cloudColor[9] = glm::vec4(1.0, 4.0, 2.0, 1.0);

  glm::vec3 vert(0.0f, 1.0f, 0.0f);
  while(!glfwWindowShouldClose(window)) {
    // Check for user input
    processInput(window);

    // Rendering commands here
    glClearColor(0.0f, 0.f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    shader.use();

    camera.front = -1.0f * glm::normalize(camera.pos);
    view = glm::lookAt(camera.pos, camera.pos + camera.front, camera.up);
    projection = glm::perspective(glm::radians(camera.fov), ASPECT_RATIO, camera.nearField, camera.farField);

    // Pass matrices to the shader program
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);

    for (int i = 0; i < nClouds; ++i) {
      shader.setVec4("cloudColor", cloudColor[i]);

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
  int maxIter = 500;
  for (int iiter = 0; iiter < maxIter; ++iiter) {
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
      srt = procrustes(current, target);
      clouds[i] = srtWarp(clouds[i], srt);
      float meanSquaredError = eigenMSE(target, current);

      std::cout << "   Loss (MSE) (cloud" << i << "): "  << meanSquaredError;
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

    int nVertices = clouds[0].rows();

    float *vertices[NCLOUDS];
    std::vector<float> vecs[NCLOUDS];

    for (int i = 0; i < NCLOUDS; ++i) {
      vecs[i] = eigenToVector(clouds[i]);
      vertices[i] = vecs[i].data();
    }

    std::future<void> drawResult(std::async(drawThread, vertices, nVertices, NCLOUDS, sizeof(float)*nVertices*3));
    std::string startGPA;
    std::cin >> startGPA;
    std::future<void> updateModelsMatrixResult(std::async(updateModelsMatrix, clouds, NCLOUDS));

    drawResult.get();
    updateModelsMatrixResult.get();

    return 0;
}
