// ----------------------------------------------------------------------------
// main.cpp
//
//  Created on: Wed Oct 21 11:32:52 2020
//      Author: Kiwon Um (originally designed by Tamy Boubekeur)
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: IGR202 - Practical - Shadow (DO NOT distribute!)
//
// http://www.opengl-tutorial.org/intermediate-tutorials/tutorial-16-shadow-mapping/
//
// Copyright 2020 Kiwon Um and Tamy Boubekeur
//
// The copyright to the computer program(s) herein is the property of Kiwon Um.
// The program(s) may be used and/or copied only with the written permission of
// Kiwon Um or in accordance with the terms and conditions stipulated in the
// agreement/contract under which the program(s) have been supplied.
// ----------------------------------------------------------------------------

#define _USE_MATH_DEFINES

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/ext.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <ios>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>
#include <exception>
#include <iomanip>
#include <array>
#include <ctime>

#include "Error.h"
#include "ShaderProgram.h"
#include "Camera.h"
#include "Mesh.h"

#include "RayTracer.h"
#include "EnvMap.h"
#include "Texture2D.h"
#include "FrogSelectAnim.h"


#include <glm/gtx/quaternion.hpp>

#include "stb_image.h"

const std::string DEFAULT_MESH_FILENAME("data/frog1.obj");

GLFWwindow *g_window = nullptr;
int g_windowWidth = 1024;
int g_windowHeight = 768;

std::shared_ptr<Camera> g_cam;

float g_meshScale = 1.0; 
bool g_rotatingP = false;
bool g_panningP = false;
bool g_zoomingP = false;
double g_baseX = 0.0, g_baseY = 0.0;
glm::vec3 g_baseTrans(0.0);
glm::vec3 g_baseRot(0.0);

float g_appTimer = 0.0;
float g_appTimerLastColckTime;
bool g_appTimerStoppedP = true;

unsigned int g_availableTextureSlot = 0;

GLuint g_skyTex = 0;
GLuint g_skyVao = 0;
std::shared_ptr<ShaderProgram> g_skyShader = nullptr;

bool g_doRayTrace = false;

static bool g_showRayTrace = false;
static GLuint g_rtTex = 0;
static int g_rtW = 800, g_rtH = 600;

static std::shared_ptr<ShaderProgram> g_rtShader;
static GLuint g_rtVao = 0;

static GLuint g_waterTex = 0; 


static FrogSelectAnim g_frogSelect;


static std::shared_ptr<Mesh> makeUvSphere(float r=1.0f, int rings=16, int sectors=32) {
  auto m = std::make_shared<Mesh>();
  auto& P  = m->vertexPositions();
  auto& N  = m->vertexNormals();
  auto& UV = m->vertexTexCoords();
  auto& I  = m->triangleIndices();

  P.clear(); N.clear(); UV.clear(); I.clear();

  for (int i = 0; i <= rings; ++i) {
    float v = float(i) / float(rings);
    float theta = v * float(M_PI); 

    for (int j = 0; j <= sectors; ++j) {
      float u = float(j) / float(sectors);
      float phi = u * 2.0f * float(M_PI); 

      glm::vec3 n(
        std::sin(theta) * std::cos(phi),
        std::cos(theta),
        std::sin(theta) * std::sin(phi)
      );

      P.push_back(r * n);
      N.push_back(glm::normalize(n));
      UV.push_back(glm::vec2(u, v));
    }
  }

  auto idx = [&](int i, int j){ return i*(sectors+1) + j; };
  for (int i = 0; i < rings; ++i) {
    for (int j = 0; j < sectors; ++j) {
      int i0 = idx(i, j);
      int i1 = idx(i+1, j);
      int i2 = idx(i+1, j+1);
      int i3 = idx(i, j+1);
      I.push_back(glm::uvec3(i0, i1, i2));
      I.push_back(glm::uvec3(i0, i2, i3));
    }
  }

  m->init();
  return m;
}

struct WaterParticle {
  glm::vec3 pos = glm::vec3(0);
  glm::vec3 vel = glm::vec3(0);
};

struct WaterEmitter {
  float rand01() const { return float(std::rand()) / float(RAND_MAX); }
  float rand11() const { return 2.0f * rand01() - 1.0f; }   

  std::shared_ptr<Mesh> sphere;

  glm::vec3 mouthLocal = glm::vec3(21.0f, 31.5f, -2.0f);

  bool running = false;
  bool looping = true;
  float groundY = -3.0f;

  float spitSpeed = 3.5f;
  float liftSpeed = 1.3f;
  float rightBias = 0.6f;
  glm::vec3 gravity = glm::vec3(0, -9.81f, 0);

  float radius = 0.025f; 

  static constexpr int kBallsPerBatch = 15;
  static constexpr int kNumBatches    = 60;
  static constexpr int kTotalBalls    = kBallsPerBatch * kNumBatches;

  std::array<WaterParticle, kTotalBalls> p;

  float batchDelay = 0.015f;
  float timeSinceStart = 0.0f;
  bool  batchActive[kNumBatches] = {false};

  int idx(int batch, int i) const { return batch * kBallsPerBatch + i; }

  glm::vec3 mouthWorld(const glm::mat4& frogMat) const {
    return glm::vec3(frogMat * glm::vec4(mouthLocal, 1.0f));
  }

  glm::vec3 computeV0(const glm::mat4& frogMat) const {
    glm::vec3 right = glm::normalize(glm::vec3(frogMat[0]));
    glm::vec3 up    = glm::normalize(glm::vec3(frogMat[1]));
    glm::vec3 fwd   = glm::normalize(glm::vec3(frogMat[2]));
    return fwd * spitSpeed + up * liftSpeed + right * rightBias;
  }

  void launchBatch(int b, const glm::mat4& frogMat) {
    glm::vec3 mouth = mouthWorld(frogMat);
    glm::vec3 right = glm::normalize(glm::vec3(frogMat[0]));
    glm::vec3 up    = glm::normalize(glm::vec3(frogMat[1]));
    glm::vec3 v0    = computeV0(frogMat);

    for (int i = 0; i < kBallsPerBatch; ++i) {
      float a = (i - 4.5f) * 0.03f;
      float c = ((i % 3) - 1) * 0.02f;

      float ja = rand11() * 0.03f;
      float jc = rand11() * 0.02f;  
      float jf = rand11() * 0.02f;  

      glm::vec3 mouth = mouthWorld(frogMat);
      glm::vec3 right = glm::normalize(glm::vec3(frogMat[0]));
      glm::vec3 up    = glm::normalize(glm::vec3(frogMat[1]));
      glm::vec3 fwd   = glm::normalize(glm::vec3(frogMat[2]));

      p[idx(b, i)].pos = mouth + right * (rand11() * 0.01f) + up * (rand11() * 0.01f);
      p[idx(b, i)].vel = v0 + right * (a + ja) + up * (c + jc) + fwd * (jf);
    }
    batchActive[b] = true;
  }

  void start(const glm::mat4& frogMat) {
    running = true;
    timeSinceStart = 0.0f;
    for (int b = 0; b < kNumBatches; ++b) batchActive[b] = false;

    glm::vec3 mouth = mouthWorld(frogMat);
    for (auto& pi : p) pi.pos = mouth;

    launchBatch(0, frogMat);
  }

  void relaunchOne(WaterParticle& pi, const glm::mat4& frogMat, int iInBatch) {
    float a = (iInBatch - 4.5f) * 0.02f;
    float c = ((iInBatch % 3) - 1) * 0.01f;

    float ja = rand11() * 0.15f; 
    float jc = rand11() * 0.02f;  
    float jf = rand11() * 0.06f;   

    glm::vec3 mouth = mouthWorld(frogMat);
    glm::vec3 right = glm::normalize(glm::vec3(frogMat[0]));
    glm::vec3 up    = glm::normalize(glm::vec3(frogMat[1]));
    glm::vec3 fwd   = glm::normalize(glm::vec3(frogMat[2]));
    glm::vec3 v0 = fwd * spitSpeed + up * liftSpeed + right * rightBias;

    pi.pos = mouth + right * (rand11() * 0.01f) + up * (rand11() * 0.01f);
    pi.vel = v0 + right * (a + ja) + up * (c + jc) + fwd * (jf);
  }


  void update(float dt, const glm::mat4& frogMat) {
    if (!running) {
      glm::vec3 mouth = mouthWorld(frogMat);
      for (auto& pi : p) pi.pos = mouth;
      return;
    }

    timeSinceStart += dt;

    for (int b = 0; b < kNumBatches; ++b) {
      if (!batchActive[b] && timeSinceStart >= b * batchDelay) {
        launchBatch(b, frogMat);
      }
    }

    for (int b = 0; b < kNumBatches; ++b) {
      if (!batchActive[b]) continue;

      for (int i = 0; i < kBallsPerBatch; ++i) {
        auto& pi = p[idx(b,i)];

        pi.vel += gravity * dt;
        pi.pos += pi.vel * dt;

        if (pi.pos.y <= groundY) {
          relaunchOne(pi, frogMat, i);
        }
      }
    }
  }


  void render(const std::shared_ptr<ShaderProgram>& sh) const {
    if (!sphere) return;

    sh->set("material.useTexture", 1);
    sh->set("material.albedo", glm::vec3(1.0f));
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, g_waterTex);
    sh->set("material.albedoTex", 0);

    for (auto& pi : p) {
      glm::mat4 M = glm::translate(glm::mat4(1.0f), pi.pos)
                  * glm::scale(glm::mat4(1.0f), glm::vec3(radius));
      sh->set("modelMat", M);
      sh->set("normMat", glm::mat3(glm::inverseTranspose(M)));
      sphere->render();
    }

    glBindTexture(GL_TEXTURE_2D, 0);
    sh->set("material.useTexture", 0);
  }
};


static WaterEmitter g_water;



GLuint loadTextureFromFileToGPU(const std::string &filename)
{
  int width, height, numComponents;
  stbi_set_flip_vertically_on_load(true);

  unsigned char *data = stbi_load(
    filename.c_str(),
    &width,
    &height,
    &numComponents, 
    0);

  GLuint texID;
  glGenTextures(1, &texID);
  glBindTexture(GL_TEXTURE_2D, texID);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexImage2D(
    GL_TEXTURE_2D,
    0,
    (numComponents == 1 ? GL_RED : numComponents == 3 ? GL_RGB : GL_RGBA), 
    width,
    height,
    0,
    (numComponents == 1 ? GL_RED : numComponents == 3 ? GL_RGB : GL_RGBA), 
    GL_UNSIGNED_BYTE,
    data);

  glGenerateMipmap(GL_TEXTURE_2D);

  stbi_image_free(data);
  glBindTexture(GL_TEXTURE_2D, 0); 
  return texID;
}

GLuint loadHDRTexture2D(const std::string& filename) {
  stbi_set_flip_vertically_on_load(true);
  int w, h, n;
  float* data = stbi_loadf(filename.c_str(), &w, &h, &n, 0);
  if(!data) {
    throw std::runtime_error(std::string("Failed to load HDR: ") + filename);
  }

  GLenum format = (n == 4) ? GL_RGBA : GL_RGB;

  GLuint tex;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, w, h, 0, format, GL_FLOAT, data);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  stbi_image_free(data);
  glBindTexture(GL_TEXTURE_2D, 0);
  return tex;
}

class FboShadowMap {
public:
  GLuint getTextureId() const { return _depthMapTexture; }

  bool allocate(unsigned int width=1024, unsigned int height=768)
  {
    glGenFramebuffers(1, &_depthMapFbo);
    glBindFramebuffer(GL_FRAMEBUFFER, _depthMapFbo);

    _depthMapTextureWidth = width;
    _depthMapTextureHeight = height;

    glGenTextures(1, &_depthMapTexture);
    glBindTexture(GL_TEXTURE_2D, _depthMapTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT16, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, _depthMapTexture, 0);

    glDrawBuffer(GL_NONE);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) {
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      return true;
    } else {
      std::cout << "PROBLEM IN FBO FboShadowMap::allocate(): FBO NOT successfully created" << std::endl;
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      return false;
    }
  }

  void bindFbo()
  {
    glViewport(0, 0, _depthMapTextureWidth, _depthMapTextureHeight);
    glBindFramebuffer(GL_FRAMEBUFFER, _depthMapFbo);
    glClear(GL_DEPTH_BUFFER_BIT);

  }

  void free() { glDeleteFramebuffers(1, &_depthMapFbo); }

  void savePpmFile(std::string const &filename)
  {
    std::ofstream output_image(filename.c_str());

    int i, j, k;
    float *pixels = new float[_depthMapTextureWidth*_depthMapTextureHeight];

    glReadBuffer(GL_COLOR_ATTACHMENT0);
    glReadPixels(0, 0, _depthMapTextureWidth, _depthMapTextureHeight, GL_DEPTH_COMPONENT , GL_FLOAT, pixels);

    output_image << "P3" << std::endl;
    output_image << _depthMapTextureWidth << " " << _depthMapTextureHeight << std::endl;
    output_image << "255" << std::endl;

    k = 0;
    for(i=0; i<_depthMapTextureWidth; ++i) {
      for(j=0; j<_depthMapTextureHeight; ++j) {
        output_image <<
          static_cast<unsigned int>(255*pixels[k]) << " " <<
          static_cast<unsigned int>(255*pixels[k]) << " " <<
          static_cast<unsigned int>(255*pixels[k]) << " ";
        k = k+1;
      }
      output_image << std::endl;
    }
    delete [] pixels;
    output_image.close();
  }

private:
  GLuint _depthMapFbo;
  GLuint _depthMapTexture;
  unsigned int _depthMapTextureWidth;
  unsigned int _depthMapTextureHeight;
};


struct Light {
  FboShadowMap shadowMap;
  glm::mat4 depthMVP;
  unsigned int shadowMapTexOnGPU;

  glm::vec3 position;
  glm::vec3 color;
  float intensity;

  void setupCameraForShadowMapping(
    std::shared_ptr<ShaderProgram> shader_shadow_map_Ptr,
    const glm::vec3 scene_center,
    const float scene_radius)
  {
  }

  void allocateShadowMapFbo(unsigned int w=800, unsigned int h=600)
  {
    shadowMap.allocate(w, h);
  }
  void bindShadowMap()
  {
    shadowMap.bindFbo();
  }
};


struct Scene {
  std::vector<Light> lights;
  std::shared_ptr<Mesh> back_rock = nullptr;
  std::shared_ptr<Mesh> stage = nullptr;
  std::shared_ptr<Mesh> rock = nullptr;
  std::shared_ptr<Mesh> frog = nullptr;

  GLuint back_rockTexture = 0;
  GLuint stageTexture = 0;


  glm::mat4 backRockMat = glm::mat4(1.0);
  glm::mat4 stageMat = glm::mat4(1.0);
  glm::mat4 frogMat = glm::mat4(1.0);

  glm::mat4 rockMat1 = glm::mat4(1.0);
  glm::mat4 rockMat2 = glm::mat4(1.0);
  glm::mat4 rockMat3 = glm::mat4(1.0);
  

  glm::vec3 scene_center = glm::vec3(0);
  float scene_radius = 1.f;

  std::shared_ptr<ShaderProgram> mainShader, shadomMapShader;

  bool saveShadowMapsPpm = false;

  void render()
  {

    if (g_showRayTrace) {
      glDisable(GL_DEPTH_TEST);

      g_rtShader->use();
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, g_rtTex);
      g_rtShader->set("rtTex", 0);

      glBindVertexArray(g_rtVao);
      glDrawArrays(GL_TRIANGLES, 0, 3);
      glBindVertexArray(0);

      glBindTexture(GL_TEXTURE_2D, 0);
      glEnable(GL_DEPTH_TEST);
      return; 
    }

    glEnable(GL_CULL_FACE);

    
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, g_windowWidth, g_windowHeight);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glCullFace(GL_BACK);

    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    g_skyShader->use();
    g_skyShader->set("invView", glm::inverse(g_cam->computeViewMatrix()));
    g_skyShader->set("invProj", glm::inverse(g_cam->computeProjectionMatrix()));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, g_skyTex);
    g_skyShader->set("skyEquirect", 0);

    glBindVertexArray(g_skyVao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);

    g_skyShader->stop();

    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);



    mainShader->use();

    mainShader->set("camPos", g_cam->getPosition());
    mainShader->set("viewMat", g_cam->computeViewMatrix());
    mainShader->set("projMat", g_cam->computeProjectionMatrix());

    Light &L = lights[0];
    mainShader->set("light.position",  L.position);
    mainShader->set("light.color",     L.color);
    mainShader->set("light.intensity", L.intensity);




    mainShader->set("material.albedo", glm::vec3(0.29, 0.51, 0.82));
    mainShader->set("modelMat", backRockMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(backRockMat)));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, back_rockTexture);
    mainShader->set("material.useTexture", 1);
    mainShader->set("material.albedoTex", 0);

    glDisable(GL_CULL_FACE);
    back_rock->render();
    glEnable(GL_CULL_FACE);
    
    mainShader->set("material.useTexture", 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    

    mainShader->set("material.albedo", glm::vec3(0.6f, 0.6f, 0.6f));
    mainShader->set("modelMat", stageMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(stageMat)));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, stageTexture);
    mainShader->set("material.useTexture", 1);
    mainShader->set("material.albedoTex", 0);

    stage->render();

    mainShader->set("material.useTexture", 0);
    glBindTexture(GL_TEXTURE_2D, 0);

    auto drawRock = [&](const glm::mat4& M){
      mainShader->set("material.albedo", glm::vec3(0.6f, 0.6f, 0.6f));
      mainShader->set("modelMat", M);
      mainShader->set("normMat", glm::mat3(glm::inverseTranspose(M)));
      rock->render();
    };
    drawRock(rockMat1);
    drawRock(rockMat2);
    drawRock(rockMat3);

    mainShader->set("material.albedo", glm::vec3(0.6f, 0.6f, 0.6f));
    mainShader->set("modelMat", frogMat);
    mainShader->set("normMat", glm::mat3(glm::inverseTranspose(frogMat)));

    frog->render();

    glDisable(GL_CULL_FACE);
    g_water.render(mainShader);
    glEnable(GL_CULL_FACE);

    mainShader->stop();
  }
};

Scene g_scene;

void printHelp()
{
  std::cout <<
    "> Help:" << std::endl <<
    "    Mouse commands:" << std::endl <<
    "    * Left button: rotate camera" << std::endl <<
    "    * Middle button: zoom" << std::endl <<
    "    * Right button: pan camera" << std::endl <<
    "    Keyboard commands:" << std::endl <<
    "    * H: print this help" << std::endl <<
    "    * R: turn ray trace on" << std::endl <<
    "    * K: turn ray tracing off" << std::endl <<
    "    * P: toggle frog animation" << std::endl <<
    "    * L: apply geometry filtering on frog" << std::endl <<
    "    * U: reset geometry filtering on frog" << std::endl <<
    "    * ESC: quit the program" << std::endl;
}

void windowSizeCallback(GLFWwindow *window, int width, int height)
{
  g_windowWidth = width;
  g_windowHeight = height;
  g_cam->setAspectRatio(static_cast<float>(width)/static_cast<float>(height));
  glViewport(0, 0, (GLint)width, (GLint)height);
}


void keyCallback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
  if(action == GLFW_PRESS && key == GLFW_KEY_H) {
    printHelp();
  } else if(action == GLFW_PRESS && key == GLFW_KEY_ESCAPE) {
    glfwSetWindowShouldClose(window, true); 
  } else if(action == GLFW_PRESS && key == GLFW_KEY_R) {
    g_doRayTrace = true;
  } else if (action == GLFW_PRESS && key == GLFW_KEY_K) {
    g_showRayTrace = false;
  } else if(action == GLFW_PRESS && key == GLFW_KEY_P) {
      glm::mat4 invV = glm::inverse(g_cam->computeViewMatrix());
      glm::vec3 right   = glm::normalize(glm::vec3(invV[0]));
      glm::vec3 up      = glm::normalize(glm::vec3(invV[1]));
      glm::vec3 forward = glm::normalize(-glm::vec3(invV[2]));

      g_frogSelect.toggleStart(g_scene.frogMat, g_cam->getPosition(), right, up, forward);
  } else if(action == GLFW_PRESS && key == GLFW_KEY_L) {
      if(g_scene.frog) {
        g_scene.frog->bilateralFilterWelded(2, 2.0f, 0.6f, 1e-6f);
        g_scene.frog->updatePositionsAndNormalsOnGPU();
      }
  } else if (action == GLFW_PRESS && key == GLFW_KEY_U) {
      if (g_scene.frog) {
        g_scene.frog->restoreState();
        g_scene.frog->updatePositionsAndNormalsOnGPU();
      }
  } else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
      if (!g_water.running) {
        g_water.start(g_scene.frogMat);
      } else {
      }
    }

}

void cursorPosCallback(GLFWwindow *window, double xpos, double ypos)
{
  int width, height;
  glfwGetWindowSize(window, &width, &height);
  const float normalizer = static_cast<float>((width + height)/2);
  const float dx = static_cast<float>((g_baseX - xpos) / normalizer);
  const float dy = static_cast<float>((ypos - g_baseY) / normalizer);
  if(g_rotatingP) {
    const glm::vec3 dRot(-dy*M_PI, dx*M_PI, 0.0);
    g_cam->setRotation(g_baseRot + dRot);
  } else if(g_panningP) {
    g_cam->setPosition(g_baseTrans + g_meshScale*glm::vec3(dx, dy, 0.0));
  } else if(g_zoomingP) {
    g_cam->setPosition(g_baseTrans + g_meshScale*glm::vec3(0.0, 0.0, dy));
  }
}

void mouseButtonCallback(GLFWwindow *window, int button, int action, int mods)
{
  if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
    if(!g_rotatingP) {
      g_rotatingP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseRot = g_cam->getRotation();
    }
  } else if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
    g_rotatingP = false;
  } else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
    if(!g_panningP) {
      g_panningP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseTrans = g_cam->getPosition();
    }
  } else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE) {
    g_panningP = false;
  } else if(button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_PRESS) {
    if(!g_zoomingP) {
      g_zoomingP = true;
      glfwGetCursorPos(window, &g_baseX, &g_baseY);
      g_baseTrans = g_cam->getPosition();
    }
  } else if(button == GLFW_MOUSE_BUTTON_MIDDLE && action == GLFW_RELEASE) {
    g_zoomingP = false;
  }
}

void initGLFW()
{
  if(!glfwInit()) {
    std::cerr << "ERROR: Failed to init GLFW" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  g_window = glfwCreateWindow(g_windowWidth, g_windowHeight, "IGR202 - Practical - Shadow", nullptr, nullptr);
  if(!g_window) {
    std::cerr << "ERROR: Failed to open window" << std::endl;
    glfwTerminate();
    std::exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(g_window);

  glfwGetFramebufferSize(g_window, &g_windowWidth, &g_windowHeight);

  glfwSetWindowSizeCallback(g_window, windowSizeCallback);
  glfwSetKeyCallback(g_window, keyCallback);
  glfwSetCursorPosCallback(g_window, cursorPosCallback);
  glfwSetMouseButtonCallback(g_window, mouseButtonCallback);
}

void clear();
void exitOnCriticalError(const std::string &message)
{
  std::cerr << "> [Critical error]" << message << std::endl;
  std::cerr << "> [Clearing resources]" << std::endl;
  clear();
  std::cerr << "> [Exit]" << std::endl;
  std::exit(EXIT_FAILURE);
}

void initOpenGL()
{
  if(!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    exitOnCriticalError("[Failed to initialize OpenGL context]");

  glEnable(GL_DEBUG_OUTPUT);   
  glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS); 
  glDebugMessageCallback(debugMessageCallback, 0);

  glCullFace(GL_BACK); 
  glEnable(GL_CULL_FACE);
  glDepthFunc(GL_LESS);   
  glEnable(GL_DEPTH_TEST);     
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f); 

  try {
    g_scene.mainShader = ShaderProgram::genBasicShaderProgram("src/vertexShader.glsl", "src/fragmentShader.glsl");

    g_skyShader = ShaderProgram::genBasicShaderProgram("src/vertexShaderSky.glsl",
                                                   "src/fragmentShaderSky.glsl");
    g_skyTex = loadHDRTexture2D("data/farmland_overcast_4k.hdr");
    glGenVertexArrays(1, &g_skyVao);



    g_scene.mainShader->stop();

  } catch(std::exception &e) {
    exitOnCriticalError(std::string("[Error loading shader program]") + e.what());
  }
}

void initScene(const std::string &meshFilename)
{
  int width, height;
  glfwGetWindowSize(g_window, &width, &height);
  g_cam = std::make_shared<Camera>();
  g_cam->setAspectRatio(static_cast<float>(width)/static_cast<float>(height));
  {
    g_scene.back_rock = std::make_shared<Mesh>();
    try {
      loadOBJ("data/rock_back.obj", g_scene.back_rock);
    } catch(std::exception &e) {
      exitOnCriticalError(std::string("[Error loading back_rock mesh]") + e.what());
    }
    g_scene.back_rock->init();

    g_scene.stage = std::make_shared<Mesh>();
    try{
      loadOBJ("data/stage.obj", g_scene.stage);
    }catch(std::exception &e){
      exitOnCriticalError(std::string("[Error loading stage mesh]") + e.what());
    }
    g_scene.stage->init();

    g_scene.rock = std::make_shared<Mesh>();
    try{
      loadOBJ("data/rock.obj", g_scene.rock);
    }catch(std::exception &e){
      exitOnCriticalError(std::string("[Error loading rock mesh]") + e.what());
    }
    g_scene.rock->init();

    g_scene.frog = std::make_shared<Mesh>();
    try{
      loadOBJ("data/frog_decimated.obj", g_scene.frog);
    }catch(std::exception &e){
      exitOnCriticalError(std::string("[Error loading frog mesh]") + e.what());
    }
    g_scene.frog->saveState();
    g_scene.frog->init();

    
    glm::vec3 Stage_Position(-0.05f, -0.55f, -5.5f);
    glm::vec3 Full_Object_Scale(0.129f);

    glm::mat4 stageTranslate = glm::translate(glm::mat4(1.0f), Stage_Position);
    glm::mat4 stageRotate = glm::rotate(glm::mat4(1.0f), glm::radians(-35.0f), glm::vec3(0,1,0));
    glm::mat4 stageScale = glm::scale(glm::mat4(1.0f), Full_Object_Scale);

    g_scene.stageMat = 
        stageTranslate *
        stageRotate *
        stageScale;

    glm::vec3 Wall_Position = Stage_Position + glm::vec3(0.0f, 0.05f, 0.0f);
    glm::mat4 wallTranslate = glm::translate(glm::mat4(1.0f), Wall_Position);

    g_scene.backRockMat =
        wallTranslate *
        stageRotate *
        stageScale;

    glm::mat4 rockScale = glm::scale(glm::mat4(1.0f), glm::vec3(0.2f));
    g_scene.rockMat1 = stageTranslate * glm::translate(glm::mat4(1.0f), glm::vec3(-0.2f, -0.7f, 1.8f)) * rockScale;
    g_scene.rockMat2 = stageTranslate * glm::translate(glm::mat4(1.0f), glm::vec3(-1.2f, -0.7f, 2.0f)) * rockScale;
    g_scene.rockMat3 = stageTranslate * glm::translate(glm::mat4(1.0f), glm::vec3(-0.2f, -10.7f, 1.5f)) * rockScale;

    glm::mat4 frogRotate = glm::rotate(glm::mat4(1.0f), glm::radians(-125.0f), glm::vec3(0,1,0));

    g_scene.frogMat = stageTranslate * glm::translate(glm::mat4(1.0f), glm::vec3(-1.2f, -0.46f, 1.95f)) * frogRotate * glm::scale(glm::mat4(1.0f), glm::vec3(0.015f));
    g_frogSelect.initFromFrogMat(g_scene.frogMat);

    g_water.sphere = makeUvSphere(1.0f, 16, 32);
    g_water.groundY = -3.0f;
    g_water.update(0.0f, g_scene.frogMat);
  }

  GLuint back_rockTex = loadTextureFromFileToGPU("data/rock_back_texture.png");
  g_scene.back_rockTexture = back_rockTex;

  GLuint stageTex = loadTextureFromFileToGPU("data/wood_table_diff_2k.jpg");
  g_scene.stageTexture = stageTex;

  g_waterTex = loadTextureFromFileToGPU("data/water.png");
  std::cout << "waterTex id = " << g_waterTex << "\n";

  g_scene.lights.clear();
  g_scene.lights.push_back(Light());
  Light &L = g_scene.lights[0];

  L.position  = glm::vec3(0.0f, 1.0f, 0.5f);
  L.color     = glm::vec3(0.6f, 0.6f, 0.6f);
  L.intensity = 3.0f;


  g_scene.scene_center = glm::vec3(0.0f);
  g_scene.scene_radius = 1.0f;
  g_meshScale = g_scene.scene_radius;


  g_cam->setPosition(g_scene.scene_center + glm::vec3(0.0, 0.0, 3.0*g_meshScale));
  g_cam->setNear(g_meshScale/100.f);
  g_cam->setFar(12.0*g_meshScale);
}

static void initRaytraceDisplay(){
  g_rtShader = ShaderProgram::genBasicShaderProgram(
    "src/vertexShaderRaytrace.glsl",
    "src/fragmentShaderRaytrace.glsl"
  );

  glGenVertexArrays(1, &g_rtVao);

  glGenTextures(1, &g_rtTex);
  glBindTexture(GL_TEXTURE_2D, g_rtTex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, g_rtW, g_rtH, 0, GL_RGB, GL_FLOAT, nullptr);
  glBindTexture(GL_TEXTURE_2D, 0);
}



static void appendMeshToRTScene(RTScene& rt, const Mesh& mesh, const glm::mat4& modelMat, int matId){
  glm::mat3 normalMat = glm::transpose(glm::inverse(glm::mat3(modelMat)));

  const auto& P = mesh.vertexPositions();
  const auto& N = mesh.vertexNormals();
  const auto& T = mesh.triangleIndices();
  const auto& UV = mesh.vertexTexCoords();


  for (size_t i = 0; i < T.size(); ++i) {
    glm::uvec3 triIdx = T[i];

    auto xformP = [&](uint32_t idx){
      return glm::vec3(modelMat * glm::vec4(P[idx], 1.f));
    };
    auto xformN = [&](uint32_t idx){
      return glm::normalize(normalMat * N[idx]);
    };

    RTTriangle tri;

    if(!UV.empty()) {
      tri.uv0 = UV[triIdx[0]];
      tri.uv1 = UV[triIdx[1]];
      tri.uv2 = UV[triIdx[2]];
    } else {
      tri.uv0 = tri.uv1 = tri.uv2 = glm::vec2(0.0f);
    }

    tri.p0 = xformP(triIdx[0]);
    tri.p1 = xformP(triIdx[1]);
    tri.p2 = xformP(triIdx[2]);
    tri.n0 = xformN(triIdx[0]);
    tri.n1 = xformN(triIdx[1]);
    tri.n2 = xformN(triIdx[2]);
    tri.matId = matId;

    rt.tris.push_back(tri);
  }

}


void init(const std::string &meshFilename)
{
  initGLFW();                 
  initOpenGL();                
  initScene(meshFilename);      
}

void clear()
{
  g_cam.reset();
  g_scene.mainShader.reset();
  g_scene.shadomMapShader.reset();
  glfwDestroyWindow(g_window);
  glfwTerminate();
}

void render()
{
  g_scene.render();
}


void update(float currentTime)
{
  static float lastTime = currentTime;   
  float dt = currentTime - lastTime;  
  lastTime = currentTime;

  if (dt > 1.0f/60.0f) dt = 1.0f/60.0f;
  if (dt < 0.0f) dt = 0.0f;

  g_frogSelect.update(dt, g_scene.frogMat);

  g_water.update(dt, g_scene.frogMat);


}


void usage(const char *command)
{
  std::cerr << "Usage : " << command << " [<file.off>]" << std::endl;
  std::exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  if(argc > 2) usage(argv[0]);
  init(argc==1 ? DEFAULT_MESH_FILENAME : argv[1]);
  initRaytraceDisplay(); 

  std::srand(std::time(nullptr));

  while(!glfwWindowShouldClose(g_window)) {
    update(static_cast<float>(glfwGetTime()));
    render();

    if (g_doRayTrace) {
      g_doRayTrace = false;
      int W = g_rtW, H = g_rtH;

      RTScene rt;

      rt.mats.clear();
      rt.textures.clear();

      int texWall = (int)rt.textures.size();
      rt.textures.push_back(Texture2D());
      rt.textures.back().load("data/rock_back_texture.png", true);

      int texStage = (int)rt.textures.size();
      rt.textures.push_back(Texture2D());
      rt.textures.back().load("data/wood_table_diff_2k.jpg", true);

      int matRock = (int)rt.mats.size();
      rt.mats.push_back(RTMaterial());
      rt.mats.back().albedo = glm::vec3(0.65f);
      rt.mats.back().shadowCatcher = false;
      rt.mats.back().useTexture = false;
      rt.mats.back().texId = -1;

      int matWall = (int)rt.mats.size();
      rt.mats.push_back(RTMaterial());
      rt.mats.back().albedo = glm::vec3(1.0f);
      rt.mats.back().shadowCatcher = false;
      rt.mats.back().useTexture = true;
      rt.mats.back().texId = texWall;

      int matStage = (int)rt.mats.size();
      rt.mats.push_back(RTMaterial());
      rt.mats.back().albedo = glm::vec3(1.0f);
      rt.mats.back().shadowCatcher = false;
      rt.mats.back().useTexture = true;
      rt.mats.back().texId = texStage;

      int matFrog = (int)rt.mats.size();
      rt.mats.push_back(RTMaterial());
      rt.mats.back().albedo = glm::vec3(0.35f, 0.95f, 0.35f);
      rt.mats.back().shadowCatcher = false;
      rt.mats.back().useTexture = false;
      rt.mats.back().texId = -1;

      int matGround = (int)rt.mats.size();
      rt.mats.push_back(RTMaterial());
      rt.mats.back().albedo = glm::vec3(1.0f);
      rt.mats.back().shadowCatcher = true;
      rt.mats.back().useTexture = false;
      rt.mats.back().texId = -1;

      appendMeshToRTScene(rt, *g_scene.back_rock, g_scene.backRockMat, matWall);
      appendMeshToRTScene(rt, *g_scene.stage, g_scene.stageMat, matStage);
      appendMeshToRTScene(rt, *g_scene.rock, g_scene.rockMat1, matRock);
      appendMeshToRTScene(rt, *g_scene.rock, g_scene.rockMat2, matRock);
      appendMeshToRTScene(rt, *g_scene.frog, g_scene.frogMat, matFrog);

      RTCamera cam;
      cam.pos = g_cam->getPosition();
      cam.invView = glm::inverse(g_cam->computeViewMatrix());
      cam.fovYDegrees = g_cam->getFov();
      cam.aspect = g_cam->getAspectRatio();

      glm::mat4 invV = cam.invView;
      glm::vec3 right = glm::vec3(invV[0]);
      glm::vec3 up = glm::vec3(invV[1]);
      glm::vec3 forward = -glm::vec3(invV[2]);

      RTLight L;
      L.position  = cam.pos + (-1.5f)*right + (3.0f)*up + (3.0f)*forward;
      L.color     = glm::vec3(1.f);
      L.intensity = 15.0f;
      
      
      EnvMap env;
      env.loadHDR("data/farmland_overcast_4k.hdr");


      RayTracer tracer(W, H);
      tracer.setEnvMap(&env);
      tracer.buildBVH(rt);
      tracer.setGround(-1.925f, matGround, 0.6f);

      auto pixels = tracer.render(rt, cam, L);
      
      glBindTexture(GL_TEXTURE_2D, g_rtTex);
      glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, W, H, GL_RGB, GL_FLOAT, pixels.data());
      glBindTexture(GL_TEXTURE_2D, 0);

      g_showRayTrace = true;
    }


    glfwSwapBuffers(g_window);
    glfwPollEvents();
  }
  clear();
  std::cout << " > Quit" << std::endl;
  return EXIT_SUCCESS;
}
