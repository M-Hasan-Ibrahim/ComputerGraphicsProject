#pragma once
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>

struct FrogSelectAnim
{
  // state
  bool  active   = false;
  bool  selected = false;
  bool  target   = false;

  float time = 0.0f;
  float dur  = 0.9f;    // seconds

  // knobs
  float scaleSelect = 0.05f;
  float spinTurns   = 0.9f;    // 1 turn = 360 deg
  float forwardExtra = 0.0f;   // move closer to camera
  glm::vec3 centerOffset = glm::vec3(0.0f, -1.2f, -2.5f); // (right, up, forward) offsets in camera basis

  // saved base pose
  glm::mat4 baseMat = glm::mat4(1.0f);
  float baseScale = 1.0f;

  // curve
  glm::vec3 fromPos{0}, toPos{0}, ctrlPos{0};

  // call once after frog is placed
  void initFromFrogMat(const glm::mat4& frogMat);

  // call on P press
  void toggleStart(const glm::mat4& currentFrogMat,
                   const glm::vec3& camPos,
                   const glm::vec3& right,
                   const glm::vec3& up,
                   const glm::vec3& forward);

  // call every frame
  void update(float dt, glm::mat4& frogMat);
};
