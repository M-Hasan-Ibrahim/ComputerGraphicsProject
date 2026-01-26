#include "FrogSelectAnim.h"
#include <algorithm>
#include <cmath>

static inline float clamp01(float x) { return std::max(0.f, std::min(1.f, x)); }
static inline float smoothstep01(float t) { t = clamp01(t); return t*t*(3.f - 2.f*t); }

static inline glm::vec3 bezier2(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2, float t) {
  float a = 1.f - t;
  return a*a*p0 + 2.f*a*t*p1 + t*t*p2;
}

void FrogSelectAnim::initFromFrogMat(const glm::mat4& frogMat)
{
  baseMat = frogMat;
  baseScale = glm::length(glm::vec3(baseMat[0])); // uniform scale assumption
}

void FrogSelectAnim::toggleStart(const glm::mat4& currentFrogMat,
                                 const glm::vec3& camPos,
                                 const glm::vec3& right,
                                 const glm::vec3& up,
                                 const glm::vec3& forward)
{
  if(active) return;

  active = true;
  time = 0.0f;

  fromPos = glm::vec3(currentFrogMat[3]);

  bool goingToCenter = !selected;
  target = goingToCenter;

  if(goingToCenter) {
    float dist = glm::dot((fromPos - camPos), forward);

    glm::vec3 centerPos =
      camPos + forward * (dist - forwardExtra)
            + right   * centerOffset.x
            + up      * centerOffset.y
            + forward * centerOffset.z;

    toPos = centerPos;

    // curve strength (edit these if you want)
    ctrlPos = 0.5f*(fromPos + toPos) + 0.35f*up + 0.25f*right;
  } else {
    toPos = glm::vec3(baseMat[3]);
    ctrlPos = 0.5f*(fromPos + toPos) + 0.35f*up - 0.25f*right;
  }
}

void FrogSelectAnim::update(float dt, glm::mat4& frogMat)
{
  if(!active) return;

  time += dt;
  float t = time / dur;
  float u = smoothstep01(t);

  glm::vec3 pos = bezier2(fromPos, ctrlPos, toPos, u);

  // scale
  float s0 = baseScale;
  float s1 = scaleSelect;
  float sc = target ? (s0 + (s1 - s0)*u) : (s1 + (s0 - s1)*u);

  // spin
  float angle = target
    ? (-glm::two_pi<float>() * spinTurns * u)
    : (-glm::two_pi<float>() * spinTurns * (1.0f - u));

  glm::quat qSpin = glm::angleAxis(angle, glm::vec3(0,1,0));

  // base rotation: remove translation + remove base scale
  glm::mat4 RS = baseMat;
  RS[3] = glm::vec4(0,0,0,1);

  float invHome = (baseScale > 1e-8f) ? (1.0f / baseScale) : 1.0f;
  glm::mat4 Ronly = glm::scale(glm::mat4(1.0f), glm::vec3(invHome)) * RS;

  glm::quat qBase = glm::quat_cast(glm::mat3(Ronly));
  glm::quat q = qSpin * qBase;

  frogMat =
    glm::translate(glm::mat4(1.0f), pos) *
    glm::mat4_cast(q) *
    glm::scale(glm::mat4(1.0f), glm::vec3(sc));

  if(t >= 1.0f) {
    active = false;
    selected = target;

    if(!selected) frogMat = baseMat; // snap home exactly
  }
}
