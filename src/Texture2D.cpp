#include "Texture2D.h"
#include <algorithm>
#include <cmath>

#include "stb_image.h"

static inline float srgbToLinear(float c) {
  return std::pow(std::max(c, 0.0f), 2.2f);
}

float Texture2D::wrap01(float x) {
  x = x - std::floor(x);
  return x;
}

float Texture2D::clamp01(float x) {
  return std::max(0.0f, std::min(1.0f, x));
}

glm::vec3 Texture2D::texel(int x, int y) const {
  if(_w == 0 || _h == 0) return glm::vec3(1,0,1);
  x = (x % _w + _w) % _w;
  y = (y % _h + _h) % _h;
  int idx = (y*_w + x)*3;
  return glm::vec3(_rgb[idx+0], _rgb[idx+1], _rgb[idx+2]);
}

bool Texture2D::load(const std::string& filename, bool flipVertically) {
  stbi_set_flip_vertically_on_load(flipVertically ? 1 : 0);

  unsigned char* data = stbi_load(filename.c_str(), &_w, &_h, &_comp, 3);
  if(!data) {
    _w = _h = _comp = 0;
    _rgb.clear();
    return false;
  }

  _rgb.resize(_w*_h*3);
  for(int i=0; i<_w*_h*3; ++i) {
    float c = data[i] / 255.0f;
    _rgb[i] = srgbToLinear(c);
  }

  stbi_image_free(data);
  return true;
}

glm::vec3 Texture2D::sample(const glm::vec2& uvIn) const {
  if(_w==0 || _h==0) return glm::vec3(0.2f,0.2f,0.2f);

  float u = wrap01(uvIn.x);
  float v = wrap01(uvIn.y);

  float fx = u * (_w - 1);
  float fy = v * (_h - 1);

  int x0 = (int)std::floor(fx);
  int y0 = (int)std::floor(fy);
  int x1 = x0 + 1;
  int y1 = y0 + 1;

  float tx = fx - x0;
  float ty = fy - y0;

  glm::vec3 c00 = texel(x0,y0), c10 = texel(x1,y0);
  glm::vec3 c01 = texel(x0,y1), c11 = texel(x1,y1);

  glm::vec3 c0 = (1-tx)*c00 + tx*c10;
  glm::vec3 c1 = (1-tx)*c01 + tx*c11;
  return (1-ty)*c0 + ty*c1;
}
