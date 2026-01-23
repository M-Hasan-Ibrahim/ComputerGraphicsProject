#pragma once
#include <vector>
#include <string>
#include <glm/glm.hpp>

class Texture2D {
public:
  bool load(const std::string& filename, bool flipVertically = true);
  glm::vec3 sample(const glm::vec2& uv) const;

  int w() const { return _w; }
  int h() const { return _h; }

private:
  int _w = 0, _h = 0, _comp = 0;
  std::vector<float> _rgb;

  static float wrap01(float x);
  static float clamp01(float x);
  glm::vec3 texel(int x, int y) const;
};
