#pragma once
#include <glm/glm.hpp>
#include <string>
#include <vector>

class EnvMap {
public:
  bool loadHDR(const std::string& filename);
  glm::vec3 sample(const glm::vec3& dir) const;

private:
  int _w=0, _h=0;
  std::vector<float> _rgb;

  glm::vec3 texel(int x, int y) const;
};
