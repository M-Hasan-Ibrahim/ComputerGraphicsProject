#version 450 core            // minimal GL version support expected from the GPU

struct LightSource {
  vec3 position;
  vec3 color;
  float intensity;
};


uniform LightSource light;
// TODO: shadow maps

struct Material {
  vec3 albedo;
  // TODO: textures
  int useTexture;        
  sampler2D albedoTex;
};

uniform Material material;

uniform vec3 camPos;

in vec3 fPositionModel;
in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexCoord;

out vec4 colorOut; // shader output: the color response attached to this fragment

float pi = 3.1415927;

// TODO: shadows
void main() {
  vec3 n = normalize(fNormal);

  vec3 wi = normalize(light.position - fPosition);
  vec3 Li = light.color * light.intensity;

  vec3 albedo = material.albedo;
  if (material.useTexture == 1) {
    albedo *= texture(material.albedoTex, fTexCoord).rgb;
  }

  vec3 radiance = Li * albedo * max(dot(n, wi), 0.0);
  colorOut = vec4(radiance, 1.0);
}
