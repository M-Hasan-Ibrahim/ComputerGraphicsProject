#version 450 core
in vec2 uv;
out vec4 colorOut;
uniform sampler2D rtTex;

void main() {
  vec3 c = texture(rtTex, vec2(uv.x, 1.0 - uv.y)).rgb;

  // optional simple tone-map / gamma (keeps it from looking blown-out)
  c = c / (c + vec3(1.0));
  c = pow(c, vec3(1.0/2.2));
  colorOut = vec4(c, 1.0);
}
