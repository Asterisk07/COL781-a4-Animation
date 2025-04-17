#include "skeletal.hpp"

// ...

#include "camera.hpp"

#include <iostream>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;
CameraControl camCtl;

int main() {
  int width = 640, height = 480;
  if (!r.initialize("Skeleton Animation", width, height))
    return EXIT_FAILURE;

  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(glm::vec3(1, 1, 2), glm::vec3(0, 0, 0),
                              glm::vec3(0, 1, 0));
  program = r.createShaderProgram(r.vsBlinnPhong(), r.fsBlinnPhong());

  Skeleton::initialize(r);

  while (!r.shouldQuit()) {
    float t = SDL_GetTicks64() * 1e-3;
    Skeleton::update(t);

    camCtl.update();
    auto &camera = camCtl.camera;

    r.clear(glm::vec4(0.4, 0.4, 0.4, 1.0));
    r.enableDepthTest();
    r.useShaderProgram(program);

    r.setUniform(program, "view", camera.getViewMatrix());
    r.setUniform(program, "projection", camera.getProjectionMatrix());
    r.setUniform(program, "lightPos", camera.position);
    r.setUniform(program, "viewPos", camera.position);
    r.setUniform(program, "lightColor", glm::vec3(1.0f));

    r.setupFilledFaces();
    r.setUniform(program, "ambientColor", 0.2f * glm::vec3(1));
    r.setUniform(program, "extdiffuseColor", 0.9f * glm::vec3(1, 0.7, 0.2));
    r.setUniform(program, "intdiffuseColor", 0.4f * glm::vec3(1, 0.7, 0.2));
    r.setUniform(program, "specularColor", 0.6f * glm::vec3(1));
    r.setUniform(program, "phongExponent", 20.0f);

    Skeleton::draw(r, program);

    r.show();
  }
}
