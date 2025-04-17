#include "camera.hpp"
#include "skeletal.hpp"
#include <iostream>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;
CameraControl camCtl;

// Pointers to your bones
Bone *root = nullptr;
std::vector<Bone *> allBones;

int main() {
  int width = 640, height = 480;
  if (!r.initialize("Skeleton Animation", width, height))
    return EXIT_FAILURE;

  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(glm::vec3(1, 1, 2), glm::vec3(0, 0, 0),
                              glm::vec3(0, 1, 0));
  program = r.createShaderProgram(r.vsBlinnPhong(), r.fsBlinnPhong());

  // ✨ BUILD BONE HIERARCHY: 1 → 2 → 2 each
  root = createBone(r, nullptr, glm::vec3(0.0f), glm::vec3(0, 0, 1));
  allBones.push_back(root);

  Bone *left =
      createBone(r, root, glm::vec3(-0.3f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  Bone *right =
      createBone(r, root, glm::vec3(0.3f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  allBones.push_back(left);
  allBones.push_back(right);

  Bone *l1 =
      createBone(r, left, glm::vec3(-0.1f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  Bone *l2 =
      createBone(r, left, glm::vec3(0.1f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  Bone *r1 =
      createBone(r, right, glm::vec3(-0.1f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  Bone *r2 =
      createBone(r, right, glm::vec3(0.1f, 1.0f, 0.0f), glm::vec3(0, 0, 1));
  allBones.insert(allBones.end(), {l1, l2, r1, r2});

  while (!r.shouldQuit()) {
    float t = SDL_GetTicks64() * 1e-3f;

    // ✨ Animate bones (e.g. simple sinusoidal waving)
    root->hinge_angle = 0.0f;
    left->hinge_angle = 0.5f * sin(t);
    right->hinge_angle = -0.5f * sin(t);
    l1->hinge_angle = 0.5f * cos(2 * t);
    l2->hinge_angle = -0.3f * cos(1.5f * t);
    r1->hinge_angle = 0.5f * cos(2 * t + 1.0f);
    r2->hinge_angle = -0.3f * cos(1.5f * t + 0.5f);

    updateBoneTransforms(root);
    updateAllMeshes(r, root);

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

    // ✨ draw all bones recursively
    std::function<void(Bone *)> drawRecursive = [&](Bone *bone) {
      r.setUniform(program, "model", bone->world_transform);
      r.drawTriangles(bone->object);
      r.drawEdges(bone->object);
      for (Bone *child : bone->children)
        drawRecursive(child);
    };
    drawRecursive(root);

    r.show();
  }
}
