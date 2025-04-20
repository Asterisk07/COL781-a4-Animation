#include "camera.hpp"
#include "skeletal.hpp"
#include <iostream>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;
CameraControl camCtl;

Bone *root = nullptr;
std::vector<Joint> joints;
std::vector<std::vector<float>> keyframes;
std::vector<float> keyframeTimes;

// Catmull-Rom interpolation with unequal spacing
float catmullRom(float p0, float p1, float p2, float p3, float t, float dt1,
                 float dt2, float dt3) {
  float t1 =
      ((p2 - p0) / (dt1 + dt2) - (p1 - p0) / dt1) * dt2 + (p2 - p1) / dt2;
  float t2 =
      ((p3 - p1) / (dt2 + dt3) - (p2 - p1) / dt2) * dt2 + (p2 - p1) / dt2;

  float a = 2 * (p1 - p2) + t1 + t2;
  float b = -3 * (p1 - p2) - 2 * t1 - t2;
  float c = t1;
  float d = p1;

  return a * t * t * t + b * t * t + c * t + d;
}

int main() {
  int width = 640, height = 480;
  if (!r.initialize("Skeleton Animation", width, height))
    return EXIT_FAILURE;

  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(glm::vec3(1, 1, 2), glm::vec3(0, 0, 0),
                              glm::vec3(0, 1, 0));
  program = r.createShaderProgram(r.vsBlinnPhong(), r.fsBlinnPhong());

  /*---- START OF USER INPUT ----*/
  // maybe the distance in y cordainte is distance from centre ofprevious to
  // start ofnew one ofc, thats why 0.5f
  // Bone *addJoint(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
  //   glm::vec3 axis, std::vector<Joint> &jointList, float length,
  //   float width, float depth)
  float l = 0.2f, b = 0.8f, h = 0.2f;
  root = addJoint(r, nullptr, glm::vec3(0.4f, -0.0f, -1.0f), glm::vec3(1, 0, 0),
                  joints, 1.0f, 1.0f, 1.0f); // shoulder (fixed)
  Bone *forearm = addJoint(r, root, glm::vec3(0.0f, 0.5f, 0.0f),
                           glm::vec3(0, 0, 1), joints, l, b, h); // elbow
  Bone *hand = addJoint(r, forearm, glm::vec3(0.0f, 0.5f, 0.0f),
                        glm::vec3(0, 0, 1), joints, l, b, h); // elbow

  keyframeTimes = {0.0f, 1.0f, 2.0f, 3.0f, 4.0f};

  keyframes = {
      {0.0f, 0.0f, 0.0f},  // frame 1: straight arm
      {0.0f, 2.0f, -0.0f}, // frame 2: elbow bent
      {0.0f, 4.0f, 0.0f},  // frame 3: elbow bent other way
      {0.0f, 1.5f, 0.0f},  // frame 4: back to straight
      {0.0f, 2.0f, 0.0f}   // frame 4: back to straight
  };
  /*---- END OF USER INPUT ----*/

  while (!r.shouldQuit()) {
    float t = SDL_GetTicks64() * 1e-3f;
    float T = fmod(t, keyframeTimes.back());

    // 🔍 Find keyframe interval [k1, k2]
    int k = 0;
    while (k < keyframeTimes.size() - 1 && T >= keyframeTimes[k + 1])
      k++;

    int k0 = std::max(0, k - 1);
    int k1 = k;
    int k2 = std::min((int)keyframeTimes.size() - 1, k + 1);
    int k3 = std::min((int)keyframeTimes.size() - 1, k + 2);

    float t0 = keyframeTimes[k0];
    float t1 = keyframeTimes[k1];
    float t2 = keyframeTimes[k2];
    float t3 = keyframeTimes[k3];

    float u = (T - t1) / (t2 - t1); // normalize t between k1 and k2

    for (int i = 0; i < joints.size(); ++i) {
      float a0 = keyframes[k0][i];
      float a1 = keyframes[k1][i];
      float a2 = keyframes[k2][i];
      float a3 = keyframes[k3][i];
      float angle = catmullRom(a0, a1, a2, a3, u, t1 - t0, t2 - t1, t3 - t2);
      setTheta(joints[i], angle);
    }

    updateJointHierarchy(root);
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

    // 🦴 Draw recursively
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
