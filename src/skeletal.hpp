
#ifndef SKELETAL_HPP // "If SKELETAL_HPP has NOT been defined yet..."
#define SKELETAL_HPP // "Then define SKELETAL_HPP"

#include "camera.hpp"
#include "hw.hpp"

#include <functional>
#include <glm/gtc/matrix_transform.hpp>

#include <iostream>
#include <vector>

namespace GL = COL781::OpenGL;

// namespace Skeleton {

struct Bone {
  Bone *parent;
  std::vector<Bone *> children;

  glm::vec3 joint_pos;
  glm::vec3 axis;
  float hinge_angle;

  glm::mat4 local_transform;
  glm::mat4 world_transform;

  std::vector<glm::vec3> mesh_vertices_local;
  std::vector<glm::vec3> mesh_normals_local;

  GL::Object object;
  GL::AttribBuf vertexBuf, normalBuf;
};

struct Joint {
  Bone *bone;
  float theta; // editable angle
};

Bone *addJoint(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
               glm::vec3 axis, std::vector<Joint> &jointList);

// Optional helper
void setTheta(Joint &joint, float theta);
void updateJointHierarchy(Bone *root);
// Character state
// std::vector<Bone> bones;
// Bone *root;
// glm::vec3 root_pos;
// glm::quat root_rot; // tk include later

// Animation keyframes
struct Keyframe {
  float time;
  std::vector<float> hinge_angles; // one per bone
};

// std::vector<Keyframe> keyframes;

// // Called once at startup
// void initialize(GL::Rasterizer &r);

// // Called every frame with time `t` (in seconds)
// void update(float t);

// // Called every frame to draw
// void draw(GL::Rasterizer &r, GL::ShaderProgram &program);
void createBoxMesh(std::vector<glm::vec3> &vertices,
                   std::vector<glm::vec3> &normals,
                   std::vector<glm::ivec3> &triangles);

// 2. Builder Function: createBone
Bone *createBone(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
                 glm::vec3 axis);

void updateBoneTransforms(Bone *bone,
                          const glm::mat4 &parent_transform = glm::mat4(1.0f));
void updateBoneMesh(GL::Rasterizer &r, Bone *bone);
void updateAllMeshes(GL::Rasterizer &r, Bone *root);

// } // namespace Skeleton

#endif
