#include "skeletal.hpp"
#include <algorithm>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>
#include <iostream>

namespace Animation {

// Bone Implementation
Bone::Bone(const std::string &name, const glm::vec3 &offset,
           const glm::vec3 &rotAxis)
    : m_name(name), m_localOffset(offset),
      m_rotationAxis(glm::normalize(rotAxis)), m_currentRotation(0.0f),
      m_parent(nullptr), m_localTransform(1.0f), m_worldTransform(1.0f) {}

void Bone::addChild(std::shared_ptr<Bone> child) {
  child->m_parent = shared_from_this();
  m_children.push_back(child);
}

void Bone::attachMesh(std::shared_ptr<BoxMesh> mesh) { m_attachedMesh = mesh; }

void Bone::setRotation(float angle) { m_currentRotation = angle; }

float Bone::getRotation() const { return m_currentRotation; }

const glm::vec3 &Bone::getLocalOffset() const { return m_localOffset; }

const glm::vec3 &Bone::getRotationAxis() const { return m_rotationAxis; }

void Bone::updateTransforms(const glm::mat4 &parentTransform) {
  // Create local rotation matrix
  glm::quat rotation = glm::angleAxis(m_currentRotation, m_rotationAxis);
  glm::mat4 rotationMatrix = glm::toMat4(rotation);

  // Apply local offset
  glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), m_localOffset);

  // Combine local transformations
  m_localTransform = translationMatrix * rotationMatrix;

  // Combine with parent transform to get world transform
  m_worldTransform = parentTransform * m_localTransform;

  // Recursively update all children
  for (auto &child : m_children) {
    child->updateTransforms(m_worldTransform);
  }
}

const glm::mat4 &Bone::getWorldTransform() const { return m_worldTransform; }

std::shared_ptr<Bone> Bone::getParent() const { return m_parent; }

const std::vector<std::shared_ptr<Bone>> &Bone::getChildren() const {
  return m_children;
}

const std::string &Bone::getName() const { return m_name; }

// BoxMesh Implementation
BoxMesh::BoxMesh(const glm::vec3 &dimensions, const glm::vec3 &color)
    : m_dimensions(dimensions), m_color(color) {
  generateGeometry();
}

void BoxMesh::generateGeometry() {
  // Calculate half-dimensions for convenience
  glm::vec3 half = m_dimensions * 0.5f;

  // Define the 8 vertices of the box
  m_vertices = {
      // Front face
      glm::vec3(-half.x, -half.y, half.z), // 0
      glm::vec3(half.x, -half.y, half.z),  // 1
      glm::vec3(half.x, half.y, half.z),   // 2
      glm::vec3(-half.x, half.y, half.z),  // 3

      // Back face
      glm::vec3(-half.x, -half.y, -half.z), // 4
      glm::vec3(half.x, -half.y, -half.z),  // 5
      glm::vec3(half.x, half.y, -half.z),   // 6
      glm::vec3(-half.x, half.y, -half.z)   // 7
  };

  // Define normals (simplified - just using axis directions)
  m_normals.resize(8);
  // Front/back normals
  for (int i = 0; i < 4; i++)
    m_normals[i] = glm::vec3(0, 0, 1);
  for (int i = 4; i < 8; i++)
    m_normals[i] = glm::vec3(0, 0, -1);

  // Define triangles (2 per face, 6 faces)
  m_triangles = {// Front face
                 glm::ivec3(0, 1, 2), glm::ivec3(0, 2, 3),
                 // Back face
                 glm::ivec3(5, 4, 7), glm::ivec3(5, 7, 6),
                 // Right face
                 glm::ivec3(1, 5, 6), glm::ivec3(1, 6, 2),
                 // Left face
                 glm::ivec3(4, 0, 3), glm::ivec3(4, 3, 7),
                 // Top face
                 glm::ivec3(3, 2, 6), glm::ivec3(3, 6, 7),
                 // Bottom face
                 glm::ivec3(4, 5, 1), glm::ivec3(4, 1, 0)};

  // Define edges for wireframe rendering
  m_edges = {
      // Front face
      glm::ivec2(0, 1), glm::ivec2(1, 2), glm::ivec2(2, 3), glm::ivec2(3, 0),
      // Back face
      glm::ivec2(4, 5), glm::ivec2(5, 6), glm::ivec2(6, 7), glm::ivec2(7, 4),
      // Connecting edges
      glm::ivec2(0, 4), glm::ivec2(1, 5), glm::ivec2(2, 6), glm::ivec2(3, 7)};
}

const std::vector<glm::vec3> &BoxMesh::getVertices() const {
  return m_vertices;
}

const std::vector<glm::vec3> &BoxMesh::getNormals() const { return m_normals; }

const std::vector<glm::ivec3> &BoxMesh::getTriangles() const {
  return m_triangles;
}

const std::vector<glm::ivec2> &BoxMesh::getEdges() const { return m_edges; }

const glm::vec3 &BoxMesh::getColor() const { return m_color; }

// Keyframe Implementation
Keyframe::Keyframe(float timestamp, const glm::vec3 &rootPos)
    : m_timestamp(timestamp), m_rootPosition(rootPos) {}

void Keyframe::setBoneRotation(const std::string &boneName, float rotation) {
  m_boneRotations[boneName] = rotation;
}

float Keyframe::getBoneRotation(const std::string &boneName) const {
  auto it = m_boneRotations.find(boneName);
  if (it != m_boneRotations.end()) {
    return it->second;
  }
  return 0.0f; // Default to no rotation if not specified
}

float Keyframe::getTimestamp() const { return m_timestamp; }

const glm::vec3 &Keyframe::getRootPosition() const { return m_rootPosition; }

const std::map<std::string, float> &Keyframe::getAllRotations() const {
  return m_boneRotations;
}

// Skeleton Implementation
Skeleton::Skeleton() : m_position(0.0f) {}

void Skeleton::setRootBone(std::shared_ptr<Bone> root) {
  m_rootBone = root;
  addBone(root);
}

std::shared_ptr<Bone> Skeleton::getRootBone() const { return m_rootBone; }

void Skeleton::addBone(std::shared_ptr<Bone> bone) {
  m_boneMap[bone->getName()] = bone;

  // Recursively add all children
  for (const auto &child : bone->getChildren()) {
    addBone(child);
  }
}

std::shared_ptr<Bone> Skeleton::getBoneByName(const std::string &name) const {
  auto it = m_boneMap.find(name);
  if (it != m_boneMap.end()) {
    return it->second;
  }
  return nullptr;
}

void Skeleton::setPosition(const glm::vec3 &position) { m_position = position; }

const glm::vec3 &Skeleton::getPosition() const { return m_position; }

void Skeleton::updateTransforms() {
  if (m_rootBone) {
    // Start with a translation matrix for the skeleton's position
    glm::mat4 rootTransform = glm::translate(glm::mat4(1.0f), m_position);
    m_rootBone->updateTransforms(rootTransform);
  }
}

// AnimationTimeline Implementation
AnimationTimeline::AnimationTimeline() : m_duration(0.0f), m_looping(true) {}

void AnimationTimeline::addKeyframe(std::shared_ptr<Keyframe> keyframe) {
  // Insert keyframe in sorted order by timestamp
  auto it = std::lower_bound(m_keyframes.begin(), m_keyframes.end(), keyframe,
                             [](const std::shared_ptr<Keyframe> &a,
                                const std::shared_ptr<Keyframe> &b) {
                               return a->getTimestamp() < b->getTimestamp();
                             });

  m_keyframes.insert(it, keyframe);

  // Update duration if this is now the last keyframe
  if (keyframe->getTimestamp() > m_duration) {
    m_duration = keyframe->getTimestamp();
  }
}

void AnimationTimeline::setLooping(bool loop) { m_looping = loop; }

std::shared_ptr<Keyframe> AnimationTimeline::getKeyframe(int index) const {
  if (index >= 0 && index < m_keyframes.size()) {
    return m_keyframes[index];
  }
  return nullptr;
}

int AnimationTimeline::getKeyframeCount() const { return m_keyframes.size(); }

float AnimationTimeline::getDuration() const { return m_duration; }

bool AnimationTimeline::isLooping() const { return m_looping; }

void AnimationTimeline::getKeyframeIndices(float t, int &idx1, int &idx2,
                                           float &factor) const {
  // Handle looping time
  if (m_looping && t > m_duration) {
    t = std::fmod(t, m_duration);
  }

  // Clamp time to animation duration
  t = std::max(0.0f, std::min(t, m_duration));

  // Find keyframes surrounding time t
  for (int i = 0; i < m_keyframes.size() - 1; i++) {
    if (t <= m_keyframes[i + 1]->getTimestamp()) {
      idx1 = i;
      idx2 = i + 1;

      // Calculate interpolation factor
      float t1 = m_keyframes[idx1]->getTimestamp();
      float t2 = m_keyframes[idx2]->getTimestamp();
      factor = (t - t1) / (t2 - t1);

      return;
    }
  }

  // If we get here, t is at or beyond the last keyframe
  idx1 = m_keyframes.size() - 1;
  idx2 = m_looping ? 0 : idx1;
  factor = 0.0f; // No interpolation needed at the end
}

float AnimationTimeline::estimateDerivative(float p0, float p1,
                                            float p2) const {
  // Simple finite difference formula for derivative estimation
  return 0.5f * (p2 - p0);
}

glm::vec3 AnimationTimeline::estimateDerivative(const glm::vec3 &p0,
                                                const glm::vec3 &p1,
                                                const glm::vec3 &p2) const {
  // Component-wise derivative estimation
  return 0.5f * (p2 - p0);
}

glm::vec3 AnimationTimeline::interpolatePosition(float t, int idx1,
                                                 int idx2) const {
  // Get the 4 control points needed for Catmull-Rom
  int idx0 =
      (idx1 > 0) ? (idx1 - 1) : (m_looping ? (m_keyframes.size() - 1) : idx1);
  int idx3 =
      (idx2 < m_keyframes.size() - 1) ? (idx2 + 1) : (m_looping ? 0 : idx2);

  const glm::vec3 &p0 = m_keyframes[idx0]->getRootPosition();
  const glm::vec3 &p1 = m_keyframes[idx1]->getRootPosition();
  const glm::vec3 &p2 = m_keyframes[idx2]->getRootPosition();
  const glm::vec3 &p3 = m_keyframes[idx3]->getRootPosition();

  // Estimate tangents
  glm::vec3 m1 = estimateDerivative(p0, p1, p2);
  glm::vec3 m2 = estimateDerivative(p1, p2, p3);

  // Catmull-Rom interpolation using Hermite basis
  float u = t;
  float u2 = u * u;
  float u3 = u2 * u;

  float c1 = 2.0f * u3 - 3.0f * u2 + 1.0f;
  float c2 = u3 - 2.0f * u2 + u;
  float c3 = -2.0f * u3 + 3.0f * u2;
  float c4 = u3 - u2;

  return c1 * p1 + c2 * m1 + c3 * p2 + c4 * m2;
}

float AnimationTimeline::interpolateRotation(
    float t, int idx1, int idx2, const std::string &boneName) const {
  // Get the 4 control points needed for Catmull-Rom
  int idx0 =
      (idx1 > 0) ? (idx1 - 1) : (m_looping ? (m_keyframes.size() - 1) : idx1);
  int idx3 =
      (idx2 < m_keyframes.size() - 1) ? (idx2 + 1) : (m_looping ? 0 : idx2);

  float r0 = m_keyframes[idx0]->getBoneRotation(boneName);
  float r1 = m_keyframes[idx1]->getBoneRotation(boneName);
  float r2 = m_keyframes[idx2]->getBoneRotation(boneName);
  float r3 = m_keyframes[idx3]->getBoneRotation(boneName);

  // Estimate tangents
  float m1 = estimateDerivative(r0, r1, r2);
  float m2 = estimateDerivative(r1, r2, r3);

  // Catmull-Rom interpolation using Hermite basis
  float u = t;
  float u2 = u * u;
  float u3 = u2 * u;

  float c1 = 2.0f * u3 - 3.0f * u2 + 1.0f;
  float c2 = u3 - 2.0f * u2 + u;
  float c3 = -2.0f * u3 + 3.0f * u2;
  float c4 = u3 - u2;

  return c1 * r1 + c2 * m1 + c3 * r2 + c4 * m2;
}

void AnimationTimeline::applyAnimation(std::shared_ptr<Skeleton> skeleton,
                                       float t) const {
  if (m_keyframes.empty())
    return;

  // Single keyframe case
  if (m_keyframes.size() == 1) {
    skeleton->setPosition(m_keyframes[0]->getRootPosition());

    // Apply rotations to all bones
    const auto &rotations = m_keyframes[0]->getAllRotations();
    for (const auto &[boneName, rotation] : rotations) {
      auto bone = skeleton->getBoneByName(boneName);
      if (bone) {
        bone->setRotation(rotation);
      }
    }

    return;
  }

  // Find keyframes to interpolate between
  int idx1, idx2;
  float interpolationFactor;
  getKeyframeIndices(t, idx1, idx2, interpolationFactor);

  // Interpolate root position
  glm::vec3 position = interpolatePosition(interpolationFactor, idx1, idx2);
  skeleton->setPosition(position);

  // Get all bone names from both keyframes
  std::set<std::string> boneNames;
  for (const auto &[name, _] : m_keyframes[idx1]->getAllRotations()) {
    boneNames.insert(name);
  }

  for (const auto &[name, _] : m_keyframes[idx2]->getAllRotations()) {
    boneNames.insert(name);
  }

  // Interpolate rotations for all bones
  for (const auto &boneName : boneNames) {
    auto bone = skeleton->getBoneByName(boneName);
    if (bone) {
      float rotation =
          interpolateRotation(interpolationFactor, idx1, idx2, boneName);
      bone->setRotation(rotation);
    }
  }
}

// Animator Implementation
Animator::Animator() : m_currentTime(0.0f) {}

void Animator::setSkeleton(std::shared_ptr<Skeleton> skeleton) {
  m_skeleton = skeleton;
}

std::shared_ptr<Skeleton> Animator::getSkeleton() const { return m_skeleton; }

void Animator::setAnimation(std::shared_ptr<AnimationTimeline> animation) {
  m_animation = animation;
  m_currentTime = 0.0f;
}

std::shared_ptr<AnimationTimeline> Animator::getAnimation() const {
  return m_animation;
}

void Animator::update(float deltaTime) {
  if (!m_skeleton || !m_animation)
    return;

  // Update time
  m_currentTime += deltaTime;

  // Apply animation to skeleton
  m_animation->applyAnimation(m_skeleton, m_currentTime);

  // Update skeleton transforms
  m_skeleton->updateTransforms();
}

float Animator::getCurrentTime() const { return m_currentTime; }

// Add the method implementation for getAttachedMesh()
std::shared_ptr<BoxMesh> Bone::getAttachedMesh() const {
  return m_attachedMesh;
}

} // namespace Animation
