#pragma once

#include "camera.hpp"
#include "hw.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace Animation {

// Forward declarations
class Skeleton;
class Bone;
class BoxMesh;
class Keyframe;
class AnimationTimeline;

// Represents a single bone in the skeletal hierarchy
class Bone : public std::enable_shared_from_this<Bone> {
private:
  std::string m_name;
  glm::vec3 m_localOffset;  // Offset relative to parent
  glm::vec3 m_rotationAxis; // Axis around which the bone rotates
  float m_currentRotation;  // Current rotation angle in radians

  std::shared_ptr<Bone> m_parent;
  std::vector<std::shared_ptr<Bone>> m_children;

  glm::mat4 m_localTransform; // Local transform matrix
  glm::mat4 m_worldTransform; // World transform matrix

  std::shared_ptr<BoxMesh> m_attachedMesh;

public:
  Bone(const std::string &name, const glm::vec3 &offset,
       const glm::vec3 &rotAxis);

  void addChild(std::shared_ptr<Bone> child);
  void attachMesh(std::shared_ptr<BoxMesh> mesh);

  void setRotation(float angle);
  float getRotation() const;

  const glm::vec3 &getLocalOffset() const;
  const glm::vec3 &getRotationAxis() const;

  void updateTransforms(const glm::mat4 &parentTransform);
  const glm::mat4 &getWorldTransform() const;

  std::shared_ptr<Bone> getParent() const;
  const std::vector<std::shared_ptr<Bone>> &getChildren() const;

  const std::string &getName() const;
  std::shared_ptr<BoxMesh> getAttachedMesh() const;
};

// A simple box mesh attached to a bone
class BoxMesh {
private:
  glm::vec3 m_dimensions;
  glm::vec3 m_color;

  std::vector<glm::vec3> m_vertices;
  std::vector<glm::vec3> m_normals;
  std::vector<glm::ivec3> m_triangles;
  std::vector<glm::ivec2> m_edges;

public:
  BoxMesh(const glm::vec3 &dimensions, const glm::vec3 &color);

  void generateGeometry();

  const std::vector<glm::vec3> &getVertices() const;
  const std::vector<glm::vec3> &getNormals() const;
  const std::vector<glm::ivec3> &getTriangles() const;
  const std::vector<glm::ivec2> &getEdges() const;

  const glm::vec3 &getColor() const;
};

// Represents a single frame in an animation
class Keyframe {
private:
  float m_timestamp;
  glm::vec3 m_rootPosition;
  std::map<std::string, float> m_boneRotations;

public:
  Keyframe(float timestamp, const glm::vec3 &rootPos);

  void setBoneRotation(const std::string &boneName, float rotation);
  float getBoneRotation(const std::string &boneName) const;

  float getTimestamp() const;
  const glm::vec3 &getRootPosition() const;

  const std::map<std::string, float> &getAllRotations() const;
};

// Complete skeleton made of multiple bones
class Skeleton {
private:
  std::shared_ptr<Bone> m_rootBone;
  std::map<std::string, std::shared_ptr<Bone>> m_boneMap;
  glm::vec3 m_position;

public:
  Skeleton();

  void setRootBone(std::shared_ptr<Bone> root);
  std::shared_ptr<Bone> getRootBone() const;

  void addBone(std::shared_ptr<Bone> bone);
  std::shared_ptr<Bone> getBoneByName(const std::string &name) const;

  void setPosition(const glm::vec3 &position);
  const glm::vec3 &getPosition() const;

  void updateTransforms();
};

// Animation timeline managing keyframes and interpolation
class AnimationTimeline {
private:
  std::vector<std::shared_ptr<Keyframe>> m_keyframes;
  float m_duration;
  bool m_looping;

  // Helper methods for Catmull-Rom interpolation
  glm::vec3 interpolatePosition(float t, int idx1, int idx2) const;
  float interpolateRotation(float t, int idx1, int idx2,
                            const std::string &boneName) const;
  float estimateDerivative(float p0, float p1, float p2) const;
  glm::vec3 estimateDerivative(const glm::vec3 &p0, const glm::vec3 &p1,
                               const glm::vec3 &p2) const;

public:
  AnimationTimeline();

  void addKeyframe(std::shared_ptr<Keyframe> keyframe);
  void setLooping(bool loop);

  std::shared_ptr<Keyframe> getKeyframe(int index) const;
  int getKeyframeCount() const;

  float getDuration() const;
  bool isLooping() const;

  // Apply interpolated animation to skeleton at time t
  void applyAnimation(std::shared_ptr<Skeleton> skeleton, float t) const;

  // Find indices and interpolation factor for time t
  void getKeyframeIndices(float t, int &idx1, int &idx2, float &factor) const;
};

// Animator class to manage skeletons and animations
class Animator {
private:
  std::shared_ptr<Skeleton> m_skeleton;
  std::shared_ptr<AnimationTimeline> m_animation;
  float m_currentTime;

public:
  Animator();

  void setSkeleton(std::shared_ptr<Skeleton> skeleton);
  std::shared_ptr<Skeleton> getSkeleton() const;

  void setAnimation(std::shared_ptr<AnimationTimeline> animation);
  std::shared_ptr<AnimationTimeline> getAnimation() const;

  void update(float deltaTime);
  float getCurrentTime() const;
};

} // namespace Animation