#include "skeletal.hpp"

namespace Skeleton {

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

// Character state
std::vector<Bone> bones;
Bone *root;
glm::vec3 root_pos;
// glm::quat root_rot; // tk include later

// Animation keyframes
struct Keyframe {
  float time;
  std::vector<float> hinge_angles; // one per bone
};

std::vector<Keyframe> keyframes;

void initialize(GL::Rasterizer &r) {
  // Create skeleton (torso → upper arm → forearm, etc.)
  // Set mesh vertices per bone
  // Create rasterizer buffers for each bone's mesh

  // Store buffer handles into each Bone object
}

void update(float t) {
  // Interpolate hinge angles from keyframes using Catmull-Rom
  // Update each bone's transform
  // Update each bone’s mesh vertex buffer

  for (Bone &bone : bones) {
    // bone.world_transform = ... (from hierarchy)
    // transform mesh vertices: world_pos = bone.world_transform *
    // mesh_vertex_local write to vertexBuf
  }
}

void draw(GL::Rasterizer &r, GL::ShaderProgram &program) {
  for (Bone &bone : bones) {
    r.setUniform(program, "model", bone.world_transform);
    r.drawTriangles(bone.object);
    r.drawEdges(bone.object); // if desired
  }
}

} // namespace Skeleton
