#include "skeletal.hpp"

// namespace Skeleton {

// Creates a simple flat rectangular mesh centered at origin
void createBoxMesh(std::vector<glm::vec3> &vertices,
                   std::vector<glm::vec3> &normals,
                   std::vector<glm::ivec3> &triangles) {
  // Define a flat rectangle (quad split into 2 triangles)
  vertices = {
      glm::vec3(-0.05f, 0.0f, 0.0f), // bottom left
      glm::vec3(0.05f, 0.0f, 0.0f),  // bottom right
      glm::vec3(0.05f, 1.0f, 0.0f),  // top right
      glm::vec3(-0.05f, 1.0f, 0.0f)  // top left
  };

  triangles = {glm::ivec3(0, 1, 2), glm::ivec3(0, 2, 3)};

  glm::vec3 normal = glm::vec3(0, 0, 1);
  normals = {normal, normal, normal, normal};
}
// 2. Builder Function: createBone
// Inside skeletal.cpp
Bone *createBone(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
                 glm::vec3 axis) {
  Bone *b = new Bone();
  b->parent = parent;
  b->joint_pos = joint_pos;
  b->axis = axis;
  if (parent)
    parent->children.push_back(b);

  std::vector<glm::vec3> verts, norms;
  std::vector<glm::ivec3> tris;
  createBoxMesh(verts, norms, tris); // mesh created

  b->mesh_vertices_local = verts;
  b->mesh_normals_local = norms;

  b->object = r.createObject();
  b->vertexBuf =
      r.createVertexAttribs(b->object, 0, verts.size(), verts.data());
  b->normalBuf =
      r.createVertexAttribs(b->object, 1, norms.size(), norms.data());
  r.createTriangleIndices(b->object, tris.size(), tris.data());

  std::vector<glm::ivec2> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  r.createEdgeIndices(b->object, edges.size(), edges.data());

  return b;
}

// Bone *createBone(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
//                  glm::vec3 axis) {
//   Bone *bone = new Bone();
//   bone->parent = parent;
//   bone->joint_pos = joint_pos;
//   bone->axis = axis;

//   if (parent)
//     parent->children.push_back(bone);

//   std::vector<glm::vec3> verts, norms;
//   std::vector<glm::ivec3> tris;
//   createBoxMesh(verts, norms, tris);

//   bone->mesh_vertices_local = verts;
//   bone->mesh_normals_local = norms;

//   bone->object = r.createObject();
//   bone->vertexBuf =
//       r.createVertexAttribs(bone->object, 0, verts.size(), verts.data());
//   bone->normalBuf =
//       r.createVertexAttribs(bone->object, 1, norms.size(), norms.data());
//   r.createTriangleIndices(bone->object, tris.size(), tris.data());

//   std::vector<glm::ivec2> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
//   r.createEdgeIndices(bone->object, edges.size(), edges.data());

//   return bone;
// }

// void initialize(GL::Rasterizer &r) {
//   // Create skeleton (torso → upper arm → forearm, etc.)
//   // Set mesh vertices per bone
//   // Create rasterizer buffers for each bone's mesh

//   // Store buffer handles into each Bone object
//   for (Bone &bone : bones) {
//     std::vector<glm::vec3> local_vertices, local_normals;
//     std::vector<glm::ivec3> triangles;

//     createBoxMesh(local_vertices, local_normals, triangles);

//     bone.mesh_vertices_local = local_vertices;
//     bone.mesh_normals_local = local_normals;

//     bone.object = r.createObject();
//     bone.vertexBuf = r.createVertexAttribs(
//         bone.object, 0, local_vertices.size(), local_vertices.data());
//     bone.normalBuf = r.createVertexAttribs(bone.object, 1,
//     local_normals.size(),
//                                            local_normals.data());
//     r.createTriangleIndices(bone.object, triangles.size(), triangles.data());

//     // Optional: wireframe
//     std::vector<glm::ivec2> edges = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
//     r.createEdgeIndices(bone.object, edges.size(), edges.data());
//   }
// }

// void update(float t) {
//   // Interpolate hinge angles from keyframes using Catmull-Rom
//   // Update each bone's transform
//   // Update each bone’s mesh vertex buffer

//   for (Bone &bone : bones) {
//     // bone.world_transform = ... (from hierarchy)
//     // transform mesh vertices: world_pos = bone.world_transform *
//     // mesh_vertex_local write to vertexBuf
//   }
// }

// void draw(GL::Rasterizer &r, GL::ShaderProgram &program) {
//   for (Bone &bone : bones) {
//     r.setUniform(program, "model", bone.world_transform);
//     r.drawTriangles(bone.object);
//     r.drawEdges(bone.object); // if desired
//   }
// }

// 3. Update Bone
void updateBoneTransforms(Bone *bone, const glm::mat4 &parent_transform) {
  glm::mat4 T = glm::translate(glm::mat4(1.0f), bone->joint_pos);
  glm::mat4 R = glm::rotate(glm::mat4(1.0f), bone->hinge_angle, bone->axis);
  bone->local_transform = T * R;
  bone->world_transform = parent_transform * bone->local_transform;

  for (Bone *child : bone->children)
    updateBoneTransforms(child, bone->world_transform);
}

// 4 . Update mesh
void updateBoneMesh(GL::Rasterizer &r, Bone *bone) {
  std::vector<glm::vec3> world_verts;
  for (auto &v : bone->mesh_vertices_local) {
    glm::vec4 pos = bone->world_transform * glm::vec4(v, 1.0f);
    world_verts.push_back(glm::vec3(pos));
  }
  r.updateVertexAttribs(bone->vertexBuf, world_verts.size(),
                        world_verts.data());
}

void updateAllMeshes(GL::Rasterizer &r, Bone *root) {
  updateBoneMesh(r, root);
  for (Bone *child : root->children)
    updateAllMeshes(r, child);
}

// } // namespace Skeleton
Bone *addJoint(GL::Rasterizer &r, Bone *parent, glm::vec3 joint_pos,
               glm::vec3 axis, std::vector<Joint> &jointList) {
  Bone *b = createBone(r, parent, joint_pos, axis);
  jointList.push_back({b, 0.0f});
  return b;
}

void setTheta(Joint &joint, float theta) {
  joint.theta = theta;
  joint.bone->hinge_angle = theta;
}

void updateJointHierarchy(Bone *root) { updateBoneTransforms(root); }
