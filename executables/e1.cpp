// #include "animator.hpp"
#include "camera.hpp"
#include "skeletal.hpp"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <memory>
// #include <set>
#include <vector>

using namespace COL781;
namespace GL = COL781::OpenGL;
CameraControl m_cameraControl;
using namespace Animation;

class AnimationRenderer {
private:
  GL::Rasterizer &m_rasterizer;
  GL::ShaderProgram m_shaderProgram;
  // CameraControl m_cameraControl;

  // Animator and scene elements
  std::shared_ptr<Animator> m_animator;

  std::vector<GL::AttribBuf> m_vertexBuffers;
  std::vector<GL::AttribBuf> m_normalBuffers;

  // Helper functions
  void createHumanoidSkeleton();
  void createKeyframes();
  // void setupCamera(int width, int height);
  void renderSkeleton();
  void updateMeshGeometry();

public:
  std::vector<GL::Object> m_meshObjects;
  AnimationRenderer(GL::Rasterizer &rasterizer, int width, int height);

  void initialize();
  void update(float deltaTime);
  void render();
};

AnimationRenderer::AnimationRenderer(GL::Rasterizer &rasterizer, int width,
                                     int height)
    : m_rasterizer(rasterizer) {

  // Initialize the camera controller
  m_cameraControl.initialize(width, height);
  // setupCamera(width, height);
  m_cameraControl.camera.setCameraView(
      glm::vec3(0.0f, 1.0f, 5.0f), // Position
      glm::vec3(0.0f, 1.0f, 0.0f), // Look-at target
      glm::vec3(0.0f, 1.0f, 0.0f)  // Up vector
  );

  // Create shader program
  m_shaderProgram = m_rasterizer.createShaderProgram(
      m_rasterizer.vsBlinnPhong(), m_rasterizer.fsBlinnPhong());

  // Create the animator
  m_animator = std::make_shared<Animator>();
}

void AnimationRenderer::createHumanoidSkeleton() {
  // Create skeleton
  auto skeleton = std::make_shared<Skeleton>();

  // Create bones
  // Torso (root)
  auto torso = std::make_shared<Bone>("torso", glm::vec3(0.0f, 1.0f, 0.0f),
                                      glm::vec3(0.0f, 1.0f, 0.0f));
  auto torsoBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.8f, 1.0f, 0.4f), glm::vec3(0.8f, 0.4f, 0.2f));
  torso->attachMesh(torsoBoundingBox);

  // Head
  auto head = std::make_shared<Bone>("head", glm::vec3(0.0f, 0.6f, 0.0f),
                                     glm::vec3(0.0f, 0.0f, 1.0f));
  auto headBoundingBox = std::make_shared<BoxMesh>(glm::vec3(0.5f, 0.5f, 0.5f),
                                                   glm::vec3(0.9f, 0.8f, 0.7f));
  head->attachMesh(headBoundingBox);

  // Upper Arms
  auto leftUpperArm =
      std::make_shared<Bone>("leftUpperArm", glm::vec3(-0.45f, 0.2f, 0.0f),
                             glm::vec3(0.0f, 0.0f, 1.0f));
  auto leftUpperArmBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.3f, 0.6f, 0.3f), glm::vec3(0.7f, 0.3f, 0.2f));
  leftUpperArm->attachMesh(leftUpperArmBoundingBox);

  auto rightUpperArm =
      std::make_shared<Bone>("rightUpperArm", glm::vec3(0.45f, 0.2f, 0.0f),
                             glm::vec3(0.0f, 0.0f, 1.0f));
  auto rightUpperArmBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.3f, 0.6f, 0.3f), glm::vec3(0.7f, 0.3f, 0.2f));
  rightUpperArm->attachMesh(rightUpperArmBoundingBox);

  // Lower Arms
  auto leftLowerArm =
      std::make_shared<Bone>("leftLowerArm", glm::vec3(0.0f, -0.6f, 0.0f),
                             glm::vec3(0.0f, 0.0f, 1.0f));
  auto leftLowerArmBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.25f, 0.6f, 0.25f), glm::vec3(0.8f, 0.6f, 0.5f));
  leftLowerArm->attachMesh(leftLowerArmBoundingBox);

  auto rightLowerArm =
      std::make_shared<Bone>("rightLowerArm", glm::vec3(0.0f, -0.6f, 0.0f),
                             glm::vec3(0.0f, 0.0f, 1.0f));
  auto rightLowerArmBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.25f, 0.6f, 0.25f), glm::vec3(0.8f, 0.6f, 0.5f));
  rightLowerArm->attachMesh(rightLowerArmBoundingBox);

  // Upper Legs
  auto leftUpperLeg =
      std::make_shared<Bone>("leftUpperLeg", glm::vec3(-0.25f, -0.75f, 0.0f),
                             glm::vec3(1.0f, 0.0f, 0.0f));
  auto leftUpperLegBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.3f, 0.7f, 0.3f), glm::vec3(0.3f, 0.3f, 0.7f));
  leftUpperLeg->attachMesh(leftUpperLegBoundingBox);

  auto rightUpperLeg =
      std::make_shared<Bone>("rightUpperLeg", glm::vec3(0.25f, -0.75f, 0.0f),
                             glm::vec3(1.0f, 0.0f, 0.0f));
  auto rightUpperLegBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.3f, 0.7f, 0.3f), glm::vec3(0.3f, 0.3f, 0.7f));
  rightUpperLeg->attachMesh(rightUpperLegBoundingBox);

  // Lower Legs
  auto leftLowerLeg =
      std::make_shared<Bone>("leftLowerLeg", glm::vec3(0.0f, -0.7f, 0.0f),
                             glm::vec3(1.0f, 0.0f, 0.0f));
  auto leftLowerLegBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.25f, 0.7f, 0.25f), glm::vec3(0.4f, 0.4f, 0.8f));
  leftLowerLeg->attachMesh(leftLowerLegBoundingBox);

  auto rightLowerLeg =
      std::make_shared<Bone>("rightLowerLeg", glm::vec3(0.0f, -0.7f, 0.0f),
                             glm::vec3(1.0f, 0.0f, 0.0f));
  auto rightLowerLegBoundingBox = std::make_shared<BoxMesh>(
      glm::vec3(0.25f, 0.7f, 0.25f), glm::vec3(0.4f, 0.4f, 0.8f));
  rightLowerLeg->attachMesh(rightLowerLegBoundingBox);

  // Set up the bone hierarchy
  torso->addChild(head);
  torso->addChild(leftUpperArm);
  torso->addChild(rightUpperArm);
  torso->addChild(leftUpperLeg);
  torso->addChild(rightUpperLeg);

  leftUpperArm->addChild(leftLowerArm);
  rightUpperArm->addChild(rightLowerArm);
  leftUpperLeg->addChild(leftLowerLeg);
  rightUpperLeg->addChild(rightLowerLeg);

  // Set up the skeleton
  skeleton->setRootBone(torso);
  skeleton->setPosition(glm::vec3(0.0f, 0.0f, 0.0f));

  // Create objects for rendering
  for (const auto &boneName :
       {"torso", "head", "leftUpperArm", "rightUpperArm", "leftLowerArm",
        "rightLowerArm", "leftUpperLeg", "rightUpperLeg", "leftLowerLeg",
        "rightLowerLeg"}) {
    auto bone = skeleton->getBoneByName(boneName);
    if (bone && bone->getAttachedMesh()) {
      GL::Object obj = m_rasterizer.createObject();
      const auto &mesh = bone->getAttachedMesh();
      const auto &vertices = mesh->getVertices();
      const auto &normals = mesh->getNormals();
      GL::AttribBuf vertexBuf = m_rasterizer.createVertexAttribs(
          obj, 0, vertices.size(), vertices.data());
      GL::AttribBuf normalBuf = m_rasterizer.createVertexAttribs(
          obj, 1, normals.size(), normals.data());

      m_vertexBuffers.push_back(vertexBuf);
      m_normalBuffers.push_back(normalBuf);

      m_rasterizer.createTriangleIndices(obj, mesh->getTriangles().size(),
                                         mesh->getTriangles().data());
      m_rasterizer.createEdgeIndices(obj, mesh->getEdges().size(),
                                     mesh->getEdges().data());

      m_meshObjects.push_back(obj); // right place
    }
  }

  // Set the skeleton to the animator
  m_animator->setSkeleton(skeleton);
}

void AnimationRenderer::createKeyframes() {
  auto timeline = std::make_shared<AnimationTimeline>();

  // Keyframe 0: Neutral pose (T-pose) at t=0
  auto keyframe0 =
      std::make_shared<Keyframe>(0.0f, glm::vec3(0.0f, 0.0f, 0.0f));
  keyframe0->setBoneRotation("torso", 0.0f);
  keyframe0->setBoneRotation("head", 0.0f);
  keyframe0->setBoneRotation("leftUpperArm", 0.0f);
  keyframe0->setBoneRotation("rightUpperArm", 0.0f);
  keyframe0->setBoneRotation("leftLowerArm", 0.0f);
  keyframe0->setBoneRotation("rightLowerArm", 0.0f);
  keyframe0->setBoneRotation("leftUpperLeg", 0.0f);
  keyframe0->setBoneRotation("rightUpperLeg", 0.0f);
  keyframe0->setBoneRotation("leftLowerLeg", 0.0f);
  keyframe0->setBoneRotation("rightLowerLeg", 0.0f);
  timeline->addKeyframe(keyframe0);

  // Keyframe 1: Right arm raised, left leg forward at t=1.0
  auto keyframe1 =
      std::make_shared<Keyframe>(1.0f, glm::vec3(0.0f, 0.0f, 0.0f));
  keyframe1->setBoneRotation("torso", 0.0f);
  keyframe1->setBoneRotation("head", 0.2f);
  keyframe1->setBoneRotation("rightUpperArm", 1.5f); // Raised up
  keyframe1->setBoneRotation("rightLowerArm", 0.5f); // Slight bend
  keyframe1->setBoneRotation("leftUpperArm", 0.3f);  // Slight back movement
  keyframe1->setBoneRotation("leftLowerArm", 0.0f);
  keyframe1->setBoneRotation("leftUpperLeg", 0.6f);   // Forward step
  keyframe1->setBoneRotation("leftLowerLeg", -0.3f);  // Slight bend
  keyframe1->setBoneRotation("rightUpperLeg", -0.3f); // Back leg
  keyframe1->setBoneRotation("rightLowerLeg", 0.5f);  // Bend for balance
  timeline->addKeyframe(keyframe1);

  auto keyframe1_5 =
      std::make_shared<Keyframe>(1.5f, glm::vec3(0.0f, -0.15f, 0.0f));
  keyframe1_5->setBoneRotation("torso", 0.15f); // Halfway between 1.0 and 2.5
  keyframe1_5->setBoneRotation("head", 0.15f);  // Halfway looking down
  keyframe1_5->setBoneRotation("rightUpperArm",
                               1.15f); // Coming down from raised
  keyframe1_5->setBoneRotation("rightLowerArm", 0.75f); // More bent
  keyframe1_5->setBoneRotation("leftUpperArm", 0.55f);  // Moving outward
  keyframe1_5->setBoneRotation("leftLowerArm", 0.5f);   // Starting to bend
  keyframe1_5->setBoneRotation("leftUpperLeg", 0.45f);  // Halfway bent
  keyframe1_5->setBoneRotation("rightUpperLeg", -0.3f); // Maintaining bend
  keyframe1_5->setBoneRotation("leftLowerLeg", 0.0f);   // Straightening
  keyframe1_5->setBoneRotation("rightLowerLeg", 0.4f);  // Bent
  timeline->addKeyframe(keyframe1_5);

  // Keyframe 2: Return to neutral at t=2.0
  auto keyframe2 =
      std::make_shared<Keyframe>(2.0f, glm::vec3(0.0f, 0.0f, 0.0f));
  keyframe2->setBoneRotation("torso", 0.0f);
  keyframe2->setBoneRotation("head", 0.0f);
  keyframe2->setBoneRotation("leftUpperArm", 0.0f);
  keyframe2->setBoneRotation("rightUpperArm", 0.0f);
  keyframe2->setBoneRotation("leftLowerArm", 0.0f);
  keyframe2->setBoneRotation("rightLowerArm", 0.0f);
  keyframe2->setBoneRotation("leftUpperLeg", 0.0f);
  keyframe2->setBoneRotation("rightUpperLeg", 0.0f);
  keyframe2->setBoneRotation("leftLowerLeg", 0.0f);
  keyframe2->setBoneRotation("rightLowerLeg", 0.0f);
  timeline->addKeyframe(keyframe2);

  // Keyframe 3: Left arm raised, right leg forward at t=3.0
  auto keyframe3 =
      std::make_shared<Keyframe>(3.0f, glm::vec3(0.0f, 0.0f, 0.0f));
  keyframe3->setBoneRotation("torso", 0.0f);
  keyframe3->setBoneRotation("head", -0.2f);
  keyframe3->setBoneRotation("leftUpperArm", -1.5f);  // Raised up
  keyframe3->setBoneRotation("leftLowerArm", 0.5f);   // Slight bend
  keyframe3->setBoneRotation("rightUpperArm", -0.3f); // Slight back movement
  keyframe3->setBoneRotation("rightLowerArm", 0.0f);
  keyframe3->setBoneRotation("rightUpperLeg", 0.6f);  // Forward step
  keyframe3->setBoneRotation("rightLowerLeg", -0.3f); // Slight bend
  keyframe3->setBoneRotation("leftUpperLeg", -0.3f);  // Back leg
  keyframe3->setBoneRotation("leftLowerLeg", 0.5f);   // Bend for balance
  timeline->addKeyframe(keyframe3);

  auto keyframe3_5 =
      std::make_shared<Keyframe>(3.5f, glm::vec3(0.0f, -0.15f, 0.0f));
  keyframe3_5->setBoneRotation("torso", 0.15f); // Halfway between 3.0 and 4.0
  keyframe3_5->setBoneRotation("head", -0.15f); // Halfway looking down
  keyframe3_5->setBoneRotation("leftUpperArm",
                               -1.15f); // Coming down from raised
  keyframe3_5->setBoneRotation("leftLowerArm", 0.75f);   // More bent
  keyframe3_5->setBoneRotation("rightUpperArm", -0.55f); // Moving outward
  keyframe3_5->setBoneRotation("rightLowerArm", 0.5f);   // Starting to bend
  keyframe3_5->setBoneRotation("rightUpperLeg", 0.45f);  // Halfway bent
  keyframe3_5->setBoneRotation("leftUpperLeg", -0.3f);   // Maintaining bend
  keyframe3_5->setBoneRotation("rightLowerLeg", 0.0f);   // Straightening
  keyframe3_5->setBoneRotation("leftLowerLeg", 0.4f);    // Bent
  timeline->addKeyframe(keyframe3_5);

  // Keyframe 4: Return to neutral at t=4.0
  auto keyframe4 =
      std::make_shared<Keyframe>(4.0f, glm::vec3(0.0f, 0.0f, 0.0f));
  keyframe4->setBoneRotation("torso", 0.0f);
  keyframe4->setBoneRotation("head", 0.0f);
  keyframe4->setBoneRotation("leftUpperArm", 0.0f);
  keyframe4->setBoneRotation("rightUpperArm", 0.0f);
  keyframe4->setBoneRotation("leftLowerArm", 0.0f);
  keyframe4->setBoneRotation("rightLowerArm", 0.0f);
  keyframe4->setBoneRotation("leftUpperLeg", 0.0f);
  keyframe4->setBoneRotation("rightUpperLeg", 0.0f);
  keyframe4->setBoneRotation("leftLowerLeg", 0.0f);
  keyframe4->setBoneRotation("rightLowerLeg", 0.0f);
  timeline->addKeyframe(keyframe4);

  // Set looping to true
  timeline->setLooping(true);

  // Set the animation to the animator
  m_animator->setAnimation(timeline);
}

void AnimationRenderer::initialize() {
  createHumanoidSkeleton();

  // tkcomment std::cout //tkcomment << "Before creating keyframes : \n";

  for (size_t i = 0; i < m_meshObjects.size(); ++i) {
    // tkcomment std::cout //tkcomment << m_meshObjects[i].nTris //tkcomment <<
    // " tris" //tkcomment << m_meshObjects[i].nEdges
    // tkcomment << "edges \n";
  }
  createKeyframes();
}
void AnimationRenderer::updateMeshGeometry() {
  // Get the skeleton
  auto skeleton = m_animator->getSkeleton();
  if (!skeleton)
    return;

  // Transform mesh geometry based on bone positions
  int meshIndex = 0;
  for (const auto &boneName :
       {"torso", "head", "leftUpperArm", "rightUpperArm", "leftLowerArm",
        "rightLowerArm", "leftUpperLeg", "rightUpperLeg", "leftLowerLeg",
        "rightLowerLeg"}) {

    auto bone = skeleton->getBoneByName(boneName);
    if (!bone || !bone->getAttachedMesh() || meshIndex >= m_meshObjects.size())
      continue;

    // Get mesh data
    const auto &mesh = bone->getAttachedMesh();
    const auto &vertices = mesh->getVertices();
    const auto &normals = mesh->getNormals();

    // Get bone transformation
    const glm::mat4 &boneTransform = bone->getWorldTransform();

    // Pre-allocate vectors with correct size
    std::vector<glm::vec3> transformedVertices(vertices.size());
    std::vector<glm::vec3> transformedNormals(normals.size());

    // Transform vertices and normals
    for (size_t i = 0; i < vertices.size(); ++i) {
      // Transform vertex
      glm::vec4 worldPos = boneTransform * glm::vec4(vertices[i], 1.0f);
      transformedVertices[i] = glm::vec3(worldPos);

      // Transform normal
      glm::mat3 normalMatrix =
          glm::transpose(glm::inverse(glm::mat3(boneTransform)));
      transformedNormals[i] = glm::normalize(normalMatrix * normals[i]);
    }

    // Update buffers with error checking
    if (!transformedVertices.empty() && !transformedNormals.empty()) {
      m_rasterizer.updateVertexAttribs(m_vertexBuffers[meshIndex],
                                       transformedVertices.size(),
                                       transformedVertices.data());
      m_rasterizer.updateVertexAttribs(m_normalBuffers[meshIndex],
                                       transformedNormals.size(),
                                       transformedNormals.data());
    }

    ++meshIndex;
  }
}
void AnimationRenderer::update(float deltaTime) {
  // Update camera
  m_cameraControl.update();

  // Update animation
  m_animator->update(deltaTime);

  // Update mesh geometry based on new bone positions
  updateMeshGeometry();
}

void AnimationRenderer::render() {
  // Clear the screen
  m_rasterizer.clear(glm::vec4(0.4f, 0.4f, 0.4f, 1.0f));

  // Enable depth testing
  m_rasterizer.enableDepthTest();

  // Use shader program
  m_rasterizer.useShaderProgram(m_shaderProgram);

  // Set common uniforms
  m_rasterizer.setUniform(m_shaderProgram, "model", glm::mat4(1.0f));
  m_rasterizer.setUniform(m_shaderProgram, "view",
                          m_cameraControl.camera.getViewMatrix());
  m_rasterizer.setUniform(m_shaderProgram, "projection",
                          m_cameraControl.camera.getProjectionMatrix());
  m_rasterizer.setUniform(m_shaderProgram, "lightPos",
                          m_cameraControl.camera.position);
  m_rasterizer.setUniform(m_shaderProgram, "viewPos",
                          m_cameraControl.camera.position);
  m_rasterizer.setUniform(m_shaderProgram, "lightColor",
                          glm::vec3(1.0f, 1.0f, 1.0f));

  // Render filled faces
  m_rasterizer.setupFilledFaces();

  glm::vec3 orange(1.0f, 0.6f, 0.2f);
  glm::vec3 white(1.0f, 1.0f, 1.0f);
  m_rasterizer.setUniform(m_shaderProgram, "ambientColor", 0.2f * white);
  m_rasterizer.setUniform(m_shaderProgram, "extdiffuseColor", 0.9f * orange);
  m_rasterizer.setUniform(m_shaderProgram, "intdiffuseColor", 0.4f * orange);
  m_rasterizer.setUniform(m_shaderProgram, "specularColor", 0.6f * white);
  m_rasterizer.setUniform(m_shaderProgram, "phongExponent", 20.0f);

  // Draw all mesh objects (render only a few objects at a time to debug)
  // tkcomment std::cout //tkcomment << "Rendering triangles for " //tkcomment
  // << m_meshObjects.size()
  // tkcomment << " objects\n";
  // tkcomment std::cout //tkcomment << "Number of vertices\n";
  for (size_t i = 0; i < m_meshObjects.size(); ++i) {
    // tkcomment std::cout //tkcomment << m_meshObjects[i].nTris //tkcomment <<
    // " tris" //tkcomment << m_meshObjects[i].nEdges
    // tkcomment << "edges \n";
    // //tkcomment std::cout //tkcomment << "Triangles: " //tkcomment <<
    // object.nTris //tkcomment << "\n"
    //       //tkcomment << "Edges: " //tkcomment << object.nEdges //tkcomment
    //       << std::endl;

    m_rasterizer.drawTriangles(m_meshObjects[i]);
  }

  // Render wireframe
  m_rasterizer.setupWireFrame();

  glm::vec3 black(0.0f, 0.0f, 0.0f);
  m_rasterizer.setUniform(m_shaderProgram, "ambientColor", black);
  m_rasterizer.setUniform(m_shaderProgram, "extdiffuseColor", black);
  m_rasterizer.setUniform(m_shaderProgram, "intdiffuseColor", black);
  m_rasterizer.setUniform(m_shaderProgram, "specularColor", black);
  m_rasterizer.setUniform(m_shaderProgram, "phongExponent", 0.0f);

  // Draw wireframe edges
  for (const auto &obj : m_meshObjects) {
    m_rasterizer.drawEdges(obj);
  }

  // Show the rendered frame
  m_rasterizer.show();
}
// void AnimationRenderer::render() {
//   //tkcomment std::cout //tkcomment << "render called\n";
//   // Clear the screen
//   m_rasterizer.clear(glm::vec4(0.4f, 0.4f, 0.4f, 1.0f));

//   // Enable depth testing
//   m_rasterizer.enableDepthTest();

//   // Use shader program
//   m_rasterizer.useShaderProgram(m_shaderProgram);

//   // Set common uniforms
//   m_rasterizer.setUniform(m_shaderProgram, "model", glm::mat4(1.0f));
//   m_rasterizer.setUniform(m_shaderProgram, "view",
//                           m_cameraControl.camera.getViewMatrix());
//   m_rasterizer.setUniform(m_shaderProgram, "projection",
//                           m_cameraControl.camera.getProjectionMatrix());
//   m_rasterizer.setUniform(m_shaderProgram, "lightPos",
//                           m_cameraControl.camera.position);
//   m_rasterizer.setUniform(m_shaderProgram, "viewPos",
//                           m_cameraControl.camera.position);
//   m_rasterizer.setUniform(m_shaderProgram, "lightColor",
//                           glm::vec3(1.0f, 1.0f, 1.0f));

//   // Render filled faces
//   m_rasterizer.setupFilledFaces();

//   glm::vec3 orange(1.0f, 0.6f, 0.2f);
//   glm::vec3 white(1.0f, 1.0f, 1.0f);
//   m_rasterizer.setUniform(m_shaderProgram, "ambientColor", 0.2f * white);
//   m_rasterizer.setUniform(m_shaderProgram, "extdiffuseColor", 0.9f * orange);
//   m_rasterizer.setUniform(m_shaderProgram, "intdiffuseColor", 0.4f * orange);
//   m_rasterizer.setUniform(m_shaderProgram, "specularColor", 0.6f * white);
//   m_rasterizer.setUniform(m_shaderProgram, "phongExponent", 20.0f);

//   // Draw all mesh objects
//   //tkcomment std::cout //tkcomment << "Number of objects : " //tkcomment <<
//   m_meshObjects.size() //tkcomment <<
//   "\n";
//   // //tkcomment std::cout //tkcomment << "Meshobjects: " //tkcomment <<
//   m_meshObjects.size() //tkcomment <<
//   "\n"; for (auto i : m_meshObjects) {
//   }
//   int i = 0;
//   // // for (const auto &obj : m_meshObjects) {
//   // //   m_rasterizer.drawTriangles(obj);
//   // // }
//   // for (size_t i = 0; i < m_meshObjects.size(); ++i) {
//   //   // if (i % 2 == 0) { // Only odd indices
//   //   if (i < 4)
//   //     m_rasterizer.drawTriangles(m_meshObjects[i]);
//   //   // }
//   // }
//   for (const auto &obj : m_meshObjects) {
//     // i++;
//     // if (i > 3)
//     //   break;
//     m_rasterizer.drawTriangles(obj);
//   }

//   // Render wireframe
//   m_rasterizer.setupWireFrame();

//   glm::vec3 black(0.0f, 0.0f, 0.0f);
//   m_rasterizer.setUniform(m_shaderProgram, "ambientColor", black);
//   m_rasterizer.setUniform(m_shaderProgram, "extdiffuseColor", black);
//   m_rasterizer.setUniform(m_shaderProgram, "intdiffuseColor", black);
//   m_rasterizer.setUniform(m_shaderProgram, "specularColor", black);
//   m_rasterizer.setUniform(m_shaderProgram, "phongExponent", 0.0f);

//   // Draw wireframe edges
//   for (const auto &obj : m_meshObjects) {
//     m_rasterizer.drawEdges(obj);
//   }

//   // Show the rendered frame
//   m_rasterizer.show();
// }

int main() {
  // Initialize the rasterizer
  int width = 800, height = 600;
  GL::Rasterizer rasterizer;

  if (!rasterizer.initialize("Skeletal Animation", width, height)) {
    std::cerr << "Failed to initialize rasterizer" << std::endl;
    // std::cerr // tkcomment << "Failed to initialize rasterizer" //tkcomment
    // <<
    //           // std::endl;
    return EXIT_FAILURE;
  }

  // Create the animation renderer
  AnimationRenderer renderer(rasterizer, width, height);

  // Initialize the scene
  renderer.initialize();

  // Main loop
  float framerate = 1e-3f;
  float lastTime = SDL_GetTicks64() * framerate;

  // while (!rasterizer.shouldQuit()) {
  //   // Calculate delta time
  //   float currentTime = SDL_GetTicks64() * framerate;
  //   float deltaTime = currentTime - lastTime;
  //   lastTime = currentTime;

  //   // Update and render
  //   renderer.update(deltaTime);
  //   renderer.render();
  // }

  // In main() function, modify the main loop:
  while (!rasterizer.shouldQuit()) {
    // Calculate delta time
    float currentTime = SDL_GetTicks64() * 1e-3f;
    float deltaTime = currentTime - lastTime;
    lastTime = currentTime;

    // Update and render
    renderer.update(deltaTime);
    renderer.render();

    // Add a small delay to limit framerate
    // SDL_Delay(200); // ~60 FPS
  }

  return EXIT_SUCCESS;
}