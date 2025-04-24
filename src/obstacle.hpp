#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include "spring.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>

namespace COL781 {

// Forward declaration
class Particle;

// Base class for obstacles
class Obstacle {
public:
  virtual ~Obstacle() {}

  // Check collision with a particle and provide collision info
  virtual bool checkCollision(Particle *particle, glm::vec3 &collisionPoint,
                              glm::vec3 &normal, float &penetrationDepth) = 0;

  // Get velocity of the obstacle at a given point (for relative velocity
  // calculation)
  virtual glm::vec3 getVelocityAtPoint(const glm::vec3 &point) const = 0;

  // Draw the obstacle
  virtual void draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                    const Camera &camera) = 0;

  // Update obstacle position/orientation
  virtual void update(float dt) = 0;
};

// Sphere obstacle
class SphereObstacle : public Obstacle {
public:
  SphereObstacle(const glm::vec3 &center, float radius,
                 const glm::vec3 &linearVel = glm::vec3(0.0f),
                 const glm::vec3 &angularVel = glm::vec3(0.0f),
                 float collisionRadius = -1.0f);

  bool checkCollision(Particle *particle, glm::vec3 &collisionPoint,
                      glm::vec3 &normal, float &penetrationDepth) override;

  glm::vec3 getVelocityAtPoint(const glm::vec3 &point) const override;

  void draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
            const Camera &camera) override;

  void update(float dt) override;

private:
  glm::vec3 center;
  float radius;
  float collisionRadius; // Radius used for collision detection (can be larger
                         // than visual radius)
  glm::vec3 linearVelocity;
  glm::vec3 angularVelocity; // Angular velocity around each axis
  glm::quat orientation;

  // For rendering
  OpenGL::Object sphereObj;
  OpenGL::AttribBuf vertexBuf, normalBuf, texBuf;
  bool renderingInitialized = false;

  void initializeRendering(OpenGL::Rasterizer &r);
  std::vector<glm::vec3> generateSphereVertices(int resolution = 32);
  std::vector<glm::vec3>
  generateSphereNormals(const std::vector<glm::vec3> &vertices);
  std::vector<glm::vec2> generateSphereTexCoords(int resolution = 32);
  std::vector<glm::ivec3> generateSphereTriangles(int resolution = 32);
};

// Plane obstacle (infinite plane)
class PlaneObstacle : public Obstacle {
public:
  PlaneObstacle(const glm::vec3 &point, const glm::vec3 &normal,
                float size = 10.0f,
                const glm::vec3 &linearVel = glm::vec3(0.0f),
                const glm::vec3 &angularVel = glm::vec3(0.0f));

  bool checkCollision(Particle *particle, glm::vec3 &collisionPoint,
                      glm::vec3 &normal, float &penetrationDepth) override;

  glm::vec3 getVelocityAtPoint(const glm::vec3 &point) const override;

  void draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
            const Camera &camera) override;

  void update(float dt) override;

private:
  glm::vec3 point;  // A point on the plane
  glm::vec3 normal; // Plane normal (normalized)
  float size;       // Size of the visual plane quad
  glm::vec3 linearVelocity;
  glm::vec3 angularVelocity;
  glm::quat orientation;

  // For rendering
  OpenGL::Object planeObj;
  OpenGL::AttribBuf vertexBuf, normalBuf, texBuf;
  bool renderingInitialized = false;

  void initializeRendering(OpenGL::Rasterizer &r);
};

} // namespace COL781

#endif // OBSTACLE_HPP