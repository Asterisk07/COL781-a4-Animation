#include "obstacle.hpp"
#include "spring.hpp"
#include <vector>

namespace COL781 {

// SphereObstacle implementation
SphereObstacle::SphereObstacle(const glm::vec3 &center, float radius,
                               const glm::vec3 &linearVel,
                               const glm::vec3 &angularVel,
                               float collisionRadius)
    : center(center), radius(radius), linearVelocity(linearVel),
      angularVelocity(angularVel),
      orientation(glm::quat(1.0f, 0.0f, 0.0f, 0.0f)) {

  // If collision radius not specified, use the visual radius
  this->collisionRadius = (collisionRadius < 0.0f) ? radius : collisionRadius;
}

bool SphereObstacle::checkCollision(Particle *particle,
                                    glm::vec3 &collisionPoint,
                                    glm::vec3 &normal,
                                    float &penetrationDepth) {
  glm::vec3 direction = particle->position - center;
  float distance = glm::length(direction);

  // Check if particle is inside the sphere
  if (distance <= collisionRadius) {
    // Normalize the direction
    normal =
        (distance > 0.0f) ? direction / distance : glm::vec3(0.0f, 1.0f, 0.0f);

    // Calculate penetration depth
    penetrationDepth = collisionRadius - distance;

    // Collision point is on the sphere surface in the direction from center to
    // particle
    collisionPoint = center + normal * collisionRadius;

    return true;
  }

  return false;
}

glm::vec3 SphereObstacle::getVelocityAtPoint(const glm::vec3 &point) const {
  // Linear velocity + (angular velocity cross (point - center))
  return linearVelocity + glm::cross(angularVelocity, point - center);
}

void SphereObstacle::update(float dt) {
  // Update position
  center += linearVelocity * dt;

  // Update orientation
  if (glm::length(angularVelocity) > 0.0f) {
    // Convert angular velocity to a quaternion change
    float angle = glm::length(angularVelocity) * dt;
    glm::vec3 axis = glm::normalize(angularVelocity);
    glm::quat rot = glm::angleAxis(angle, axis);

    // Apply rotation
    orientation = rot * orientation;
    orientation = glm::normalize(orientation);
  }
}

void SphereObstacle::initializeRendering(OpenGL::Rasterizer &r) {
  auto vertices = generateSphereVertices();
  auto normals = generateSphereNormals(vertices);
  auto texCoords = generateSphereTexCoords();
  auto triangles = generateSphereTriangles();

  sphereObj = r.createObject();

  vertexBuf =
      r.createVertexAttribs(sphereObj, 0, vertices.size(), vertices.data());
  normalBuf =
      r.createVertexAttribs(sphereObj, 1, normals.size(), normals.data());
  texBuf =
      r.createVertexAttribs(sphereObj, 2, texCoords.size(), texCoords.data());

  r.createTriangleIndices(sphereObj, triangles.size(), triangles.data());

  renderingInitialized = true;
}

std::vector<glm::vec3> SphereObstacle::generateSphereVertices(int resolution) {
  std::vector<glm::vec3> vertices;

  // Generate vertices for a unit sphere
  for (int j = 0; j <= resolution; ++j) {
    float theta = j * glm::pi<float>() / resolution;
    float sinTheta = sin(theta);
    float cosTheta = cos(theta);

    for (int i = 0; i <= resolution; ++i) {
      float phi = i * 2.0f * glm::pi<float>() / resolution;
      float sinPhi = sin(phi);
      float cosPhi = cos(phi);

      float x = cosPhi * sinTheta;
      float y = cosTheta;
      float z = sinPhi * sinTheta;

      vertices.push_back(glm::vec3(x, y, z));
    }
  }

  return vertices;
}

std::vector<glm::vec3>
SphereObstacle::generateSphereNormals(const std::vector<glm::vec3> &vertices) {
  // For a sphere centered at origin, normals are just normalized vertices
  std::vector<glm::vec3> normals;
  for (const auto &v : vertices) {
    normals.push_back(glm::normalize(v));
  }
  return normals;
}

std::vector<glm::vec2> SphereObstacle::generateSphereTexCoords(int resolution) {
  std::vector<glm::vec2> texCoords;

  for (int j = 0; j <= resolution; ++j) {
    float v = (float)j / resolution;

    for (int i = 0; i <= resolution; ++i) {
      float u = (float)i / resolution;
      texCoords.push_back(glm::vec2(u, v));
    }
  }

  return texCoords;
}

std::vector<glm::ivec3>
SphereObstacle::generateSphereTriangles(int resolution) {
  std::vector<glm::ivec3> triangles;

  for (int j = 0; j < resolution; ++j) {
    for (int i = 0; i < resolution; ++i) {
      int p1 = j * (resolution + 1) + i;
      int p2 = p1 + 1;
      int p3 = (j + 1) * (resolution + 1) + i;
      int p4 = p3 + 1;

      triangles.push_back(glm::ivec3(p1, p2, p3));
      triangles.push_back(glm::ivec3(p2, p4, p3));
    }
  }

  return triangles;
}

void SphereObstacle::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                          const Camera &camera) {
  if (!renderingInitialized) {
    initializeRendering(r);
  }

  // Create model matrix for sphere
  glm::mat4 model = glm::mat4(1.0f);
  model = glm::translate(model, center);
  model = model * glm::toMat4(orientation);
  model = glm::scale(model, glm::vec3(radius));

  // Set uniforms and draw
  r.setupFilledFaces();
  r.setUniform(program, "model", model);

  // Set lighting parameters
  glm::vec3 sphereColor(0.7f, 0.2f, 0.2f); // Red sphere
  glm::vec3 white(1.0f, 1.0f, 1.0f);
  r.setUniform(program, "ambientColor", 0.2f * white);
  r.setUniform(program, "extdiffuseColor", 0.8f * sphereColor);
  r.setUniform(program, "intdiffuseColor", 0.4f * sphereColor);
  r.setUniform(program, "specularColor", 0.6f * white);
  r.setUniform(program, "phongExponent", 20.0f);

  r.drawTriangles(sphereObj);
}

// PlaneObstacle implementation
PlaneObstacle::PlaneObstacle(const glm::vec3 &point, const glm::vec3 &normal,
                             float size, const glm::vec3 &linearVel,
                             const glm::vec3 &angularVel)
    : point(point), normal(glm::normalize(normal)), size(size),
      linearVelocity(linearVel), angularVelocity(angularVel),
      orientation(glm::quat(1.0f, 0.0f, 0.0f, 0.0f)) {}

bool PlaneObstacle::checkCollision(Particle *particle,
                                   glm::vec3 &collisionPoint, glm::vec3 &normal,
                                   float &penetrationDepth) {
  // Calculate signed distance from point to plane
  float signedDistance = glm::dot(particle->position - point, this->normal);

  // Check if particle is on the wrong side of the plane
  if (signedDistance <= 0.0f) {
    normal = this->normal;
    penetrationDepth = -signedDistance;
    collisionPoint = particle->position - normal * penetrationDepth;
    return true;
  }

  return false;
}

glm::vec3 PlaneObstacle::getVelocityAtPoint(const glm::vec3 &point) const {
  // For a plane, we need to consider both linear and angular velocity
  // The angular velocity creates a rotational effect around the plane's point
  return linearVelocity + glm::cross(angularVelocity, point - this->point);
}

void PlaneObstacle::update(float dt) {
  // Update position
  point += linearVelocity * dt;

  // Update orientation if there's angular velocity
  if (glm::length(angularVelocity) > 0.0f) {
    float angle = glm::length(angularVelocity) * dt;
    glm::vec3 axis = glm::normalize(angularVelocity);
    glm::quat rot = glm::angleAxis(angle, axis);

    // Apply rotation to orientation and normal
    orientation = rot * orientation;
    orientation = glm::normalize(orientation);

    // Rotate the normal
    normal = glm::rotate(orientation, glm::vec3(0.0f, 1.0f, 0.0f));
    normal = glm::normalize(normal);
  }
}

void PlaneObstacle::initializeRendering(OpenGL::Rasterizer &r) {
  // Create a quad for the plane
  std::vector<glm::vec3> vertices = {glm::vec3(-size / 2, 0.0f, -size / 2),
                                     glm::vec3(size / 2, 0.0f, -size / 2),
                                     glm::vec3(size / 2, 0.0f, size / 2),
                                     glm::vec3(-size / 2, 0.0f, size / 2)};

  std::vector<glm::vec3> normals(4, glm::vec3(0.0f, 1.0f, 0.0f));

  std::vector<glm::vec2> texCoords = {
      glm::vec2(0.0f, 0.0f), glm::vec2(1.0f, 0.0f), glm::vec2(1.0f, 1.0f),
      glm::vec2(0.0f, 1.0f)};

  std::vector<glm::ivec3> triangles = {glm::ivec3(0, 1, 2),
                                       glm::ivec3(0, 2, 3)};

  planeObj = r.createObject();

  vertexBuf =
      r.createVertexAttribs(planeObj, 0, vertices.size(), vertices.data());
  normalBuf =
      r.createVertexAttribs(planeObj, 1, normals.size(), normals.data());
  texBuf =
      r.createVertexAttribs(planeObj, 2, texCoords.size(), texCoords.data());

  r.createTriangleIndices(planeObj, triangles.size(), triangles.data());

  renderingInitialized = true;
}

void PlaneObstacle::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                         const Camera &camera) {
  if (!renderingInitialized) {
    initializeRendering(r);
  }

  // Find the rotation from (0,1,0) to the plane normal
  glm::vec3 defaultNormal(0.0f, 1.0f, 0.0f);
  glm::vec3 rotationAxis = glm::cross(defaultNormal, normal);
  float rotationAngle = acos(glm::dot(defaultNormal, normal));

  // Create model matrix for the plane
  glm::mat4 model = glm::mat4(1.0f);
  model = glm::translate(model, point);

  // If we need to rotate the plane
  if (glm::length(rotationAxis) > 0.001f) {
    model = model * glm::rotate(glm::mat4(1.0f), rotationAngle,
                                glm::normalize(rotationAxis));
  }

  // Set uniforms and draw
  r.setupFilledFaces();
  r.setUniform(program, "model", model);

  // Set lighting parameters
  glm::vec3 planeColor(0.5f, 0.5f, 0.5f); // Gray plane
  glm::vec3 white(1.0f, 1.0f, 1.0f);
  r.setUniform(program, "ambientColor", 0.2f * white);
  r.setUniform(program, "extdiffuseColor", 0.8f * planeColor);
  r.setUniform(program, "intdiffuseColor", 0.4f * planeColor);
  r.setUniform(program, "specularColor", 0.6f * white);
  r.setUniform(program, "phongExponent", 20.0f);

  r.drawTriangles(planeObj);
}

} // namespace COL781