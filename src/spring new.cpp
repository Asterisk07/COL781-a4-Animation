#include "spring.hpp"
#include "obstacle.hpp"
#include <glm/glm.hpp>
#include <iostream>

namespace COL781 {

// Particle implementation
Particle::Particle(const glm::vec3 &pos, float m, bool fixed)
    : position(pos), velocity(glm::vec3(0.0f)), force(glm::vec3(0.0f)), mass(m),
      isFixed(fixed) {}

void Particle::applyForce(const glm::vec3 &f) { force += f; }

void Particle::update(float dt) {
  if (isFixed)
    return;

  // Semi-implicit Euler integration
  velocity += force / mass * dt;
  position += velocity * dt;

  // Reset force for next step
  force = glm::vec3(0.0f);
}

void Particle::reset(const glm::vec3 &pos) {
  position = pos;
  velocity = glm::vec3(0.0f);
  force = glm::vec3(0.0f);
}

// Spring implementation
Spring::Spring(Particle *particle1, Particle *particle2, float rest_length,
               float spring_const, float damping_const)
    : p1(particle1), p2(particle2), restLength(rest_length), ks(spring_const),
      kd(damping_const) {}

void Spring::applyForce() {
  // Vector from p1 to p2
  glm::vec3 direction = p2->position - p1->position;
  float length = glm::length(direction);

  if (length > 0.00001f) { // Avoid division by zero
    direction /= length;   // Normalize the direction

    // Compute relative velocity
    glm::vec3 velDiff = p2->velocity - p1->velocity;

    // Spring force (Hooke's law) + Damping force
    float springForce = ks * (length - restLength);
    float dampingForce = kd * glm::dot(velDiff, direction);

    glm::vec3 force = direction * (springForce + dampingForce);

    // Apply forces to particles
    p1->applyForce(force);
    p2->applyForce(-force); // Equal and opposite force
  }
}

// ClothSystem implementation
ClothSystem::ClothSystem(const ClothConfig &config) : config(config) {
  // Set default values for collision parameters if not provided
  if (config.restitution <= 0.0f) {
    this->config.restitution = 0.5f; // Default restitution
  }
  if (config.friction <= 0.0f) {
    this->config.friction = 0.3f; // Default friction
  }

  initialize();
}

ClothSystem::~ClothSystem() {
  // Clean up particles and springs
  for (auto p : particles) {
    delete p;
  }
  for (auto s : springs) {
    delete s;
  }
  // Note: We don't delete obstacles as they are owned by the main program
}

void ClothSystem::initialize() {
  // Create particles
  createParticles();

  // Create structural springs
  createSprings(SpringType::STRUCTURAL);

  // Create shear springs
  createSprings(SpringType::SHEAR);

  // Create bending springs
  createSprings(SpringType::BENDING);

  // Initialize triangle indices for rendering
  triangles.clear();
  for (int j = 0; j < config.resolutionY - 1; ++j) {
    for (int i = 0; i < config.resolutionX - 1; ++i) {
      int p1 = getParticleIndex(i, j);
      int p2 = getParticleIndex(i + 1, j);
      int p3 = getParticleIndex(i, j + 1);
      int p4 = getParticleIndex(i + 1, j + 1);

      // Two triangles per grid cell
      triangles.push_back(glm::ivec3(p1, p2, p3));
      triangles.push_back(glm::ivec3(p2, p4, p3));
    }
  }

  // Initialize edges for debug rendering
  edges.clear();
  for (auto spring : springs) {
    if (spring) {
      int p1Index = -1, p2Index = -1;
      for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i] == spring->p1)
          p1Index = i;
        if (particles[i] == spring->p2)
          p2Index = i;
      }
      if (p1Index >= 0 && p2Index >= 0) {
        edges.push_back(glm::ivec2(p1Index, p2Index));
      }
    }
  }
}

void ClothSystem::reset() {
  // Recreate all particles and springs
  for (auto p : particles) {
    delete p;
  }
  for (auto s : springs) {
    delete s;
  }
  particles.clear();
  springs.clear();

  // Reinitialize
  initialize();
}

void ClothSystem::update(float dt) {
  // Apply gravity to all particles
  for (auto p : particles) {
    p->applyForce(glm::vec3(0.0f, -config.gravity * p->mass, 0.0f));
  }

  // Apply spring forces
  for (auto s : springs) {
    s->applyForce();
  }

  // Handle collisions with obstacles
  handleCollisions();

  // Update particles
  for (auto p : particles) {
    p->update(dt);
  }

  // Update normals for rendering
  updateNormals();
}

int ClothSystem::getParticleIndex(int i, int j) const {
  return j * config.resolutionX + i;
}

void ClothSystem::createParticles() {
  particles.clear();

  float dx = config.width / (config.resolutionX - 1);
  float dy = config.height / (config.resolutionY - 1);

  for (int j = 0; j < config.resolutionY; ++j) {
    for (int i = 0; i < config.resolutionX; ++i) {
      // Initial position (x, y, z)
      glm::vec3 pos(i * dx, 0.0f, j * dy);

      // Adjust position to center the cloth
      pos.x -= config.width / 2.0f;
      pos.z -= config.height / 2.0f;

      // Check if this is a fixed corner
      bool isFixed = false;
      int cornerIndex = -1;

      if (i == 0 && j == 0)
        cornerIndex = 0; // Top-left
      else if (i == config.resolutionX - 1 && j == 0)
        cornerIndex = 1; // Top-right
      else if (i == 0 && j == config.resolutionY - 1)
        cornerIndex = 2; // Bottom-left
      else if (i == config.resolutionX - 1 && j == config.resolutionY - 1)
        cornerIndex = 3; // Bottom-right

      if (cornerIndex >= 0) {
        for (int corner : config.fixedCorners) {
          if (corner == cornerIndex) {
            isFixed = true;
            break;
          }
        }
      }

      // Create and add the particle
      Particle *p = new Particle(pos, config.particleMass, isFixed);
      particles.push_back(p);
    }
  }
}

void ClothSystem::createSprings(SpringType type) {
  switch (type) {
  case SpringType::STRUCTURAL: {
    // Create structural springs (connecting adjacent particles horizontally and
    // vertically)
    for (int j = 0; j < config.resolutionY; ++j) {
      for (int i = 0; i < config.resolutionX; ++i) {
        // Get current particle
        int idx = getParticleIndex(i, j);
        Particle *p = particles[idx];

        // Connect to right neighbor
        if (i < config.resolutionX - 1) {
          int rightIdx = getParticleIndex(i + 1, j);
          Particle *rightP = particles[rightIdx];
          float restLength = glm::distance(p->position, rightP->position);
          Spring *s = new Spring(p, rightP, restLength, config.kStructural,
                                 config.dStructural);
          springs.push_back(s);
        }

        // Connect to bottom neighbor
        if (j < config.resolutionY - 1) {
          int bottomIdx = getParticleIndex(i, j + 1);
          Particle *bottomP = particles[bottomIdx];
          float restLength = glm::distance(p->position, bottomP->position);
          Spring *s = new Spring(p, bottomP, restLength, config.kStructural,
                                 config.dStructural);
          springs.push_back(s);
        }
      }
    }
    break;
  }

  case SpringType::SHEAR: {
    // Create shear springs (connecting diagonal particles)
    for (int j = 0; j < config.resolutionY - 1; ++j) {
      for (int i = 0; i < config.resolutionX - 1; ++i) {
        // Get top-left particle
        int tlIdx = getParticleIndex(i, j);
        Particle *tlP = particles[tlIdx];

        // Get bottom-right particle
        int brIdx = getParticleIndex(i + 1, j + 1);
        Particle *brP = particles[brIdx];

        // Connect top-left to bottom-right
        float restLength = glm::distance(tlP->position, brP->position);
        Spring *s1 =
            new Spring(tlP, brP, restLength, config.kShear, config.dShear);
        springs.push_back(s1);

        // Get top-right particle
        int trIdx = getParticleIndex(i + 1, j);
        Particle *trP = particles[trIdx];

        // Get bottom-left particle
        int blIdx = getParticleIndex(i, j + 1);
        Particle *blP = particles[blIdx];

        // Connect top-right to bottom-left
        restLength = glm::distance(trP->position, blP->position);
        Spring *s2 =
            new Spring(trP, blP, restLength, config.kShear, config.dShear);
        springs.push_back(s2);
      }
    }
    break;
  }

  case SpringType::BENDING: {
    // Create bending springs (connecting particles two steps away)
    for (int j = 0; j < config.resolutionY; ++j) {
      for (int i = 0; i < config.resolutionX; ++i) {
        // Get current particle
        int idx = getParticleIndex(i, j);
        Particle *p = particles[idx];

        // Connect to particle two steps to the right
        if (i < config.resolutionX - 2) {
          int rightIdx = getParticleIndex(i + 2, j);
          Particle *rightP = particles[rightIdx];
          float restLength = glm::distance(p->position, rightP->position);
          Spring *s = new Spring(p, rightP, restLength, config.kBending,
                                 config.dBending);
          springs.push_back(s);
        }

        // Connect to particle two steps down
        if (j < config.resolutionY - 2) {
          int bottomIdx = getParticleIndex(i, j + 2);
          Particle *bottomP = particles[bottomIdx];
          float restLength = glm::distance(p->position, bottomP->position);
          Spring *s = new Spring(p, bottomP, restLength, config.kBending,
                                 config.dBending);
          springs.push_back(s);
        }
      }
    }
    break;
  }
  }
}

void ClothSystem::updateNormals() {
  // Initialize vertex and normal buffers
  vertices.resize(particles.size());
  normals.resize(particles.size(), glm::vec3(0.0f));

  // Copy particle positions to vertices
  for (size_t i = 0; i < particles.size(); ++i) {
    vertices[i] = particles[i]->position;
  }

  // Calculate normals for each triangle and accumulate at vertices
  for (const auto &tri : triangles) {
    glm::vec3 v1 = vertices[tri.y] - vertices[tri.x];
    glm::vec3 v2 = vertices[tri.z] - vertices[tri.x];
    glm::vec3 normal = glm::cross(v1, v2);

    // Accumulate normals (will be normalized later)
    normals[tri.x] += normal;
    normals[tri.y] += normal;
    normals[tri.z] += normal;
  }

  // Normalize all normals
  for (auto &normal : normals) {
    normal = glm::normalize(normal);
  }
}

void ClothSystem::initializeRendering(OpenGL::Rasterizer &r) {
  // Create object for rendering
  object = r.createObject();

  // Create vertex and normal buffers
  vertexBuf =
      r.createVertexAttribs(object, 0, vertices.size(), vertices.data());
  normalBuf = r.createVertexAttribs(object, 1, normals.size(), normals.data());

  // Create triangle indices
  r.createTriangleIndices(object, triangles.size(), triangles.data());
}

void ClothSystem::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                       const Camera &camera) {
  // Initialize rendering if not done yet
  if (object.tri_vao == 0) {
    initializeRendering(r);
  }

  // Update vertex and normal buffers
  r.updateVertexAttribs(vertexBuf, vertices.size(), vertices.data());
  r.updateVertexAttribs(normalBuf, normals.size(), normals.data());

  // Draw cloth
  r.setupFilledFaces();
  glm::vec3 clothColor(0.7f, 0.7f, 0.9f); // Cloth color
  glm::vec3 white(1.0f, 1.0f, 1.0f);
  r.setUniform(program, "ambientColor", 0.2f * white);
  r.setUniform(program, "extdiffuseColor", 0.8f * clothColor);
  r.setUniform(program, "intdiffuseColor", 0.4f * clothColor);
  r.setUniform(program, "specularColor", 0.6f * white);
  r.setUniform(program, "phongExponent", 20.0f);
  r.drawTriangles(object);

  // Draw all obstacles
  for (auto obstacle : obstacles) {
    obstacle->draw(r, program, camera);
  }
}

// void ClothSystem::initializeRendering(OpenGL::Rasterizer &r) {
//   // Create vertex and normal buffers
//   vertexBuf = r.createAttribBuf(vertices);
//   normalBuf = r.createAttribBuf(normals);

//   // Create object and set attributes
//   object = r.createObject();
//   r.setVertexAttribs(object, vertexBuf);
//   r.setVertexAttribs(object, normalBuf, 1);
//   r.setTriangleIndices(object, triangles);
// }

// void ClothSystem::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
//                        const Camera &camera) {
//   // Initialize rendering if not done yet
//   if (object == 0) {
//     initializeRendering(r);
//   }

//   // Update vertex and normal buffers
//   r.updateAttribBuf(vertexBuf, vertices);
//   r.updateAttribBuf(normalBuf, normals);

//   // Draw cloth
//   r.setUniform(program, "objectColor",
//                glm::vec3(0.7f, 0.7f, 0.9f)); // Cloth color
//   r.drawObject(object);

//   // Draw all obstacles
//   for (auto obstacle : obstacles) {
//     obstacle->draw(r, program, camera);
//   }
// }

void ClothSystem::addObstacle(Obstacle *obstacle) {
  if (obstacle) {
    obstacles.push_back(obstacle);
  }
}

void ClothSystem::handleCollisions() {
  // Check each particle against each obstacle
  for (auto p : particles) {
    if (!p->isFixed) { // Skip fixed particles
      for (auto obstacle : obstacles) {
        handleCollision(p, obstacle);
      }
    }
  }
}

void ClothSystem::handleCollision(Particle *particle, Obstacle *obstacle) {
  glm::vec3 collisionPoint, normal;
  float penetrationDepth;

  // Check if there is a collision
  if (obstacle->checkCollision(particle, collisionPoint, normal,
                               penetrationDepth)) {
    // Get obstacle velocity at the collision point
    glm::vec3 obstacleVelocity = obstacle->getVelocityAtPoint(collisionPoint);

    // Calculate relative velocity
    glm::vec3 relativeVelocity = particle->velocity - obstacleVelocity;

    // Normal component of relative velocity
    float normalVelocityMagnitude = glm::dot(relativeVelocity, normal);

    // Only respond to collision if we're moving into the obstacle
    if (normalVelocityMagnitude < 0.0f) {
      // 1. Position correction to move particle out of obstacle
      particle->position += normal * penetrationDepth;

      // 2. Velocity response: normal component with restitution
      glm::vec3 normalVelocity = normal * normalVelocityMagnitude;
      glm::vec3 tangentialVelocity = relativeVelocity - normalVelocity;

      // Apply restitution (coefficient of restitution = epsilon)
      particle->velocity -= (1.0f + config.restitution) * normalVelocity;

      // 3. Apply friction (tangential resistance)
      float tangentialSpeed = glm::length(tangentialVelocity);
      if (tangentialSpeed > 0.0001f) { // Avoid division by zero
        glm::vec3 frictionDir = tangentialVelocity / tangentialSpeed;
        float frictionMagnitude =
            config.friction * std::abs(normalVelocityMagnitude);

        // Apply friction force (limited by tangential velocity magnitude)
        float maxFriction = tangentialSpeed;
        float appliedFriction = std::min(frictionMagnitude, maxFriction);

        // Apply friction impulse
        particle->velocity -= frictionDir * appliedFriction;
      }

      // Add obstacle velocity back to particle
      particle->velocity += obstacleVelocity;
    }
  }
}

} // namespace COL781