#include "spring.hpp"
#include <iostream>

namespace COL781 {

// Particle Implementation
Particle::Particle(const glm::vec3 &pos, float m, bool fixed)
    : position(pos), velocity(glm::vec3(0.0f)), force(glm::vec3(0.0f)), mass(m),
      isFixed(fixed) {}

void Particle::applyForce(const glm::vec3 &f) {
  if (!isFixed) {
    force += f;
  }
}

void Particle::update(float dt) {
  if (!isFixed) {
    // Semi-implicit Euler integration
    glm::vec3 acceleration = force / mass;
    velocity += acceleration * dt;
    position += velocity * dt;

    // Reset force for next frame
    force = glm::vec3(0.0f);
  }
}

void Particle::reset(const glm::vec3 &pos) {
  position = pos;
  velocity = glm::vec3(0.0f);
  force = glm::vec3(0.0f);
}

// Spring Implementation
Spring::Spring(Particle *particle1, Particle *particle2, float rest_length,
               float spring_const, float damping_const)
    : p1(particle1), p2(particle2), restLength(rest_length), ks(spring_const),
      kd(damping_const) {}

void Spring::applyForce() {
  // Calculate spring direction vector
  glm::vec3 direction = p2->position - p1->position;
  float distance = glm::length(direction);

  // Avoid division by zero
  if (distance < 0.00001f) {
    return;
  }

  // Normalize direction
  glm::vec3 directionNorm = direction / distance;

  // Calculate spring force using Hooke's law (F = -k * (|x| - L0) * x_hat)
  float displacement = distance - restLength;
  glm::vec3 springForce = ks * displacement * directionNorm;

  // Calculate damping force (F_d = -kd * (v_rel · x_hat) * x_hat)
  glm::vec3 relativeVelocity = p2->velocity - p1->velocity;
  float dampingMagnitude = kd * glm::dot(relativeVelocity, directionNorm);
  glm::vec3 dampingForce = dampingMagnitude * directionNorm;

  // Total force
  glm::vec3 totalForce = springForce + dampingForce;

  // Apply forces to particles (equal and opposite)
  p1->applyForce(totalForce);
  p2->applyForce(-totalForce);
}

// ClothSystem Implementation
ClothSystem::ClothSystem(const ClothConfig &cfg) : config(cfg) { initialize(); }

ClothSystem::~ClothSystem() {
  // Clean up memory
  for (auto &p : particles) {
    delete p;
  }
  for (auto &s : springs) {
    delete s;
  }
  particles.clear();
  springs.clear();
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
}

void ClothSystem::reset() {
  // Delete all particles and springs
  for (auto &p : particles) {
    delete p;
  }
  for (auto &s : springs) {
    delete s;
  }
  particles.clear();
  springs.clear();

  // Reinitialize
  initialize();
}

void ClothSystem::update(float dt) {
  // Apply gravity
  for (auto &p : particles) {
    // Apply gravity force: F = m * g
    p->applyForce(glm::vec3(0.0f, -config.gravity * p->mass, 0.0f));
  }

  // Apply spring forces
  for (auto &s : springs) {
    s->applyForce();
  }

  // Update particles
  for (auto &p : particles) {
    p->update(dt);
  }
}

int ClothSystem::getParticleIndex(int i, int j) const {
  return i * config.resolutionX + j;
}

void ClothSystem::createParticles() {
  float dx = config.width / (config.resolutionX - 1);
  float dz = config.height / (config.resolutionY - 1);

  // Create particles in a grid
  for (int i = 0; i < config.resolutionY; ++i) {
    for (int j = 0; j < config.resolutionX; ++j) {
      // Calculate position (initially horizontal sheet in XZ plane)
      float x = j * dx;
      float y = 0.0f; // Initially flat horizontal cloth
      float z = i * dz;

      // Check if this particle should be fixed
      bool isFixed = false;

      // Check corners: 0=top-left, 1=top-right, 2=bottom-left, 3=bottom-right
      for (int corner : config.fixedCorners) {
        if ((corner == 0 && i == 0 && j == 0) || // Top-left
            (corner == 1 && i == 0 &&
             j == config.resolutionX - 1) || // Top-right
            (corner == 2 && i == config.resolutionY - 1 &&
             j == 0) || // Bottom-left
            (corner == 3 && i == config.resolutionY - 1 &&
             j == config.resolutionX - 1)) { // Bottom-right
          isFixed = true;
          break;
        }
      }

      // Create particle and add to array
      Particle *p =
          new Particle(glm::vec3(x, y, z), config.particleMass, isFixed);
      particles.push_back(p);
    }
  }
}

void ClothSystem::createSprings(SpringType type) {
  float dx = config.width / (config.resolutionX - 1);
  float dz = config.height / (config.resolutionY - 1);

  // Use appropriate spring constants and damping constants based on spring type
  float springConstant = 0.0f;
  float dampingConstant = 0.0f;

  switch (type) {
  case SpringType::STRUCTURAL:
    springConstant = config.kStructural;
    dampingConstant = config.dStructural;
    break;
  case SpringType::SHEAR:
    springConstant = config.kShear;
    dampingConstant = config.dShear;
    break;
  case SpringType::BENDING:
    springConstant = config.kBending;
    dampingConstant = config.dBending;
    break;
  }

  // Create springs based on type
  for (int i = 0; i < config.resolutionY; ++i) {
    for (int j = 0; j < config.resolutionX; ++j) {
      int idx = getParticleIndex(i, j);

      // Structural springs (connect adjacent particles)
      if (type == SpringType::STRUCTURAL) {
        // Horizontal springs
        if (j < config.resolutionX - 1) {
          int rightIdx = getParticleIndex(i, j + 1);
          Spring *s = new Spring(particles[idx], particles[rightIdx], dx,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }

        // Vertical springs
        if (i < config.resolutionY - 1) {
          int downIdx = getParticleIndex(i + 1, j);
          Spring *s = new Spring(particles[idx], particles[downIdx], dz,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }
      }

      // Shear springs (connect diagonal particles)
      else if (type == SpringType::SHEAR) {
        // Diagonal springs
        if (i < config.resolutionY - 1 && j < config.resolutionX - 1) {
          int diagIdx = getParticleIndex(i + 1, j + 1);
          float diagLength = std::sqrt(dx * dx + dz * dz);
          Spring *s = new Spring(particles[idx], particles[diagIdx], diagLength,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }

        // Other diagonal springs
        if (i < config.resolutionY - 1 && j > 0) {
          int diagIdx = getParticleIndex(i + 1, j - 1);
          float diagLength = std::sqrt(dx * dx + dz * dz);
          Spring *s = new Spring(particles[idx], particles[diagIdx], diagLength,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }
      }

      // Bending springs (connect particles two apart)
      else if (type == SpringType::BENDING) {
        // Horizontal bending
        if (j < config.resolutionX - 2) {
          int bendIdx = getParticleIndex(i, j + 2);
          Spring *s = new Spring(particles[idx], particles[bendIdx], 2.0f * dx,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }

        // Vertical bending
        if (i < config.resolutionY - 2) {
          int bendIdx = getParticleIndex(i + 2, j);
          Spring *s = new Spring(particles[idx], particles[bendIdx], 2.0f * dz,
                                 springConstant, dampingConstant);
          springs.push_back(s);
        }
      }
    }
  }
}

void ClothSystem::updateNormals() {
  // Reset all normals
  for (int i = 0; i < normals.size(); ++i) {
    normals[i] = glm::vec3(0.0f);
  }

  // Calculate face normals and add to vertex normals
  for (const auto &tri : triangles) {
    // Get vertices
    const glm::vec3 &p0 = vertices[tri[0]];
    const glm::vec3 &p1 = vertices[tri[1]];
    const glm::vec3 &p2 = vertices[tri[2]];

    // Calculate face normal
    glm::vec3 e1 = p1 - p0;
    glm::vec3 e2 = p2 - p0;
    glm::vec3 normal = glm::normalize(glm::cross(e1, e2));

    // Add to vertex normals
    normals[tri[0]] += normal;
    normals[tri[1]] += normal;
    normals[tri[2]] += normal;
  }

  // Normalize all normals
  for (auto &normal : normals) {
    normal = glm::normalize(normal);
  }
}

void ClothSystem::initializeRendering(OpenGL::Rasterizer &r) {
  // Create object for rendering
  object = r.createObject();

  // Initialize vertices and normals
  vertices.resize(particles.size());
  normals.resize(particles.size(), glm::vec3(0.0f, 1.0f, 0.0f));

  // Update vertices from particle positions
  for (int i = 0; i < particles.size(); ++i) {
    vertices[i] = particles[i]->position;
  }

  // Create vertex and normal buffers
  vertexBuf =
      r.createVertexAttribs(object, 0, vertices.size(), vertices.data());
  normalBuf = r.createVertexAttribs(object, 1, normals.size(), normals.data());

  // Create triangles
  triangles.clear();
  for (int i = 0; i < config.resolutionY - 1; ++i) {
    for (int j = 0; j < config.resolutionX - 1; ++j) {
      int idx = getParticleIndex(i, j);
      int rightIdx = getParticleIndex(i, j + 1);
      int downIdx = getParticleIndex(i + 1, j);
      int diagIdx = getParticleIndex(i + 1, j + 1);

      // First triangle
      triangles.push_back(glm::ivec3(idx, rightIdx, diagIdx));

      // Second triangle
      triangles.push_back(glm::ivec3(idx, diagIdx, downIdx));
    }
  }

  // Create triangle indices
  r.createTriangleIndices(object, triangles.size(), triangles.data());

  // Create edges
  edges.clear();

  // Horizontal edges
  for (int i = 0; i < config.resolutionY; ++i) {
    for (int j = 0; j < config.resolutionX - 1; ++j) {
      int idx = getParticleIndex(i, j);
      int rightIdx = getParticleIndex(i, j + 1);
      edges.push_back(glm::ivec2(idx, rightIdx));
    }
  }

  // Vertical edges
  for (int i = 0; i < config.resolutionY - 1; ++i) {
    for (int j = 0; j < config.resolutionX; ++j) {
      int idx = getParticleIndex(i, j);
      int downIdx = getParticleIndex(i + 1, j);
      edges.push_back(glm::ivec2(idx, downIdx));
    }
  }

  // Create edge indices
  r.createEdgeIndices(object, edges.size(), edges.data());
}

void ClothSystem::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                       const Camera &camera) {
  // Initialize rendering data if not already done
  if (object.tri_vao == 0) {
    initializeRendering(r);
  }

  // Update vertex positions from particles
  for (int i = 0; i < particles.size(); ++i) {
    vertices[i] = particles[i]->position;
  }

  // Update vertex buffers
  r.updateVertexAttribs(vertexBuf, vertices.size(), vertices.data());

  // Update normals
  updateNormals();
  r.updateVertexAttribs(normalBuf, normals.size(), normals.data());

  // Draw triangles
  r.setupFilledFaces();
  glm::vec3 blue(0.2f, 0.4f, 0.8f);
  glm::vec3 white(1.0f, 1.0f, 1.0f);
  r.setUniform(program, "ambientColor", 0.2f * white);
  r.setUniform(program, "extdiffuseColor", 0.8f * blue);
  r.setUniform(program, "intdiffuseColor", 0.4f * blue);
  r.setUniform(program, "specularColor", 0.6f * white);
  r.setUniform(program, "phongExponent", 20.0f);
  r.drawTriangles(object);

  // Draw edges
  r.setupWireFrame();
  glm::vec3 black(0.0f, 0.0f, 0.0f);
  r.setUniform(program, "ambientColor", black);
  r.setUniform(program, "extdiffuseColor", black);
  r.setUniform(program, "intdiffuseColor", black);
  r.setUniform(program, "specularColor", black);
  r.setUniform(program, "phongExponent", 0.0f);
  r.drawEdges(object);
}

} // namespace COL781