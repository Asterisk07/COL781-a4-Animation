#ifndef SPRING_HPP
#define SPRING_HPP

#include "camera.hpp"
#include "hw.hpp"
#include <memory>
#include <vector>

namespace COL781 {

// Forward declarations
class Particle;
class Spring;
class ClothSystem;
class Obstacle;

// Class to represent a particle in the system
class Particle {
public:
  glm::vec3 position; // Current position
  glm::vec3 velocity; // Current velocity
  glm::vec3 force;    // Accumulated force
  float mass;         // Particle mass
  bool isFixed;       // Whether particle is fixed in space

  Particle(const glm::vec3 &pos, float m, bool fixed = false);
  void applyForce(const glm::vec3 &f);
  void update(float dt); // Update position using semi-implicit Euler
  void reset(const glm::vec3 &pos);
};

// Class to represent a spring connecting two particles
class Spring {
public:
  Particle *p1;     // First particle
  Particle *p2;     // Second particle
  float restLength; // Rest length of the spring
  float ks;         // Spring constant
  float kd;         // Damping constant

  Spring(Particle *particle1, Particle *particle2, float rest_length,
         float spring_const, float damping_const);
  void applyForce(); // Calculate and apply spring forces to particles
};

// Enum to define types of springs
enum class SpringType { STRUCTURAL, SHEAR, BENDING };

// Class to represent a cloth as a mass-spring system
class ClothSystem {
public:
  // Configuration parameters for the cloth
  struct ClothConfig {
    float width;        // Physical width of cloth (meters)
    float height;       // Physical height of cloth (meters)
    int resolutionX;    // Number of particles in X direction
    int resolutionY;    // Number of particles in Y direction
    float particleMass; // Mass of each particle (kg)
    float kStructural;  // Spring constant for structural springs
    float kShear;       // Spring constant for shear springs
    float kBending;     // Spring constant for bending springs
    float dStructural;  // Damping constant for structural springs
    float dShear;       // Damping constant for shear springs
    float dBending;     // Damping constant for bending springs
    float gravity;      // Gravity acceleration (m/s^2)
    float restitution;  // Coefficient of restitution (epsilon) for collisions
    float friction;     // Coefficient of friction (mu) for collisions
    std::vector<int>
        fixedCorners; // Indices of fixed corners (0=top-left, 1=top-right,
                      // 2=bottom-left, 3=bottom-right)
  };

  ClothSystem(const ClothConfig &config);
  ~ClothSystem();

  void initialize();     // Initialize the cloth system
  void reset();          // Reset the cloth to its initial state
  void update(float dt); // Update the system by one time step
  void draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
            const Camera &camera); // Draw the cloth

  // Add obstacle to the system
  void addObstacle(Obstacle *obstacle);

  // Handle collisions with all obstacles
  void handleCollisions();

  // Getters
  int getNumParticles() const { return particles.size(); }
  int getNumSprings() const { return springs.size(); }
  const std::vector<Particle *> &getParticles() const { return particles; }
  const std::vector<Spring *> &getSprings() const { return springs; }
  const ClothConfig &getConfig() const { return config; }

private:
  ClothConfig config;
  std::vector<Particle *> particles;
  std::vector<Spring *> springs;
  std::vector<Obstacle *> obstacles;

  // Helper variables for rendering
  OpenGL::Object object;
  OpenGL::AttribBuf vertexBuf, normalBuf;
  std::vector<glm::vec3> vertices;
  std::vector<glm::vec3> normals;
  std::vector<glm::ivec3> triangles;
  std::vector<glm::ivec2> edges;

  // Helper functions
  int getParticleIndex(int i, int j) const;
  void createParticles();
  void createSprings(SpringType type);
  void updateNormals();
  void initializeRendering(OpenGL::Rasterizer &r);

  // Handle collision between a particle and an obstacle
  void handleCollision(Particle *particle, Obstacle *obstacle);
};

} // namespace COL781

#endif // SPRING_HPP