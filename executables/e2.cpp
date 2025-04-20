#include "spring.hpp"
#include <iostream>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;

GL::Rasterizer r;
GL::ShaderProgram program;
CameraControl camCtl;

// Global variables for simulation
ClothSystem *cloth = nullptr;
float simulationTime = 0.0f;
float timeStep = 0.0005f;
bool paused = false;

void initializeScene() {
  // Setup the cloth simulation parameters
  ClothSystem::ClothConfig config;
  config.width = 1.0f;     // 1m width
  config.height = 1.0f;    // 1m height
  config.resolutionX = 21; // 21x21 grid (20 segments)
  config.resolutionY = 21;
  config.particleMass = 0.01f; // 10g per particle

  // Spring constants - follow the relationship: kBending << kShear <
  // kStructural
  config.kStructural = 500.0f; // Strong structural springs
  config.kShear = 100.0f;      // Medium shear springs
  config.kBending = 10.0f;     // Weak bending springs

  // Damping constants - much less than the corresponding spring constants
  config.dStructural = 5.0f;
  config.dShear = 1.0f;
  config.dBending = 0.1f;

  config.gravity = 9.8f;        // Standard gravity
  config.fixedCorners = {0, 1}; // Fix top-left and top-right corners

  // Create cloth system
  cloth = new ClothSystem(config);
}

void updateScene(float t) {
  if (!paused) {
    // Perform multiple mini-steps for better stability
    for (int i = 0; i < 10; ++i) {
      cloth->update(timeStep);
    }
    simulationTime += timeStep * 10.0f;
  }
}

void handleInput() {
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    if (event.type == SDL_KEYDOWN) {
      switch (event.key.keysym.sym) {
      case SDLK_SPACE:
        // Toggle pause
        paused = !paused;
        break;
      case SDLK_r:
        // Reset simulation
        cloth->reset();
        simulationTime = 0.0f;
        break;
      case SDLK_ESCAPE:
        SDL_Quit();
        exit(0);
        break;
      }
    }
  }
}

int main() {
  int width = 800, height = 600;
  if (!r.initialize("Cloth Simulation", width, height)) {
    return EXIT_FAILURE;
  }

  // Initialize camera
  camCtl.initialize(width, height);
  camCtl.camera.setCameraView(vec3(0.5f, 0.5f, 2.0f), vec3(0.5f, 0.0f, 0.5f),
                              vec3(0.0f, 1.0f, 0.0f));

  // Create shaders
  program = r.createShaderProgram(r.vsBlinnPhong(), r.fsBlinnPhong());

  // Initialize scene
  initializeScene();

  // Main loop
  while (!r.shouldQuit()) {
    // Handle user input
    handleInput();

    // Update simulation
    float t = SDL_GetTicks64() * 1e-3;
    updateScene(t);

    // Update camera
    camCtl.update();
    Camera &camera = camCtl.camera;

    // Clear screen
    r.clear(vec4(0.4f, 0.4f, 0.4f, 1.0f));
    r.enableDepthTest();
    r.useShaderProgram(program);

    // Set common shader uniforms
    r.setUniform(program, "model", mat4(1.0f));
    r.setUniform(program, "view", camera.getViewMatrix());
    r.setUniform(program, "projection", camera.getProjectionMatrix());
    r.setUniform(program, "lightPos", camera.position);
    r.setUniform(program, "viewPos", camera.position);
    r.setUniform(program, "lightColor", vec3(1.0f, 1.0f, 1.0f));

    // Draw the cloth
    cloth->draw(r, program, camera);

    // Display simulation info
    std::cout << "\rSimulation time: " << simulationTime << " seconds"
              << (paused ? " (PAUSED)" : "") << std::flush;

    // Show the rendered frame
    r.show();
  }

  // Clean up
  delete cloth;

  return 0;
}