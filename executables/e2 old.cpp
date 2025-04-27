#include "spring.hpp"
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp> // JSON for Modern C++
#include <string>

using namespace COL781;
namespace GL = COL781::OpenGL;
using namespace glm;
using json = nlohmann::json;
#include <stdexcept>

// Config file path
const std::string CONFIG_FILE = "config/config_current.json";

GL::Rasterizer r;
GL::ShaderProgram program;
CameraControl camCtl;

// Global variables for simulation
ClothSystem *cloth = nullptr;
float simulationTime = 0.0f;
float timeStep = 10.0f; // Default value, will be overridden
bool paused = false;

// Load cloth configuration from a JSON file
ClothSystem::ClothConfig loadClothConfig(const std::string &configFile) {
  ClothSystem::ClothConfig config;

  try {
    // Read JSON file
    std::ifstream file(configFile);
    if (!file.is_open()) {
      std::cerr << "Error: Could not open config file: " << configFile
                << std::endl;
      std::exit(EXIT_FAILURE);
      // Return default config if file can't be opened
      //   return config;
    }

    json j;
    file >> j;

    try {
      // Load timestep
      timeStep = j.at("timeStep").get<float>();

      // Load cloth dimensions
      config.width = j.at("width").get<float>();
      config.height = j.at("height").get<float>();
      config.resolutionX = j.at("resolutionX").get<int>();
      config.resolutionY = j.at("resolutionY").get<int>();
      config.particleMass = j.at("particleMass").get<float>();

      // Load spring constants
      config.kStructural = j.at("kStructural").get<float>();
      config.kShear = j.at("kShear").get<float>();
      config.kBending = j.at("kBending").get<float>();

      // Load damping constants
      config.dStructural = j.at("dStructural").get<float>();
      config.dShear = j.at("dShear").get<float>();
      config.dBending = j.at("dBending").get<float>();

      // Load gravity
      config.gravity = j.at("gravity").get<float>();

      // Load fixed corners
      config.fixedCorners = j.at("fixedCorners").get<std::vector<int>>();

    } catch (const nlohmann::json::exception &e) {
      std::cerr << "FATAL CONFIG ERROR: Missing required field - " << e.what()
                << std::endl;
      std::exit(EXIT_FAILURE); // Crash immediately
    }
    std::cerr << "Loaded file " << std::endl;
    {
    }
    {
    }
    // Load fixed corners
    if (j.contains("fixedCorners")) {
      config.fixedCorners.clear();
      for (auto &corner : j["fixedCorners"]) {
        config.fixedCorners.push_back(corner);
      }
    } else {
      config.fixedCorners = {0, 1}; // Default to top-left and top-right corners
    }

    std::cout << "Successfully loaded config from: " << configFile << std::endl;
  } catch (json::exception &e) {
    std::cerr << "JSON parsing error: " << e.what() << std::endl;
    std::cerr << "Using default configuration." << std::endl;
  } catch (std::exception &e) {
    std::cerr << "Error loading config: " << e.what() << std::endl;
    std::cerr << "Using default configuration." << std::endl;
  }

  return config;
}

void initializeScene() {
  // Load the cloth simulation parameters from config file
  ClothSystem::ClothConfig config = loadClothConfig(CONFIG_FILE);

  // Create cloth system with loaded configuration
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