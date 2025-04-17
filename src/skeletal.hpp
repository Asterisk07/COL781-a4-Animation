
#ifndef SKELETAL_HPP // "If SKELETAL_HPP has NOT been defined yet..."
#define SKELETAL_HPP // "Then define SKELETAL_HPP"

#include "camera.hpp"
#include "hw.hpp"

#include <iostream>
#include <vector>

namespace GL = COL781::OpenGL;

namespace Skeleton {

// Called once at startup
void initialize(GL::Rasterizer &r);

// Called every frame with time `t` (in seconds)
void update(float t);

// Called every frame to draw
void draw(GL::Rasterizer &r, GL::ShaderProgram &program);

} // namespace Skeleton

#endif
