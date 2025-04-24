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
