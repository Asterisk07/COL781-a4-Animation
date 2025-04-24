void ClothSystem::initializeRendering(OpenGL::Rasterizer &r) {
  // Create vertex and normal buffers
  vertexBuf = r.createAttribBuf(vertices);
  normalBuf = r.createAttribBuf(normals);

  // Create object and set attributes
  object = r.createObject();
  r.setVertexAttribs(object, vertexBuf);
  r.setVertexAttribs(object, normalBuf, 1);
  r.setTriangleIndices(object, triangles);
}

void ClothSystem::draw(OpenGL::Rasterizer &r, OpenGL::ShaderProgram &program,
                       const Camera &camera) {
  // Initialize rendering if not done yet
  if (object == 0) {
    initializeRendering(r);
  }

  // Update vertex and normal buffers
  r.updateAttribBuf(vertexBuf, vertices);
  r.updateAttribBuf(normalBuf, normals);

  // Draw cloth
  r.setUniform(program, "objectColor",
               glm::vec3(0.7f, 0.7f, 0.9f)); // Cloth color
  r.drawObject(object);

  // Draw all obstacles
  for (auto obstacle : obstacles) {
    obstacle->draw(r, program, camera);
  }
}
