# COL781 Assignment 4: Keyframing and Simulation

Demonstration video : https://www.youtube.com/watch?v=C-d_979t52o

Some snapshots:
<img width="600" alt="cloth" src="https://github.com/user-attachments/assets/b6193e41-2def-4746-888e-50475615a11b" />
<img width="600" alt="collision 1" src="https://github.com/user-attachments/assets/627d2c41-53b5-4627-8721-0e7da4076ae1" />
<img width="600" alt="collision 2" src="https://github.com/user-attachments/assets/effc491a-f43e-4123-ae4d-9fbafa8f8f86" />

Make sure that [glm](https://github.com/g-truc/glm) and [SDL2](https://www.libsdl.org/) are installed. Ideally, these should be installed by your package manager rather than manually (at least, if you are on Linux or Mac).

Then compile the code using the standard CMake procedure:

- The first time, run `cmake -B build` from the project root to create a `build/` directory and initialize a build system there.
- Then, every time you want to compile the code, run `cmake --build build` (again from the project root). Then the example programs will be created under `build/`.
