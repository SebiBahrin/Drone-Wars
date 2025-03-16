# Shooter Drone Game

## Overview
This project is a real-time 3D simulation featuring a controllable drone, enemy drones, obstacles, projectiles, and an interactive terrain. The simulation demonstrates collision detection, physics-based movement, and simple enemy behavior within an OpenGL environment.

## Features
### **Drone Control**
- The player can move the drone in a 3D environment using keyboard controls.
- The drone’s orientation and movement are updated in real time.

### **Enemy Drones**
- Enemy drones are randomly spawned on the map.
- They navigate toward designated corners of the terrain.
- When in range, they fire projectiles at the player.

### **Obstacle Generation**
- The environment includes randomly placed obstacles like trees and buildings.
- Objects use different meshes (e.g., cylinders, cones, cubes, and parallelepipeds) with custom dimensions.

### **Projectile System**
- Both the player and enemy drones can shoot projectiles.
- Projectiles are rendered as small spheres.
- They are removed when exiting world boundaries or colliding with an object.

### **Collision Detection**
- Uses Axis-Aligned Bounding Box (AABB) and sphere collision detection.
- Handles interactions between drones, obstacles, and projectiles.

### **Camera Modes**
- Supports first-person and third-person camera views that follow the drone.

## **Setup and Compilation**
### **Setup**
1. **Extract the Framework**
   - Download and extract the framework from `gfx-framework.zip`
   (https://github.com/UPB-Graphics/gfx-framework)
2. **Install CMake**
   - Use your package manager to install CMake if it's not already installed.

### **Generate the Build Files**
Open a terminal and navigate to the root folder of the project (the folder containing `CMakeLists.txt`). Then, run:
```bash
mkdir build
cd build
cmake ..
```

### **Build the Project**
#### **Windows**
Run:
```bash
cmake --build .
```
Or open the generated `.sln` file in Visual Studio and press **Ctrl+Shift+B** to build.

#### **Linux and macOS**
Run:
```bash
cmake --build .
# or
make
```

### **Modify the Main Entry Point**
Navigate to the `src` folder, open `main.c`, and locate:
```cpp
// Create a new 3D world and start running it
World *world = new gfxc::SimpleScene();
```
Replace it with:
```cpp
// Create a new 3D world and start running it
World *world = new m1::Lab7();
```

### **Replace the Lab Folder**
- Open the `lab_m1` folder and replace the `lab7` folder with the version included in this project.

### **Run the Game**
From the root folder, navigate to the build directory and run:
```bash
cd build
make
cd bin/Debug/
./GFXFramework
```

## **Controls**
### **Movement**
- `W / S` - Move forward/backward
- `A / D` - Move left/right
- `UP / DOWN` - Move upward/downward

### **Rotation**
- `LEFT / RIGHT` - Rotate the drone around the Y-axis

### **Actions**
- `SPACE` - Shoot a projectile
- `C` - Toggle between first-person and third-person camera views
- `Q` - Affects collision behavior with the ground (hold to modify behavior)

---

Enjoy playing my Shooter Drone Game!
