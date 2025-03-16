#include "lab_m1/lab7/lab7.h"

#include <vector>
#include <string>
#include <iostream>

using namespace std;
using namespace m1;


/*
 *  To find out more about FrameStart, Update, FrameEnd
 *  and the order in which they are called, see world.cpp.
 */


Lab7::Lab7()
{
}


Lab7::~Lab7()
{
}

Mesh* Lab7::createParallelepiped(const std::string& name, glm::vec3 dimensions, const glm::vec3& color)
{
    float dx = dimensions.x / 2.0f;
    float dy = dimensions.y / 2.0f;
    float dz = dimensions.z / 2.0f;

    Mesh* parallelepipedMesh = new Mesh(name);

    std::vector<VertexFormat> vertices = {
        VertexFormat(glm::vec3(-dx, -dy,  dz), color),
        VertexFormat(glm::vec3( dx, -dy,  dz), color),
        VertexFormat(glm::vec3( dx,  dy,  dz), color),
        VertexFormat(glm::vec3(-dx,  dy,  dz), color),

        VertexFormat(glm::vec3(-dx, -dy, -dz), color),
        VertexFormat(glm::vec3( dx, -dy, -dz), color),
        VertexFormat(glm::vec3( dx,  dy, -dz), color),
        VertexFormat(glm::vec3(-dx,  dy, -dz), color)
    };

    std::vector<unsigned int> indices = {
        0, 1, 2, 0, 2, 3,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        2, 3, 7, 2, 7, 6,
        0, 3, 7, 0, 7, 4,
        1, 2, 6, 1, 6, 5
    };

    parallelepipedMesh->InitFromData(vertices, indices);
    return parallelepipedMesh;
}

Mesh* Lab7::createCube(const std::string& name, float size, const glm::vec3& color)
{
    float halfSize = size / 2.0f;

    Mesh* cubeMesh = new Mesh(name);

    std::vector<VertexFormat> vertices = {
        VertexFormat(glm::vec3(-halfSize, -halfSize,  halfSize), color),
        VertexFormat(glm::vec3( halfSize, -halfSize,  halfSize), color),
        VertexFormat(glm::vec3( halfSize,  halfSize,  halfSize), color),
        VertexFormat(glm::vec3(-halfSize,  halfSize,  halfSize), color),

        VertexFormat(glm::vec3(-halfSize, -halfSize, -halfSize), color),
        VertexFormat(glm::vec3( halfSize, -halfSize, -halfSize), color),
        VertexFormat(glm::vec3( halfSize,  halfSize, -halfSize), color),
        VertexFormat(glm::vec3(-halfSize,  halfSize, -halfSize), color)
    };

    std::vector<unsigned int> indices = {
        0, 1, 2, 0, 2, 3,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        2, 3, 7, 2, 7, 6,
        0, 3, 7, 0, 7, 4,
        1, 2, 6, 1, 6, 5
    };

    cubeMesh->InitFromData(vertices, indices);
    return cubeMesh;
}

void Lab7::createDrone(const glm::vec3& initialPosition)
{
    dronePosition = initialPosition;
    glm::vec3 baseColor = glm::vec3(0.5f, 0.5f, 0.5f);

    Mesh* body1 = createParallelepiped("body1", glm::vec3(2.0f, 0.2f, 0.2f), baseColor);
    meshes[body1->GetMeshID()] = body1;

    Mesh* body2 = createParallelepiped("body2", glm::vec3(2.0f, 0.2f, 0.2f), baseColor);
    meshes[body2->GetMeshID()] = body2;

    Mesh* cube1 = createCube("cube1", 0.2f, baseColor); 
    meshes[cube1->GetMeshID()] = cube1;

    Mesh* cube2 = createCube("cube2", 0.2f, baseColor); 
    meshes[cube2->GetMeshID()] = cube2;

    Mesh* cube3 = createCube("cube3", 0.2f, baseColor); 
    meshes[cube3->GetMeshID()] = cube3;

    Mesh* cube4 = createCube("cube4", 0.2f, baseColor); 
    meshes[cube4->GetMeshID()] = cube4;

    Mesh* propeller1 = createParallelepiped("propeller1", glm::vec3(0.2f, 0.05f, 0.7f), glm::vec3(0.0f, 0.0f, 0.0f));
    meshes[propeller1->GetMeshID()] = propeller1;

    Mesh* propeller2 = createParallelepiped("propeller2", glm::vec3(0.2f, 0.05f, 0.7f), glm::vec3(0.0f, 0.0f, 0.0f));
    meshes[propeller2->GetMeshID()] = propeller2;

    Mesh* propeller3 = createParallelepiped("propeller3", glm::vec3(0.2f, 0.05f, 0.7f), glm::vec3(0.0f, 0.0f, 0.0f));
    meshes[propeller3->GetMeshID()] = propeller3;

    Mesh* propeller4 = createParallelepiped("propeller4", glm::vec3(0.2f, 0.05f, 0.7f), glm::vec3(0.0f, 0.0f, 0.0f));
    meshes[propeller4->GetMeshID()] = propeller4;
}

void Lab7::updateDrone()
{
    glm::mat4 droneModelMatrix = glm::mat4(1);
    droneModelMatrix = glm::translate(droneModelMatrix, dronePosition);
    droneModelMatrix = glm::rotate(droneModelMatrix, glm::radians(droneRotationY), glm::vec3(0, 1, 0));

    RenderMesh(meshes["arrow"], shaders["Simple"], glm::scale(droneModelMatrix, glm::vec3(1, 0.1f, 0.1f)));
    RenderMesh(meshes["arrow"], shaders["Simple"], glm::scale(glm::rotate(droneModelMatrix, glm::radians(90.0f), glm::vec3(0, 0, 1)), glm::vec3(0.1f, 1, 0.1f)));
    RenderMesh(meshes["arrow"], shaders["Simple"], glm::scale(glm::rotate(droneModelMatrix, glm::radians(90.0f), glm::vec3(1, 0, 0)), glm::vec3(0.1f, 0.1f, 1)));

    glm::vec3 baseColor = glm::vec3(0.5f, 0.5f, 0.5f);

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["body1"], shaders["LabShader"], modelMatrix, baseColor);
    }

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(-45.0f), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["body2"], shaders["LabShader"], modelMatrix, baseColor);
    }

    float cubeSize = 0.2f;
    float propellerHeight = 0.05f;

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, 0.2f, 0.9f));
        RenderSimpleMesh(meshes["cube1"], shaders["LabShader"], modelMatrix, baseColor);

        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
        modelMatrix = glm::rotate(modelMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["propeller1"], shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.0f, 0.0f));
    }

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, 0.2f, -0.9f));
        RenderSimpleMesh(meshes["cube2"], shaders["LabShader"], modelMatrix, baseColor);

        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
        modelMatrix = glm::rotate(modelMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["propeller2"], shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.0f, 0.0f));
    }

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
        modelMatrix = glm::translate(modelMatrix, glm::vec3(-0.9f, 0.2f, 0.0f));
        RenderSimpleMesh(meshes["cube3"], shaders["LabShader"], modelMatrix, baseColor);

        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
        modelMatrix = glm::rotate(modelMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["propeller3"], shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.0f, 0.0f));
    }

    {
        glm::mat4 modelMatrix = droneModelMatrix;
        modelMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.9f, 0.2f, 0.0f));
        RenderSimpleMesh(meshes["cube4"], shaders["LabShader"], modelMatrix, baseColor);

        modelMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
        modelMatrix = glm::rotate(modelMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
        RenderSimpleMesh(meshes["propeller4"], shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.0f, 0.0f));
    }
}

Mesh* Lab7::createTerrain(const std::string& name, int m, int n, const glm::vec3& color, const glm::vec3& gridOffset)
{
    std::vector<VertexFormat> vertices;
    std::vector<unsigned int> indices;

    float cellWidth = 1.0f;
    float cellHeight = 1.0f;

    // Generăm vertecșii cu un offset
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            glm::vec3 position = glm::vec3(j * cellWidth, 0.0f, i * cellHeight) + gridOffset;
            vertices.emplace_back(position, color);
        }
    }

    // Generăm indicii pentru triunghiuri
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int topLeft = i * (n + 1) + j;
            int topRight = topLeft + 1;
            int bottomLeft = topLeft + (n + 1);
            int bottomRight = bottomLeft + 1;

            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    // Creăm mesh-ul
    Mesh* terrainMesh = new Mesh(name);
    terrainMesh->InitFromData(vertices, indices);
    return terrainMesh;
}

void Lab7::createTree(const glm::vec3& position)
{
    // Trunchiul copacului
    Mesh* trunk = createCylinder("trunk", 20, 1.0f, 0.2f, glm::vec3(0.55f, 0.27f, 0.07f));
    glm::mat4 modelMatrix = glm::mat4(1);
    modelMatrix = glm::translate(modelMatrix, position);
    RenderSimpleMesh(trunk, shaders["LabShader"], modelMatrix, glm::vec3(0.55f, 0.27f, 0.07f)); 

    // Conurile frunzelor
    Mesh* leaves = createCone("leaves", 20, 1.0f, 0.6f, glm::vec3(0.0f, 0.5f, 0.0f)); 
    modelMatrix = glm::translate(glm::mat4(1), position + glm::vec3(0.0f, 0.7f, 0.0f));
    modelMatrix = glm::scale(modelMatrix, glm::vec3(1.4f));
    RenderSimpleMesh(leaves, shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.5f, 0.0f));
    modelMatrix = glm::translate(glm::mat4(1), position + glm::vec3(0.0f, 1.4f, 0.0f));
    RenderSimpleMesh(leaves, shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.5f, 0.0f)); 
}

Mesh* Lab7::createCone(const std::string& name, int segments, float height, float radius, const glm::vec3& color)
{
    std::vector<VertexFormat> vertices;
    std::vector<unsigned int> indices;

    vertices.emplace_back(glm::vec3(0.0f, height, 0.0f), color); 

    for (int i = 0; i < segments; ++i) {
        float angle = 2.0f * M_PI * i / segments;
        float x = radius * cos(angle);
        float z = radius * sin(angle);
        vertices.emplace_back(glm::vec3(x, 0.0f, z), color);
    }

    for (int i = 1; i <= segments; ++i) {
        indices.push_back(0);
        indices.push_back(i); 
        indices.push_back(i % segments + 1); 
    }

    vertices.emplace_back(glm::vec3(0.0f, 0.0f, 0.0f), color); 
    int centerIndex = static_cast<int>(vertices.size()) - 1;

    for (int i = 1; i <= segments; ++i) {
        indices.push_back(centerIndex);
        indices.push_back(i % segments + 1);
        indices.push_back(i); 
    }

    Mesh* coneMesh = new Mesh(name);
    coneMesh->InitFromData(vertices, indices);
    return coneMesh;
}

Mesh* Lab7::createCylinder(const std::string& name, int segments, float height, float radius, const glm::vec3& color)
{
    std::vector<VertexFormat> vertices;
    std::vector<unsigned int> indices;

    for (int i = 0; i < segments; ++i) {
        float angle = 2.0f * M_PI * i / segments;
        float x = radius * cos(angle);
        float z = radius * sin(angle);
        vertices.emplace_back(glm::vec3(x, 0.0f, z), color);
    }

    for (int i = 0; i < segments; ++i) {
        float angle = 2.0f * M_PI * i / segments;
        float x = radius * cos(angle);
        float z = radius * sin(angle);
        vertices.emplace_back(glm::vec3(x, height, z), color);
    }

    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        indices.push_back(i);
        indices.push_back(next);
        indices.push_back(i + segments);
        indices.push_back(next);
        indices.push_back(next + segments);
        indices.push_back(i + segments);
    }

    vertices.emplace_back(glm::vec3(0.0f, 0.0f, 0.0f), color);
    int bottomCenter = static_cast<int>(vertices.size()) - 1;

    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        indices.push_back(bottomCenter);
        indices.push_back(next);
        indices.push_back(i);
    }

    vertices.emplace_back(glm::vec3(0.0f, height, 0.0f), color);
    int topCenter = static_cast<int>(vertices.size()) - 1;

    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        indices.push_back(topCenter);
        indices.push_back(i + segments);
        indices.push_back(next + segments);
    }

    Mesh* cylinderMesh = new Mesh(name);
    cylinderMesh->InitFromData(vertices, indices);
    return cylinderMesh;
}

void Lab7::createBuilding(const glm::vec3& position, const glm::vec3& dimensions, const glm::vec3& color)
{
    Mesh* building = createParallelepiped("building", dimensions, color);
    glm::mat4 modelMatrix = glm::mat4(1);
    modelMatrix = glm::translate(modelMatrix, position);
    RenderSimpleMesh(building, shaders["LabShader"], modelMatrix, color);
}

// Generates a random position on the XZ plane at a fixed Y level.
glm::vec3 Lab7::generateRandomPosition(float xmin, float xmax, float zmin, float zmax, float y)
{
    // Generate a random x coordinate within [xmin, xmax]
    float x = xmin + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (xmax - xmin)));
    // Generate a random z coordinate within [zmin, zmax]
    float z = zmin + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (zmax - zmin)));
    return glm::vec3(x, y, z);
}

// Generates a random position in 3D space for a drone, with random y between ymin and ymax.
glm::vec3 Lab7::generateRandomPositionDrone(float xmin, float xmax, float zmin, float zmax, float ymin, float ymax)
{
    // Generate a random x coordinate within [xmin, xmax]
    float x = xmin + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (xmax - xmin)));
    // Generate a random z coordinate within [zmin, zmax]
    float z = zmin + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (zmax - zmin)));
    // Generate a random y coordinate within [ymin, ymax]
    float y = ymin + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (ymax - ymin)));
    return glm::vec3(x, y, z);
}

// Checks if the given position is valid compared to existing positions,
// ensuring that the new position is at least minDistance away from all existing positions.
bool Lab7::isPositionValid(const glm::vec3& position, const std::vector<glm::vec3>& existingPositions, float minDistance)
{
    for (const auto& existing : existingPositions) {
        // If the distance between the new position and an existing position is less than minDistance, it's invalid.
        if (glm::distance(position, existing) < minDistance) {
            return false;
        }
    }
    return true;
}

// Generates decorative obstacles (trees or buildings) with given parameters.
void Lab7::generateDecor(int numObstacles, float xmin, float xmax, float zmin, float zmax, float minDistance)
{
    // Clear any previous obstacle data.
    obstaclePositions.clear();
    obstacleTypes.clear();
    obstacleDimensions.clear();

    // Loop to generate the specified number of obstacles.
    for (int i = 0; i < numObstacles; ++i) {
        glm::vec3 position;
        // Keep generating a random position until it is valid (i.e., does not collide with existing obstacles).
        do {
            position = generateRandomPosition(xmin, xmax, zmin, zmax, 0.0f);
        } while (!isPositionValid(position, obstaclePositions, minDistance));

        // Randomly decide whether to create a tree or a building.
        if (rand() % 2 == 0) {
            // Create a tree obstacle.
            obstacleTypes.push_back("tree");
            obstaclePositions.push_back(position);
            // Randomly set dimensions for the tree (trunk height, trunk radius, foliage height, foliage radius).
            float trunkHeight = 1.0f + static_cast<float>(rand()) / RAND_MAX * 2.0f;
            float trunkRadius = 0.1f + static_cast<float>(rand()) / RAND_MAX * 0.2f;
            float foliageHeight = 1.0f + static_cast<float>(rand()) / RAND_MAX * 1.0f;
            float foliageRadius = 0.5f + static_cast<float>(rand()) / RAND_MAX * 1.0f;
            // Store only trunkHeight, trunkRadius, and foliageHeight in dimensions.
            obstacleDimensions.push_back(glm::vec3(trunkHeight, trunkRadius, foliageHeight));
        } else {
            // Create a building obstacle.
            // Generate random dimensions for the building.
            glm::vec3 dimensions = glm::vec3(2.0f + rand() % 3, 5.0f + rand() % 5, 2.0f + rand() % 3);
            glm::vec3 adjustedPosition = position;
            // Clamp the position to ensure the building fits within the boundaries.
            adjustedPosition.x = glm::clamp(adjustedPosition.x, xmin + dimensions.x / 2.0f, xmax - dimensions.x / 2.0f);
            adjustedPosition.z = glm::clamp(adjustedPosition.z, zmin + dimensions.z / 2.0f, zmax - dimensions.z / 2.0f);
            obstacleTypes.push_back("building");
            obstaclePositions.push_back(adjustedPosition);
            obstacleDimensions.push_back(dimensions);
        }
    }
}

// Renders and updates all decorative obstacles.
void Lab7::updateDecor()
{
    // Iterate through each obstacle.
    for (size_t i = 0; i < obstaclePositions.size(); ++i) {
        if (obstacleTypes[i] == "tree") {
            // For trees, retrieve dimensions for trunk and foliage.
            glm::vec3 dimensions = obstacleDimensions[i];
            float trunkHeight = dimensions.x;
            float trunkRadius = dimensions.y;
            float foliageHeight = dimensions.z;

            // Render the trunk of the tree.
            glm::mat4 modelMatrix = glm::mat4(1);
            modelMatrix = glm::translate(modelMatrix, obstaclePositions[i]);
            modelMatrix = glm::scale(modelMatrix, glm::vec3(trunkRadius, trunkHeight, trunkRadius));
            RenderSimpleMesh(createCylinder("trunk", 20, 1.0f, 1.0f, glm::vec3(0.55f, 0.27f, 0.07f)), shaders["LabShader"], modelMatrix, glm::vec3(0.55f, 0.27f, 0.07f));

            // Render the first foliage cone.
            Mesh* foliage = createCone("foliage1", 20, 1.0f, 1.0f, glm::vec3(0.0f, 0.5f, 0.0f));
            modelMatrix = glm::translate(glm::mat4(1), obstaclePositions[i] + glm::vec3(0.0f, trunkHeight, 0.0f));
            modelMatrix = glm::scale(modelMatrix, glm::vec3(foliageHeight, foliageHeight, foliageHeight));
            RenderSimpleMesh(foliage, shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.5f, 0.0f));

            // Render the second foliage cone (positioned higher for layered effect).
            Mesh* foliage2 = createCone("foliage2", 20, 1.0f, 1.0f, glm::vec3(0.0f, 0.5f, 0.0f));
            modelMatrix = glm::translate(glm::mat4(1), obstaclePositions[i] + glm::vec3(0.0f, trunkHeight + foliageHeight * 0.8f, 0.0f));
            modelMatrix = glm::scale(modelMatrix, glm::vec3(foliageHeight * 0.8f, foliageHeight * 0.8f, foliageHeight * 0.8f));
            RenderSimpleMesh(foliage2, shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 0.5f, 0.0f));
        } else if (obstacleTypes[i] == "building") {
            // For buildings, calculate the center position based on dimensions.
            glm::vec3 dimensions = obstacleDimensions[i];
            glm::vec3 buildingPosition = obstaclePositions[i] + glm::vec3(0.0f, dimensions.y / 2.0f, 0.0f);
            // Create the building mesh.
            Mesh* buildingMesh = createParallelepiped("building_" + std::to_string(i), dimensions, glm::vec3(0.57f, 0.58f, 0.59f));
            glm::mat4 modelMatrix = glm::mat4(1);
            modelMatrix = glm::translate(modelMatrix, buildingPosition);
            RenderSimpleMesh(buildingMesh, shaders["LabShader"], modelMatrix, glm::vec3(0.59f, 0.59f, 0.59f));
        }
    }
}

// Checks if two axis-aligned bounding boxes (AABB) intersect.
bool Lab7::checkAABBCollision(const glm::vec3& min1, const glm::vec3& max1, const glm::vec3& min2, const glm::vec3& max2)
{
    return (min1.x <= max2.x && max1.x >= min2.x) &&
           (min1.y <= max2.y && max1.y >= min2.y) &&
           (min1.z <= max2.z && max1.z >= min2.z);
}

// Calculates the drone's axis-aligned bounding box based on its current position.
void Lab7::getDroneAABB(glm::vec3& min, glm::vec3& max)
{
    // Half sizes for the drone approximated dimensions.
    glm::vec3 halfSize = glm::vec3(1.0f, 0.6f, 1.0f);

    // Calculate min and max points of the AABB.
    min = dronePosition - halfSize;
    max = dronePosition + halfSize;
}

// Calculates the axis-aligned bounding box (AABB) for an obstacle (tree or building).
void Lab7::getObstacleAABB(const glm::vec3& position, const glm::vec3& dimensions, glm::vec3& min, glm::vec3& max, bool isTree)
{
    if (isTree) {
        // For trees, use trunk and foliage dimensions.
        float trunkHeight = dimensions.x;    // Trunk height
        float trunkRadius = dimensions.y;      // Trunk radius
        float foliageHeight = dimensions.z;    // Foliage height

        // Expand the AABB to include the foliage.
        min = position - glm::vec3(trunkRadius * 2, 0, trunkRadius * 2); // Expand on X and Z axes
        max = position + glm::vec3(trunkRadius * 2, trunkHeight + foliageHeight * 2, trunkRadius * 2); // Include trunk and extra foliage height
    } else {
        // For buildings, calculate the AABB based on dimensions.
        min = position - dimensions / 2.0f;
        max = position + dimensions / 2.0f;

        // Adjust Y coordinates to reflect the building's height.
        min.y = position.y - dimensions.y;
        max.y = position.y + dimensions.y;
    }
}

// Checks for collisions between the drone and obstacles (or ground) and handles the consequences.
void Lab7::checkCollisions()
{
    // Calculate the drone's AABB.
    glm::vec3 droneMin, droneMax;
    getDroneAABB(droneMin, droneMax);

    bool collisionDetected = false;
    bool groundCollision = false;

    // Check collision with the ground.
    if (droneMin.y <= 0.0f) {
        // If the Q key is not held down, register the collision.
        if (!window->KeyHold(GLFW_KEY_Q)) {
            groundCollision = true;
            if (!isColliding) {
                std::cout << "Drona a atins solul!" << std::endl;
                // Decrease the drone's health.
                droneHealth -= 10;
                std::cout << "Viata ramasa a dronei: " << droneHealth << std::endl;
            }
            // Ensure the drone doesn't go below a minimum y value.
            if (dronePosition.y < 0.01f) {
                dronePosition.y = 0.01f;
            }
        }
    }

    // Check for collisions with each obstacle.
    for (size_t i = 0; i < obstaclePositions.size(); ++i) {
        glm::vec3 obstacleMin, obstacleMax;

        if (obstacleTypes[i] == "tree") {
            // Calculate AABB for tree obstacles.
            getObstacleAABB(obstaclePositions[i], obstacleDimensions[i], obstacleMin, obstacleMax, true);
        } else if (obstacleTypes[i] == "building") {
            // Calculate AABB for building obstacles.
            getObstacleAABB(obstaclePositions[i], obstacleDimensions[i], obstacleMin, obstacleMax, false);
        }

        // Check if the drone's AABB intersects with the obstacle's AABB.
        if (checkAABBCollision(droneMin, droneMax, obstacleMin, obstacleMax)) {
            collisionDetected = true;
            if (!isColliding) {
                if (obstacleTypes[i] == "tree") {
                    std::cout << "Coliziune detectată cu un copac!" << std::endl;
                } else {
                    std::cout << "Coliziune detectată cu obstacolul " << i << "!" << std::endl;
                }

                // Decrease the drone's health upon collision.
                droneHealth -= 10;
                std::cout << "Viata ramasa a dronei: " << droneHealth << std::endl;
            }
            break;
        }
    }

    // If a collision was detected, revert the drone's position to its previous state.
    if (collisionDetected) {
        dronePosition = previousDronePosition;
    }

    // Update the collision state flag.
    isColliding = collisionDetected || groundCollision;
    
    // Check if drone health has depleted.
    if (droneHealth <= 0) {
        if (counter == 0) {
            std::cout << "Game Over! Drona a fost distrusa." << std::endl;
            counter++;
        }
        droneHealth = 0;
    }
}


bool Lab7::checkProjectileCollision(const glm::vec3& projectilePosition, const glm::vec3& obstaclePosition, const glm::vec3& dimensions, bool isCone)
{
    if (!isCone) {
        // AABB collision check
        glm::vec3 minObstacle = obstaclePosition - glm::vec3(dimensions.x / 2.0f, dimensions.y / 2.0f, dimensions.z / 2.0f);
        glm::vec3 maxObstacle = obstaclePosition + glm::vec3(dimensions.x / 2.0f, dimensions.y / 2.0f, dimensions.z / 2.0f);

        return (projectilePosition.x >= minObstacle.x && projectilePosition.x <= maxObstacle.x) &&
               (projectilePosition.y >= minObstacle.y && projectilePosition.y <= maxObstacle.y) &&
               (projectilePosition.z >= minObstacle.z && projectilePosition.z <= maxObstacle.z);
    } else {
        // Cone collision check (frunze)
        glm::vec3 relativePos = projectilePosition - obstaclePosition;

        // Check if the projectile is within the height of the cone
        float height = dimensions.z;  // `z` represents cone height
        float radius = dimensions.y;  // `y` represents cone base radius
        if (relativePos.y < 0.0f || relativePos.y > height) {
            return false;
        }

        // Calculate the cone radius at the projectile's height
        float coneRadiusAtHeight = radius * (1.0f - relativePos.y / height);

        // Check if the projectile is within the cone's radius
        float distanceSquared = relativePos.x * relativePos.x + relativePos.z * relativePos.z;
        return distanceSquared <= coneRadiusAtHeight * coneRadiusAtHeight;
    }
}

Mesh* Lab7::createSphere(const std::string& name, int stacks, int slices, const glm::vec3& color)
{
    std::vector<VertexFormat> vertices;
    std::vector<unsigned int> indices;

    for (int stack = 0; stack <= stacks; ++stack) {
        float phi = glm::pi<float>() * float(stack) / stacks;
        for (int slice = 0; slice <= slices; ++slice) {
            float theta = 2.0f * glm::pi<float>() * float(slice) / slices;

            // Spherical coordinates
            float x = sin(phi) * cos(theta);
            float y = cos(phi);
            float z = sin(phi) * sin(theta);

            vertices.emplace_back(glm::vec3(x, y, z), color);
        }
    }

    for (int stack = 0; stack < stacks; ++stack) {
        for (int slice = 0; slice < slices; ++slice) {
            int first = stack * (slices + 1) + slice;
            int second = first + slices + 1;

            indices.push_back(first);
            indices.push_back(second);
            indices.push_back(first + 1);

            indices.push_back(second);
            indices.push_back(second + 1);
            indices.push_back(first + 1);
        }
    }

    Mesh* sphere = new Mesh(name);
    sphere->InitFromData(vertices, indices);
    return sphere;
}

void Lab7::shootProjectile()
{
    // Calculate the direction the drone is facing
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(droneRotationY), glm::vec3(0, 1, 0));
    glm::vec3 forward = glm::normalize(glm::vec3(rotationMatrix[2])); 

    // Create a new projectile
    Projectile projectile;
    projectile.position = dronePosition; // Start at the drone's position
    projectile.direction = -forward;
    projectile.speed = projectileSpeed;

    // Add the projectile to the list
    projectiles.push_back(projectile);
}

void Lab7::updateProjectiles(float deltaTime)
{
    for (size_t i = 0; i < projectiles.size();) {
        Projectile& projectile = projectiles[i];

        projectile.position += projectile.direction * projectile.speed * deltaTime;

        glm::mat4 modelMatrix = glm::translate(glm::mat4(1), projectile.position);
        modelMatrix = glm::scale(modelMatrix, glm::vec3(0.1f));
        RenderSimpleMesh(meshes["sphere"], shaders["LabShader"], modelMatrix, glm::vec3(1.0f, 0.0f, 0.0f));

        bool collided = false;
        for (size_t j = 0; j < enemies.size(); ++j) {
            EnemyDrone& enemy = enemies[j];

            // Calculăm distanța între proiectil și drona inamică
            float distance = glm::distance(projectile.position, enemy.position);

            if (distance < 0.5f) {
                std::cout << "Drona inamică distrusă!" << std::endl;

                // Animația de distrugere
                enemy.isBeingDestroyed = true;
                enemy.destructionProgress = 0.0f;

                playerScore += 10;
                std::cout << "Scor: " << playerScore << std::endl;

                projectiles.erase(projectiles.begin() + i);
                collided = true;
                break;
            }
        }

        if (!collided) {
            if (projectile.position.x < -25.0f || projectile.position.x > 25.0f ||
                projectile.position.z < -25.0f || projectile.position.z > 25.0f ||
                projectile.position.y < 0.0f || projectile.position.y > 50.0f) {
                projectiles.erase(projectiles.begin() + i);
            } else {
                ++i;
            }
        }
    }
}

void Lab7::createEnemyDrones()
{
    /*
    // Pozițiile de spawn ale inamicilor (margini)
    std::vector<glm::vec3> enemyPositions = {
        glm::vec3(-25.0f, 10.0f, -25.0f),
        glm::vec3(-25.0f, 10.0f, 25.0f),
        glm::vec3(25.0f, 10.0f, -25.0f),
        glm::vec3(25.0f, 10.0f, 25.0f)
    }; 

    // Spawn fiecare inamic
    for (const glm::vec3& position : enemyPositions) {
        EnemyDrone enemy = { position, 0.0f }; // Poziția inițială și rotația
        enemies.push_back(enemy);
    } */

    for (int i = 0; i < 4; ++i) {
        glm::vec3 randomPosition = generateRandomPositionDrone(-25.0f, 25.0f, -25.0f, 25.0f, 1.0f, 10.0f);
        EnemyDrone enemy = { randomPosition, 0.0f, 0.0f };
        enemies.push_back(enemy);
    }


}

void Lab7::moveEnemyDrones(float deltaTime)
{
    // Colțuri harta
    std::vector<glm::vec3> cornerPositions = {
        glm::vec3(-25.0f, 10.0f, -25.0f),
        glm::vec3(-25.0f, 10.0f, 25.0f),
        glm::vec3(25.0f, 10.0f, -25.0f),
        glm::vec3(25.0f, 10.0f, 25.0f)
    };

    glm::vec3 playerPosition = dronePosition;

    for (EnemyDrone& enemy : enemies) {
        // Actualizăm cooldown-ul pentru lansarea proiectilului
        enemy.fireCooldown -= deltaTime;

        // Găsim cel mai apropiat colț pentru deplasare
        glm::vec3 closestCorner;
        float minDistance = std::numeric_limits<float>::max();

        for (const glm::vec3& corner : cornerPositions) {
            float distance = glm::distance(enemy.position, corner);
            if (distance < minDistance) {
                minDistance = distance;
                closestCorner = corner;
            }
        }

        // Calculăm direcția spre colț și actualizăm poziția
        glm::vec3 directionToCorner = glm::normalize(closestCorner - enemy.position);
        enemy.position += directionToCorner * enemySpeed * deltaTime;

        // Calculăm rotația dronei spre jucător
        glm::vec3 directionToPlayer = glm::normalize(glm::vec3(
            playerPosition.x - enemy.position.x,
            0.0f,
            playerPosition.z - enemy.position.z
        ));
        enemy.rotationY = atan2(-directionToPlayer.x, -directionToPlayer.z);

        if (minDistance < 1.0f) {
            enemy.position = generateRandomPosition(-25.0f, 25.0f, -25.0f, 25.0f, 10.0f);
        }

        if (enemy.fireCooldown <= 0.0f) {
            fireProjectile(enemy.position, playerPosition);
            enemy.fireCooldown = 3.0f;
        }
    }
}

void Lab7::updateEnemyProjectiles(float deltaTime)
{
    for (size_t i = 0; i < enemyProjectiles.size(); ++i) {
        EnemyProjectile& projectile = enemyProjectiles[i];

       projectile.position += projectile.direction * projectile.speed * deltaTime;

        glm::mat4 modelMatrix = glm::mat4(1);
        modelMatrix = glm::translate(modelMatrix, projectile.position);
        modelMatrix = glm::scale(modelMatrix, glm::vec3(0.2f));
        RenderSimpleMesh(meshes["sphere"], shaders["LabShader"], modelMatrix, glm::vec3(1.0f, 0.0f, 0.0f));

        if (checkSpherePlayerCollision(projectile.position, 0.2f)) {
            std::cout << "Coliziune cu proiectil! Viata scade cu 10." << std::endl;
            droneHealth -= 10;
            std::cout << "Viata ramasa: " << droneHealth << std::endl;

            enemyProjectiles.erase(enemyProjectiles.begin() + i);
            --i;

            if (droneHealth <= 0) {
                std::cout << "Game Over! Jucatorul a fost distrus." << std::endl;
            }
        } else {
            if (projectile.position.x < -25.0f || projectile.position.x > 25.0f ||
                projectile.position.z < -25.0f || projectile.position.z > 25.0f ||
                projectile.position.y < 0.0f || projectile.position.y > 50.0f) {
                enemyProjectiles.erase(enemyProjectiles.begin() + i);
            }
        }
    }
}

bool Lab7::checkSpherePlayerCollision(const glm::vec3& spherePosition, float sphereRadius)
{
    glm::vec3 playerMin, playerMax;
    getDroneAABB(playerMin, playerMax);

    glm::vec3 playerCenter = (playerMin + playerMax) / 2.0f;

    float distance = glm::distance(spherePosition, playerCenter);

    sphereRadius = 0.5f;

    return distance <= sphereRadius;
}

void Lab7::updateEnemyDrones(float deltaTime)
{
    glm::vec3 enemyColor = glm::vec3(0.5f, 0.0f, 0.5f); // Changed to purple
    glm::vec3 frontCubeColor = glm::vec3(0.5f, 0.5f, 0.5f);
    glm::vec3 propellerColor = glm::vec3(0.0f, 0.0f, 0.0f); 

    for (size_t i = 0; i < enemies.size(); ++i) {
        EnemyDrone& enemy = enemies[i];

        if (enemy.isBeingDestroyed) {
            enemy.destructionProgress += deltaTime;

            float scale = 1.0f - enemy.destructionProgress;
            float dropAmount = enemy.destructionProgress * 2.0f;
            float rotationAmount = enemy.destructionProgress * 10.0f;

            if (scale <= 0.0f) {
                glm::vec3 newPosition = generateRandomPositionDrone(-25.0f, 25.0f, -25.0f, 25.0f, 1.0f, 10.0f);
                enemies[i] = EnemyDrone(newPosition, 0.0f, 3.0f);
                continue;
            }

            glm::mat4 modelMatrix = glm::mat4(1);
            modelMatrix = glm::translate(modelMatrix, enemy.position - glm::vec3(0.0f, dropAmount, 0.0f));
            modelMatrix = glm::rotate(modelMatrix, glm::radians(rotationAmount), glm::vec3(0, 1, 0));
            modelMatrix = glm::scale(modelMatrix, glm::vec3(scale));

            glm::mat4 bodyMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["body1"], shaders["LabShader"], bodyMatrix, enemyColor);

            bodyMatrix = glm::rotate(modelMatrix, glm::radians(-45.0f), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["body2"], shaders["LabShader"], bodyMatrix, enemyColor);

            glm::mat4 cubeMatrix;

            cubeMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, 0.2f, 0.9f));
            RenderSimpleMesh(meshes["cube1"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller1"], shaders["LabShader"], cubeMatrix, propellerColor);

            cubeMatrix = glm::translate(modelMatrix, glm::vec3(0.0f, 0.2f, -0.9f));
            RenderSimpleMesh(meshes["cube2"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller2"], shaders["LabShader"], cubeMatrix, propellerColor);

            cubeMatrix = glm::translate(modelMatrix, glm::vec3(-0.9f, 0.2f, 0.0f));
            RenderSimpleMesh(meshes["cube3"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller3"], shaders["LabShader"], cubeMatrix, propellerColor);

            cubeMatrix = glm::translate(modelMatrix, glm::vec3(0.9f, 0.2f, 0.0f));
            RenderSimpleMesh(meshes["cube4"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller4"], shaders["LabShader"], cubeMatrix, propellerColor);

        } else {
            glm::mat4 modelMatrix = glm::mat4(1);
            modelMatrix = glm::translate(modelMatrix, enemy.position);
            modelMatrix = glm::rotate(modelMatrix, enemy.rotationY, glm::vec3(0, 1, 0)); 

            glm::mat4 bodyMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["body1"], shaders["LabShader"], bodyMatrix, enemyColor);

            bodyMatrix = glm::rotate(modelMatrix, glm::radians(-45.0f), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["body2"], shaders["LabShader"], bodyMatrix, enemyColor);

            glm::mat4 cubeMatrix;
            float cubeSize = 0.2f;
            float propellerHeight = 0.05f;

            cubeMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, 0.2f, 0.9f));
            RenderSimpleMesh(meshes["cube1"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller1"], shaders["LabShader"], cubeMatrix, glm::vec3(0.0f, 0.0f, 0.0f));

            cubeMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, 0.2f, -0.9f));
            RenderSimpleMesh(meshes["cube2"], shaders["LabShader"], cubeMatrix, frontCubeColor);

            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller2"], shaders["LabShader"], cubeMatrix, glm::vec3(0.0f, 0.0f, 0.0f)); 

            cubeMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(-0.9f, 0.2f, 0.0f));
            RenderSimpleMesh(meshes["cube3"], shaders["LabShader"], cubeMatrix, enemyColor);

            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller3"], shaders["LabShader"], cubeMatrix, glm::vec3(0.0f, 0.0f, 0.0f));

            cubeMatrix = glm::rotate(modelMatrix, glm::radians(45.0f), glm::vec3(0, 1, 0));
            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.9f, 0.2f, 0.0f));
            RenderSimpleMesh(meshes["cube4"], shaders["LabShader"], cubeMatrix, frontCubeColor);

            cubeMatrix = glm::translate(cubeMatrix, glm::vec3(0.0f, cubeSize / 2.0f + propellerHeight / 2.0f, 0.0f));
            cubeMatrix = glm::rotate(cubeMatrix, glm::radians(propellerRotationAngle), glm::vec3(0, 1, 0));
            RenderSimpleMesh(meshes["propeller4"], shaders["LabShader"], cubeMatrix, glm::vec3(0.0f, 0.0f, 0.0f));
        }
    }
}

void Lab7::fireProjectile(const glm::vec3& startPosition, const glm::vec3& targetPosition)
{
    glm::vec3 direction = glm::normalize(targetPosition - startPosition);

    EnemyProjectile projectile = {
        startPosition,
        direction,
        10.0f //viteza
    };

    enemyProjectiles.push_back(projectile);
}

void Lab7::createArrow()
{
    Mesh* leftArm = createParallelepiped("leftArm", glm::vec3(0.25f, 0.1f, 0.5f), glm::vec3(1.0f, 1.0f, 0.0f));  
    meshes[leftArm->GetMeshID()] = leftArm;

    Mesh* rightArm = createParallelepiped("rightArm", glm::vec3(0.25f, 0.1f, 0.5f), glm::vec3(1.0f, 1.0f, 0.0f)); 
    meshes[rightArm->GetMeshID()] = rightArm;

    Mesh* connector = createParallelepiped("connector", glm::vec3(0.25f, 0.1f, 1.0f), glm::vec3(1.0f, 1.0f, 0.0f));  
    meshes[connector->GetMeshID()] = connector;
}

void Lab7::updateArrowDirection()
{

    // Găsim cea mai apropiată dronă inamică
    float minDistance = std::numeric_limits<float>::max();
    glm::vec3 closestEnemyPosition;

    for (const EnemyDrone& enemy : enemies) {
        float distance = glm::distance(dronePosition, enemy.position);
        if (distance < minDistance) {
            minDistance = distance;
            closestEnemyPosition = enemy.position;
        }
    }

    arrowDirection = glm::normalize(closestEnemyPosition - dronePosition);
}

void Lab7::renderArrow()
{
    glm::vec3 arrowPosition = dronePosition + glm::normalize(arrowDirection) * 2.0f;
    arrowPosition.y += 1.5f;

    float angle = atan2(arrowDirection.x, arrowDirection.z);

    glm::mat4 modelMatrix = glm::mat4(1);
    modelMatrix = glm::translate(modelMatrix, arrowPosition);
    modelMatrix = glm::rotate(modelMatrix, angle, glm::vec3(0, 1, 0));

    glm::mat4 leftArmMatrix = modelMatrix;
    leftArmMatrix = glm::translate(leftArmMatrix, glm::vec3(0.0f, 0.15f, 0.9f));
    leftArmMatrix = glm::rotate(leftArmMatrix, glm::radians(45.0f), glm::vec3(1, 0, 0));
    RenderSimpleMesh(meshes["leftArm"], shaders["LabShader"], leftArmMatrix, glm::vec3(1.0f, 1.0f, 0.0f));

    glm::mat4 rightArmMatrix = modelMatrix;
    rightArmMatrix = glm::translate(rightArmMatrix, glm::vec3(0.0f, -0.15f, 0.9f));
    rightArmMatrix = glm::rotate(rightArmMatrix, glm::radians(-45.0f), glm::vec3(1, 0, 0));
    RenderSimpleMesh(meshes["rightArm"], shaders["LabShader"], rightArmMatrix, glm::vec3(1.0f, 1.0f, 0.0f));

    glm::mat4 connectorMatrix = modelMatrix;
    connectorMatrix = glm::translate(connectorMatrix, glm::vec3(0.0f, 0.0f, 0.5f));
    RenderSimpleMesh(meshes["connector"], shaders["LabShader"], connectorMatrix, glm::vec3(1.0f, 1.0f, 0.0f));
}

void Lab7::Init()
{
    
    {
        glm::vec3 gridOffset = glm::vec3(-25.0f, 0.0f, -25.0f);
        Mesh* mesh = createTerrain("terrain", 50, 50, glm::vec3(0.0f, 1.0f, 0.0f), gridOffset); // 50x50 grid
        meshes[mesh->GetMeshID()] = mesh;
    }

    Mesh* sphereMesh = createSphere("sphere", 20, 20, glm::vec3(1.0f, 0.0f, 0.0f));
    meshes[sphereMesh->GetMeshID()] = sphereMesh;

    createArrow();


    createDrone(glm::vec3(0.0f, 5.0f, 0.0f));
    generateDecor(10, -25.0f, 25.0f, -25.0f, 25.0f, 10.0f);
    createEnemyDrones(); 

    {
        Shader* shader = new Shader("LabShader");
        shader->AddShader(PATH_JOIN(window->props.selfDir, SOURCE_PATH::M1, "lab7", "shaders", "VertexShader.glsl"), GL_VERTEX_SHADER);
        shader->AddShader(PATH_JOIN(window->props.selfDir, SOURCE_PATH::M1, "lab7", "shaders", "FragmentShader.glsl"), GL_FRAGMENT_SHADER);
        shader->CreateAndLink();
        shaders[shader->GetName()] = shader;
    }

    {
        lightPosition = glm::vec3(0, 1, 1);
        materialShininess = 30;
        materialKd = 0.5;
        materialKs = 0.5;
    }
}

void Lab7::FrameStart()
{
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glm::ivec2 resolution = window->GetResolution();
    glViewport(0, 0, resolution.x, resolution.y);
}

void Lab7::Update(float deltaTimeSeconds)
{
    if (droneHealth > 0) {
    propellerRotationAngle += 360.0f * deltaTimeSeconds; 

    if (propellerRotationAngle >= 360.0f)
        propellerRotationAngle -= 360.0f;

    // Actualizăm camera în funcție de modul selectat
    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(droneRotationY), glm::vec3(0, 1, 0));
    glm::vec3 forward = glm::normalize(glm::vec3(rotationMatrix[2])); 
    glm::vec3 up = glm::vec3(0, 1, 0);

    if (isFirstPerson) {
        glm::vec3 cameraPosition = dronePosition - forward * cameraOffsetFirstPerson.z + up * cameraOffsetFirstPerson.y;
        GetSceneCamera()->SetPosition(cameraPosition);
        GetSceneCamera()->SetRotation(glm::vec3(0, glm::radians(droneRotationY), 0));
    } else {
        glm::vec3 cameraPosition = dronePosition - forward * cameraOffsetThirdPerson.z + up * cameraOffsetThirdPerson.y;
        GetSceneCamera()->SetPosition(cameraPosition);
        GetSceneCamera()->SetRotation(glm::vec3(0, glm::radians(droneRotationY), 0));
    }

   {
        glm::mat4 modelMatrix = glm::mat4(1);
        RenderSimpleMesh(meshes["terrain"], shaders["LabShader"], modelMatrix, glm::vec3(0.0f, 1.0f, 0.0f)); 
    }

    updateDrone();

    updateDecor();

    checkCollisions();

    updateProjectiles(deltaTimeSeconds);

    updateEnemyDrones(deltaTimeSeconds);

    moveEnemyDrones(deltaTimeSeconds);

    updateEnemyDrones(deltaTimeSeconds);

    updateEnemyProjectiles(deltaTimeSeconds);

    updateArrowDirection();

    //renderArrow();

    }
}

void Lab7::FrameEnd()
{
   // DrawCoordinateSystem();
}

void Lab7::RenderSimpleMesh(Mesh* mesh, Shader* shader, const glm::mat4& modelMatrix, const glm::vec3& color)
{
    if (!mesh || !shader || !shader->GetProgramID())
        return;

    glUseProgram(shader->program);

    glUniform3fv(glGetUniformLocation(shader->program, "light_position"), 1, glm::value_ptr(lightPosition));
    glm::vec3 eyePosition = GetSceneCamera()->m_transform->GetWorldPosition();
    glUniform3fv(glGetUniformLocation(shader->program, "eye_position"), 1, glm::value_ptr(eyePosition));

    glUniform1i(glGetUniformLocation(shader->program, "material_shininess"), materialShininess);
    glUniform1f(glGetUniformLocation(shader->program, "material_kd"), materialKd);
    glUniform1f(glGetUniformLocation(shader->program, "material_ks"), materialKs);

    glUniform3fv(glGetUniformLocation(shader->program, "object_color"), 1, glm::value_ptr(color));

 

    GLint loc_model_matrix = glGetUniformLocation(shader->program, "Model");
    glUniformMatrix4fv(loc_model_matrix, 1, GL_FALSE, glm::value_ptr(modelMatrix));

    glm::mat4 viewMatrix = GetSceneCamera()->GetViewMatrix();
    int loc_view_matrix = glGetUniformLocation(shader->program, "View");
    glUniformMatrix4fv(loc_view_matrix, 1, GL_FALSE, glm::value_ptr(viewMatrix));

    glm::mat4 projectionMatrix = GetSceneCamera()->GetProjectionMatrix();
    int loc_projection_matrix = glGetUniformLocation(shader->program, "Projection");
    glUniformMatrix4fv(loc_projection_matrix, 1, GL_FALSE, glm::value_ptr(projectionMatrix));

    glBindVertexArray(mesh->GetBuffers()->m_VAO);
    glDrawElements(mesh->GetDrawMode(), static_cast<int>(mesh->indices.size()), GL_UNSIGNED_INT, 0);
}


/*
 *  These are callback functions. To find more about callbacks and
 *  how they behave, see input_controller.h.
 */


void Lab7::OnInputUpdate(float deltaTime, int mods)
{
    float speed = 5.0f;

    previousDronePosition = dronePosition;

    glm::mat4 rotationMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(droneRotationY), glm::vec3(0, 1, 0));

    glm::vec3 forward = glm::normalize(glm::vec3(rotationMatrix[2])); 
    glm::vec3 right = glm::normalize(glm::vec3(rotationMatrix[0]));   
    glm::vec3 up = glm::vec3(0, 1, 0);                                

    if (window->KeyHold(GLFW_KEY_W)) dronePosition -= forward * speed * deltaTime; 
    if (window->KeyHold(GLFW_KEY_S)) dronePosition += forward * speed * deltaTime; 
    if (window->KeyHold(GLFW_KEY_A)) dronePosition -= right * speed * deltaTime;   
    if (window->KeyHold(GLFW_KEY_D)) dronePosition += right * speed * deltaTime;   

    if (window->KeyHold(GLFW_KEY_UP)) dronePosition += up * speed * deltaTime;
    if (window->KeyHold(GLFW_KEY_DOWN)) dronePosition -= up * speed * deltaTime;
    if (window->KeyHold(GLFW_KEY_LEFT)) droneRotationY += 60.0f * deltaTime;
    if (window->KeyHold(GLFW_KEY_RIGHT)) droneRotationY -= 60.0f * deltaTime;
}

void Lab7::OnKeyPress(int key, int mods)
{
    if (key == GLFW_KEY_C) {
        isFirstPerson = !isFirstPerson; 
    }
    if (key == GLFW_KEY_SPACE) {
        shootProjectile();
    }
}


void Lab7::OnKeyRelease(int key, int mods)
{
    // Add key release event
}


void Lab7::OnMouseMove(int mouseX, int mouseY, int deltaX, int deltaY)
{
    // Add mouse move event
}


void Lab7::OnMouseBtnPress(int mouseX, int mouseY, int button, int mods)
{
    // Add mouse button press event
}


void Lab7::OnMouseBtnRelease(int mouseX, int mouseY, int button, int mods)
{
    // Add mouse button release event
}


void Lab7::OnMouseScroll(int mouseX, int mouseY, int offsetX, int offsetY)
{
}


void Lab7::OnWindowResize(int width, int height)
{
}