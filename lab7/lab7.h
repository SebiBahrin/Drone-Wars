#pragma once

#include "components/simple_scene.h"
#include "components/transform.h"

namespace m1
{
struct Projectile {
    glm::vec3 position;
    glm::vec3 direction;
    float speed;
};

struct EnemyDrone {
    glm::vec3 position;
    float rotationY;
    float fireCooldown;
    bool isBeingDestroyed = false;
    float destructionProgress = 0.0f;

    EnemyDrone(const glm::vec3& pos, float rotY, float cooldown)
        : position(pos), rotationY(rotY), fireCooldown(cooldown) {}
};

struct EnemyProjectile {
    glm::vec3 position;
    glm::vec3 direction;
    float speed;
};

class Lab7 : public gfxc::SimpleScene
{
public:
    Lab7();
    ~Lab7();

    void Init() override;

private:
    void FrameStart() override;
    void Update(float deltaTimeSeconds) override;
    void FrameEnd() override;

    void RenderSimpleMesh(Mesh* mesh, Shader* shader, const glm::mat4& modelMatrix, const glm::vec3& color = glm::vec3(1));

    void OnInputUpdate(float deltaTime, int mods) override;
    void OnKeyPress(int key, int mods) override;
    void OnKeyRelease(int key, int mods) override;
    void OnMouseMove(int mouseX, int mouseY, int deltaX, int deltaY) override;
    void OnMouseBtnPress(int mouseX, int mouseY, int button, int mods) override;
    void OnMouseBtnRelease(int mouseX, int mouseY, int button, int mods) override;
    void OnMouseScroll(int mouseX, int mouseY, int offsetX, int offsetY) override;
    void OnWindowResize(int width, int height) override;

    glm::vec3 lightPosition;
    unsigned int materialShininess;
    float materialKd;
    float materialKs;

    Mesh* createCube(const std::string& name, float size, const glm::vec3& color);
    Mesh* createParallelepiped(const std::string& name, glm::vec3 dimensions, const glm::vec3& color);
    void createDrone(const glm::vec3& initialPosition);
    void updateDrone();

    float rotationAngle = 90;
    float propellerRotationAngle = 0.0f;
    glm::vec3 dronePosition = glm::vec3(0.0f, 0.5f, 0.0f);
    float droneRotationY = 0.0f;

    bool isFirstPerson = false;
    glm::vec3 cameraOffsetThirdPerson = glm::vec3(0.0f, 2.0f, -5.0f);
    glm::vec3 cameraOffsetFirstPerson = glm::vec3(0.0f, 0.5f, 0.5f);

    Mesh* createTerrain(const std::string& name, int m, int n, const glm::vec3& color, const glm::vec3& gridOffset);

    void createTree(const glm::vec3& position);
    Mesh* createCone(const std::string& name, int segments, float height, float radius, const glm::vec3& color);
    Mesh* createCylinder(const std::string& name, int segments, float height, float radius, const glm::vec3& color);

    void generateDecor(int numObstacles, float xmin, float xmax, float zmin, float zmax, float minDistance);
    glm::vec3 generateRandomPosition(float xmin, float xmax, float zmin, float zmax, float y);
    bool isPositionValid(const glm::vec3& position, const std::vector<glm::vec3>& existingPositions, float minDistance);
    void createBuilding(const glm::vec3& position, const glm::vec3& dimensions, const glm::vec3& color);

    std::vector<glm::vec3> obstaclePositions;
    std::vector<std::string> obstacleTypes;
    std::vector<glm::vec3> obstacleDimensions;

    void updateDecor();

    bool checkAABBCollision(const glm::vec3& min1, const glm::vec3& max1, const glm::vec3& min2, const glm::vec3& max2);

    void getDroneAABB(glm::vec3& min, glm::vec3& max);
    void getObstacleAABB(const glm::vec3& position, const glm::vec3& dimensions, glm::vec3& min, glm::vec3& max, bool isTree);

    bool isColliding = false;
    void checkCollisions();
    glm::vec3 previousDronePosition;

    std::vector<Projectile> projectiles;
    float projectileSpeed = 10.0f;
    float projectileLifetime = 5.0f;
    float elapsedTime = 0.0f;

    void shootProjectile();
    void updateProjectiles(float deltaTime);

    Mesh* createSphere(const std::string& name, int stacks, int slices, const glm::vec3& color);
    bool checkProjectileCollision(const glm::vec3& projectilePosition, const glm::vec3& obstaclePosition, const glm::vec3& dimensions, bool isCone);

    int droneHealth = 100;
    int counter = 0;

    std::vector<EnemyDrone> enemies;
    void createEnemyDrones();
    void updateEnemyDrones(float deltaTime);

    float enemySpeed = 2.0f;
    void moveEnemyDrones(float deltaTime);

    std::vector<EnemyProjectile> enemyProjectiles;
    void fireProjectile(const glm::vec3& startPosition, const glm::vec3& targetPosition);

    void updateEnemyProjectiles(float deltaTime);
    bool checkSpherePlayerCollision(const glm::vec3& spherePosition, float sphereRadius);
    int playerScore = 0;

    glm::vec3 generateRandomPositionDrone(float xmin, float xmax, float zmin, float zmax, float ymin, float ymax);
    void updateArrowDirection();
    void renderArrow();
    glm::vec3 arrowDirection;
    void createArrow();
};
}