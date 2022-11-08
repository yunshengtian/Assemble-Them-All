#pragma once

#include "Simulation.h"

namespace redmax {

class SimEnvGenerator {
public:
    static Simulation* createSinglePendulumTest(std::string integrator = "BDF2");
    static Simulation* createSinglePendulumObjTest(std::string integrator = "BDF2");
    static Simulation* createMultiPendulumTest(int num_links, std::string integrator = "BDF2");
    static Simulation* createPrismaticTest(std::string integrator = "BDF2");
    static Simulation* createCableTest(std::string integrator = "BDF2");
    static Simulation* createFree2DTest(std::string integrator = "BDF2");
    static Simulation* createFree3DEulerTest(std::string integrator = "BDF2");
    static Simulation* createFree3DExpTest(std::string integrator = "BDF2");
    static Simulation* createGroundContactTest(std::string integrator = "BDF2");
    static Simulation* createBoxContactTest(std::string integrator = "BDF2");
    static Simulation* createTorqueFingerDemo(std::string integrator = "BDF2");
    static Simulation* createTorqueFingerFlickDemo(std::string integrator = "BDF2");
    static Simulation* createSphereGroundContactTest(std::string integrator = "BDF2");
};

}