#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters();

    const static int F_GRAVITY = 1;

    const static int R_SPHERE  = 0;
    const static int R_2BY4    = 1;
    const static int R_BUNNY   = 2;
    const static int R_CUSTOM  = 3;
    const static int R_PLANE   = 4;

    enum UseCage {C_NEVER, C_ALWAYS, C_BROADPHASE};

    bool simRunning;
    double timeStep;
    double NewtonTolerance;
    int NewtonMaxIters;
    double penaltyStiffness;

    int activeForces;
    double gravityG;

    double bodyDensity;
    int launchBody;
    double launchVel;
    bool randomLaunchOrientation;

    bool randomLaunchAngVel;
    double randomLaunchVelMagnitude;

    UseCage useCage;
};

#endif // SIMPARAMETERS_H
