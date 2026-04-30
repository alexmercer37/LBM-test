#include "../inc/main.hpp"
#include "../inc/LBM.hpp"

int main()
{

    int nx = 200;
    int ny = 30;
    int nstep = 500;
    double dt = 1.0;

    PipePhaseChangeLBM sim(nx, ny, nstep, dt);
    sim.Run();

    return 0;
}