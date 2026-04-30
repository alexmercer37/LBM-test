#include "../inc/LBM.hpp"

PipePhaseChangeLBM::PipePhaseChangeLBM(int nx_, int ny_, int nstep_, double dt_)
    : nx(nx_), ny(ny_), nstep(nstep_), dt(dt_)
{

    tau_f = 1.0;
    tau_t = 1.2;
    tau_phi = 1.2;

    u_in = 0.01;
    T_in = 1.7;

    rho_l = 5.0;
    rho_g = 0.4;

    Tsat = 0.5;
    Kphase = 0.1;

    latent = 2.0;
    cp = 1.0;

    rho.assign(nx, vector<double>(ny, rho_l));
    ux.assign(nx, vector<double>(ny, 0.0));
    uy.assign(nx, vector<double>(ny, 0.0));
    T.assign(nx, vector<double>(ny, T_in));
    phi.assign(nx, vector<double>(ny, 1.0));
    mdot.assign(nx, vector<double>(ny, 0.0));

    f.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));
    f_new.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));

    g.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));
    g_new.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));

    h.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));
    h_new.assign(9, vector<vector<double>>(nx, vector<double>(ny, 0.0)));

    Twall = 2.2;
    liquid_mass = 0.0;
    vapor_mass = 0.0;
}

double PipePhaseChangeLBM::feq(int i, double rho_, double ux_, double uy_)
{
    double eu = ex[i] * ux_ + ey[i] * uy_;
    double u2 = ux_ * ux_ + uy_ * uy_;
    return w[i] * rho_ * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
}

double PipePhaseChangeLBM::geq(int i, double T_, double ux_, double uy_)
{
    double eu = ex[i] * ux_ + ey[i] * uy_;
    double u2 = ux_ * ux_ + uy_ * uy_;
    return w[i] * T_ * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
}

double PipePhaseChangeLBM::heq(int i, double phi_, double ux_, double uy_)
{
    double eu = ex[i] * ux_ + ey[i] * uy_;
    double u2 = ux_ * ux_ + uy_ * uy_;
    return w[i] * phi_ * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
}

void PipePhaseChangeLBM::Initialize()
{
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            ux[x][y] = 0.0;
            uy[x][y] = 0.0;

            phi[x][y] = 1.0;
            T[x][y] = T_in;

            rho[x][y] = phi[x][y] * rho_l + (1.0 - phi[x][y]) * rho_g;

            for (int i = 0; i < 9; i++)
            {
                f[i][x][y] = feq(i, rho[x][y], ux[x][y], uy[x][y]);
                g[i][x][y] = geq(i, T[x][y], ux[x][y], uy[x][y]);
                h[i][x][y] = heq(i, phi[x][y], ux[x][y], uy[x][y]);
            }
        }
    }
}

void PipePhaseChangeLBM::ComputeMdot()
{
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            double dT = T[i][j] - Tsat;

            if (dT > 0.0 && phi[i][j] > 0.01)
            {
                mdot[i][j] = Kphase * dT * phi[i][j];

                // double max_evap = 0.8 * phi[i][j] / dt;
                // if (mdot[i][j] > max_evap)
                // {
                //     mdot[i][j] = max_evap;
                // }
            }

            else
            {
                mdot[i][j] = 0.0;
            }
        }
    }
}

void PipePhaseChangeLBM::CollisionF()
{
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            for (int i = 0; i < 9; i++)
            {
                double feqv = feq(i, rho[x][y], ux[x][y], uy[x][y]);
                f_new[i][x][y] = f[i][x][y] - (dt / tau_f) * (f[i][x][y] - feqv);
            }
        }
    }
}

void PipePhaseChangeLBM::CollisionT()
{
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {

            double sourceT = -(latent / cp) * mdot[x][y];

            for (int i = 0; i < 9; i++)
            {
                double geqv = geq(i, T[x][y], ux[x][y], uy[x][y]);
                g_new[i][x][y] = g[i][x][y] - (dt / tau_t) * (g[i][x][y] - geqv) + dt * w[i] * sourceT;
            }
        }
    }
}

void PipePhaseChangeLBM::CollisionPhi()
{
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {

            double sourcePhi = -mdot[x][y];

            for (int i = 0; i < 9; i++)
            {
                double heqv = heq(i, phi[x][y], ux[x][y], uy[x][y]);
                h_new[i][x][y] = h[i][x][y] - (dt / tau_phi) * (h[i][x][y] - heqv) + dt * w[i] * sourcePhi;
            }
        }
    }
}

void PipePhaseChangeLBM::Streaming()
{

    for (int i = 0; i < 9; i++)
        for (int x = 0; x < nx; x++)
            for (int y = 0; y < ny; y++)
            {
                f[i][x][y] = 0.0;
                g[i][x][y] = 0.0;
                h[i][x][y] = 0.0;
            }

    for (int i = 0; i < 9; i++)
    {
        for (int x = 0; x < nx; x++)
        {
            for (int y = 0; y < ny; y++)
            {
                int xn = x + ex[i];
                int yn = y + ey[i];

                if (xn >= 0 && xn < nx && yn >= 0 && yn < ny)
                {
                    f[i][xn][yn] += f_new[i][x][y];
                    g[i][xn][yn] += g_new[i][x][y];
                    h[i][xn][yn] += h_new[i][x][y];
                }
            }
        }
    }
}

void PipePhaseChangeLBM::Macroscopic()
{
    for (int x = 0; x < nx; x++)
    {
        for (int y = 0; y < ny; y++)
        {
            double jx = 0.0;
            double jy = 0.0;

            double temp = 0.0;
            double ph = 0.0;

            for (int i = 0; i < 9; i++)
            {
                jx += f[i][x][y] * ex[i];
                jy += f[i][x][y] * ey[i];

                temp += g[i][x][y];
                ph += h[i][x][y];
            }

            phi[x][y] = std::max(0.0, std::min(1.0, ph));

            rho[x][y] = phi[x][y] * rho_l + (1.0 - phi[x][y]) * rho_g;
            if (rho[x][y] < 1e-8)
                rho[x][y] = 1e-8;

            ux[x][y] = jx / rho[x][y];
            uy[x][y] = jy / rho[x][y];

            T[x][y] = temp;
        }
    }
}

void PipePhaseChangeLBM::ApplyBoundary()
{
    int x = 0;

    for (int y = 0; y < ny; y++)
    {
        ux[x][y] = u_in;
        uy[x][y] = 0.0;
        T[x][y] = T_in;
        phi[x][y] = 1.0;

        rho[x][y] = phi[x][y] * rho_l + (1.0 - phi[x][y]) * rho_g;

        for (int i = 0; i < 9; i++)
        {
            f[i][x][y] = feq(i, rho[x][y], ux[x][y], uy[x][y]);
            g[i][x][y] = geq(i, T[x][y], ux[x][y], uy[x][y]);
            h[i][x][y] = heq(i, phi[x][y], ux[x][y], uy[x][y]);
        }
    }

    x = nx - 1;

    for (int y = 0; y < ny; y++)
    {
        ux[x][y] = ux[x - 1][y];
        uy[x][y] = uy[x - 1][y];
        T[x][y] = T[x - 1][y];
        phi[x][y] = phi[x - 1][y];
        rho[x][y] = rho[x - 1][y];

        for (int i = 0; i < 9; i++)
        {
            double feq_out = feq(i, rho[x][y], ux[x][y], uy[x][y]);
            double feq_in = feq(i, rho[x - 1][y], ux[x - 1][y], uy[x - 1][y]);
            f[i][x][y] = feq_out + (f[i][x - 1][y] - feq_in); // 非平衡外推

            double geq_out = geq(i, T[x][y], ux[x][y], uy[x][y]);
            double geq_in = geq(i, T[x - 1][y], ux[x - 1][y], uy[x - 1][y]);
            g[i][x][y] = geq_out + (g[i][x - 1][y] - geq_in);

            double heq_out = heq(i, phi[x][y], ux[x][y], uy[x][y]);
            double heq_in = heq(i, phi[x - 1][y], ux[x - 1][y], uy[x - 1][y]);
            h[i][x][y] = heq_out + (h[i][x - 1][y] - heq_in);
        }
    }

    for (int x = 0; x < nx; x++)
    {
        int yb = 0;
        int yt = ny - 1;

        std::swap(f[2][x][yb], f[4][x][yb]);
        std::swap(f[5][x][yb], f[7][x][yb]);
        std::swap(f[6][x][yb], f[8][x][yb]);

        std::swap(f[2][x][yt], f[4][x][yt]);
        std::swap(f[5][x][yt], f[7][x][yt]);
        std::swap(f[6][x][yt], f[8][x][yt]);

        T[x][yb] = Twall;
        T[x][yt] = Twall;

        rho[x][yb] = phi[x][yb] * rho_l + (1.0 - phi[x][yb]) * rho_g;
        rho[x][yt] = phi[x][yt] * rho_l + (1.0 - phi[x][yt]) * rho_g;

        ux[x][yb] = uy[x][yb] = 0.0;
        ux[x][yt] = uy[x][yt] = 0.0;

        for (int i = 0; i < 9; i++)
        {
            g[i][x][yb] = geq(i, T[x][yb], ux[x][yb], uy[x][yb]);
            g[i][x][yt] = geq(i, T[x][yt], ux[x][yt], uy[x][yt]);
            h[i][x][yb] = heq(i, phi[x][yb], ux[x][yb], uy[x][yb]);
            h[i][x][yt] = heq(i, phi[x][yt], ux[x][yt], uy[x][yt]);
        }
    }
}

void PipePhaseChangeLBM::export_results(const string &filename)
{
    std::ofstream file(filename);
    // file << "x,y,rho,ux,uy,T,phi,mdot,liquid_mass,vapor_mass,outlet_vapor_mass" << std::endl;
    file << "x,y,phi,mdot" << std::endl;

    double liquid_mass = 0.0;
    double vapor_mass = 0.0;
    double outlet_vapor_mass = 0.0;

    for (int x = 0; x < nx; ++x)
    {
        for (int y = 0; y < ny; ++y)
        {
            liquid_mass += phi[x][y] * rho_l;
            vapor_mass += (1.0 - phi[x][y]) * rho_g;
        }
    }

    int x = nx - 1;
    for (int y = 0; y < ny; ++y)
    {
        outlet_vapor_mass += (1.0 - phi[x][y]) * rho_g * ux[x][y];
    }

    for (int j = 0; j < nx; ++j)
    {
        for (int k = 0; k < ny; ++k)
        {
            file << j << "," << k << ","
                 << std::fixed
                 << rho[j][k] << ","
                 //  << ux[j][k] << ","
                 //  << uy[j][k] << ","
                 //  << T[j][k] << ","
                 << phi[j][k] << ","
                 << mdot[j][k] << endl;
            //  << liquid_mass << ","
            //  << vapor_mass << ","
            //  << outlet_vapor_mass << std::endl;
        }
    }

    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

void PipePhaseChangeLBM::Run()
{
    Initialize();

    for (int step = 0; step < nstep; step++)
    {
        ComputeMdot();

        CollisionF();
        CollisionT();
        CollisionPhi();

        Streaming();
        Macroscopic();
        ApplyBoundary();

        if (step % 10 == 0)
        {
            string folder = "data/";
            string filename = folder + "results_" + to_string(step) + ".csv";
            export_results(filename);
        }
    }
}