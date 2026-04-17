#include "../inc/D1Q3.hpp"

void D1Q3::export_results(const string &filename, const std::vector<double> &omega)
{
    ofstream outFile(filename);
    outFile << "x,T,Flux\n";

    for (int i = 0; i < m; i++)
    {
        double x_pos = i * dx;
        double flux = (f[1][i] - f[2][i]) * (1.0 - 0.5 * omega[i]);

        outFile << x_pos << "," << T[i] << "," << flux << "\n";
    }

    outFile.close();
    std::cout << "Successfully saved to: " << filename << std::endl;
}

void D1Q3::export_results(const string &filename)
{
    ofstream outFile(filename);
    outFile << "x,T,Flux\n";

    for (int i = 0; i < m; i++)
    {
        double x_pos = i * dx;
        double flux = (f[1][i] - f[2][i]) * (1.0 - 0.5 * omega);

        outFile << x_pos << "," << T[i] << "," << flux << "\n";
    }

    outFile.close();
    std::cout << "Successfully saved to: " << filename << std::endl;
}

#ifdef D1Q3_T

void D1Q3::D1Q3_t()
{

    reset();

    for (int k1 = 1; k1 <= nstep; ++k1)
    {
        for (int i = 0; i < m; ++i)
        {
            double current_T = T[i];

            for (int q = 0; q < 3; ++q)
            {
                f[q][i] = (1.0 - omega) * f[q][i] + omega * w[q] * current_T;
            }
        }

        for (int i = m - 1; i > 0; --i)
        {
            f[1][i] = f[1][i - 1];
        }

        for (int i = 0; i < m - 1; ++i)
        {
            f[2][i] = f[2][i + 1];
        }

        f[1][0] = Tw_left - f[2][0] - f[0][0];
        f[2][m - 1] = Tw_right - f[1][m - 1] - f[0][m - 1];

        for (int i = 0; i < m; i++)
        {
            T[i] = f[0][i] + f[1][i] + f[2][i];
        }
    }

    D1Q3::export_results("T_D1Q3.csv");
}

#endif

#ifdef D1Q3_F

void D1Q3::D1Q3_f()
{

    reset();

    T.assign(m, 0.0);

    for (int q = 0; q < 3; ++q)
    {

        f[q].assign(m, 0.0);
    }

    for (int k1 = 1; k1 <= nstep; ++k1)
    {
        for (int i = 0; i < m; ++i)
        {
            double current_T = T[i];

            for (int q = 0; q < 3; ++q)
            {
                f[q][i] = (1.0 - omega) * f[q][i] + omega * w[q] * current_T;
            }
        }

        for (int i = m - 1; i > 0; --i)
        {
            f[1][i] = f[1][i - 1];
        }

        for (int i = 0; i < m - 1; ++i)
        {
            f[2][i] = f[2][i + 1];
        }

        f[1][0] = f[2][0] + Tw_left * dx / (3 * omega);
        f[2][m - 1] = Tw_right - f[1][m - 1] - f[0][m - 1];

        for (int i = 0; i < m; i++)
        {
            T[i] = f[0][i] + f[1][i] + f[2][i];
        }
    }

    D1Q3::export_results("F_D1Q3.csv");
}

#endif

#ifdef D1Q3_Alpha

void D1Q3::D1Q3_alpha()
{

    reset();

    double alpha_start = 0.1;
    double alpha_end = 1.0;

    vector<double> alpha(m);
    vector<double> omega(m);

    for (int q = 0; q < 3; ++q)
    {
        for (int i = 0; i < m; ++i)
        {
            f[q][i] = w[q] * T[i];
        }
    }

    for (int i = 0; i < m; ++i)
    {
        alpha[i] = alpha_start + i * (alpha_end - alpha_start) / (m - 1);
        omega[i] = 1.0 / (3 * alpha[i] * dt / (dx * dx) + 0.5);

    } // 设置平衡态

    for (int k1 = 1; k1 <= nstep; ++k1)
    {
        for (int i = 0; i < m; ++i)
        {
            double current_T = T[i];
            double curOmega = omega[i];

            for (int q = 0; q < 3; ++q)
            {
                f[q][i] = (1.0 - curOmega) * f[q][i] + curOmega * w[q] * current_T;
            }
        }

        for (int i = m - 1; i > 0; --i)
        {
            f[1][i] = f[1][i - 1];
        }

        for (int i = 0; i < m - 1; ++i)
        {
            f[2][i] = f[2][i + 1];
        }

        f[1][0] = Tw_left - f[2][0] - f[0][0];
        f[2][m - 1] = Tw_right - f[1][m - 1] - f[0][m - 1];

        for (int i = 0; i < m; i++)
        {
            T[i] = f[0][i] + f[1][i] + f[2][i];
        }
    }

    export_results("A_D1Q3.csv", omega);
}

#endif

#ifdef D1Q3_S

void D1Q3::D1Q3_s()
{
    reset();

    T.assign(m, 0);
    f.assign(3, vector<double>(m, 0.0));

    for (int k1 = 0; k1 < nstep; ++k1)
    {
        for (int i = 0; i < m; ++i)
        {
            double current_T = T[i];

            for (int q = 0; q < 3; ++q)
            {
                f[q][i] = (1.0 - omega) * f[q][i] + omega * w[q] * current_T + dt * w[q] * S;
            }
        }

        for (int i = m - 1; i > 0; --i)
        {
            f[1][i] = f[1][i - 1];
        }

        for (int i = 0; i < m - 1; ++i)
        {
            f[2][i] = f[2][i + 1];
        }

        f[1][0] = Tw_left - f[2][0] - f[0][0];
        f[2][m - 1] = Tw_right - f[1][m - 1] - f[0][m - 1];

        for (int i = 0; i < m; i++)
        {
            T[i] = f[0][i] + f[1][i] + f[2][i];
        }
    }

    D1Q3::export_results("S_D1Q3.csv");
}
#endif

#ifdef D1Q3_V

void D1Q3::D1Q3_v()
{
    reset();
    f.assign(3, vector<double>(m, 0));

    for (int k1 = 0; k1 < nstep; ++k1)
    {

        for (int i = 0; i < m; ++i)
        {
            double current_T = T[i];

            for (int q = 0; q < 3; ++q)
            {
                f[q][i] = (1 - omega) * f[q][i] + omega * w[q] * current_T;
            }
        }

        for (int i = m - 1; i > 0; --i)
        {
            f[1][i] = f[1][i - 1];
        }

        for (int i = 0; i < m - 1; ++i)
        {
            f[2][i] = f[2][i + 1];
        }

        f[1][0] = Tw_left - f[2][0] - f[0][0];
        f[2][m - 1] = Tw_right - f[1][m - 1] - f[0][m - 1];

        for (int i = 0; i < m; i++)
        {
            T[i] = f[0][i] + f[1][i] + f[2][i];
        }
    }

    D1Q3::export_results("V_D1Q3.csv");
}

#endif

D1Q3::~D1Q3()
{

    std::cout << "Object destroyed: D1Q3 simulation resources have been released." << std::endl;
}