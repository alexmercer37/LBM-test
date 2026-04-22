#include "../inc/D2Q9_P.hpp"

double D2Q9_T::calculate_dpdT(double rho_val, double T_val)
{
    double C_omega = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;

    double sqrt_Tr = sqrt(Tc / T_val);
    double alpha_sqrt = 1.0 + C_omega * (1.0 - sqrt_Tr);

    double alpha_prime = -(C_omega * alpha_sqrt) / sqrt(T_val);

    double term1_den = 1.0 - b * rho_val;
    double term2_den = 1.0 + 2.0 * b * rho_val - b * b * rho_val * rho_val;

    double dpdT = (rho_val * R) / term1_den - (a * alpha_prime * rho_val * rho_val) / term2_den;

    return dpdT;
}

double D2Q9_T::get_dTdt(int j, int k, const std::vector<std::vector<double>> &current_T)
{
    double local_u = this->u[j][k];
    double local_v = this->v[j][k];
    double local_rho = this->rho[j][k]; // 提取密度标量
    double T_val = current_T[j][k];     // 提取温度标量

    double lap_T = current_T[(j + 1) % m][k] + current_T[(j - 1 + m) % m][k] +
                   current_T[j][k + 1] + current_T[j][k - 1] - 4.0 * current_T[j][k];

    double dTdx = (current_T[(j + 1) % m][k] - current_T[(j - 1 + m) % m][k]) / 2.0;
    double dTdy = (current_T[j][k + 1] - current_T[j][k - 1]) / 2.0;
    double advection = -(local_u * dTdx + local_v * dTdy);

    double dpdT = D2Q9_T::calculate_dpdT(local_rho, T_val);

    double div_u = (u[(j + 1) % m][k] - u[(j - 1 + m) % m][k]) / 2.0 +
                   (v[j][k + 1] - v[j][k - 1]) / 2.0;

    double source = -(T_val * dpdT * div_u) / (local_rho * Cv + 1e-10);

    return alpha * lap_T + source + advection;
}

void D2Q9_T::rk4_step()
{

    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            k1[j][k] = get_dTdt(j, k, current_T);

    // --- Step 2: Compute k2 ---
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            T_tmp[j][k] = current_T[j][k] + 0.5 * dt * k1[j][k];
    // 应用壁面边界条件(如 T_wall) 到 T_tmp
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            k2[j][k] = get_dTdt(j, k, T_tmp);

    // --- Step 3: Compute k3 ---
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            T_tmp[j][k] = current_T[j][k] + 0.5 * dt * k2[j][k];
    // 应用边界条件
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            k3[j][k] = get_dTdt(j, k, T_tmp);

    // --- Step 4: Compute k4 ---
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            T_tmp[j][k] = current_T[j][k] + dt * k3[j][k];
    // 应用边界条件
    for (int j = 0; j < m; ++j)
        for (int k = 1; k < n - 1; ++k)
            k4[j][k] = get_dTdt(j, k, T_tmp);

    // --- Final Update ---
    for (int j = 0; j < m; ++j)
    {
        for (int k = 1; k < n - 1; ++k)
        {
            current_T[j][k] += (dt / 6.0) * (k1[j][k] + 2.0 * k2[j][k] + 2.0 * k3[j][k] + k4[j][k]);
        }
    }
}

void D2Q9_T::compute_rho_u()
{
    for (int j = 0; j < m; ++j)
    {
        for (int k = 0; k < n; ++k)
        {

            double r = 0.0, jx = 0.0, jy = 0.0;

            for (int i = 0; i < 9; ++i)
            {
                r += f[i][j][k];
                jx += f[i][j][k] * ex[i];
                jy += f[i][j][k] * ey[i];
            }

            rho[j][k] = r;
            // 计算密度 rho

            u[j][k] = (jx + 0.5 * dt * Fx[j][k]) / r;
            v[j][k] = (jy + 0.5 * dt * Fy[j][k]) / r;
            // 加上力项修正 (Guenst-Li 方案核心)

            psi[j][k] = r;
        }
    }
}

void D2Q9_T::compute_i_foeces()
{
    for (int j = 0; j < m; ++j)
    {
        for (int k = 0; k < n; ++k)
        {
            double force_x = 0.0;
            double force_y = 0.0;

            for (int i = 0; i < 9; ++i)
            {
                int next_j = (j + ex[i] + m) % m;
                int next_k = k + ey[i];

                double psi_neighbor;

                if (next_k < 0 || next_k >= n)
                {
                    psi_neighbor = psi_wall;
                    // 壁面润湿性虚拟密度
                }

                else
                {
                    psi_neighbor = psi[next_j][next_k];
                }

                force_x += w_f[i] * psi_neighbor * ex[i];
                force_y += w_f[i] * psi_neighbor * ey[i];
            }
            Fx[j][k] = -G * psi[j][k] * force_x;
            Fy[j][k] = -G * psi[j][k] * force_y;
        }
    }
}

void D2Q9_P::D2Q9_p()
{
    D2Q9_T solver(101, 101, 50000, 1.0, 0.1);

    for (int step = 0; step < nstep; ++step)
    {

        solver.compute_i_foeces();

        solver.compute_rho_u();

        solver.rk4_step();

        for (int j = 0; j < m; ++j)
        {
            for (int k = 0; k < n; ++k)
            {

                for (int i = 0; i < 9; ++i)
                {
                    f_local[i] = f[i][j][k];
                }

                for (int i = 0; i < 9; ++i)
                {
                    m_local[0] += f_local[i]; // 密度
                }

                m_local[1] = -4.0 * f_local[0] - f_local[1] - f_local[2] - f_local[3] - f_local[4] + 2.0 * (f_local[5] + f_local[6] + f_local[7] + f_local[8]);
                // 能量

                m_local[2] = 4.0 * f_local[0] - 2.0 * (f_local[1] + f_local[2] + f_local[3] + f_local[4]) + f_local[5] + f_local[6] + f_local[7] + f_local[8];
                // 能量平方

                m_local[3] = f_local[1] - f_local[3] + f_local[5] - f_local[6] - f_local[7] + f_local[8];
                // x方向动量 (jx)

                m_local[4] = -2.0 * (f_local[1] - f_local[3]) + f_local[5] - f_local[6] - f_local[7] + f_local[8];
                // x方向热流 (qx)

                m_local[5] = f_local[2] - f_local[4] + f_local[5] + f_local[6] - f_local[7] - f_local[8];
                // y方向动量 (jy)

                m_local[6] = -2.0 * (f_local[2] - f_local[4]) + f_local[5] + f_local[6] - f_local[7] - f_local[8];
                // y方向热流 (qy)

                m_local[7] = f_local[1] - f_local[2] + f_local[3] - f_local[4];
                // 应力张量分量 (pxx)

                m_local[8] = f_local[5] - f_local[6] + f_local[7] - f_local[8];
                // 应力张量分量 (pxy)

                double r = rho[j][k];
                double ux = u[j][k];
                double uy = v[j][k];
                double u2 = ux * ux + uy * uy;

                meq = {r,
                       r * (-2.0 + 3.0 * u2),
                       r * (1.0 - 3.0 * u2),
                       r * ux,
                       -r * ux,
                       r * uy,
                       -r * uy,
                       r * (ux * ux - uy * uy),
                       r * ux * uy};
                // 计算平衡态矩阵m_ep

                double FX = Fx[j][k]; // 包含 Fm, Fb, Fw
                double FY = Fy[j][k];
                double psi_val = psi[j][k];

                double term_sigma = (sigma * (FX * FX + FY * FY)) / (psi_val * psi_val * dt * (tau_e - 0.5));

                S_local = {
                    0.0,
                    6.0 * (ux * FX + uy * FY) + term_sigma,
                    -6.0 * (ux * FX + uy * FY) - term_sigma,
                    FX,
                    -FX,
                    FY,
                    -FY,
                    2.0 * (ux * FX - uy * FY),
                    ux * FY + uy * FX};

                for (int i = 0; i < 9; ++i)
                {
                    m_local[i] = m_local[i] - s_diag[i] * (m_local[i] - meq[i]) + dt * (1.0 - 0.5 * s_diag[i]) * S_local[i];
                } // 碰撞

                f_post = {
                    (1.0 / 9.0) * (m_local[0] - m_local[1] + m_local[2]),
                    (1.0 / 18.0) * (2.0 * m_local[0] - m_local[1] - 2.0 * m_local[2] + 3.0 * m_local[3] - 3.0 * m_local[4] + 9.0 * m_local[7]),
                    (1.0 / 18.0) * (2.0 * m_local[0] - m_local[1] - 2.0 * m_local[2] + 3.0 * m_local[5] - 3.0 * m_local[6] - 9.0 * m_local[7]),
                    (1.0 / 18.0) * (2.0 * m_local[0] - m_local[1] - 2.0 * m_local[2] - 3.0 * m_local[3] + 3.0 * m_local[4] + 9.0 * m_local[7]),
                    (1.0 / 18.0) * (2.0 * m_local[0] - m_local[1] - 2.0 * m_local[2] - 3.0 * m_local[5] + 3.0 * m_local[6] - 9.0 * m_local[7]),
                    (1.0 / 36.0) * (4.0 * m_local[0] + 2.0 * m_local[1] + m_local[2] + 6.0 * m_local[3] + 3.0 * m_local[4] + 6.0 * m_local[5] + 3.0 * m_local[6] + 9.0 * m_local[8]),
                    (1.0 / 36.0) * (4.0 * m_local[0] + 2.0 * m_local[1] + m_local[2] - 6.0 * m_local[3] - 3.0 * m_local[4] + 6.0 * m_local[5] + 3.0 * m_local[6] - 9.0 * m_local[8]),
                    (1.0 / 36.0) * (4.0 * m_local[0] + 2.0 * m_local[1] + m_local[2] - 6.0 * m_local[3] - 3.0 * m_local[4] - 6.0 * m_local[5] - 3.0 * m_local[6] + 9.0 * m_local[8]),
                    (1.0 / 36.0) * (4.0 * m_local[0] + 2.0 * m_local[1] + m_local[2] + 6.0 * m_local[3] + 3.0 * m_local[4] - 6.0 * m_local[5] - 3.0 * m_local[6] - 9.0 * m_local[8])};
                // 转换回分布函数

                for (int i = 0; i < 9; ++i)
                {
                    int next_j = (j + ex[i] + m) % m;
                    int next_k = k + ey[i];

                    if (next_k >= 0 && next_k < n)
                    {
                        f_new[i][next_j][next_k] = f_post[i];
                    }
                }
                // 迁移
            }
        }

        for (int j = 0; j < m; ++j)
        {
            int k = 0;
            // 无滑移

            double r_w = f_new[0][j][k] + f_new[1][j][k] + f_new[3][j][k] + 2.0 * (f_new[4][j][k] + f_new[7][j][k] + f_new[8][j][k]);
            // 求密度

            f_new[2][j][k] = f_new[4][j][k];
            f_new[5][j][k] = f_new[7][j][k] - 0.5 * (f_new[1][j][k] - f_new[3][j][k]);
            f_new[6][j][k] = f_new[8][j][k] + 0.5 * (f_new[1][j][k] - f_new[3][j][k]);
            // 下壁面

            f_new[5][j][k] -= 0.25 * dt * Fx[j][k];
            f_new[6][j][k] += 0.25 * dt * Fx[j][k];
            // 考虑外力对边界的影响
        }

        for (int j = 0; j < m; ++j)
        {
            int k = n - 1;
            // 上壁面无滑移

            double r_w = f_new[0][j][k] + f_new[1][j][k] + f_new[3][j][k] + 2.0 * (f_new[2][j][k] + f_new[5][j][k] + f_new[6][j][k]);
            // 计算密度

            f_new[4][j][k] = f_new[2][j][k];
            f_new[7][j][k] = f_new[5][j][k] + 0.5 * (f_new[1][j][k] - f_new[3][j][k]);
            f_new[8][j][k] = f_new[6][j][k] - 0.5 * (f_new[1][j][k] - f_new[3][j][k]);
            // 下壁面

            f_new[7][j][k] += 0.25 * dt * Fx[j][k];
            f_new[8][j][k] -= 0.25 * dt * Fx[j][k];
            // 考虑外力对边界的影响
        }

        for (int k = 1; k < n - 1; ++k)
        {
            int j = 0;
            double r_in = (f_new[0][j][k] + f_new[2][j][k] + f_new[4][j][k] + 2.0 * (f_new[3][j][k] + f_new[6][j][k] + f_new[7][j][k])) / (1.0 - u_in);

            f_new[1][j][k] = f_new[3][j][k] + (2.0 / 3.0) * r_in * u_in;
            f_new[5][j][k] = f_new[7][j][k] - 0.5 * (f_new[2][j][k] - f_new[4][j][k]) + (1.0 / 6.0) * r_in * u_in;
            f_new[8][j][k] = f_new[6][j][k] + 0.5 * (f_new[2][j][k] - f_new[4][j][k]) + (1.0 / 6.0) * r_in * u_in;
        }

        for (int k = 1; k < n - 1; ++k)
        {
            int j = m - 1;
            double u_out = -1.0 + (f_new[0][j][k] + f_new[2][j][k] + f_new[4][j][k] + 2.0 * (f_new[1][j][k] + f_new[5][j][k] + f_new[8][j][k])) / rho_out;

            f_new[3][j][k] = f_new[1][j][k] - (2.0 / 3.0) * rho_out * u_out;
            f_new[7][j][k] = f_new[5][j][k] + 0.5 * (f_new[2][j][k] - f_new[4][j][k]) - (1.0 / 6.0) * rho_out * u_out;
            f_new[6][j][k] = f_new[8][j][k] - 0.5 * (f_new[2][j][k] - f_new[4][j][k]) - (1.0 / 6.0) * rho_out * u_out;
        }

        for (int i = 0; i < 9; ++i)
        {

            for (int j = 0; j < m; ++j)
            {

                for (int k = 0; k < n; ++k)
                {
                    f[i][j][k] = f_new[i][j][k];
                }
            }
        }

        for (int j = 0; j < m; ++j)
        {
            for (int k = 0; k < n; ++k)
            {
                double r_temp = 0.0;
                double ux_temp = 0.0;
                double uy_temp = 0.0;

                for (int i = 0; i < 9; ++i)
                {
                    r_temp += f[i][j][k];
                    ux_temp += f[i][j][k] * ex[i];
                    uy_temp += f[i][j][k] * ey[i];
                }

                rho[j][k] = r_temp;

                u[j][k] = ux_temp / r_temp + 0.5 * Fx[j][k] / r_temp;
                v[j][k] = uy_temp / r_temp + 0.5 * Fy[j][k] / r_temp;

                psi[j][k] = 1.0 - exp(-rho[j][k]);
            }
        }

        if (step % 1000 == 0)
        {
            // write_vtk(step);
        }
    }
}
D2Q9_P::~D2Q9_P()
{

    std::cout << "Object destroyed: D1Q3 simulation resources have been released." << std::endl;
}