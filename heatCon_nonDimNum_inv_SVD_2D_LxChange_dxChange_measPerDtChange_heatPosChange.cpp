#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#define _USE_MATH_DEFINES

using namespace std;
using namespace Eigen;

#include "saveData.h"
#include "inputData.h"
#include "kronDelta.h"
#include "kronDeltaConst.h"
#include "kronDeltaConstInv.h"
#include "kronDeltaSoC.h"
#include "kronDeltaPulse.h"
#include "removeRowCol.h"

//--------Caution--------
// dx and dy is assumed to be equal.
// if dx and dy is not equal, we need to multiple dy or dx with W and M.
//--------Caution--------

int main()
{
    // Initialize inverse variables
    int i, k, j, l, n, q, measFreq, mCols, wRows, Nx_inv, Ny_inv, Nxy_inv, Nx_inv_INIT, Ny_inv_INIT, p_xNum, p_yNum, kappa_x_inv_skipNum, count_p_x, count_p_y, count_pn;
    int con, filter, rankMMt, heatSeek_iterNum, heatSeek_iterIndex, heatSeek_iterNumMax, heatSeek_iterIndexMax, iterStart, W_initialRows, W_finalRows, W_namaRows, MMtRows, iterBound, startIndex, startNum, iterJump;
    int nonDim_dx_iterIndex, nonDim_dx_iterIndexMax, filter_low, centreTemp;
    int measDx, measDy, measDt;
    int nonDim_tau_iterIndexMax, nonDim_tau_iterIndex;
    int measFreq_start;
    int gammaRow;
    int measPerDt_final, measPerDt_minus, measPerDt_max, measPerDt_iterIndex, measPerDt_iterIndexMax;
    int errTempMatrix_extract_iterIndex, errTempMatrix_xLen, initialLoop;
    int nonDim_Lx_iterIndex, nonDim_Lx_iterIndexMax;
    int inv_tempStart_idx;
    double dx_inv, dy_inv, dt_inv, forwardDt, startTime, measEnd, kappa_x_inv, kappa_y_inv, alpha, Lx, Ly, lambda, cp, rho;
    double bound_x0, bound_y0, bound_Nxm1, bound_Nym1, inputBound;
    double maxErr, rankPerc;
    double tau;
    double inv_tempStart, initialT_in;

    // Initialize FDS variables
    //*******Caution********
    // Heating position is different for even and odd Nx & Ny
    // When Nx & Ny is even, heating position can be in the middle
    // When Nx & Ny is odd, heating position cannot be in the middle
    //*******Caution********
    int Nx_FDS_woBound, Ny_FDS_woBound, Nx_FDS_wBound, Ny_FDS_wBound, Nxy_FDS_woBound, flatTimeIndex, flatXIndex, flatXEndIndex, timeCountIndex;
    int d_xaxis, d_yaxis, d_delay, d_SoC;
    int cond;
    int Nx_FDS, Ny_FDS, boundHosei_FDS, W_frontBackHosei, pointMulti;
    int heatPos_x_index, heatPos_y_index;
    double iterTime, iterPos_x, iterPos_y, endTime, endTimeIndex;
    double condTime, condSoC, heatPos_x, heatPos_y, Nx_SoC, R, kappa_x_FDS, kappa_y_FDS;
    double Ea_e, A_e, W_e, H_e, Ea_sei, A_sei, W_sei, H_sei, Ea_n, A_n, W_n, H_n, Ea_p, A_p, W_p, H_p, Con_SoC, A_SoC, Ea_SoC, R_b, V_SoC, C_SoC, eff_SoC, H_SoC, vol_SoC;
    double qdot_e, qdot_sei, qdot_n, qdot_p, qdot_total, qdot_SoC;
    double dx_FDS, dy_FDS, dt_FDS;

    // FDS Variables 1
    dt_FDS = 0.001;

    cout << "--------dx & tau Change Variables--------" << endl;
    cout << "dx change total iteration = " << endl;
    cin >> nonDim_dx_iterIndexMax;
    cout << "--------FDS Variables--------" << endl;
    cout << "Set iteration end time [s] = " << endl;
    cin >> endTime;
    cout << "Set condition time that heating starts [s] = " << endl;
    cin >> condTime;
    pointMulti = 1;
    cout << "--------Inverse Variables--------" << endl;
    cout << "Set 1st division number of x-axis space (bound to bound) [lowest is Nx=4] [ascending order] = " << endl;
    cin >> Nx_inv_INIT;
    cout << "Set 1st division number of y-axis space (bound to bound) [lowest is Ny=4] [ascending order] = " << endl;
    cin >> Ny_inv_INIT;

    // Iteration index settings
    nonDim_dx_iterIndex = 0;
    nonDim_tau_iterIndex = 0;
    n = 0;
    measPerDt_iterIndex = 0;
    measPerDt_minus = 0;
    errTempMatrix_xLen = 16;
    initialLoop = 0;

    // Other settings
    Nx_inv = Nx_inv_INIT;
    Ny_inv = Ny_inv_INIT;
    inv_tempStart = 313;
    initialT_in = 300;

    // Initialize Lx & Ly
    MatrixXd Lx_vec(7, 1);
    Lx_vec << 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0;
    MatrixXd Ly_vec(7, 1);
    Ly_vec = Lx_vec;

    MatrixXd Nx_SoC_vec(7, 7);
    Nx_SoC_vec << 
    155, 155, 155, 160, 160, 160, 160, 
    210, 210, 210, 210, 210, 210, 210,
    235, 235, 235, 235, 235, 235, 240,
    290, 290, 290, 290, 290, 290, 290,
    315, 315, 315, 315, 315, 315, 315,
    370, 370, 370, 370, 370, 370, 370,
    395, 395, 395, 395, 395, 395, 395;

    inputBound = 300; //[K]

    MatrixXi measFreq_vec(1,1);
    measFreq_vec << 24;

    nonDim_Lx_iterIndex = 0;
    nonDim_Lx_iterIndexMax = Lx_vec.rows();
    nonDim_tau_iterIndex = 0;
    nonDim_tau_iterIndexMax = measFreq_vec.rows();
    measFreq_start = measFreq_vec(0,0);

    while (nonDim_Lx_iterIndex < nonDim_Lx_iterIndexMax)
    {
        Lx = Lx_vec(nonDim_Lx_iterIndex, 0);
        Ly = Ly_vec(nonDim_Lx_iterIndex, 0);

        while (nonDim_dx_iterIndex < nonDim_dx_iterIndexMax)
        {
            Nx_SoC = Nx_SoC_vec(nonDim_Lx_iterIndex, nonDim_dx_iterIndex);

            // Inverse Variables
            heatPos_x_index = int(Nx_inv / 2);
            heatPos_y_index = int(Ny_inv / 2);
            measFreq = measFreq_start;
            dx_inv = Lx / Nx_inv;
            dy_inv = Ly / Ny_inv;
            dt_inv = 0.01;
            lambda = 398; // 45;
            cp = 389;     // 445;
            rho = 8960;   // 7850;
            alpha = lambda / (rho * cp);
            Nxy_inv = (Nx_inv - 1) * (Ny_inv - 1);
            kappa_x_inv_skipNum = Ny_inv - 1;
            p_xNum = 2 * (Ny_inv - 1);
            p_yNum = 2 * ((Nx_inv - 1) - 2);
            measPerDt_max = 2 * (Ny_inv - 1) + 2 * (Nx_inv - 3);
            kappa_x_inv = alpha * dt_inv / (dx_inv * dx_inv);
            kappa_y_inv = alpha * dt_inv / (dy_inv * dy_inv);
            measPerDt_iterIndexMax = measPerDt_max - 1;

            cout << "p_xNum = " << p_xNum << "\np_yNum = " << p_yNum << endl;

            // Set boundary temperature
            bound_x0 = inputBound;
            bound_Nxm1 = inputBound;
            bound_y0 = inputBound;
            bound_Nym1 = inputBound;

            // Forward variables 2
            Nx_FDS = pointMulti * Nx_inv;
            Ny_FDS = pointMulti * Ny_inv;
            dx_FDS = Lx / Nx_FDS;
            dy_FDS = Ly / Ny_FDS;

            // Inverse Variables 2
            measDt = dt_inv / dt_FDS;
            measDx = dx_inv / dx_FDS;
            measDy = dy_inv / dy_FDS;

// FDS starts
#include "heatCon_nonDimNum_FDS_2D_youkaihou.h"

            // Initialize A 1
            MatrixXd A = MatrixXd::Zero(Nxy_inv, Nxy_inv);
            MatrixXd An = MatrixXd::Zero(Nxy_inv, Nxy_inv);
            MatrixXd I_A = MatrixXd::Identity(Nxy_inv, Nxy_inv);

            // Initialize p
            MatrixXd p = MatrixXd::Zero(measPerDt_max, Nxy_inv);
            MatrixXd pn = MatrixXd::Zero(measPerDt_max, Nxy_inv);
            MatrixXd I_p = MatrixXd::Identity(Nxy_inv, measPerDt_max);
            MatrixXd p_initial = MatrixXd::Zero(measPerDt_max, Nxy_inv);

            // Initialize errTempMatrix & errTempMatrix_extract
            MatrixXd errTempMatrix = MatrixXd::Zero(measPerDt_iterIndexMax, errTempMatrix_xLen);
            MatrixXd errTempMatrix_extract = MatrixXd::Zero(measPerDt_iterIndexMax * nonDim_tau_iterIndexMax, errTempMatrix_xLen);
            errTempMatrix_extract_iterIndex = 0;

            // Measured vector W 1
            MatrixXd W_nama = load_csv<MatrixXd>("measForward_FDS_expHeating_Tflat_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_bound" + std::to_string(int(T_spaceMatrix(0, 0))) + "_Nx" + std::to_string(Nx_FDS) + "_Ny" + std::to_string(Ny_FDS) + "_NxSoC" + std::to_string(int(Nx_SoC)) + "_dt" + std::to_string(int(dt_FDS * 1000)) + "m.csv");
            W_namaRows = W_nama.rows();

            while (nonDim_tau_iterIndex < nonDim_tau_iterIndexMax)
            {
                //Set measFreq
                measFreq = measFreq_vec(nonDim_tau_iterIndex, 0);

                // Initialize W_initial
                VectorXd W_initial = VectorXd::Zero(measPerDt_max * measFreq);
                W_initialRows = W_initial.rows();

                while (measPerDt_iterIndex < measPerDt_iterIndexMax)
                {

                    // Set measPerDt_final
                    measPerDt_final = measPerDt_max - measPerDt_minus;

                    // Initialize W_final
                    VectorXd W_final = VectorXd::Zero(measPerDt_final * measFreq);
                    W_finalRows = W_final.rows();

                    // Iteration settings
                    if ( nonDim_tau_iterIndex == 0 && measPerDt_iterIndex == 0 )
                    {
                        startTime = inv_tempStart_idx*dt_FDS;
                        startNum = (startTime * measDt) / dt_inv;
                        startIndex = startNum;
                        iterBound = (W_namaRows - startNum) / measDt;
                        heatSeek_iterNum = 1;
                        iterJump = 1;
                        heatSeek_iterNumMax = startNum + heatSeek_iterNum * iterJump * measDt;
                        heatSeek_iterIndexMax = heatSeek_iterNumMax - 1;
                        heatSeek_iterIndex = startIndex;
                    }
                    else if (nonDim_tau_iterIndex > 0 || measPerDt_iterIndex > 0 )
                    {
                    }
                    else
                    {
                        cout << "Error at iteration settings" << endl;
                        return 0;
                    }

                    if (nonDim_tau_iterIndex == 0 && measPerDt_iterIndex == 0)
                    {
                        // A 2
                        for (i = 0; i < Nxy_inv; i++)
                        {
                            A(i, i) = 1 - 2 * kappa_x_inv - 2 * kappa_y_inv;
                        }
                        for (i = (Nx_inv - 1); i < Nxy_inv; i++)
                        {
                            A(i - (Nx_inv - 1), i) = kappa_y_inv;
                            A(i, i - (Nx_inv - 1)) = kappa_y_inv;
                        }
                        j = 1;
                        for (i = 1; i < Nxy_inv; i++)
                        {
                            if (j % kappa_x_inv_skipNum != 0)
                            {
                                A(i - 1, i) = kappa_x_inv;
                                A(i, i - 1) = kappa_x_inv;
                                j++;
                            }
                            else if (j % kappa_x_inv_skipNum == 0)
                            {
                                j++;
                            }
                            else
                            {
                                cout << "error at kappa_x subtitution" << endl;
                                return 0;
                            }
                        }
                        An = A * I_A;
                        
                        j = 0;
                        for (i = 0; i < p_xNum; i++)
                        {
                            // xNum in p
                            p(i, j) = 1;
                            i++;
                            p(i, j + (p_yNum / 2) + 1) = -1;
                            j += (p_yNum / 2) + 2;

                        }
                        cout << "\nCompleted p_x" << endl;

                        i = 0;
                        for (j = p_xNum; j < p_xNum + p_yNum; j++)
                        {
                            // yNum in p
                            p(j, i + 1) = 1;
                            p(j + 1, i + (Nxy_inv - 1) - (p_yNum / 2)) = -1;
                            j++;
                            i++;

                        }
                        cout << "\nCompleted p_y" << endl;

                        for (j = 0; j < measPerDt_max; j++)
                        {
                            for (i = 0; i < Nxy_inv; i++)
                            {
                                p_initial(j, i) = p(j, i);
                            }
                        }
                    }
                    else if (nonDim_tau_iterIndex > 0 || measPerDt_iterIndex > 0)
                    {
                        cout << "At A p loop, p rows = " << p.rows() << " p cols = " << p.cols() << endl;
                    }
                    else
                    {
                        cout << "Error at A & p iteration" << endl;
                        return 0;
                    }

                    // M
                    MatrixXd M = MatrixXd::Zero(measPerDt_final * measFreq, Nxy_inv);
                    MatrixXd Mt = MatrixXd::Zero(Nxy_inv, measPerDt_final * measFreq);
                    MatrixXd MMt = MatrixXd::Zero(measPerDt_final * measFreq, measPerDt_final * measFreq);

                    for (i = 0; i < measFreq; i++)
                    {

                        // Resizing p & pn
                        if (initialLoop > 0 && measPerDt_iterIndex == 0)
                        {
                            p.conservativeResize(measPerDt_max, Nxy_inv);
                            pn.conservativeResize(measPerDt_max, Nxy_inv);

                            // Reinitialize p & pn

                            for (j = 0; j < measPerDt_max; j++)
                            {
                                for (k = 0; k < Nxy_inv; k++)
                                {
                                    p(j, k) = p_initial(j, k);
                                }
                            }
                        }

                        pn = (p * An);
                        k = 0;

                        // Making M
                        for (j = (i * measPerDt_final); j < ((i + 1) * measPerDt_final); j++)
                        {
                            for (l = 0; l < Nxy_inv; l++)
                            {
                                M(j, l) = pn(k, l);
                            }
                            k++;
                        }
                        An = An * A;
                    }
                    An = A * I_A;

                    // remove row from p
                    removeRow(p, measPerDt_final - 1); // measPerDt_final-1 converts measPerDt to index
                    removeRow(pn, measPerDt_final - 1);

                    // M*Mt
                    cout << "Finding M and its characteristics..." << endl;
                    Mt = M.transpose();
                    MMt = M * Mt;
                    MMtRows = MMt.rows();

                    // Rank of MMt
                    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(MMt);
                    rankMMt = lu_decomp.rank();
                    cout << "\nRank of MMt = " << rankMMt << endl;

                    // MMt eigenvalues & eigenvectors
                    EigenSolver<MatrixXd> es(MMt);
                    MatrixXd gamma = es.pseudoEigenvalueMatrix();

                    // Checking if measFreq is suitable for size of MMt
                    gammaRow = gamma.rows();
                    if (gammaRow >= Nxy_inv)
                    {
                    }
                    else if (gammaRow < Nxy_inv)
                    {
                        nonDim_tau_iterIndex++;
                        continue;
                    }
                    else
                    {
                        cout << "Error at checking if measFreq is suitable for size of MMt" << endl;
                        return 0;
                    }

                    VectorXd gamma_show = VectorXd::Zero(Nxy_inv);
                    for (i = 0; i < Nxy_inv; i++)
                    {
                        gamma_show(i) = gamma(i, i);
                    }

                    // Find w of W*Gamma*Vt
                    MatrixXd w = es.pseudoEigenvectors();
                    wRows = w.rows();

                    // Auto calc filter
                    filter = rankMMt;
                    filter_low = rankMMt;
                    cout << "Calculating filter..." << endl;
                    if (filter == Nxy_inv)
                    {
                    }
                    else if (filter < Nxy_inv)
                    {
                        for (i = 0; i < (Nxy_inv - rankMMt); i++)
                        {
                            if (gamma(filter_low + i, filter_low + i) > 0)
                            {
                                filter++;
                            }
                            else if (gamma(filter_low + i, filter_low + i) < 0)
                            {
                                break;
                            }
                            else
                            {
                                cout << "error at Auto calc filter inner if" << endl;
                                return 0;
                            }
                        }
                    }
                    else
                    {
                        cout << "error at Auto calc filter outer if" << endl;
                        return 0;
                    }
                    cout << "Done!\nfilter = " << filter << endl;

                    // Find vt of W*Gamma*Vt
                    MatrixXd vt = MatrixXd::Zero(Nxy_inv, Nxy_inv);

                    for (i = 0; i < filter; i++)
                    {
                        MatrixXd vRow = (w.col(i).transpose() * M) / pow(gamma(i, i), 0.5);
                        for (j = 0; j < Nxy_inv; j++)
                        {
                            vt(i, j) = vRow(0, j);
                        }
                    }

                    MatrixXd VFinal = MatrixXd::Zero(Nxy_inv, heatSeek_iterNum);
                    MatrixXd VFinal_nama = MatrixXd::Zero(Nxy_inv, heatSeek_iterNum);
                    MatrixXd W_extract = MatrixXd::Zero(measPerDt_max * measFreq, heatSeek_iterNum * measPerDt_max);
                    VectorXd initialT = VectorXd::Zero(Nxy_inv);
                    VectorXd W_bound = VectorXd::Zero(measPerDt_max * measFreq);

                    // Record initial Temp (not supposed to have in inverse analysis)
                    k = 0;
                    l = 0;
                    for (i = 0; i < Ny_inv - 1; i++)
                    {
                        l += (measDx - 1) + Nx_FDS_woBound * (measDy - 1);
                        for (j = 0; j < Nx_inv - 1; j++)
                        {
                            initialT(k) = W_nama(startIndex, l);
                            l += measDx;
                            k++;
                        }
                    }
                    
                    while (heatSeek_iterIndex < heatSeek_iterIndexMax)
                    {

                        cout << "measurement data W starts at heatSeek_iterIndex of " << heatSeek_iterIndex + 1 << endl;
                        cout << "iteration --------> [" << heatSeek_iterIndex + 1 << "/" << heatSeek_iterIndexMax << "]" << endl;
                        cout << "loop n[" << n << "]" << endl;

                        if (measPerDt_iterIndex == 0)
                        {
                            // Insert measurement data into w
                            //--------Caution--------
                            // This is only for 2D space with boundary not equal to 0.
                            // And gradient of measurement data in corners is taken as the gradient over the x axis not y axis
                            //--------Caution--------
                            k = startIndex;
                            W_frontBackHosei = Ny_FDS_woBound * measDy * (Ny_inv - 2);
                            for (i = 0; i < measPerDt_final * measFreq; i++)
                            {
                                l = measDx - 1; // change position at x-axis

                                // All x-axis measurements
                                for (j = 0; j < p_xNum; j++)
                                {
                                    l += Ny_FDS_woBound * (measDy - 1); // change position at y-axis
                                    W_initial(i) = W_nama(k, l) - bound_x0; /// dx;
                                    i++;
                                    l += measDx * ((p_yNum / 2) + 1);
                                    W_initial(i) = bound_Nxm1 - W_nama(k, l); /// dx;
                                    i++;
                                    j++;
                                    l += 2 * (measDx - 1) + 1; // change postion at x-axis 2
                                }

                                // All y-axix measurements w/o corners
                                l = (measDx - 1);                   // change position at x-axis
                                l += Ny_FDS_woBound * (measDy - 1); // change position at y-axis
                                for (j = 0; j < p_yNum; j++)
                                {
                                    l += measDx; // change position to yNum
                                    W_initial(i) = W_nama(k, l) - bound_y0; /// dy;
                                    i++;
                                    W_initial(i) = bound_Nym1 - W_nama(k, l + W_frontBackHosei); /// dy;
                                    i++;
                                    j++;
                                }
                                i -= 1;
                                k += measDt * iterJump;
                            }
                            cout << "Completed extracting measurement data into W_initial.\nNow copying into W_final." << endl;

                            cout << "W_final rows = " << W_finalRows << " W_initial rows = " << W_initialRows << endl;
                            for (i = 0; i < (measPerDt_max * measFreq); i++)
                            {
                                W_final(i) = W_initial(i);
                            }
                        }
                        else if (measPerDt_iterIndex > 0)
                        {

                            cout << "Extracting data to fit measPerDt from W_initial to W_final" << endl;
                            q = 0;
                            for (i = 0; i < measFreq; i++)
                            {
                                for (j = 0; j < measPerDt_final; j++)
                                {
                                    W_final(q) = W_initial(j + i * measPerDt_max);
                                    q++;
                                }
                            }
                        }
                        else
                        {
                            cout << "Error at extracting data if loop" << endl;
                            return 0;
                        }

                        // Check dimensions of MMt.rows and W_final.rows
                        if (W_finalRows == MMtRows)
                        {
                        }
                        else if (W_finalRows != MMtRows)
                        {
                            cout << "MMt rows and W_final rows are not equal" << endl;
                            cout << "MMt rows = " << MMtRows << " W_final rows = " << W_finalRows << endl;
                            return 0;
                        }
                        else
                        {
                            cout << "Error at checking dimensions of MMt and W_final." << endl;
                            return 0;
                        }

                        // Copying W_final to W_extract
                        for (i = 0; i < measPerDt_final * measFreq; i++)
                        {
                            W_extract(i, (n * measPerDt_max) + measPerDt_iterIndex) = W_final(i);
                        }

                        // Find coefficient b of w
                        MatrixXd b = (w.inverse()) * W_final;

                        VectorXd d = VectorXd::Zero(Nxy_inv);

                        for (i = 0; i < filter; i++)
                        {
                            d(i) = b(i) / pow(gamma(i, i), 0.5);
                        }

                        // Find final V

                        MatrixXd v = vt.transpose();
                        MatrixXd V = d(0) * v.col(0);

                        for (i = 0; i < filter - 1; i++)
                        {
                            V += d(i + 1) * v.col(i + 1);
                        }

                        for (i = 0; i < Nxy_inv; i++)
                        {
                            VFinal_nama(i, n) = V(i, 0);
                        }

                        for (i = 0; i < Nxy_inv; i++)
                        {
                            VFinal(i, n) = VFinal_nama(i, n) + inputBound;
                        }

                        heatSeek_iterIndex += measDt * iterJump;
                        cout << "iterMax = " << heatSeek_iterNumMax << endl;
                        cout << "heatSeek_iterIndex = " << heatSeek_iterIndex << endl;
                        n++;
                    }

                    heatSeek_iterIndex = startIndex;
                    n = 0;
                    cout << "Final V with start time at " << startTime << " [s] and iter of " << heatSeek_iterNum << "." << endl;

                    saveData("inv_VFinal_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) +  "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", VFinal);
                    saveData("inv_Wextract_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv * 1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", W_extract);
                    saveData("inv_initialTemp_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", initialT);

                    // Error Calc
                    // errMatrix { dx1, err1; dx2, err2; .... }
                    VectorXd errVector;
                    if (heatSeek_iterNum == 1)
                    {
                        errVector = VFinal - initialT;
                    }
                    else if (heatSeek_iterNum > 1)
                    {
                        errVector = VFinal.row(n) - initialT;
                    }
                    else
                    {
                        cout << "Error at errVector" << endl;
                        return 0;
                    }

                    errVector = errVector.cwiseAbs();

                    maxErr = errVector.maxCoeff();
                    rankPerc = 1.0 * rankMMt / Nxy_inv;
                    tau = measFreq * dt_inv;

                    errTempMatrix(measPerDt_iterIndex, 0) = dt_inv;
                    errTempMatrix(measPerDt_iterIndex, 1) = dx_inv;
                    errTempMatrix(measPerDt_iterIndex, 2) = Nx_inv;
                    errTempMatrix(measPerDt_iterIndex, 3) = alpha;
                    errTempMatrix(measPerDt_iterIndex, 4) = measFreq;
                    errTempMatrix(measPerDt_iterIndex, 5) = tau;
                    errTempMatrix(measPerDt_iterIndex, 6) = measPerDt_max;
                    errTempMatrix(measPerDt_iterIndex, 7) = measPerDt_final;
                    errTempMatrix(measPerDt_iterIndex, 8) = measPerDt_final * measFreq;
                    errTempMatrix(measPerDt_iterIndex, 9) = rankMMt;
                    errTempMatrix(measPerDt_iterIndex, 10) = rankPerc;
                    errTempMatrix(measPerDt_iterIndex, 11) = maxErr;
                    errTempMatrix(measPerDt_iterIndex, 12) = heatPos_x_index;
                    errTempMatrix(measPerDt_iterIndex, 13) = Nx_inv - 1;
                    errTempMatrix(measPerDt_iterIndex, 14) = startTime;
                    errTempMatrix(measPerDt_iterIndex, 15) = startNum;

                    cout << "tauIndex [" << nonDim_tau_iterIndex << "/" << nonDim_tau_iterIndexMax << "]" << endl;
                    cout << "measPerDt_iterIndex [" << measPerDt_iterIndex << "/" << measPerDt_iterIndexMax << "]" << endl;
                    cout << "measFreq for Nx[" << Nx_inv << "] = " << measFreq << endl;
                    cout << "tau for for Nx[" << Nx_inv << "] = " << tau << endl;
                    cout << "measPerDt for Nx[" << Nx_inv << "] = " << measPerDt_final << "/" << measPerDt_max << endl;
                    cout << "\nrankMMt% = " << rankPerc << endl;
                    cout << "Max error for Nx[" << Nx_inv << "] = " << maxErr << endl;

                    measPerDt_iterIndex++;
                    measPerDt_minus++;
                    initialLoop++;

                    // check whether measPerDt_final fits
                    if ((measPerDt_max - measPerDt_minus) * measFreq > Nxy_inv)
                    {
                    }
                    else if ((measPerDt_max - measPerDt_minus) * measFreq <= Nxy_inv)
                    {
                        cout << "Breaking measPerDt loop, going next loop for tau" << endl;
                        break;
                    }
                    else
                    {
                        cout << "Error at checking measPerDt_final fittings" << endl;
                        return 0;
                    }
                }

                cout << "Importing into errTempMatrix_extract for tauIndex [" << nonDim_tau_iterIndex << "/" << nonDim_tau_iterIndexMax << "]..." << endl;
                
                for (j = 0; j < measPerDt_iterIndexMax; j++)
                {
                    for (i = 0; i < errTempMatrix_xLen; i++)
                    {
                        errTempMatrix_extract(errTempMatrix_extract_iterIndex, i) = errTempMatrix(j, i);
                    }
                    errTempMatrix_extract_iterIndex++;
                }
                cout << "Importing ends!" << endl;

                nonDim_tau_iterIndex++;
                measPerDt_iterIndex = 0;
                measPerDt_final = measPerDt_max;
                measPerDt_minus = 0;
            }
            cout << "errTempMatrix for dx[" << dx_inv << "] Nx[" << Nx_inv << "] tau[" << tau << "measFreq[" << measFreq << "] = \n"
                 << errTempMatrix_extract << endl;
            saveData("inv_errTempMatrix_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measPerDtMax" + std::to_string(measPerDt_max) + "_dt" + std::to_string(int(dt_inv * 1000)) + "m_ptM" + std::to_string(pointMulti) + "_dxChange.csv", errTempMatrix_extract);

            nonDim_dx_iterIndex++;
            Nx_inv++;
            Ny_inv++;
            nonDim_tau_iterIndex = 0;
        }
        nonDim_dx_iterIndex = 0;
        Nx_inv = Nx_inv_INIT;
        Ny_inv = Ny_inv_INIT;
        nonDim_Lx_iterIndex++;
    }

    return 0;
}
