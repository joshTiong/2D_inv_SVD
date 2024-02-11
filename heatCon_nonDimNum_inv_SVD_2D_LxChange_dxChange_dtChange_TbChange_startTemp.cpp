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
    int nonDim_tau_iterIndexMax, nonDim_tau_iterIndex, nonDim_dt_iterIndexMax, nonDim_dt_iterIndex;
    int gammaRow;
    int measPerDt_final, measPerDt_minus, measPerDt_max, measPerDt_iterIndex, measPerDt_iterIndexMax;
    int errTempMatrix_extract_iterIndex, errTempMatrix_xLen, initialLoop;
    int nonDim_Lx_iterIndex, nonDim_Lx_iterIndexMax, nonDim_inputBound_iterIndex, nonDim_inputBound_iterIndexMax;
    int inv_tempStart_idx;
    double dx_inv, dy_inv, forwardDt, startTime, measEnd, kappa_x_inv, kappa_y_inv, alpha, Lx, Ly, lambda, cp, rho;
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

    cout << "--------Common Variables--------" << endl;
    // cout << "Set length of 1st x-axis space [m] = " << endl;
    // cin >> Lx;
    // cout << "Set length of 1st y-axis space [m] = " << endl;
    // cin >> Ly;
    // cout << "Set boundary temperature = " << endl;
    // cin >> inputBound;
    cout << "--------dx & tau Change Variables--------" << endl;
    cout << "dx change total iteration = " << endl;
    cin >> nonDim_dx_iterIndexMax;
    //cout << "tau change total iteration = " << endl;
    //cin >> nonDim_tau_iterIndexMax;
    cout << "--------FDS Variables--------" << endl;
    //cout << "Set iteration end time [s] = " << endl;
    //cin >> endTime;
    endTime = 20;
    //cout << "Set condition time that heating starts [s] = " << endl;
    //cin >> condTime;
    condTime = 0;
    // cout << "Set short circuit volume divider of " << Lx << "[m3] = " << endl;
    // cin >> Nx_SoC;
    //cout << "Set point multiplier = " << endl;
    //cin >> pointMulti;
    pointMulti = 1;
    cout << "--------Inverse Variables--------" << endl;
    //cout << "Set 1st measurement frequency for one iteration of inverse analysis [ascreasing order] = " << endl;
    //cin >> measFreq_start;
    cout << "Set 1st division number of x-axis space (bound to bound) [lowest is Nx=4] [ascending order] = " << endl;
    cin >> Nx_inv_INIT;
    cout << "Set 1st division number of y-axis space (bound to bound) [lowest is Ny=4] [ascending order] = " << endl;
    cin >> Ny_inv_INIT;

    // Iteration index settings
    nonDim_dx_iterIndex = 0;
    n = 0;
    measPerDt_iterIndex = 0;
    measPerDt_minus = 0;
    errTempMatrix_xLen = 16;
    initialLoop = 0;

    // Other settings
    Nx_inv = Nx_inv_INIT;
    Ny_inv = Ny_inv_INIT;
    inv_tempStart = 288;
    initialT_in = 275;

    // Initialize Lx & Ly
    MatrixXd Lx_vec(6, 1);
    Lx_vec << 0.08, 0.1, 0.12, 0.14, 0.16, 0.18;
    //0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0;
    //1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0
    MatrixXd Ly_vec(6, 1);
    Ly_vec = Lx_vec;

    //MatrixXd Nx_vec(3,1);
    //Nx_vec << 8, 9, 10;
    //MatrixXd Ny_vec(3,1);
    //Ny_vec = Nx_vec;

    //MatrixXd Nx_SoC_vec(7, 7);
    //Nx_SoC_vec << 
    //155, 155, 155, 160, 160, 160, 160, 
    //210, 210, 210, 210, 210, 210, 210,
    //235, 235, 235, 235, 235, 235, 240,
    //290, 290, 290, 290, 290, 290, 290,
    //315, 315, 315, 315, 315, 315, 315,
    //370, 370, 370, 370, 370, 370, 370,
    //395, 395, 395, 395, 395, 395, 395;

    //MatrixXd Nx_SoC_vec(15, 3);
    //Nx_SoC_vec <<  48, 50, 50,
    //89, 90, 90, 
    //130, 130, 130, 
    //170, 170, 170,
    //210, 200, 200,
    //250, 250, 250,
    //300, 300, 300,
    //350, 350, 350,
    //400, 390, 400,
    //440, 440, 450,
    //460, 460, 460,
    //500, 500, 500,
    //550, 550, 550,
    //600, 600, 600,
    //640, 640, 640;

    //MatrixXd Nx_SoC_vec(6, 3);
    //Nx_SoC_vec << 26, 28, 28,
    //30, 31, 32,
    //33, 35, 36,
    //40, 40, 40,
    //42, 42, 43, 
    //44, 44, 46;

    MatrixXd Nx_SoC_vec(6, 3);
    Nx_SoC_vec << 73, 78, 82,
    80, 86, 90,
    96, 96, 98,
    100, 100, 104,
    102, 106, 112, 
    106, 112, 118;

    MatrixXd inputBound_vec(1, 1);
    inputBound_vec << 275; //[K]

    MatrixXd dt_inv(2, 1);
    dt_inv << 0.01, 0.02;

    MatrixXi measFreq_vec(5,1);
    measFreq_vec << 12, 15, 18, 21, 24;

    nonDim_Lx_iterIndex = 0;
    nonDim_Lx_iterIndexMax = Lx_vec.rows();
    nonDim_dt_iterIndex = 0;
    nonDim_dt_iterIndexMax = dt_inv.rows();
    nonDim_inputBound_iterIndex = 0;
    nonDim_inputBound_iterIndexMax = inputBound_vec.rows();
    nonDim_tau_iterIndex = 0;
    nonDim_tau_iterIndexMax = measFreq_vec.rows();
    //nonDim_dx_iterIndexMax = Nx_vec.rows();

    // Permission to continue start
    //     cout << "Lx_vec = \n" << Lx_vec << endl;
    //     cout << "Ly_vec = \n" << Ly_vec << endl;
    //     cout << "Continue? [1]Yes [2]End" << endl;
    //     cin >> con;
    //
    //     if ( con == 1){
    //         cout << "Continuing..." << endl;
    //     }
    //     else if (con == 2){
    //         cout << "Ending..." << endl;
    //         return 0;
    //     }
    //     else {
    //         cout << "Error at choice before M" << endl;
    //         return 0;
    //     }
    // Permission to continue end

    while (nonDim_Lx_iterIndex < nonDim_Lx_iterIndexMax)
    {
        Lx = Lx_vec(nonDim_Lx_iterIndex, 0);
        Ly = Ly_vec(nonDim_Lx_iterIndex, 0);
        // heatPos_x_even = Lx / 2;
        // heatPos_y_even = Ly / 2;

        while (nonDim_dx_iterIndex < nonDim_dx_iterIndexMax)
        {
            heatPos_x_index = int(Nx_inv / 2);
            heatPos_y_index = int(Ny_inv / 2);
            Nx_SoC = Nx_SoC_vec(nonDim_Lx_iterIndex, nonDim_dx_iterIndex);

            while (nonDim_dt_iterIndex < nonDim_dt_iterIndexMax)
            {
                while (nonDim_inputBound_iterIndex < nonDim_inputBound_iterIndexMax)
                {
                    // Inverse Variable
                    dx_inv = Lx / Nx_inv;
                    dy_inv = Ly / Ny_inv;
                    // dt_inv = 0.01;
                    lambda = 398; // 45;
                    cp = 389;     // 445;
                    rho = 8960;   // 7850;
                    alpha = lambda / (rho * cp);
                    Nxy_inv = (Nx_inv - 1) * (Ny_inv - 1);
                    kappa_x_inv_skipNum = Ny_inv - 1;
                    p_xNum = 2 * (Ny_inv - 1);
                    p_yNum = 2 * ((Nx_inv - 1) - 2);
                    measPerDt_max = 2 * (Ny_inv - 1) + 2 * (Nx_inv - 3);
                    kappa_x_inv = alpha * dt_inv(nonDim_dt_iterIndex, 0) / (dx_inv * dx_inv);
                    kappa_y_inv = alpha * dt_inv(nonDim_dt_iterIndex, 0) / (dy_inv * dy_inv);
                    measPerDt_iterIndexMax = measPerDt_max - 1;

                    cout << "p_xNum = " << p_xNum << "\np_yNum = " << p_yNum << endl;

                    // Set boundary temperature
                    inputBound = inputBound_vec(nonDim_inputBound_iterIndex, 0);
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
                    measDt = dt_inv(nonDim_dt_iterIndex, 0) / dt_FDS;
                    measDx = dx_inv / dx_FDS;
                    measDy = dy_inv / dy_FDS;

                    // FDS starts
                    if (nonDim_dt_iterIndex == 0)
                    {
#include "heatCon_nonDimNum_FDS_2D_youkaihou.h"
                    }
                    else if (nonDim_dt_iterIndex > 0)
                    {
                        // skip FDS
                    }
                    else
                    {
                        cout << "Error at FDS" << endl;
                        return 0;
                    }

                    // Initialize A 1
                    MatrixXd A = MatrixXd::Zero(Nxy_inv, Nxy_inv);
                    MatrixXd An = MatrixXd::Zero(Nxy_inv, Nxy_inv);
                    MatrixXd I_A = MatrixXd::Identity(Nxy_inv, Nxy_inv);

                    // Initialize p
                    MatrixXd p = MatrixXd::Zero(measPerDt_max, Nxy_inv);
                    MatrixXd pn = MatrixXd::Zero(measPerDt_max, Nxy_inv);
                    MatrixXd I_p = MatrixXd::Identity(Nxy_inv, measPerDt_max);
                    MatrixXd p_initial = MatrixXd::Zero(measPerDt_max, Nxy_inv);
                    // cout << "p rows = " << p.rows() << " cols = "<< p.cols() << endl;
                    // cout << "I_p rows = " << I_p.rows() << " cols = " << I_p.cols() << endl;

                    // Initialize errTempMatrix & errTempMatrix_extract
                    MatrixXd errTempMatrix = MatrixXd::Zero(measPerDt_iterIndexMax, errTempMatrix_xLen);
                    MatrixXd errTempMatrix_extract = MatrixXd::Zero(measPerDt_iterIndexMax * nonDim_tau_iterIndexMax, errTempMatrix_xLen);
                    errTempMatrix_extract_iterIndex = 0;

                    // Measured vector W 1
                    MatrixXd W_nama = load_csv<MatrixXd>("measForward_FDS_expHeating_Tflat_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_bound" + std::to_string(int(inputBound)) + "_Nx" + std::to_string(Nx_FDS) + "_Ny" + std::to_string(Ny_FDS) + "_NxSoC" + std::to_string(int(Nx_SoC)) + "_dt" + std::to_string(int(dt_FDS * 1000)) + "m.csv");
                    W_namaRows = W_nama.rows();

                    while (nonDim_tau_iterIndex < nonDim_tau_iterIndexMax)
                    {
                        //Set measFreq
                        measFreq = measFreq_vec(nonDim_tau_iterIndex, 0);

                        // Initialize W_initial
                        VectorXd W_initial = VectorXd::Zero(measPerDt_max * measFreq);
                        W_initialRows = W_initial.rows();

                        // while(measPerDt_iterIndex < measPerDt_iterIndexMax){

                        // Set measPerDt_final
                        measPerDt_final = measPerDt_max;

                        // Initialize W_final
                        VectorXd W_final = VectorXd::Zero(measPerDt_final * measFreq);
                        W_finalRows = W_final.rows();

                        // Iteration settings
                        if (nonDim_tau_iterIndex == 0 && measPerDt_iterIndex == 0 && nonDim_dt_iterIndex == 0)
                        {
                            //cout << "Inverse analysis start time [s] larger than " << dt_inv(nonDim_dt_iterIndex, 0) << "s = " << endl;
                            //cin >> startTime;
                            startTime = inv_tempStart_idx*dt_FDS;
                            startNum = (startTime * measDt) / dt_inv(nonDim_dt_iterIndex, 0);
                            startIndex = startNum;// - 1;
                            iterBound = (W_namaRows - startNum) / measDt;
                            //cout << "Iteration number count within " << iterBound << " = " << endl;
                            //cin >> heatSeek_iterNum;
                            heatSeek_iterNum = 1;
                            //cout << "Iteration jump = " << endl;
                            //cin >> iterJump;
                            iterJump = 1;
                            heatSeek_iterNumMax = startNum + heatSeek_iterNum * iterJump * measDt;
                            heatSeek_iterIndexMax = heatSeek_iterNumMax - 1;
                            heatSeek_iterIndex = startIndex;
                        }
                        else if (nonDim_tau_iterIndex > 0 || measPerDt_iterIndex > 0 || nonDim_dt_iterIndex > 0)
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
                            //saveData("inv_A_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0) * 1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", A);

                            // count_p_x = 0;
                            // count_p_y = 0;
                            j = 0;
                            for (i = 0; i < p_xNum; i++)
                            {
                                // xNum in p
                                p(i, j) = 1;
                                i++;
                                p(i, j + (p_yNum / 2) + 1) = -1;
                                j += (p_yNum / 2) + 2;

                                // Loop count
                                // if(count_p_x%1 != 0){
                                // }
                                // else if(count_p_x%1 == 0){
                                //         cout << "Loop at p_x[" << count_p_x << "]" << endl;
                                // }
                                // else{
                                //         cout << "error at loop p_x" << endl;
                                //         return 0;
                                // }

                                // count_p_x++;
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

                                // Loop count
                                // if(count_p_y%1 != 0){
                                // }
                                // else if(count_p_y%1 == 0){
                                //         cout << "Loop at p_y[" << count_p_y << "]" << endl;
                                // }
                                // else{
                                //         cout << "error at loop p_y" << endl;
                                //	  return 0;
                                // }

                                // count_p_y++;
                            }
                            cout << "\nCompleted p_y" << endl;

                            // saveData("inv_p_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", p);

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

                                // cout << "Resizing p & pn to (" << Nxy_inv << "," << measPerDt_max << ")" << endl;
                                p.conservativeResize(measPerDt_max, Nxy_inv);
                                pn.conservativeResize(measPerDt_max, Nxy_inv);

                                // Reinitialize p & pn
                                // MatrixXd p = MatrixXd::Zero(measPerDt_max, Nxy_inv);
                                // MatrixXd pn = MatrixXd::Zero(measPerDt_max, Nxy_inv);

                                for (j = 0; j < measPerDt_max; j++)
                                {
                                    for (k = 0; k < Nxy_inv; k++)
                                    {
                                        // cout << "j = " << j << " k = " << k << endl;
                                        p(j, k) = p_initial(j, k);
                                    }
                                }
                            }

                            pn = (p * An);
                            k = 0;

                            // if(i == 0){
                            //	saveData("inv_pn_i[" + std::to_string(i) + "]_" + std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt001.csv", pn);
                            //	saveData("inv_p_i[" + std::to_string(i) + "]_" + std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt001.csv", p);
                            //	saveData("inv_p_initial_i[" + std::to_string(i) + "]_" + std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt001.csv", p_initial);
                            // }

                            // cout << "Buidling M for timestep [" << i + 1 << "/" << measFreq << "]..." << endl;
                            // cout << "with p (rows,cols) = " << "(" << p.rows() << "," << p.cols() << ")" << endl;
                            // Making M
                            for (j = (i * measPerDt_final); j < ((i + 1) * measPerDt_final); j++)
                            {
                                for (l = 0; l < Nxy_inv; l++)
                                {
                                    M(j, l) = pn(k, l);
                                }
                                k++;
                            }
                            // Loop count
                            // if(count_p_y%1 != 0){
                            // }
                            // else if(count_p_y%1 == 0){
                            //         cout << "Loop at pn[" << count_pn << "]" << endl;
                            // }
                            // else{
                            //         cout << "error at loop pn" << endl;
                            //         return 0;
                            // }

                            // count_pn++;
                            An = An * A;
                        }
                        An = A * I_A;
                        // saveData("inv_M_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt001.csv", M);

                        // remove row from p
                        // removeRow(p, measPerDt_final-1);//measPerDt_final-1 converts measPerDt to index
                        // removeRow(pn, measPerDt_final-1);

                        // M*Mt
                        cout << "Finding M and its characteristics..." << endl;
                        Mt = M.transpose();
                        MMt = M * Mt;
                        MMtRows = MMt.rows();
                        // saveData("inv_MMt_2D_"+ std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_iterNum" + std::to_string(heatSeek_iterNum) + "_Nx" + std::to_string(Nx) + "_Ny" + std::to_string(Ny) + "_measDt" + std::to_string(int(measDt)) + "_measPerDt" + std::to_string(measPerDt_final) + "_measFreq" + std::to_string(measFreq) + "_dt001.csv", MMt);

                        // cout << "MMt rows = " << MMt.rows() << endl;
                        // cout << "MMt cols = " << MMt.cols() << endl;

                        // Rank of MMt
                        Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(MMt);
                        rankMMt = lu_decomp.rank();
                        cout << "\nRank of MMt = " << rankMMt << endl;

                        // MMt eigenvalues & eigenvectors
                        EigenSolver<MatrixXd> es(MMt);
                        MatrixXd gamma = es.pseudoEigenvalueMatrix();

                        // saveData("inv_gamma_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", gamma);

                        // Checking if measFreq is suitable for size of MMt
                        gammaRow = gamma.rows();
                        if (gammaRow >= Nxy_inv)
                        {
                        }
                        else if (gammaRow < Nxy_inv)
                        {
                            nonDim_tau_iterIndex++;
                            measFreq++;
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
                        // cout << "Singular Value for MMt = \n" << gamma_show << endl;

                        // Find w of W*Gamma*Vt
                        MatrixXd w = es.pseudoEigenvectors();
                        wRows = w.rows();

                        // saveData("inv_w_2D_"+ std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_iterNum" + std::to_string(heatSeek_iterNum) + "_Nx" + std::to_string(Nx) + "_Ny" + std::to_string(Ny) + "_measDt" + std::to_string(int(measDt)) + "_measPerDt" + std::to_string(measPerDt_final) + "_measFreq" + std::to_string(measFreq) + ".csv", w);

                        // cout << "w rows = " << w.rows() << endl;
                        // cout << "w cols = " << w.cols() << endl;

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

                        // cout << "vt = \n" << vt << endl;
                        cout << "Here!" << endl;

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
                                // cout << "i = " << i << " j = " << j << " k = " << k << " l = " << l << endl;
                                initialT(k) = W_nama(startIndex, l);
                                l += measDx;
                                k++;
                            }
                        }
                        // cout << "Initial temp recording ends" << endl;
                        // saveData("inv_initialTemp_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", initialT);

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
                                // W_frontBackHosei = (measDx-1) + (measDx*3) + measDx;
                                W_frontBackHosei = Ny_FDS_woBound * measDy * (Ny_inv - 2);
                                for (i = 0; i < measPerDt_final * measFreq; i++)
                                {
                                    l = measDx - 1; // change position at x-axis

                                    // All x-axis measurements
                                    for (j = 0; j < p_xNum; j++)
                                    {
                                        l += Ny_FDS_woBound * (measDy - 1); // change position at y-axis
                                        // cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << endl;
                                        W_initial(i) = W_nama(k, l) - bound_x0; /// dx;
                                        i++;
                                        l += measDx * ((p_yNum / 2) + 1);
                                        // cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << endl;
                                        W_initial(i) = bound_Nxm1 - W_nama(k, l); /// dx;
                                        i++;
                                        j++;
                                        l += 2 * (measDx - 1) + 1; // change postion at x-axis 2
                                        // saveData("inv_Wx[" + std::to_string(i) + "][" + std::to_string(j) +"]_2D.csv", W);
                                    }

                                    // All y-axix measurements w/o corners
                                    l = (measDx - 1);                   // change position at x-axis
                                    l += Ny_FDS_woBound * (measDy - 1); // change position at y-axis
                                    for (j = 0; j < p_yNum; j++)
                                    {
                                        l += measDx; // change position to yNum
                                        // cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l << endl;
                                        W_initial(i) = W_nama(k, l) - bound_y0; /// dy;
                                        i++;
                                        // cout << "i = " << i << ", j = " << j << ", k = " << k << ", l = " << l+W_frontBackHosei << endl;
                                        W_initial(i) = bound_Nym1 - W_nama(k, l + W_frontBackHosei); /// dy;
                                        i++;
                                        j++;
                                        // saveData("inv_Wy[" + std::to_string(i) + "][" + std::to_string(j) + "]_2D.csv", W);
                                    }
                                    i -= 1;
                                    k += measDt * iterJump;
                                    // cout << "i = " << i << "/" << measPerDt*measFreq << endl;
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

                            // Correction due to boundary data
                            // for(i=0; i<measFreq; i++){
                            //	for(j=0; j<p_xNum; j++){
                            //		W_bound(j + i*measPerDt) = -bound_x0*kappa_x_inv;
                            //		j++;
                            //		W_bound(j + i*measPerDt) = bound_Nxm1*kappa_x_inv;
                            //	}
                            //
                            //	for(j=0; j<p_yNum; j++){
                            //		W_bound(j + p_xNum + i*measPerDt) = -bound_y0*kappa_y_inv;
                            //		j++;
                            //		W_bound(j + p_xNum + i*measPerDt) = bound_Nym1*kappa_y_inv;			}
                            // }

                            // saveData("inv_WBound_2D_"+ std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_iterNum" + std::to_string(heatSeek_iterNum) + "_Nx" + std::to_string(Nx) + "_Ny" + std::to_string(Ny) + "_measDt" + std::to_string(int(measDt)) + "_measPerDt" + std::to_string(measPerDt) + "_measFreq" + std::to_string(measFreq) + "_dt001.csv", W_bound);

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
                            // cout << "b = \n" << b << endl;

                            VectorXd d = VectorXd::Zero(Nxy_inv);

                            for (i = 0; i < filter; i++)
                            {
                                d(i) = b(i) / pow(gamma(i, i), 0.5);
                            }
                            // cout << "d = \n" << d << endl;

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

                            // cout << "V at n[" << n << "] = \n" << VFinal << endl;

                            heatSeek_iterIndex += measDt * iterJump;
                            cout << "iterMax = " << heatSeek_iterNumMax << endl;
                            cout << "heatSeek_iterIndex = " << heatSeek_iterIndex << endl;
                            n++;
                        }

                        heatSeek_iterIndex = startIndex;
                        n = 0;
                        cout << "Final V with start time at " << startTime << " [s] and iter of " << heatSeek_iterNum << "." << endl;

                        saveData("inv_VFinal_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0))  + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0)*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", VFinal);
                        saveData("inv_Wextract_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0) * 1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", W_extract);
                        saveData("inv_initialTemp_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0)*1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", initialT);

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

                        saveData("inv_errVectorAbs_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0))  + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measFreq" + std::to_string(measFreq) + "_measPerDt" + std::to_string(measPerDt_final) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0) * 1000)) + "m_ptM" + std::to_string(pointMulti) + ".csv", errVector);

                        maxErr = errVector.maxCoeff();
                        // centreTemp = VFinal((Nxy_inv+1)/2,n);
                        rankPerc = 1.0 * rankMMt / Nxy_inv;
                        tau = measFreq * dt_inv(nonDim_dt_iterIndex, 0);

                        errTempMatrix(measPerDt_iterIndex, 0) = dt_inv(nonDim_dt_iterIndex, 0);
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

                        cout << "LxIndex [" << nonDim_Lx_iterIndex << "/" << nonDim_Lx_iterIndexMax << "]" << endl;
                        cout << "tauIndex [" << nonDim_tau_iterIndex << "/" << nonDim_tau_iterIndexMax << "]" << endl;
                        cout << "measPerDt_iterIndex [" << measPerDt_iterIndex << "/" << measPerDt_iterIndexMax << "]" << endl;
                        cout << "measFreq for Nx[" << Nx_inv << "] = " << measFreq << endl;
                        cout << "tau for for Nx[" << Nx_inv << "] = " << tau << endl;
                        cout << "measPerDt for Nx[" << Nx_inv << "] = " << measPerDt_final << "/" << measPerDt_max << endl;
                        cout << "\nrankMMt% = " << rankPerc << endl;
                        cout << "Max error for Nx[" << Nx_inv << "] = " << maxErr << endl;
                        cout << "Time step size [" << dt_inv(nonDim_dt_iterIndex, 0) << "]" << endl;
                        // cout << "Centre temp for dx[" << dx_inv << "] = " << centreTemp << endl;

                        // measPerDt_iterIndex++;
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
                        //}

                        cout << "Importing into errTempMatrix_extract for tauIndex [" << nonDim_tau_iterIndex << "/" << nonDim_tau_iterIndexMax << "]..." << endl;
                        // cout << "errTempMatrix_extract rows[" << errTempMatrix.rows() << "] cols[" << errTempMatrix.cols() << "]" << endl;
                        // cout << "errTempMatrx rows[" << errTempMatrix.rows() << "] cols[" << errTempMatrix.cols() << "]" << endl;
                        for (j = 0; j < measPerDt_iterIndexMax; j++)
                        {
                            for (i = 0; i < errTempMatrix_xLen; i++)
                            {
                                // cout << "j = " << j << " i = " << i << "." << endl;
                                errTempMatrix_extract(errTempMatrix_extract_iterIndex, i) = errTempMatrix(j, i);
                            }
                            errTempMatrix_extract_iterIndex++;
                        }
                        cout << "Importing ends!" << endl;

                        // cout << "errTempMatrix for dx[" << dx_inv << "] Nx[" << Nx_inv << "] tau[" << tau << "measFreq[" << measFreq << "] = \n" << errTempMatrix << endl;
                        // saveData("inv_errTempMatrix_2D_"+ std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"_copper_exp_bound" + std::to_string(int(bound_x0)) + "_start" + std::to_string(startNum) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measPerDtMax" + std::to_string(measPerDt_max) + "_measFreq" + std::to_string(measFreq) + "_dt" + std::to_string(dt_inv) + "_ptM" + std::to_string(pointMulti) + "_dxChange.csv", errTempMatrix);

                        nonDim_tau_iterIndex++;
                        //measFreq++;
                        measPerDt_iterIndex = 0;
                        measPerDt_final = measPerDt_max;
                        measPerDt_minus = 0;
                    }
                    cout << "errTempMatrix for dx[" << dx_inv << "] Nx[" << Nx_inv << "] tau[" << tau << "measFreq[" << measFreq << "] = \n"
                         << errTempMatrix_extract << endl;
                    saveData("inv_errTempMatrix_2D_" + std::to_string(int(Lx * 1000)) + "x" + std::to_string(int(Ly * 1000)) + "mm_copper_exp_bound" + std::to_string(int(bound_x0)) + "_Nx" + std::to_string(Nx_inv) + "_Ny" + std::to_string(Ny_inv) + "_measPerDtMax" + std::to_string(measPerDt_max) + "_dt" + std::to_string(int(dt_inv(nonDim_dt_iterIndex, 0) * 1000)) + "m_ptM" + std::to_string(pointMulti) + "_dxChange.csv", errTempMatrix_extract);

                    nonDim_tau_iterIndex = 0;
                    nonDim_inputBound_iterIndex++;

                    // Permission to continue start
                    //         if(nonDim_dx_iterIndex > 0 && nonDim_dx_iterIndex != nonDim_dx_iterIndexMax){
                    //                 cout << "Continue with loop[" << nonDim_dx_iterIndex << "] for dx[" << dx << "] ? [1]Yes [2]End" << endl;
                    //                 cin >> con;
                    //
                    //                 if ( con == 1){
                    //                         cout << "Continuing..." << endl;
                    //                 }
                    //                 else if (con == 2){
                    //                         cout << "Ending & saving..." << endl;
                    //                         break;
                    //                 }
                    //                 else {
                    //                         cout << "Error at choice before M" << endl;
                    //                         return 0;
                    //                 }
                    //         }
                    //         else if(nonDim_dx_iterIndex == nonDim_dx_iterIndexMax){
                    //
                    //         }
                    //         else{
                    //                 cout << "Error at permission of nonDim_dx_iterIndex" << endl;
                    //         }
                    // Permission to continue end
                }
                nonDim_inputBound_iterIndex = 0;
                nonDim_dt_iterIndex++;
            }
            nonDim_dt_iterIndex = 0;
            Nx_inv++;
            Ny_inv++;
            nonDim_dx_iterIndex++;
        }
        nonDim_dx_iterIndex = 0;
        Nx_inv = Nx_inv_INIT;
        Ny_inv = Ny_inv_INIT;
        nonDim_Lx_iterIndex++;
    }
    nonDim_Lx_iterIndex = 0;

    return 0;
}
