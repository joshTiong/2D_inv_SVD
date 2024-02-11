//FDS Variables 2
//heatPosIndex = heatPos/dx;
endTimeIndex = endTime/dt_FDS;
boundHosei_FDS = 1;
Nx_FDS_woBound = Nx_FDS - boundHosei_FDS;
Ny_FDS_woBound = Ny_FDS - boundHosei_FDS;
Nxy_FDS_woBound = (Nx_FDS_woBound)*(Ny_FDS_woBound);
Nx_FDS_wBound = Nx_FDS + boundHosei_FDS;
Ny_FDS_wBound = Ny_FDS + boundHosei_FDS;
R = 8.3144;
condSoC = 275;
flatTimeIndex = 0;
flatXIndex = 0;
flatXEndIndex = (Nx_FDS_woBound*Ny_FDS_woBound)-1;
timeCountIndex = 0;
kappa_x_FDS = alpha*dt_FDS/(dx_FDS*dx_FDS);
kappa_y_FDS = alpha*dt_FDS/(dy_FDS*dy_FDS);

//#include "heatPos_oddEven.h"

//Electrolyte
Ea_e = 3.014e+05;
A_e = 5.140e+25;
W_e = 510;
H_e = 3.410e+05;
//SEI
Ea_sei = 1.490e+05;
A_sei = 1.667e+15;
W_sei = 385;
H_sei = 5.780e+05;
//Cathode
Ea_n = 1.418e+05;
A_n = 2.500e+13;
W_n = 385;
H_n = 3.428e+06;
//Anode
Ea_p = 0.963e+05;
A_p = 2.000e+08;
W_p = 615;
H_p = 2.414e+05;//J/kg
//SoC (Short)
Con_SoC = 1; //-
A_SoC = 3.37e+12; //s^-1
Ea_SoC = 1.58e-19; //J
R_b = 1.380649e-23; //J/K
V_SoC = 4.2; //V
C_SoC = 2.4; //Ah
eff_SoC = 0.28;//-
H_SoC = V_SoC*C_SoC*3600*eff_SoC;//J
vol_SoC = pow((Lx/Nx_SoC),3);

//Build Space
MatrixXd T_spaceMatrix = MatrixXd::Zero(Ny_FDS_wBound, Nx_FDS_wBound);
MatrixXd T_spaceMatrix_np1 = MatrixXd::Zero(Ny_FDS_wBound, Nx_FDS_wBound);
MatrixXd I_spaceMatrix = MatrixXd::Identity(Ny_FDS_wBound, Nx_FDS_wBound);
MatrixXd T_flat = MatrixXd::Zero(endTimeIndex + 1, Nxy_FDS_woBound);
MatrixXd qdotM = MatrixXd::Zero(endTimeIndex, 11);
VectorXd d_yaxis_vec = VectorXd::Zero(Ny_FDS);

//Initialize T_spaceM
//x direction boundary
for(i=0; i<Nx_FDS_wBound; i++){
    T_spaceMatrix(0,i) = inputBound;
    T_spaceMatrix(Ny_FDS_wBound-1,i) = inputBound;
}

//y direction boundary
for(j=1; j<Ny_FDS; j++){
    T_spaceMatrix(j,0) = inputBound;
    T_spaceMatrix(j, Nx_FDS_wBound-1) = inputBound;
}

//initialize inner temp
for(j=1; j<Ny_FDS; j++){
    for(i=1; i<Nx_FDS; i++){
        T_spaceMatrix(j,i) = initialT_in;
    }
}

//saveData("measForward_FDS_expHeating_initialT_2D_"+ std::to_string(int(Lx*100)) +"x"+ std::to_string(int(Ly*100)) +"cm_copper_bound0_Nx" + std::to_string(Nx) + "_Ny" + std::to_string(Ny) + "_dt0001.csv", T_spaceMatrix);

//Make T_spaceMatrix_np1 as T_spaceMatrix
T_spaceMatrix_np1 = T_spaceMatrix*I_spaceMatrix;

//Save to T_flat
for(j=1; j<Ny_FDS; j++){
    for(i=1; i<Nx_FDS; i++){
        T_flat(flatTimeIndex, flatXIndex) = T_spaceMatrix(j,i);
        flatXIndex++;
    }
}
flatTimeIndex++;
flatXIndex = 0;

//saveData("measForward_FDS_expHeating_Tflat.csv", T_flat);

//Permission to continue start
//	cout << "initialT_spaceMatrix =\n" << T_spaceMatrix << endl;
//        cout << "Continue with FDS calc? [1]Yes [2]End" << endl;
//        cin >> cond;
//
//        if ( cond == 1){
//                cout << "Continuing..." << endl;
//        }
//        else if (cond == 2){
//                cout << "Ending..." << endl;
//                return 0;
//        }
//        else {
//                cout << "Error at choice before A" << endl;
//                return 0;
//        }
//Permission to continue end


//Youkaihou calc
while(timeCountIndex < endTimeIndex){
    //delta functions
    iterTime = (timeCountIndex+1)*dt_FDS;
    kronDeltaConstFunc(d_delay, iterTime, condTime);
    
    for(j=1; j<Ny_FDS; j++){
        //Heat Position delta function for y-axis
        iterPos_y = j*dy_FDS;
        kronDeltaFunc(d_yaxis, j, heatPos_y_index);
        
        d_yaxis_vec(j) = d_yaxis;

        
        for(i=1; i<Nx_FDS; i++){
            //Heat Position delta function for x-axis
            iterPos_x = i*dx_FDS;
            kronDeltaFunc(d_xaxis, i, heatPos_x_index);
            if(i==5 && j==5){
                //cout << "iterPos_x = " << iterPos_x << " heatPos_x = " << heatPos_x << " d_xaxis = " << d_xaxis << endl;
                //cout << "iterPos_y = " << iterPos_y << " heatPos_y = " << heatPos_y << " d_yaxis = " << d_yaxis << endl;
            } 
            
            
            
            
            //Permission to continue start
            //cout << "d_SoC = \n" << d_SoC << endl;
            //cout << "T_spaceMatrix = \n" << T_spaceMatrix << endl;
            //cout << "condSoC = \n" << condSoC << endl;
            //cout << "Continue? [1]Yes [2]End" << endl;
            //cin >> con;
            
            //if ( con == 1){
            //    cout << "Continuing..." << endl;
            //}
            //else if (con == 2){
            //    cout << "Ending..." << endl;
            //    return 0;
            //}
            //else {
            //    cout << "Error at choice before M" << endl;
            //    return 0;
            //}
            //Permission to continue end
            kronDeltaSoCFunc(d_SoC, T_spaceMatrix(j,i), condSoC);
            
            qdot_e = H_e*W_e*A_e*exp(-Ea_e/(R*T_spaceMatrix(j,i)));
            qdot_sei = H_sei*W_sei*A_sei*exp(-Ea_sei/(R*T_spaceMatrix(j,i)));
            qdot_n = H_n*W_n*A_n*exp(-Ea_n/(R*T_spaceMatrix(j,i)));
            qdot_p = H_p*W_p*A_p*exp(-Ea_p/(R*T_spaceMatrix(j,i)));
            qdot_total = qdot_e + qdot_sei + qdot_n + qdot_p;
            qdot_SoC = H_SoC*Con_SoC*A_SoC*(1/vol_SoC)*(exp(-Ea_SoC/(R_b*T_spaceMatrix(j,i))));

            
            T_spaceMatrix_np1(j,i) = kappa_x_FDS*T_spaceMatrix(j,i-1) + kappa_x_FDS*T_spaceMatrix(j,i+1) + kappa_y_FDS*T_spaceMatrix(j-1,i) + kappa_y_FDS*T_spaceMatrix(j+1,i) + (1 - 2*kappa_x_FDS - 2*kappa_y_FDS)*T_spaceMatrix(j,i)+ d_delay*d_xaxis*d_yaxis*d_SoC*(dt_FDS/(rho*cp))*qdot_SoC + d_delay*d_xaxis*d_yaxis*(dt_FDS/(rho*cp))*qdot_total;
            
            //*******Caution********
            //Save to qdotM is only for single point heating -------->
            //*******Caution********
            if(d_xaxis==1 && d_yaxis==1){
                qdotM(timeCountIndex,1) = qdot_e;
                qdotM(timeCountIndex,2) = qdot_sei;
                qdotM(timeCountIndex,3) = qdot_n;
                qdotM(timeCountIndex,4) = qdot_p;
                qdotM(timeCountIndex,5) = qdot_total;
                qdotM(timeCountIndex,6) = qdot_SoC;
                qdotM(timeCountIndex,7) = exp(-Ea_SoC/(R_b*T_spaceMatrix(j,i)));
                qdotM(timeCountIndex,8) = (-Ea_SoC/(R_b*T_spaceMatrix(j,i)));
                qdotM(timeCountIndex,9) = T_spaceMatrix(j,i);
                qdotM(timeCountIndex,10) = d_SoC;
            }
        }
    }
    T_spaceMatrix = T_spaceMatrix_np1*I_spaceMatrix;
    
    //Insert into T_flat
    for(j=1; j<Ny_FDS; j++){
        for(i=1; i<Nx_FDS; i++){
            T_flat(flatTimeIndex, flatXIndex) = T_spaceMatrix(j,i);
            flatXIndex++;
        }
    }
    
    //Loop count
    if((timeCountIndex+1)%10000 != 0){
    }
    else if((timeCountIndex+1)%10000 == 0){
        cout << "Loop progress --------> [" << timeCountIndex+1 << "/" << endTimeIndex << "]" << endl;
    }
    else{
        cout << "error at initial loop count" << endl;
        return 0;
    }
    
    timeCountIndex++;
    flatTimeIndex++;
    flatXIndex = 0;
}

for (i=0; i<endTimeIndex; i++){
    if(qdotM(i,9) > 333) {
        qdotM(1,0) = i;
        qdotM(2,0) = qdotM(1,0) - qdotM(0,0);
        break;
    }
    else if (qdotM(i,9) > inv_tempStart && qdotM(i,9) <= 333 && qdotM(0,0) == 0) {
        qdotM(0,0) = i;
        inv_tempStart_idx = i;
    }
    else if (qdotM(i,9) <= 333){
        //through
    } 
    else {
        cout << "error at finding 353K" << endl;
        cout << "Nx = " << Nx_FDS << ", Lx = " << Lx << endl;
        return 0;
    }
}

saveData("measForward_FDS_expHeating_Tflat_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_bound"+ std::to_string(int(inputBound)) +"_Nx" + std::to_string(Nx_FDS) + "_Ny" + std::to_string(Ny_FDS) + "_NxSoC" + std::to_string(int(Nx_SoC)) + "_dt" + std::to_string(int(dt_FDS*1000)) + "m.csv", T_flat);
saveData("measForward_FDS_expHeating_qdotM_2D_"+ std::to_string(int(Lx*1000)) +"x"+ std::to_string(int(Ly*1000)) +"mm_copper_bound"+ std::to_string(int(inputBound)) +"_Nx" + std::to_string(Nx_FDS) + "_Ny" + std::to_string(Ny_FDS) + "_NxSoC" + std::to_string(int(Nx_SoC)) + "_dt" + std::to_string(int(dt_FDS*1000)) + "m.csv", qdotM);

