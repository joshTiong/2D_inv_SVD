        if(Nx_FDS%2 == 0 && Ny_FDS%2 == 0){
        	heatPos_x = heatPos_x_even;
		heatPos_y = heatPos_y_even;
	}
        else if(Nx_FDS%2 != 0 && Ny_FDS%2 != 0){
		//heatPos_x_odd = ((Nx_FDS+measDx)/2)*dx_FDS;
		//heatPos_y_odd = ((Ny_FDS+measDy)/2)*dy_FDS;
		//heatPos_x = heatPos_x_odd;
		//heatPos_y = heatPos_y_odd;
                heatPos_x_odd_index = int(Nx_FDS/2);
                heatPos_x_odd_index = int(Ny_FDS/2);

        }
        else{
                cout << "Error at odd x-axis & y-axis space correction" << endl;
                return 0;
        }
