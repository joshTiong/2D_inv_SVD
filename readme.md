# Inverse Analysis Algorithm
## based on SVD of linear system

### Definitions
- SVD : Singular Value Decomposition
- linear system : a set of multiple linear algebraic equations that forms linear equality Ax=b
- FDS : Foward Direct Scheme 前進差分法

### main.cpp files
This algorithm contains three main cpp files based on their different simulation conditions. 
- heatCon_nonDimNum_inv_SVD_2D_LxChange_dxChange_dtChange_TbChange_startTemp.cpp
Inverse Analysis starts at a defined temperature. Space size (Lx, Ly), time step (dt) and boundary temperature (Tb) can be changed.
- heatCon_nonDimNum_inv_SVD_2D_LxChange_dxChange_dtChange_TbChange_startTime.cpp
Inverse Analysis starts at a defined time. Space size (Lx, Ly), time step (dt) and boundary temperature (Tb) can be changed.
- heatCon_nonDimNum_inv_SVD_2D_LxChange_dxChange_measPerDtChange_heatPosChange.cpp
Measurement frequency at time instant n_meas can be reduced from n_max until least possible n.  Space size (Lx, Ly) and time step (dt) can be changed. Position of heat source node in the space also can be edited.

### header files
This algorithm contains header files incorporated into the main cpp files. These header files contain functions that is called in the main cpp files
- heatCon_nonDimNum_FDS_2D_youkaihou.h : calculate the temperature propogation using FDS. The result is used as the validation data and temperature data from measurement nodes is used as measurement data.
- heatPos_oddEven.h : calculate the position of heat source node in even and odd configuration of nodes.
- inputData.h : read .csv into matrices or vectors.
- kronDelta.h : Kronecker delta for x- and y-axis position of heat source node. Heating terms is not zero only at defined heat source node position.
- kronDeltaConst.h : Kronecker delta for heating delay. A delay is set after heating start time.
- kronDeltaConstInv.h : Kronecker delta but not used.
- kronDeltaPulse.h : Kronecker delta for pulse type heat source node. Heating becomes a pulse where heating peaks and returns to zero.
- kronDeltaSoC.h : Kronecker delta for constant type heat source node. Heating peaks at startTime and maintain the value till the end.
- Lx_Ly_vec.h : not used
- removeRowCol.h : remove row and column of matrix function
- saveData.h : save matrix as .csv file.
