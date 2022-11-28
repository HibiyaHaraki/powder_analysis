# powder_analysis

created by Hibiya Haraki, 2022
All risks of running this script is always with you.

## Purpose
powder_analysis includes some MATLAB scripts for analyzing the data that is exported by [PIVlab](https://jp.mathworks.com/matlabcentral/fileexchange/27659-pivlab-particle-image-velocimetry-piv-tool-with-gui).

## Analyzing Process
Step 2 and 3 are the mandatory.

### Step 1. Combine

If you want to create one PIV data file including different pixel size data, please run following script.

* **combinePIVdata.m**

Above script needs following function files, so please put following files at same diretory.

* logging_func.m

### Step 2. Check truck velocity and specify pixel size

If you want to check truck velocity and specify the pixel size for analysing PIV data, please run following script.

* **check_truck_velocity.m**

Above script needs following function files, so please put following files at same diretory.

* determinePixel.m
* logging_func.m

### Step 3. Create mean-map and mean-value

If you want to create a map data including mean of multi experiment data, please run following script.

* **analyze_PIV_data.m**

Above script needs following function files, so please put following files at same diretory.

* get_PIV_XAxis.m
* get_PIV_YAxis.m
* compute_meanMap.m
* compute_meanMap_without_typevector.m
* get_PIV_Data_meanMap.m
* get_PIV_Data_meanMap_without_typevector.m
* logging_func.m

### Step 4. Analyze a map (Analyze specific time data)

If you want to analyze specific time-step, please run following scripts.

* **analyze_a_map.m** (For analyzing fluid velocity)
* **compute_pressure.m** (For analyzing pressure)
※**compute_pressure.m** have assumptions. Please check the script.

Above script needs following function files, so please put following files at same diretory.

* compute_dudx.m
* compute_dvdy.m
* solve_Poisson_rectangle.m
* logging_func.m

### Step 5. Analyze a pixel (Analyze specific pixel data)

If you want to analyze specific a pixel, please run following scripts.

* **get_reference_velocity_ratio.m**
* **analyze_a_pixel.m**

Above script needs following function files, so please put following files at same diretory.

* compute_characteristicVelocity.m
* compute_velocity_u_fluctuation.m
* compute_velocity_v_fluctuation.m
* visualize_velocity_ratio.m
* logging_func.m

### Step 6. Analyze others

There are some scripts as follows.

* **analyze_Xline.m** (Analyze specific x-line fluid velocity)
* **analyze_Yline.m** (Analyze specific y-line fluid velocity)
* **analyze_multiple_Xlines.m** (Analyze specific x-line fluid velocity with several experiment result)
* **visualize_continuity.m** (Analyze continuity of 2-d map)

Above script needs following function files, so please put following files at same diretory.

* compute_Xline.m
* compute_continuity.m
* logging_func.m

### Step 7. Create Animation

If you want to create animation, there are some scripts.

* **create_animation.m** (Create fluid velocity animation)
* **animate_continuity.m** (Create continuity animation)
* **animate_pressure.m** (Create pressure animation)

Above script needs following function files, so please put following files at same diretory.

* compute_dudx.m
* compute_dvdy.m
* solve_Poisson_rectangle.m
* compute_continuity.m
* logging_func.m
