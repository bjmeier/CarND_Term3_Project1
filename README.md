[//]: # (Image References)
[image1]: ./PID_jerk_acc.png
[image2]: ./PID_speed.png

# Udacity Self-Driving Car Nanodegree Program 
# Term 3 Project 1

## Path Planning

### Summary
This report is generated to satisfy the first project in Term 3 of the Udacity Self-Driving Car Nanodegree program, path planning.  The goal is to cretae and implement an algorithm in c++ to safely navigate a car through traffic on a highway with a speed limit of 50 mph. The Udacity simulator feeds car position and car location.  Starter functions and map data were provided to located the car on the map and convert between catesian and Frenet coordinates. 
### Rubric Points
- Code compiles correctly
Code compiles with cmake and make
- Car must complte 4.32 miles in the simulator without incident.
This is demonstrated in video linked below.
- Car Drives according to the speed limit
As demonstrated in the video linked below, the car drives just under the speed limit unless its path is obstructed.
- Car must make changes when it makes sense to do so
Demonstrated in video linked below.
- Acceleration may not exceed 10 m/s^2
This is demonstrated in video linked below.
- Jerk may not exceed 10 m/s^3
This is demonstrated in video linked below.
- Car must not spend more than three seconds between lanes.
This is demonstrated in video linked below.
- Car must change lanes when it makes sense to do so.
The logic behind lane changes is given below. The performance is demonstrated in the video linked below.


[![IMAGE ALT TEXT](https://i.ytimg.com/vi/vYhxZa9cnVI/0.jpg?time=1509368669565)](https://youtu.be/vYhxZa9cnVI "Project Video")
*Click above image to view the project video*

### Description of Code
#### Lines 1-162:  Starter code functions
#### Lines 164-267: Load data and initialize variables
#### Lines 270-308: Determine postion and speed of nearest leading and trailing car for each lne
1. Assign each car to a lane.
2. Determine distance from each other car to the controlled, ego, car.
3. Determine if a car in a lane is the closest to the ego car.  If it is, it stores this information for each lane.
#### Lines 311-339: Determine lane condition
1.  For each lane, convert distances and velocites into times

Term | Description | Units
 ------ | --------------- | -----
 `vn1`  | Speed of ego car | m/s
 `va`| Speed of nearest leading car | m/s
 `vb` | Speed of trailing car | m/s
 `ta` | Time for ego car to reach the spot occupied the nearest leading car | s
 `tb` | Time for nearset trailing car to reach the spot occupied by the ego car | s
 `tac`| Time until ego car collides with nearest leading car | s
 `tbc` | Time until nearest trailing car collides with ego car | s
 
2. Use times to assign lane condition,`lc`

Value | Description 
----- | -----
0 | Closed
1 | Wide open
2 | Somewhat open
3 | Almost closed

3.  Assign conditions based on times `ta`, `tb`, `tac`, `tbc`.

Each lane condition is first set to closed. If the wide open test is passed, the lane condition is set to wide open. If the wide open test is failed and the somewhat open test is passed, the lane condition is set to somewhat open. If the wide open and somewhat open tests are failed and the almost closed test is passed, the lane condition is set to almost closed. Definitions of each test are given in the below table.

Test | Description
----- | -------
Wide open | `ta > 2.5 && tb > 2.0 && tac > 10.0 && tbc > 10.0`
Somewhat open | `ta > 1.6 && tb > 1.0 && tac > 10.0 && tbc > 10.0`
Almost closed | `ta > 0.5 && tb > 0.5 && tac > 5.0 && tbc > 5.0`

#### Lines 341-442: Select lane
In make any lane change, it is necessary, but not sufficint that the car be within 0.5 m of the lane centerline and that the car speed be either within 5 mph of the nearest leading car or at least 45 mph.

Lane numbers and descriptions are given in the below table.

Lane Number | Lane Description
----------- | ---------
0 | Left
1 | Middle
2 | Right

Additional conditions for making a lane change are given below.

 Description | Test
------ | ------
Left to center | `lc[1] == 1 && lc[2] > 0`
Left to center | `lc[1] == 2 && lc[2] ==1`
Center to right | `lc[2] == 1`
Right to middle | `ta[2] < 2.1 && va[2] < 45.0 mph && lc[1] == 1 && lc[0] > 0`
Right to middle | `ta[2] < 2.1 && va[2] < 45.0 mph && lc[1] == 2 && lc[0] == 1`
Middle to left | `ta[1] < 2.1 && va[1] < 45.0 mph && lc[0] == 1`

Upon making a lane change the safety time, `Ts` is set to 2.0 s if the lane condition of the new lane is wide open and 1.5 s if the lane condition of the new lane is somewhat open. The above logic includes a right-lane bias as is custom in the United States.  It also used the two second "safety distance" this author was taught in drivers education school. This distance is slightly comprimised to get from a congeted lane to a free lane.  The benefits of decreased travel time and increased safety resulting from the expected reduced car denisty justifey the tradeoff.

#### Lines 424-477: Set speed
A PD controller is used to set the speed.
The error is by 
`err = sa[lane] - vn1 * Ts`
The change in the error is given by 
`derrdt = va[lane] - max(0.1, vn1)`

Constants, values and justifications are given in the below table.

Constant | Term | Value | Units | Justification
----------- | ----  | -------- | ------- | ---------
Proportional value | P | 1.0 | m/s^2 / m | Strong but not overwhelming value. Refer to the below step response plot.
Derivitive vale | D | 0.02 | m/2^2 / m/s | Light damping.  Refer to the below step resonse plot.
Jerk limit | j_max | 9.0 | m/s^3 | Stay under the given 10.0 m/s^3 limit
Accleration limit | a_max | 9.0 | m/s^2 | Stay under the given 10.0 m/s^2 limit
Time step | dt | 0.02 | s | Time between simulation steps
Number of lag steps | N_lag | # |int | Number of time steps between calculations
Acceleration | a | double | m/s^2 | acceleration
Acceleration at previous decision| ap | double | m/s^2 | acceleration at previous decision
Reference velocity | ref_vel | double | mph | speed

![alt text][image1]
![alt text][image2]
*Step response of PD controller*

The acceleration is set by the logic given below.

	a = 9.0;
	a = P * err + D * derrdt;

	a = min(a,  a_lim); // max acc constraint
	a = max(a, -a_lim); // min acc constraint
	
	// max speed constraint
	a = min(a, (49.5*0.44704 - vn1)/(N_lag*dt));
	if(a > a_lim * dt * N_lag && (49.5*0.44704 - vn1) < (ap*ap/(2*j_lim)))	
	    {
             a = ap - a_lim * dt * N_lag;
	    }

    a = min(a, ap + j_lim * dt * N_lag); // max jerk constraint
    a = max(a, ap - j_lim * dt * N_lag); // min jerk constraint

    a = min(a, (49.9 * 0.44704 - max(0.1, vn1))/(N_lag*dt)); // speed limit constraint
    a = max(a, ( 0.1 * 0.44704 - max(0.1, vn1))/(N_lag*dt)); // speed limit constraint
    ref_vel = (vn1 + a * (N_lag*dt))/0.44704;
    `

#### Lines 478-487: Set reference points
The car coordinate and bearing, `car_x` and `car_y` and `car_yaw` are set in m, m, and degrees respectively.
#### Lines 489-519: Select starting points
As is decribed below, a 10 step long (0.2 s) buffer is generated.  The starting point for setting the points.  It is the ast two poins if the buffer as long as the buffer length is at least two.  If it is not, the starting points are the car position and a point directly behind the position of the car. 
#### Lines 520-535: Add points to be used to fit the path
First, the distance to the road centerline, `s_ahead` is selected.  This is the minimum of 20 m and the velocity in m/s * 2.0s. Assuming one lane change of 4 m, this limits the during change jerk to 1 m/s^3 and the peak normal velocity to 6.0 m/s^2. This balances the risk with being betwen lanes with ride comfort. 

Then, points at `s_ahead`, `2 * s_ahead` and `3 * s_ahead` are added to the two starting points to form the set of coordinates to be fit to a spline, `ptsx` and `ptsy`.
#### Lines 538-554:  Transform points and create a spline to fit future points.
`ptsx` and `ptsy` are transformed from the map reference frame to the car coordinate reference frame and are then used to create a spline, `s`, using the `spine.h` file.  The car reference frame is one in which the car is centered at the origin and the bearing is the x-axis.
#### LInes 556-597:  Create set of points to feed back to the simulator.
1. The previous buffer points are loaded into the vectors `next_x_vals` and `next_y_vals`. 
2. The number of new points requried, `N-lag` is set to the difference of ten and the number of points remaining in the buffer.
3. The change in x-direction is estimated. It is assumed that the car will return to the track centerine in 20.0 m along a straight line. The `x_point` values are set so that the quotent of the distance along the hypotenuse and the time step, `dt` yields the reference velocity, `ref_vel`.
4. The `x_point` vaulues are fed into the spline `s` to generate new `y_points`.
5. The `x_point` and `y_point` values are transformed from the car reference frame to the map reference frame.
6. These transformed points are added to the `next_x_vals` and `next_y_vals` vectors.
#### Lines 601-602: `next_x_vals` and `next_y_vals` vectors are fed back to the simulator.
### Conclusion
A path planning algorithm has been developed and implemented in the Udacity Term 3 simulator.  
