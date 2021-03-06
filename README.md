[//]: # (Image References)
[image1]: ./PID_jerk_acc.png
[image2]: ./PID_speed.png

# [Udacity Self-Driving Car Nanodegree Program](https://www.udacity.com/course/intro-to-self-driving-cars--nd113?gclid=EAIaIQobChMIoLbZ8OeY1wIViLftCh01PgTuEAAYASAAEgK45vD_BwE)
# Term 3 Project 1

## Path Planning

### Summary
This report is generated to satisfy the first project in Term 3 of the Udacity Self-Driving Car Nanodegree program, path planning.  The purpose of the project is to create and implement an algorithm in c++ to safely navigate a car through traffic on a highway simulator. Inputs are car position and bearing data from the simulator.  Outputs are car path points. 
### Rubric Points
- Code compiles correctly.
Code compiles with cmake and make
- Car must complete 4.32 miles in the simulator without incident.
This is demonstrated in video linked below.
- Car Drives according to the speed limit.
As demonstrated in the video linked below, the car drives just under the speed limit unless its path is obstructed.
- Car must make changes when it makes sense to do so.
This is demonstrated in video linked below.
- Acceleration may not exceed 10 m/s^2.
This is demonstrated in video linked below.
- Jerk may not exceed 10 m/s^3.
This is demonstrated in video linked below.
- Car must not spend more than three seconds between lanes.
This is demonstrated in video linked below.
- Car must change lanes when it makes sense to do so.
The logic behind lane changes is given below. The performance is demonstrated in the video linked below.

[![IMAGE ALT TEXT](https://i.ytimg.com/vi/2Xphrj5-EkU/0.jpg?time=1509464417261)](https://youtu.be/2Xphrj5-EkU "Project Video")

*Click above image to view the project video*

### Description of Code
#### Lines 1-162:  Udacity provided starter code functions
#### Lines 164-267: Load data and initialize variables
#### Lines 270-308: Determine position and speed of nearest leading and trailing car for each lane
1. Assign each car to a lane.
2. Determine distance from each other car to the controlled, ego, car.
3. Determine if a car in a lane is the closest to the ego car in that lane.  If it is, it stores the car information.
#### Lines 311-339: Determine lane condition
1.  For each lane, convert distances and velocities into times.

Term | Description | Units
 ------ | --------------- | -----
 `vn1`  | Speed of ego car | m/s
 `va`| Speed of nearest leading car | m/s
 `vb` | Speed of nearest trailing car | m/s
 `ta` | Time for ego car to reach the spot occupied the nearest leading car | s
 `tb` | Time for nearest trailing car to reach the spot occupied by the ego car | s
 `tac`| Time until ego car collides with nearest leading car | s
 `tbc` | Time until nearest trailing car collides with ego car | s
 
2. Use times to assign lane condition,`lc`.

Values and description of each lane condition, `lc`, are given in the below table.

Value | Description 
----- | -----
0 | Closed
1 | Wide open
2 | Somewhat open
3 | Almost closed

3.  Assign conditions based on times `ta`, `tb`, `tac`, `tbc`.

Each lane condition is first set to `Closed`. If the `Wide open` test is passed, the lane condition is set to `Wide open`. If the `Wide open` test is failed and the `Somewhat open` test is passed, the lane condition is set to `Somewhat open`. If the `Wide open` and `Somewhat open` tests are failed and the `Almost closed` test is passed, the lane condition is set to `Almost closed`. Definitions of each test are given in the below table.

Test | Description
----- | -------
Wide open | `ta > 2.5 && tb > 2.0 && tac > 10.0 && tbc > 10.0`
Somewhat open | `ta > 1.6 && tb > 1.0 && tac > 10.0 && tbc > 10.0`
Almost closed | `ta > 0.5 && tb > 0.5 && tac > 5.0 && tbc > 5.0`

#### Lines 341-442: Select lane
To make a lane change, it is necessary, but not sufficient that the car be within 0.5 m of the lane centerline and that the car speed be either within 5 mph of the nearest leading car or at least 45 mph.

Lane numbers and descriptions are given in the below table.

Lane Number | Lane Description
----------- | ---------
0 | Left
1 | Middle
2 | Right

Lane specific conditions for making a lane change are given below.

 Description | Test
------ | ------
Left to center | `lc[1] == 1 && lc[2] > 0`
Left to center | `lc[1] == 2 && lc[2] ==1`
Center to right | `lc[2] == 1`
Right to middle | `ta[2] < 2.1 && va[2] < 45.0 mph && lc[1] == 1 && lc[0] > 0`
Right to middle | `ta[2] < 2.1 && va[2] < 45.0 mph && lc[1] == 2 && lc[0] == 1`
Middle to left | `ta[1] < 2.1 && va[1] < 45.0 mph && lc[0] == 1`

Upon making a lane change the safety time, `Ts` is set to 2.0 s if the lane condition of the new lane is `Wide open` and 1.5 s if the lane condition of the new lane is `Somewhat open`. The above logic includes a right-lane bias as is custom in the United States.  It also used the two second safety distance this author was taught in drivers education school. This distance is slightly compromised to get from a congested lane to a free lane.  The benefits of decreased travel time and increased safety resulting from the expected reduced car density justify the tradeoff.

#### Lines 424-477: Set speed
A PD controller is used to set the speed.
The error is by 
`err = sa[lane] - vn1 * Ts`
The change in the error is given by 
`derrdt = va[lane] - max(0.1, vn1)`

Constants, values and justifications are given in the below table.

Constant | Term | Value | Units | Justification/ Description
----------- | ----  | -------- | ------- | ---------
Proportional value | `P` | 1.0 | m/s^2 / m | Strong but not overwhelming value. Refer to the below step response plot.
Derivative vale | `D` | 0.02 | m/2^2 / m/s | Light damping.  Refer to the below step response plot.
Jerk limit | `j_max` | 9.0 | m/s^3 | Stay under the given 10.0 m/s^3 limit
Acceleration limit | `a_max` | 9.0 | m/s^2 | Stay under the given 10.0 m/s^2 limit
Time step | `dt` | 0.02 | s | Time between simulation steps
Number of lag steps | `N_lag` | # |int | Number of time steps between calculations
Acceleration | `a` | double | m/s^2 | acceleration
Acceleration at previous decision| `ap` | double | m/s^2 | acceleration at previous decision
Reference velocity | `ref_vel` | double | mph | speed

![alt text][image1]
![alt text][image2]
*Step response of PD controller*

The acceleration is set by the calculations listed below.

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
As described below, a buffer of size `N_points` is generated.  The buffer size is adjusted at line 611 so that about eight points remain after processing by the simulator.  On the authors computer, the number of points tended to stay between 9 and 12 points which allows the car to be projected between 0.18 and 0.24 seconds into the future.  The starting point for setting the points is the last two points if the buffer length is at least two.  If it is not, the starting points are the car position and a point directly behind the position of the car. 
#### Lines 520-535: Add points to be used to fit the path
First, the distance to the road centerline, `s_ahead` is selected.  This is the minimum of 20 m and the velocity in m/s multiplied by 2.0 s. Assuming one lane change of 4 m, this limits the during change jerk to 1 m/s^3 and the peak normal velocity to 6.0 m/s^2. This balances the risk with being between lanes with ride comfort. 

Then, points at `s_ahead`, `2 * s_ahead` and `3 * s_ahead` are added to the two starting points to form the set of coordinates to be fit to a spline, `ptsx` and `ptsy`.
#### Lines 538-554:  Transform points and create a spline to fit future points
`ptsx` and `ptsy` are transformed from the map reference frame to the car coordinate reference frame and are then used to create a spline, `s`, using functions contained in the `spine.h` file which can be found [here](http://kluge.in-chemnitz.de/opensource/spline/).  The car reference frame is one in which the car is centered at the origin and the bearing is the x-axis.
#### Lines 556-598:  Create set of points to feed back to the simulator
1. The previous buffer points are loaded into the vectors `next_x_vals` and `next_y_vals`. 
2. The number of new points required, `N-lag`, is set to the difference of ten and the number of points remaining in the buffer.
3. The change in x-direction is estimated. It is assumed that the car will return to the track centerline in 20.0 m along a straight line of distance, `dist = sqrt[(y(x = 20 m) - y(x = 0 m))^2 + (20 m) ^ 2]`The `x_point` values are set so that `dist / dt = ref_vel`.  
4. The `x_point` values are fed into the spline `s` to generate new `y_points`.
5. The `x_point` and `y_point` values are transformed from the car reference frame to the map reference frame.
6. These transformed points are added to the `next_x_vals` and `next_y_vals` vectors.
#### Lines 602-608: `next_x_vals` and `next_y_vals` vectors are fed back to the simulator.
#### Lines 609-610:  Push back speed values
#### Line 611:  Adjust the number of point remaining.
Line 611 is: '`if(prev_size < 8){N_points += 1;};  if(prev_size > 8){N_points -= 1;}`
This keeps the number of points remaining in the buffer near eight.  This helps prevent emptying the buffer on a slow machine without sacrificing responsiveness on a fast machine.    
### Conclusion
A c++ program to safely navigate a car through traffic on a highway on the Udacity simulator was successfully developed and implemented. 
