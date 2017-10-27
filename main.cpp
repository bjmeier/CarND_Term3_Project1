#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "Eigen-3.3/Eigen/Dense"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }


// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

vector<double> getAjmt(double y, double yp, double ypp, double yf, double ypf, double yppf, double xf)
// finds constants for a quintic polynomial y = f(x)
// y is the y value at x = 0
// yp = dy/dx at x = 0
// ypp = d2y/dx2 at x = 0
// xf is the final value for x
// y is the y value at x = xf
// ypf = dy/dx at x = xf
// yppf = d2y/dx2 at x = xf
// y = dy/dx = d2y/dx2 = 0 at x = xf
	{
		double a0 = y;
		double a1 = yp;
		double a2 = ypp/2;
		double b0 = yf - a0 - a1 * xf - a2 * xf * xf;
		double b1 = ypf - a1 - 2 * a2 * xf;
		double b2 = yppf - 2 * a2;
		MatrixXd b(3,1);
		b << b0, b1, b2;
		MatrixXd A(3,3);
		A <<    xf*xf*xf, xf*xf*xf*xf, xf*xf*xf*xf,
				3*xf*xf , 4*xf*xf*xf , 5*xf*xf,
				6*xf    , 12*xf      , 20*xf;
		MatrixXd alpha(3,1);
		alpha = A.inverse()*b;
		return{a0, a1, a2, alpha(0,0), alpha(1,0), alpha(2,0)};
	}

double getYjmt(vector<double> a, double x)
// a is a vector of length 6
// and contains quintic polynomial constants a0, a1, .. a5
// returns a0 + a1*x + a2*x^2+...a5*x^5
{
	double result = 0.0;
	for (int i = 0; i < 6; i++)
	{
		result += a[i]*pow(x,i);
	}
	return(result);
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);
	return {x,y};
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  double ref_vel = 0.0; // reference speed in mph
  double vn1 = 2.0; 	// speed at t - dt in m/s
  double vn2 = vn1 - 9.0/50; 	// speed at t - 2 * dt in mps
  int lane = 1;			// lane occupied by car
  double ep = 0;

  h.onMessage([&ep, &vn1, &vn2, &ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
	double dt = 1.0/50;	// time between steps in s
	double Ts = 2.0; 	// safety distance between cars in m
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {


      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON objectls

          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
          	cout << "lane = " << lane << " car_s = " << car_s << " car d = " << car_d << " car speed = " << car_speed << endl;

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.

          	auto sensor_fusion = j[1]["sensor_fusion"];

          	int prev_size = previous_path_x.size();

          	if(prev_size>0)
          	{
          		car_s = end_path_s;
          	}

          	double a = 9.0;		// acceleration in m/s^2
          	double atemp = 0.0; // temporary acceleration in m/s^2
          	double ds = 1000.0; // distance between car and car in front in m.  Initially set 1 km.
          	double dstemp = 1000.0;  // temporary distance between cars Initially set 1 km.
          	vector<double> sa(3, 1000.0); // distance between car and nearest car in front for each lane in m
          	vector<double> sb(3, -1000.0);// distance between car and nearest car behind for each lane in m
          	vector<double> va(3);// speed of vehicle agead
          	vector<double> vb(3, 0.1); // initialize vehicles behind to be near zero speed
          	va[0] = 50 * 0.44704;	// initialize speed ahead to 50 mph
          	va[1] = 50 * 0.44704;
          	va[2] = 50 * 0.44704;
          	vector<int> lc(3, 0); // initialze lane clear elements to 0 (not clear)

          	// determine closest car behind and in front for each lane
          	for(int i = 0; i<sensor_fusion.size(); i++)
          	{
          		float d = sensor_fusion[i][6];
          		double vx = sensor_fusion[i][3];
          		double vy = sensor_fusion[i][4];
          		double check_speed = sqrt(vx*vx + vy*vy);
          		double check_car_s = sensor_fusion[i][5];
          		check_car_s += ((double)prev_size * dt * check_speed);
          		double dstemp = check_car_s - car_s;
          		int cl = 0;
          		// determine and set lane of each car
          		for(int j = 0; j < 3; j++)
          		{
          			if(d <( 4 * j + 4) && d >= (4 * j))
          			{
          				cl = j;
          			}
          		}
          		// if nearest in front, set sa, va
          		if (dstemp >= 0 && dstemp < sa[cl])
          		{
          			sa[cl] = dstemp;
          			va[cl] = check_speed;
          		}
          		// if nearest behind, set sb, vb
          		if (dstemp < 0 && dstemp > sb[cl])
          		{
          			sb[cl] = dstemp;
          			vb[cl] = check_speed;
          		}
          	}


          	for(int i = 0; i < 3; i++)
          	{
          		// times below assume cars maintain constant speeds
          		double tb = -sb[i]/vb[i]; // time until car behind reaches spot occupied by car in s
          		double ta =  sa[i]/max(0.1, vn1); // time until car reaches spot occupied by car in front in s
          		double tac = sa[i]/max(0.01, (max(0.1, vn1) - va[i])); // time until car collides with car in front in s
          		double tbc =-sb[i]/max(0.01, (vb[i] - max(0.1, vn1))); // time until car is rear ended  in s

          		// determine if lane is wide open
          		if(ta > 2.5 && tb > 2.0 && tac > 10.0 && tbc > 10.0)
          		{
          			lc[i] = 1;
          		}
          		// determine if lane is open
          		else if (ta > 1.6 && tb > 1.0 && tac > 10.0 && tbc > 10.0)
          		{
          			lc[i] = 2;
          		}
          		// determine if lane is almost closed
          		else if (ta > 0.5 && tb > 0.5 && tac > 5.0 && tbc > 5.0)
          		{
          			lc[i] = 3;
          		}

          		cout << "i = " << i << " lc = " << lc[i] <<
          				" ta = " << ta << " tb = " << tb <<
						" va = " << va[i] << " vb = " << vb[i] <<
						" tac = " << tac << " tbc = " << tbc << endl;
          	}

          	// select lane
          	double vahead = va[lane]; // determine speed of car ahead in lane
          	int temp_lane = lane; // set temp_lane to current lane
          	double ta = sa[lane]/max(vn1, 0.1); // time until car reaches spot held by car ahead.
          	double derr = 2.0  + 4.0 * lane - car_d; // make sure car is in middle of lane.  Used to prevent double lane changes
          	// speed test is used to make sure car is moving close to the speed limit or the speed of the car ahead prior to changing lanes
          	int speed_test = 0;
          	if(vn1 >= 45.0*0.44704 || ((vahead - vn1) / 0.44704 <= 5.0))
          			{
          				speed_test = 1;
          			}
          	// car moves from left to center lane if right lane isn't closed and middle lane is wide open.
          	// biaed to move car to right lane.
          	if(	lane == 0 &&
          			lc[1] == 1 &&
					abs(derr) < 0.5 &&
					lc[2] > 0 &&
					speed_test == 1)
                  	{
                  		temp_lane = 1;
                  		Ts = 2.0;
                  	}
          	// car moves from left lane to center lane if right lane is wide open and center lane is mostly open
          	if(		lane == 0 &&
           			lc[1] == 2 &&
        			abs(derr) < 0.5 &&
        			lc[2] == 1 &&
        			speed_test == 1)
          	{
          		temp_lane = 1;
          		Ts = 1.5;
          	}

          	// car moves from center lane to right lane if right lane is wide open
          	if(lane == 1 &&
          			lc[2] == 1 &&
					abs(derr) < 0.5 &&
					speed_test == 1)
          	{
          		temp_lane = 2;
          		Ts = 2.0;
          	}

          	// car moves from center lane to left lane if left lane is wide open, right lane isn't wide open and car in front is moving
          	// below the speed limit
          	if(lane == 1 &&
          			lc[0] == 1 &&
					ta < 2.1 &&	va[lane] < (45.0 * 0.44704) &&
					abs(derr) < 0.5 &&
					speed_test == 1)
          	{
          		temp_lane = 0;
          		Ts = 2.0;
          	}

          	// car moves from right lane to middle lane if middle lane is wide open, left lane isn't closed and car in front is moving
          	// below the speed limit
          	if(lane == 2 &&
					lc[1] == 1 &&
					ta < 2.1 &&	va[lane] < (45.0 * 0.44704) &&
					abs(derr) < 0.5 &&
					lc[0] > 0 &&
          			speed_test == 1)
             {
                    temp_lane = 1;
                    Ts = 2.0;
             }

          	// car moves from right lane to middle lane if middle lane is open, the left lane is wide open and the car in front is moving
          	// below the speed limit
            if(lane == 2 &&
          			lc[1] == 2 &&
          			ta < 2.1 &&	va[lane] < (45.0 * 0.44704) &&
          			abs(derr) < 0.5 &&
          			lc[0] == 1 &&
          			speed_test == 1)
          	{
          		temp_lane = 1;
          		Ts = 1.5;
          	}

          	lane = temp_lane;

          	// set speed
          	// use PD control with restricitons
          	// strong, but not overwhelming P
          	// small D
          	// max acc is set to 9.0 m/s^s; max jerk is set to 9.0 m/s^3
          	// min acc is set to 9.0 m/s^s; min jerk is set to 9.0 m/s^3
			ds = sa[lane];
			cout << "ds = " << ds << " v = " << vn1 << " t = " << ds/vn1 << endl;

            double err = ds - max(0.1, vn1)*Ts;
            double derrdt = (err - ep)/dt;
            ep = err;
            double P = 1.0;
            double D = 0.02;
            double ap = (max(0.1, vn1) - max(0.05, vn2))/dt;
            atemp = P * err + D * derrdt;
            atemp = min(atemp, 9.0); // max acc constraint
            atemp = max(atemp, -9.0); // min acc constraint
            // max speed conatraint
            atemp = min(atemp, (49.5*0.44704 - vn1)/dt);
            if(atemp > 0.9*dt && (49.5*0.44704 - vn1) < (ap*ap/18.0))
            {
            	atemp = ap - 9.0 * dt;
            }

            atemp = min(atemp, ap + 9.0*dt); // max jerk constraint
            atemp = max(atemp, ap - 9.0*dt); // min jerk constraint
            atemp = min(atemp, (49.9 * 0.44704 - max(0.1, vn1))/dt); // speed limit constraint
            atemp = max(atemp, (0.1 * 0.44704 - max(0.1, vn1))/dt); // speed limit constraint
            if (atemp<a)
            {
            	a = atemp;
            }

          	cout << " j = " << (atemp - ap)/dt << " a = " << a << " ds = " << ds <<
          			" vn2 = " << vn2 << " vn1 = " << vn1 << " ap = " << ap <<
          			endl;
            cout << endl;
          	ref_vel = (vn1 + a * dt)/0.44704;

          	json msgJson;

          	vector<double> ptsx;
          	vector<double> ptsy;


          	// ref x, y, yaw states
          	double ref_x = car_x;
          	double ref_y = car_y;
          	double ref_yaw = deg2rad(car_yaw);

          	cout << "PREVIOUS POINTS" << endl;
          	for(int i = 0; i < previous_path_x.size(); i++)
          	{
          		cout << "i = " << i << " x = " << previous_path_x[i] <<
          				" y = " << previous_path_y[i] << endl;
          	}

          	// if previous size is almost empty, use the car as starting reference
          	if(prev_size< 2)
          	{
          		// use two points that make the path tangent to the car
          		double prev_car_x = car_x - cos(car_yaw);
          		double prev_car_y = car_y - sin(car_yaw);

          		ptsx.push_back(prev_car_x);
          		ptsx.push_back(car_x);

          		ptsy.push_back(prev_car_y);
          		ptsy.push_back(car_y);
          	}
          	// use the previous path's end points as the starting reference
          	else
          	{
          		//Redefine reference state as previous path end point
          		ref_x = previous_path_x[prev_size - 1];
          		ref_y = previous_path_y[prev_size - 1];

          		double ref_x_prev = previous_path_x[prev_size - 2];
          		double ref_y_prev = previous_path_y[prev_size - 2];
          		ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

          		// use two points that make the path tangent to the previous path's end point
          		ptsx.push_back(ref_x_prev);
          		ptsx.push_back(ref_x);

          		ptsy.push_back(ref_y_prev);
          		ptsy.push_back(ref_y);
          	}

          	cout << "FROM PREVIOUS POINTS" << endl;
          	for(int i = 0; i < ptsx.size(); i++)
          	{
          		cout << "i = " << i << " x = " << ptsx[i] << " y = " << ptsy[i] << endl;
          	}

           // In Frenet add evenly 30m spaced points ahead of the starting reference
          	vector<double> next_wp0 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp1 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          	vector<double> next_wp2 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);


          	ptsx.push_back(next_wp0[0]);
          	ptsx.push_back(next_wp1[0]);
          	ptsx.push_back(next_wp2[0]);

          	ptsy.push_back(next_wp0[1]);
          	ptsy.push_back(next_wp1[1]);
          	ptsy.push_back(next_wp2[1]);

          	cout << "CAR" << endl;
          	cout <<"s = " << car_s << " d = " << car_d << " x = " << car_x << " y " << car_y <<
          			" yaw = " << car_yaw << endl;;

         	cout << "PREVIOUS POINTS + WAYPOINTS" << endl;
          	for(int i = 0; i < ptsx.size(); i++)
          	{
          		cout << "i = " << i << " x = " << ptsx[i] << " y = " << ptsy[i] << endl;
          	}

      		double vmps = ref_vel * 0.44704;
      	    double dist_inc = vmps * dt;
      		//cout<<"vmps = " << vmps << " dt = " << dt << " vmph = " << " dist_inc = " << dist_inc << endl;
          	for(int i = 0; i < ptsx.size(); i++)
          	{
          		// shift car ref angle to 0
          		double shift_x = ptsx[i] - ref_x;
          		double shift_y = ptsy[i] - ref_y;

          		ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
          		ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
          	}

         	cout << "TRANSFORMED PREVIOUS POINTS + WAYPOINTS" << endl;
          	for(int i = 0; i < ptsx.size(); i++)
          	{
          		cout << "i = " << i << " x = " << ptsx[i] << " y = " << ptsy[i] << endl;
          	}

          	// create a spline
          	tk::spline s;

          	// set (x, y) points to the spline
          	s.set_points(ptsx, ptsy);

          	// Define the actual (x, y) points used for the planner
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// Start with all of the previos path points
          	for(int i = 0; i < previous_path_x.size(); i++)
          	{
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}

          	cout<<"POINTS FROM PREV" << endl;
          	for(int i = 0; i < next_x_vals.size(); i++)
          	{
          		cout << "i = " << i << " x = " << next_x_vals[i] << " y = " << next_y_vals[i] << endl;
          	}

          	double target_x = max(20.0, vmps*3.0);
          	double target_y = s(target_x);
          	double target_dist = sqrt((target_x) * (target_x) + (target_y) * (target_y));

          	double x_add_on = 0;

          	cout<<"TARGET"<< endl;
          	cout<<"x = " << target_x << " y = " << target_y << " dist = " << target_dist << endl;

          	cout << "NEW POINTS UNROTATED" << endl;
          	// fill up the rest of our path planner after filling it with previous points
          	for(int i = 1; i <= 10 - previous_path_x.size(); i++)
          	{
          		double N = (target_dist/(dt * vmps));
          		double x_point = x_add_on + (target_x) / N;
          		double y_point = s(x_point);

          		cout << "i = " << i << " x = " << x_point << " y = " << y_point << endl;
          		x_add_on = x_point;

          		double x_ref = x_point;
          		double y_ref = y_point;

          		// rotate back to normal
          		x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
          		y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

          		x_point += ref_x;
          		y_point += ref_y;

          		next_x_vals.push_back(x_point);
          		next_y_vals.push_back(y_point);
          	}

          	cout<<"NEXT VALS"<< endl;
         	cout << "PREVIOUS POINTS + WAYPOINTS" << endl;
          	for(int i = 0; i < next_x_vals.size(); i++)
          	{
          		cout << "i = " << i << " x = " << next_x_vals[i] << " y = " << next_y_vals[i] << endl;
          	}


          	// END
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          	vn2 = vn1;
          	vn1 = ref_vel * 0.44704;
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
