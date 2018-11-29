#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <time.h>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "GNaiveBayes.h"
#include "BehaviourPlanning.h"
#include "TrajectoryGeneration.h"


using namespace std;

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

	double angle = theta - heading;
	
	if (angle > pi()){angle -= 2*pi();}
	if (angle < -pi()){angle += 2*pi();}	

	if(abs(angle) > pi()/2)
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

  /***
  Waypoints, constants and containers 
  ***/
  
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
	map_waypoints_x.push_back(x + 6 * d_x); // (+ 6 * dx) I center reference way points on lane #1, for a better resolution
  	map_waypoints_y.push_back(y + 6 * d_y); // (+ 6 * dx) I center reference way points on lane #1, for a better resolution
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  // Interpolation of waypoints according to s
  tk::spline spline_waypoints_x;
  tk::spline spline_waypoints_y;
  tk::spline spline_waypoints_dtheta;
  
  // I add the last point
  map_waypoints_s.push_back(max_s);
  map_waypoints_x.push_back(map_waypoints_x[0]);
  map_waypoints_y.push_back(map_waypoints_y[0]);
  map_waypoints_dx.push_back(map_waypoints_dx[0]);
  map_waypoints_dy.push_back(map_waypoints_dy[0]);
  
  // dtheta will be the direction of d's axe according to s
  vector<double> map_waypoints_dtheta;
  for (int i = 0; i < map_waypoints_dx.size(); i++)
  {
	  double dtheta = atan2(map_waypoints_dy[i], map_waypoints_dx[i]);
	  if(i > 0)
	  {
		  while(dtheta - map_waypoints_dtheta.back() > M_PI){dtheta -= 2*M_PI;}
		  while(dtheta - map_waypoints_dtheta.back() < -M_PI){dtheta += 2*M_PI;}
	  }
	  map_waypoints_dtheta.push_back(dtheta);
  }
  
  spline_waypoints_x.set_points(map_waypoints_s, map_waypoints_x);
  spline_waypoints_y.set_points(map_waypoints_s, map_waypoints_y);
  spline_waypoints_dtheta.set_points(map_waypoints_s, map_waypoints_dtheta);
  
  // Delta time
  const double dt = 0.02; //20 ms
  
  // I add 2 vector container to memorize complete state (position, velocity, acceleration) of end of previous path
  vector<double> advanced_end_path_s;
  vector<double> advanced_end_path_d;
  
  //Vehicle object for Behaviour planning
  BPVehicle ego_b(1, 0, 0, 0, 2);
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
			  &dt, &max_s, &advanced_end_path_s, &advanced_end_path_d,
			  &spline_waypoints_x, &spline_waypoints_y, &spline_waypoints_dtheta,
			  &ego_b] 
			   (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"]; // degree
          	double car_speed = j[1]["speed"]; //mph

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];
			

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	
			int prev_size = previous_path_x.size();
			int tr_length = 120; // total number of points sent in trajectory
			int points_to_add = tr_length - prev_size;
			
			
			/*** 
			Prediction
			***/
			
			// Prediction with 1s step for behaviour planning, and predictions with 0.2s step for trajectory generation
			map<int,vector<vector<int>>> predictions_behaviourP;   // trajectory by second with trajectory {lane, s}
			map<int,vector<vector<double>>> predictions_trajectoryG;     // trajectory by 200 ms with trajectory {s, d}
			
			
			// erase too far sensor_fusion detection
			vector<int> to_erase;
			for (int i = 0; i < sensor_fusion.size(); i++)
			{
				double s_other = sensor_fusion[i][5];
				double d_other = sensor_fusion[i][6];
				if ((abs(s_other - car_s) > 200 && max_s - abs(s_other - car_s) > 200) 
					|| d_other < 0 || d_other > 12)
					{to_erase.push_back(i);}
			}
			for (int i = 0; i < to_erase.size(); i++){sensor_fusion.erase(sensor_fusion.begin()+to_erase[i]);}
			
			
			// Sort sensor_fusion id by descending s order and calculate vs_other and vd_other
			vector<vector<double>> sorted_sensor_fusion;
			vector<vector<double>>::iterator it_sensor;
			for (int i = 0; i < sensor_fusion.size(); i++){
				vector<double> sensor_sample = sensor_fusion[i];
				
				double dtheta = spline_waypoints_dtheta(sensor_sample[5]);   // angle for d axes
				double vs_other = -1 * sensor_sample[3] * sin(dtheta) + sensor_sample[4] * cos(dtheta);  // vs = - vx . sin(dtheta) + vy . cos(dtheta) 
				if(vs_other > 25){vs_other = 25;}
				sensor_sample.push_back(vs_other);
				double vd_other = sensor_sample[3] * cos(dtheta) + sensor_sample[4] * sin(dtheta);  // vd = vx . cos(dtheta) + vy . sin(dtheta)
				if(vd_other > 2){vd_other = 2;}
				if(vd_other < -2){vd_other = -2;}
				sensor_sample.push_back(vd_other);
				
				if(car_s > 6000 && sensor_sample[5] < 1000){sensor_sample[5] += max_s;}
				if(car_s < 1000 && sensor_sample[5] > 6000){sensor_sample[5] -= max_s;}
				if(sorted_sensor_fusion.size() == 0){sorted_sensor_fusion.push_back(sensor_sample);}
				else{
					int pos = 0;
					for (int n = 0; n < sorted_sensor_fusion.size(); n++){
						if(sensor_sample[5] < sorted_sensor_fusion[n][5]){pos = n+1;}
					}
					it_sensor = sorted_sensor_fusion.begin() + pos;
					sorted_sensor_fusion.insert(it_sensor, sensor_sample);
				}
			}

			
			map<int,vector<vector<int>>>::iterator it;
			// Prediction for Behaviour Planning
			bool ego_inserted = false;
			for (int i = 0; i < sorted_sensor_fusion.size(); i++){
				
				// We insert ego vehicle after ahead vehicles and before behind vehicle
				if (ego_inserted == false && (sorted_sensor_fusion[i][5] < car_s || i == sorted_sensor_fusion.size() - 1)){
					int lane_ego;
					if (car_d < 4){lane_ego = 0;}
					else if (car_d < 8){lane_ego = 1;}
					else{lane_ego = 2;}
					BPVehicle other(lane_ego, car_s, car_speed/2.24, 0, lane_ego);
					vector<vector<int>> trajectory_ego = other.trajectory_generation(predictions_behaviourP, "CS");
					trajectory_ego.push_back({trajectory_ego[5][0], 2*trajectory_ego[5][1]-trajectory_ego[4][1]});
					predictions_behaviourP.insert(pair<int,vector<vector<int>>>(-1, trajectory_ego));
					ego_inserted = true;
				}
				
				// state prediction of other vehicle
				vector<double> sample = {sorted_sensor_fusion[i][5], sorted_sensor_fusion[i][6],
										sorted_sensor_fusion[i][7], sorted_sensor_fusion[i][8]}; //{s_other, d_other, vs_other, vd_other}
				string state_predict = predict(sample);
				
				// lane association of other vehicules
				int lane_other;
				double d_other = sorted_sensor_fusion[i][6]; // a garder
				int id_other = sorted_sensor_fusion[i][0];
				
				if (d_other < 4){lane_other = 0;}
				else if (d_other < 8){lane_other = 1;}
				else{lane_other = 2;}
				
				// trajectories prediction by seconds for Behaviour Planning
				BPVehicle other(lane_other, sorted_sensor_fusion[i][5], sorted_sensor_fusion[i][7], 0, lane_other);
				vector<vector<int>> trajectory_bp;
				if ((state_predict == "right") && (lane_other == 0 || lane_other == 1) && (fmod(d_other, 4) > 3)){
					trajectory_bp = other.trajectory_generation(predictions_behaviourP, "LCR"); 
				}
				else if ((state_predict == "left") && (lane_other == 1 || lane_other == 2) && (fmod(d_other, 4) < 1)){
					trajectory_bp = other.trajectory_generation(predictions_behaviourP, "LCL");
				}
				else {
					trajectory_bp = other.trajectory_generation(predictions_behaviourP, "CS");
				}
				trajectory_bp.push_back({trajectory_bp[5][0], 2*trajectory_bp[5][1]-trajectory_bp[4][1]});	
				predictions_behaviourP.insert(pair<int,vector<vector<int>>>(id_other, trajectory_bp));
				
			}
			

			// we delete the ego vehicle
			predictions_behaviourP.erase(-1);
			
			
			// Prediction for Trajectory Generation - steptime = 200ms
			// It will be an interpolation of Prediction for behaviour palnning
			int ind = 0;
			for (it=predictions_behaviourP.begin(); it!=predictions_behaviourP.end(); ++it) 
			{	
				//use tk::spline for s
				tk::spline spline_s;
				vector<double> s_time_vector = {0, 1, 2, 3, 4, 5, 6}; //in second
				vector<double> s_vector;
				for(int i =0; i < it->second.size(); i++){
					s_vector.push_back(it->second[i][1]);
				}
				spline_s.set_points(s_time_vector, s_vector);
				
				//use tk::spline for d
				tk::spline spline_d;
				vector<double> d_time_vector = {0, 1, 4, 5}; //in second
				vector<double> d_vector = {sorted_sensor_fusion[ind][6], (double)(2+4*it->second[1][0]), 
										  (double)(2+4*it->second[4][0]), (double)(2+4*it->second[5][0])};
				ind++;
				spline_d.set_points(d_time_vector, d_vector);
				
				// From time of the end of previous trajectory, to 5s later
				vector<vector<double>> trajectory_tg;
				double step_time = 0.2;
				
				for (int t = 0; t <= 25; t++)
				{
					double T = prev_size * dt + t * step_time;
					double st = spline_s(T);
					double vst = (spline_s(T+step_time) - st) / step_time;
					double ast = ((spline_s(T+2*step_time) - spline_s(T+step_time)) / step_time - vst) / step_time;
					double dt = spline_d(T);
					double vdt = (spline_d(T+step_time) - dt) / step_time;
					double adt = ((spline_d(T+2*step_time) - spline_d(T+step_time)) / step_time - vdt) / step_time;
					trajectory_tg.push_back({st, vst, ast, dt, vdt, adt});
				}	
				predictions_trajectoryG.insert(pair<int,vector<vector<double>>>(it->first, trajectory_tg));
			}
			
			/***
			Behaviour Planning
			***/
			// lane association of ego vehicule
			double lane;
			if (car_d < 4){lane = 0;}
			else if (car_d < 8){lane = 1;}
			else{lane = 2;}
			
			ego_b.update_state(predictions_behaviourP, (int)lane, (int)car_s, (int)(car_speed / 2.24));
			string state = ego_b.state;
			
			/***
			Trajectory Generation
			***/
			
			//From here, ref corespond to end of previous trajectory - Time translation to end of previous path
			double ref_x;
			double ref_y;
			double ref_x_prev;
			double ref_y_prev;
			double ref_s;
			
			double ref_d;
			double ref_s_dot;
			double ref_d_dot;
			double ref_s_ddot;
			double ref_d_ddot;
			
			//For the first path (no previous path), i have to configure the vector
			if (prev_size == 0)
			{
				advanced_end_path_s = {car_s, 0, 0};
				advanced_end_path_d = {car_d, 0, 0};
				ref_s = car_s;
			}
			else
			{
				ref_s = end_path_s;
			}
			
			
			
			// Search and find target vehicle according to the planned state
			int target_vehicle = -1; // -1 if no target vehicle
			double distance_target = 50;
			if(state == "PLCL" || state == "PLCR"){distance_target = 20;}
			double s_target;
			double d_target;
			bool ahead_target = false;
	        
			for (map<int,vector<vector<double>>>::iterator it=predictions_trajectoryG.begin(); it!=predictions_trajectoryG.end(); ++it) 
			{	
				s_target = it->second[0][0];
				d_target = it->second[20][3];
				
				if ((state == "KL") && (d_target > 4*lane) && (d_target < 4*(lane+1)))
				{
					if((s_target - ref_s > 0) && (s_target - ref_s < distance_target))
					{
						distance_target = s_target - ref_s;
						target_vehicle = it->first;
					}
				}
				
				if ((state == "LCL") && (d_target > 4*(lane-1)) && (d_target < 4*lane))
				{
					if((s_target - ref_s > 0) && (s_target - ref_s < distance_target))
					{
						distance_target = s_target - ref_s;
						target_vehicle = it->first;
					}
				}
				
				if ((state == "LCR") && (d_target > 4*(lane+1)) & (d_target < 4*(lane+2)))
				{
					if((s_target - ref_s > 0) && (s_target - ref_s < distance_target))
					{
						distance_target = s_target - ref_s;
						target_vehicle = it->first;
					}
				}
				
				if ((state == "PLCL") && (d_target > 4*(lane-1)) && (d_target < 4*lane))
				{
					if(abs(s_target - ref_s) < distance_target)
					{
						distance_target = abs(s_target - ref_s);
						target_vehicle = it->first;
					}
					if((s_target - ref_s < 50) && (s_target - ref_s > 20))
					{
						ahead_target = true;
					}
					
				}
				
				if ((state == "PLCR") && (d_target > 4*(lane+1)) && (d_target < 4*(lane+2)))
				{
					if(abs(s_target - ref_s) < distance_target)
					{
						distance_target = abs(s_target - ref_s);
						target_vehicle = it->first;
					}
					if((s_target - ref_s < 50) && (s_target - ref_s > 20))
					{
						ahead_target = true;
					}
				}
			}
			
			// if there are no ahead vehicle for PLCR & PLCL, put target_vehicle to -1 in order to try to pass ahead 
			if (((state == "PLCL") || (state == "PLCR")) && (ahead_target == false)){target_vehicle = -1;}
			
			
			// Setting of delta vector for implementing a goal state for trajectory generation according to planned state
			vector<double> delta;
			
			if(state == "KL"){delta = {-5, 0, 0, 4*lane+2, 0, 0};}
			if(state == "LCL"){delta = {-5, 0, 0, 4*lane-2, 0, 0};}
			if(state == "LCR"){delta = {-5, 0, 0, 4*lane+6, 0, 0};}
			if(state == "PLCL"){delta = {-10, -1, 0, 4*lane+2, 0, 0};}
			if(state == "PLCR"){delta = {-10, -1, 0, 4*lane+2, 0, 0};}
			
			// Trajectory Generation function
			vector<double> ref_vect_s = advanced_end_path_s;
			vector<double> ref_vect_d = advanced_end_path_d;
			
			vector<vector<double>> coefficients = PTG(ref_vect_s, ref_vect_d,
			target_vehicle, delta, 3, predictions_trajectoryG);
			vector<double> coefficients_s = coefficients[0];
			vector<double> coefficients_d = coefficients[1];
			double time = coefficients[2][0];
			
			/***
			To finalize
			***/
			
			// Define the actual points we will use for the plane
			vector<double> next_x_vals;
          	vector<double> next_y_vals;
			
			// Start with all the previous path points from last time
			
			if (prev_size > 0)
			{
				for (int i = 0; i < prev_size; i++)
				{
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
			}
			
			// Add the new trajectory
			for(int i=1; i <= points_to_add; i++)
			{
				double s = state_at(coefficients_s, i*dt)[0];
				if (s > max_s){s -= max_s;} // to manage the s round trip 
				double d = state_at(coefficients_d, i*dt)[0];
				double x = spline_waypoints_x(s) + (d-6)*cos(spline_waypoints_dtheta(s)); // center of waypoints is lane #1
				double y = spline_waypoints_y(s) + (d-6)*sin(spline_waypoints_dtheta(s)); // center of waypoints is lane #1
				next_x_vals.push_back(x);
				next_y_vals.push_back(y);
			}
				
			//memorize the last point of path and s,d position, velocity and acceleration	
			advanced_end_path_s = state_at(coefficients_s, points_to_add*dt);
			advanced_end_path_s.pop_back();
			if (advanced_end_path_s[0] > max_s){advanced_end_path_s[0] -= max_s;}
			advanced_end_path_d = state_at(coefficients_d, points_to_add*dt);
			advanced_end_path_d.pop_back();
			
			// END
			
			
			msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
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
