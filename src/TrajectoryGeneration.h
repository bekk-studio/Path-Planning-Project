#ifndef TRAJECTORYGENERATION_H
#define TRAJECTORYGENERATION_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include "Eigen-3.3/Eigen/Dense"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Const

const int N_SAMPLES = 4;
const vector<double> SIGMA_S =  {5.0, 0.3, 0.05}; // s, s_dot, s_double_dot
const vector<double> SIGMA_D = {0.4, 0.1, 0.01};
const double SIGMA_T = 2.0;

const vector<double> MAX_JERK = {10-1, 2-.2}; // m/s/s/s
const vector<double> MAX_ACCEL= {10-1, 2-.2}; // m/s/s

const double EXPECTED_JERK_IN_ONE_SEC = 2; // m/s/s
const double EXPECTED_ACC_IN_ONE_SEC = 1; // m/s

const double SPEED_LIMIT = 48.5; //MPH
const double MAX_LAT_VELOCITY_CAPABILITY = 2; // m/s
const double VEHICLE_RADIUS = 1.5; // model vehicle as circle to simplify collision detection

double traj_step_time = 0.2; // Delta time for trajectory Generation
double target_velocity = 48/2.23694;


// tweak weights to existing cost functions

const double TIME_COST_WEIGHT = 100;
const double S_DIFF_COST_WEIGHT = 10;
const double D_DIFF_COST_WEIGHT = 50;
const double EFFICIENCY_COST_WEIGHT = 100;
const double COLLISION_COST_WEIGHT = 10000;
const double D_BUFFER_COST_WEIGHT = 50;
const double S_BUFFER_COST_WEIGHT = 100;
const double STAYS_ON_ROAD_COST_WEIGHT = 10000;
const double GO_ON_CENTER_LANE_COST_WEIGHT = 0;  
const double EXCEEDS_SPEED_LIMIT_COST_WEIGHT = 20000;
const double MIN_S_VELOCITY_COST_WEIGHT = 20000;
const double MAX_D_VELOCITY_COST_WEIGHT = 10000;


const double LON_MAX_JERK_COST_WEIGHT = 1000;
const double LON_TOTAL_JERK_COST_WEIGHT = 20;
const double LON_MAX_ACCEL_COST_WEIGHT = 100000;
const double LON_TOTAL_ACCEL_COST_WEIGHT = 20;
const double LAT_MAX_JERK_COST_WEIGHT = 1000;
const double LAT_TOTAL_JERK_COST_WEIGHT = 50;
const double LAT_MAX_ACCEL_COST_WEIGHT = 100000;
const double LAT_TOTAL_ACCEL_COST_WEIGHT = 10;


// functions
vector<vector<double>> PTG(vector<double> start_s, vector<double> start_d, int target_vehicle, vector<double> delta, 
double T, map<int,vector<vector<double>>> predictions);
double calculate_cost(vector<vector<double>> trajectory, double goal_t, 
map<int,vector<vector<double>>> predictions, vector<double> goal_s, vector<double> goal_d);
vector<double> JMT(vector<double> start, vector <double> end, double T);
double logistic(double x);
vector<double> state_at(vector<double> coefficients, double t);
double time_cost(vector<vector<double>> trajectory, double goal_t);
double s_diff_cost(vector<vector<double>>trajectory, vector<double> goal_s);
double d_diff_cost(vector<vector<double>>trajectory, vector<double> goal_d);
double collision_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions);
double d_buffer_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions);
double s_buffer_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions);
double stays_on_road_cost(vector<vector<double>>trajectory);
double exceeds_speed_limit_cost(vector<vector<double>>trajectory);
double min_s_velocity_cost(vector<vector<double>>trajectory);
double max_d_velocity_cost(vector<vector<double>>trajectory);
double efficiency_cost(vector<vector<double>>trajectory, vector<double> goal_s);
double max_accel_cost(vector<vector<double>>trajectory, int axe);
double max_jerk_cost(vector<vector<double>>trajectory, int axe);
double total_jerk_cost(vector<vector<double>>trajectory, int axe);


vector<vector<double>> PTG(vector<double> start_s, vector<double> start_d, int target_vehicle, vector<double> delta, 
double T, map<int,vector<vector<double>>> predictions)
{
	/***
    Finds the best trajectory according to WEIGHTED_COST_FUNCTIONS (global).

    arguments:
     start_s - [s, s_dot, s_ddot]

     start_d - [d, d_dot, d_ddot]

     target_vehicle - id of leading vehicle (int) which can be used to retrieve
       that vehicle from the "predictions" dictionary. This is the vehicle that 
       we are setting our trajectory relative to.

     delta - a length 6 array indicating the offset we are aiming for between us
       and the target_vehicle. So if at time 5 the target vehicle will be at 
       [100, 10, 0, 0, 0, 0] and delta is [-10, 0, 0, 4, 0, 0], then our goal 
       state for t = 5 will be [90, 10, 0, 4, 0, 0]. This would correspond to a 
       goal of "follow 10 meters behind and 4 meters to the right of target vehicle"

     T - the desired time at which we will be at the goal (relative to now as t=0)

     predictions - dictionary of {v_id : vehicle }. 

    return:
     (best_s, best_d, best_t) where best_s are the 6 coefficients representing s(t)
     best_d gives coefficients for d(t) and best_t gives duration associated w/ 
     this trajectory.
    ***/
	
	int indT = T/traj_step_time; //index of T in the prediction list
	
	// goal velocity definition
	double goal_velocity;
	double goal_position;
	if(target_vehicle == -1) //If there are no target vehicle - 
	// target position will be about start position + 70m 
	{goal_velocity = target_velocity;
	 goal_position = start_s[0] + 70;}
	else // if there are one target vehicle
	// target position will be target vehicle position at time T and velocity
	{goal_velocity = predictions[target_vehicle][indT][1];
	 if(goal_velocity > target_velocity){goal_velocity = target_velocity;} // check the outlier
	 if(goal_velocity < 0){goal_velocity =0;} // check the outlier
	 goal_position = predictions[target_vehicle][indT][0];
	 if(goal_position < start_s[0] + 20){goal_position = start_s[0] + 20;}}	// check the outlier
	
	
	vector<double> goal_s = {goal_position + delta[0], goal_velocity + delta[1] , delta[2]};
	vector<double> goal_d = {delta[3], delta[4] , delta[5]};
	
	// generate alternative goals with some errors defined by normal distribution
	// And generate trajectory with different time (step = 200ms) to reach the goal 
	vector<vector<vector<double>>> all_goals;
    double timestep = 0.2;
    double t = T - 3*timestep;
    
	unsigned seed = (unsigned)start_s[0];
	default_random_engine generator (seed);
	
	while (t <= T + 10*timestep){
        vector<vector<double>> goal = {goal_s, goal_d, {t}};
		all_goals.push_back(goal);
		
		normal_distribution<double> distribution_s(goal_s[0], SIGMA_S[0]);
		normal_distribution<double> distribution_s_dot(goal_s[1], SIGMA_S[1]);
		normal_distribution<double> distribution_s_ddot(goal_s[2], SIGMA_S[2]);
		normal_distribution<double> distribution_d(goal_d[0], SIGMA_D[0]);
		normal_distribution<double> distribution_d_dot(goal_d[1], SIGMA_D[1]);
		normal_distribution<double> distribution_d_ddot(goal_d[2], SIGMA_D[2]);
		
		
		for (int i = 0; i < N_SAMPLES; i++){
			vector<vector<double>> goal = 
			{{distribution_s(generator), distribution_s_dot(generator), distribution_s_ddot(generator)},
			 {distribution_d(generator), distribution_d_dot(generator), distribution_d_ddot(generator)},
			 {t}};
			all_goals.push_back(goal);
		}
		
        t += timestep;
	}
	
	// find best trajectory according to different cost functions
    vector<vector<double>> best;
	double min_cost = 1000000000;

	for (int i = 0; i < all_goals.size(); i++)
	{    
		vector<double> alt_s_goal = all_goals[i][0];
		vector<double> alt_d_goal = all_goals[i][1];
		double alt_t = all_goals[i][2][0];
        vector<double> s_coefficients = JMT(start_s, alt_s_goal, alt_t);
        vector<double> d_coefficients = JMT(start_d, alt_d_goal, alt_t);
        vector<vector<double>> trajectory = {s_coefficients, d_coefficients, {alt_t}};
		
		double cost = calculate_cost(trajectory, T, predictions, goal_s, goal_d);
		if (cost < min_cost){
			min_cost = cost;
			best = trajectory;
		}	
	}
	
	return best;	
}

double calculate_cost(vector<vector<double>> trajectory, double goal_t, 
map<int,vector<vector<double>>> predictions, vector<double> goal_s, vector<double> goal_d)
{
	double cost = 0;
	
	cost += TIME_COST_WEIGHT * time_cost(trajectory, goal_t);
	cost += S_DIFF_COST_WEIGHT * s_diff_cost(trajectory, goal_s);
	cost += D_DIFF_COST_WEIGHT * d_diff_cost(trajectory, goal_d);
    cost += EFFICIENCY_COST_WEIGHT * efficiency_cost(trajectory, goal_s);
	cost += COLLISION_COST_WEIGHT * collision_cost(trajectory, predictions);
	cost += D_BUFFER_COST_WEIGHT * d_buffer_cost(trajectory, predictions);
	cost += S_BUFFER_COST_WEIGHT * s_buffer_cost(trajectory, predictions);
	cost += STAYS_ON_ROAD_COST_WEIGHT * stays_on_road_cost(trajectory);
	cost += EXCEEDS_SPEED_LIMIT_COST_WEIGHT * exceeds_speed_limit_cost(trajectory);
	cost += MIN_S_VELOCITY_COST_WEIGHT * min_s_velocity_cost(trajectory);
	cost += MAX_D_VELOCITY_COST_WEIGHT * max_d_velocity_cost(trajectory);
	
	cost += LON_MAX_JERK_COST_WEIGHT * max_jerk_cost(trajectory,0);
	cost += LON_TOTAL_JERK_COST_WEIGHT * total_jerk_cost(trajectory,0);
	cost += LON_MAX_ACCEL_COST_WEIGHT * max_accel_cost(trajectory,0);
	cost += LAT_MAX_JERK_COST_WEIGHT * max_jerk_cost(trajectory,1);
	cost += LAT_TOTAL_JERK_COST_WEIGHT * total_jerk_cost(trajectory,1);
	cost += LAT_MAX_ACCEL_COST_WEIGHT * max_accel_cost(trajectory,1);
	
    return cost;
}

vector<double> JMT(vector<double> start, vector <double> end, double T)
{
    /***
    Calculate the Jerk Minimizing Trajectory that connects the initial state
    to the final state in time T.

    INPUTS

    start - the vehicles start location given as a length three array
        corresponding to initial values of [s, s_dot, s_double_dot]

    end   - the desired end state for vehicle. Like "start" this is a
        length three array.

    T     - The duration, in seconds, over which this maneuver should occur.

    OUTPUT 
    an array of length 6, each value corresponding to a coefficent in the polynomial 
    s(t) = a_0 + a_1 * t + a_2 * t**2 + a_3 * t**3 + a_4 * t**4 + a_5 * t**5

    EXAMPLE

    > JMT( [0, 10, 0], [10, 10, 0], 1)
    [0.0, 10.0, 0.0, 0.0, 0.0, 0.0]
    ***/
    MatrixXd A = MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
			    3*T*T, 4*T*T*T,5*T*T*T*T,
			    6*T, 12*T*T, 20*T*T*T;
		
	MatrixXd B = MatrixXd(3, 1);	    
	B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
			    end[1]-(start[1]+start[2]*T),
			    end[2]-start[2];
			    
	MatrixXd Ai = A.inverse();
	
	MatrixXd C = Ai*B;
	
	vector <double> result = {start[0], start[1], .5*start[2]};
	for(int i = 0; i < C.size(); i++)
	{
	    result.push_back(C.data()[i]);
	}
	
    return result; 
}

double logistic(double x)
{
    /***
    A function that returns a value between 0 and 1 for x in the 
    range [0, infinity] and -1 to 1 for x in the range [-infinity, infinity].

    Useful for cost functions.
    ***/
    return 2.0 / (1 + exp(-x)) - 1.0;
}

vector<double> state_at(vector<double> coefficients, double t)
{
	/***
	From polynomial coefficients, return state (s, s_dot, s_ddot) or (d, d_dot, d_ddot) at t time.
	***/
	double x = coefficients[0] + coefficients[1]*t + coefficients[2]*t*t 
	+ coefficients[3]*t*t*t + coefficients[4]*t*t*t*t + coefficients[5]*t*t*t*t*t;
	double x_dot = coefficients[1] + 2*coefficients[2]*t + 3*coefficients[3]*t*t 
	+ 4*coefficients[4]*t*t*t + 5*coefficients[5]*t*t*t*t;
	double x_ddot = 2*coefficients[2] + 2*3*coefficients[3]*t + 3*4*coefficients[4]*t*t 
	+ 4*5*coefficients[5]*t*t*t;
	double x_dddot = 2*3*coefficients[3] + 2*3*4*coefficients[4]*t + 3*4*5*coefficients[5]*t*t;
	
	return {x, x_dot, x_ddot, x_dddot};
}

// COST FUNCTIONS
double time_cost(vector<vector<double>> trajectory, double goal_t)
{
    /***
    Penalizes trajectories that span a duration which is longer or 
    shorter than the duration requested.
    ***/
    double t = trajectory[2][0];
    //return logistic(abs(t-goal_t) / goal_t);
	return logistic((t-goal_t) / goal_t);
}

double s_diff_cost(vector<vector<double>>trajectory, vector<double> goal_s)
{
    /***
    Penalizes trajectories whose s coordinate (and derivatives) 
    differ from the goal.
    ***/

	vector<double> coefficients_s = trajectory[0];
    double t = trajectory[2][0];

    vector<double> S = state_at(coefficients_s, t);
    double cost = 0;
    for (int i = 0; i < 3; i++)
	{
        cost += logistic(abs(S[i] - goal_s[i])/SIGMA_S[i]);
	}
    return cost;
}

double d_diff_cost(vector<vector<double>>trajectory, vector<double> goal_d)
{
    /***
    Penalizes trajectories whose d coordinate (and derivatives) 
    differ from the goal.
    ***/

	vector<double> coefficients_d = trajectory[1];
    double t = trajectory[2][0];

    vector<double> D = state_at(coefficients_d, t);
    double cost = 0;
    for (int i = 0; i < 3; i++)
	{
        cost += logistic(abs(D[i] - goal_d[i])/SIGMA_D[i]);
	}
    return cost;
}

double collision_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions)
{
    /***
    Binary cost function which penalizes collisions.
    ***/
    double t = trajectory[2][0];
	int indt = t / traj_step_time;
	double nearest = 999999;
	for (map<int,vector<vector<double>>>::iterator it=predictions.begin(); it!=predictions.end(); ++it)
	{
		vector<vector<double>> trajectory_of_vehicles = it->second;
		for (int i = 0; i <= indt; i++)
		{
			double dist = sqrt((trajectory_of_vehicles[i][0] - state_at(trajectory[0], double(i)*traj_step_time)[0])
							  *(trajectory_of_vehicles[i][0] - state_at(trajectory[0], double(i)*traj_step_time)[0])
							 + (trajectory_of_vehicles[i][3] - state_at(trajectory[1], double(i)*traj_step_time)[0])
							  *(trajectory_of_vehicles[i][3] - state_at(trajectory[1], double(i)*traj_step_time)[0]));
			
			if (dist < nearest){nearest = dist;}
		}
		
	}
    if (nearest < 2*VEHICLE_RADIUS){return 1.0;}
    else {return 0.0;}
}

double d_buffer_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions)
{
    /***
    Penalizes getting close to other vehicles.
    ***/
    double t = trajectory[2][0];
	int indt = t / traj_step_time;
	double nearest = 2*VEHICLE_RADIUS + 1;
	for (map<int,vector<vector<double>>>::iterator it=predictions.begin(); it!=predictions.end(); ++it)
	{
		vector<vector<double>> trajectory_of_vehicles = it->second;
		for (int i = 0; i <= indt; i+=5)
		{
			if (abs(trajectory_of_vehicles[i][0] - state_at(trajectory[0], double(i)*traj_step_time)[0]) < 2*VEHICLE_RADIUS)
			{
				double dist = abs(trajectory_of_vehicles[i][3] - state_at(trajectory[1], double(i)*traj_step_time)[0]);		
				if (dist < nearest){nearest = dist;}
			}
		}
	}
    return logistic((2*VEHICLE_RADIUS + 1 - nearest) / (2*VEHICLE_RADIUS + 1));
}

double s_buffer_cost(vector<vector<double>>trajectory, map<int,vector<vector<double>>> predictions)
{
    
	/***
    Penalizes getting close to other vehicle if in the same lane.
	buffer security is fixed to 20m
    ***/
	double t = trajectory[2][0];
	int indt = t / traj_step_time;
	double nearest = 20;
	for (map<int,vector<vector<double>>>::iterator it=predictions.begin(); it!=predictions.end(); ++it)
	{
		vector<vector<double>> trajectory_of_vehicles = it->second;
		for (int i = 0; i <= indt; i+=5)
		{
			if (abs(trajectory_of_vehicles[i][3] - state_at(trajectory[1], double(i)*traj_step_time)[0]) < 2*VEHICLE_RADIUS)
			{
				double dist = abs(trajectory_of_vehicles[i][0] - state_at(trajectory[0], double(i)*traj_step_time)[0]);		
				if (dist < nearest){nearest = dist;}
			}
		}
	}
    return logistic((20 - nearest) / 20);
}
    
double stays_on_road_cost(vector<vector<double>>trajectory)
{
    /***
    Penalizes trajectories too close to the limit of the road
    ***/
	
	double result = 0;
	double t = trajectory[2][0];
	int indt = t / traj_step_time;
	for (int i = 0; i <= indt; i+=5)
		{
			if ((state_at(trajectory[1], double(i)*traj_step_time)[0] < VEHICLE_RADIUS + 0.2) 
			  | (state_at(trajectory[1], double(i)*traj_step_time)[0] > 12 - (VEHICLE_RADIUS + 0.5))) 
			{result = 1;}
		}
	return result;
}


double exceeds_speed_limit_cost(vector<vector<double>>trajectory)
{
    /***
    Penalizes trajectories which exceeds the speed limit
    ***/
	
	double result = 0;
	double t = trajectory[2][0];
	int indt = t / traj_step_time;
	for (int i = 0; i <= indt; i++)
		{
			double speed = sqrt(state_at(trajectory[0], double(i)*traj_step_time)[1]
			                  * state_at(trajectory[0], double(i)*traj_step_time)[1]
			                  + state_at(trajectory[1], double(i)*traj_step_time)[1]
							  * state_at(trajectory[1], double(i)*traj_step_time)[1]);
			
			if (speed >= (SPEED_LIMIT)/2.23694){result = 1;}
		}
	return result;
}

double min_s_velocity_cost(vector<vector<double>>trajectory)
{
	/***
    Penalizes trajectories with with negative velocity according to s axe 
    ***/
	
	double result = 0;
	double t = trajectory[2][0];
	int indt = t / traj_step_time;
	for (int i = 0; i <= indt; i+=5)
		{
			double speed_s = state_at(trajectory[0], double(i)*traj_step_time)[1];
			if (speed_s < 0){result = 1;}
		}
	return result;
}

double max_d_velocity_cost(vector<vector<double>>trajectory)
{
	/***
    Penalizes trajectories which exceeds lateral velocity according to d axe 
    ***/
	double result = 0;
	double t = trajectory[2][0];
	int indt = t / traj_step_time;
	for (int i = 0; i <= indt; i+=5)
		{
			double speed_d = state_at(trajectory[1], double(i)*traj_step_time)[1];
			if (speed_d > 2){result = 1;}
		}
	return result;
}
	
double efficiency_cost(vector<vector<double>>trajectory, vector<double> goal_s)
{
    /***
    Rewards high average speeds.
    ***/
	vector<double> coefficients_s = trajectory[0];
	double ref_t = trajectory[2][0];
    double t = 0;
	double cost = 0;
	double targ_v = goal_s[1];
	while(t < ref_t)
	{
		double car_v = state_at(coefficients_s, t)[1];
		cost +=  abs(targ_v - car_v) / (ref_t / (5 * traj_step_time));
		t += 5 * traj_step_time;
	}
    return logistic(cost / targ_v);
}

    
double max_accel_cost(vector<vector<double>>trajectory, int axe)
{
    /***
    Penalizes trajectories which exceeds max acceleration according to one axe (s or d)
    ***/
	
	vector<double> coefficients = trajectory[axe];
	double t = trajectory[2][0];
    int indt = t / traj_step_time;
	double max_acc = 0;
    for (int i = 0; i <= indt; i++)
		{
			double acc = state_at(coefficients, double(i)*traj_step_time)[2];
			if (abs(acc) > max_acc)
			{
				max_acc = abs(acc);
            }				
		}
	if (max_acc > MAX_ACCEL[axe]){return 1;}
    else{return 0;}
}
    

double max_jerk_cost(vector<vector<double>>trajectory, int axe)
{
    /***
    Penalizes trajectories which exceeds max jerk according to one axe (s or d)
    ***/
	
	vector<double> coefficients = trajectory[axe];
	double t = trajectory[2][0];
    int indt = t / traj_step_time;
	double max_jerk = 0;
    for (int i = 0; i <= indt; i++)
		{
			double jerk = state_at(coefficients, double(i)*traj_step_time)[3];
			if (abs(jerk) > max_jerk)
			{
				max_jerk = abs(jerk);
            }				
		}
	if (max_jerk > MAX_JERK[axe]){return 1;}
    else{return 0;}
}

double total_jerk_cost(vector<vector<double>>trajectory, int axe)
{
    /***
    Penalizes trajectories which exceeds a jerk threshold according to one axe (s or d)
    ***/
	
	vector<double> coefficients = trajectory[axe];
	double t = trajectory[2][0];
    int indt = t / traj_step_time;
	double total_jerk = 0;
    for (int i = 0; i <= indt; i++)
		{
			total_jerk += abs(state_at(coefficients, double(i)*traj_step_time)[3]); 
		}
    total_jerk /= double(indt+1);
    return logistic(total_jerk / EXPECTED_JERK_IN_ONE_SEC - 2);
}

#endif