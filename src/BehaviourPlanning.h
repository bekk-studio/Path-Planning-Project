#ifndef BEHAVIOURPLANNING_H
#define BEHAVIOURPLANNING_H
#include <iostream>
#include <random>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <string>
#include <iterator>

using namespace std;

class BPVehicle {
public:

  int L = 1;

  int preferred_buffer = 20; // impacts "keep lane" behavior.

  int lane;
  
  int target_lane;

  int s;

  int v;

  int a;

  int target_speed = 50 / 2.24;

  int lanes_available = 3;

  int max_acceleration = 10;

  string state;

  /**
  * Constructor
  */
  BPVehicle(int lane, int s, int v, int a, int target_lane);

  /**
  * Destructor
  */
  virtual ~BPVehicle();

  void update_state(map<int,vector < vector<int> > > predictions, int lane, int s, int v);
  
  vector<vector<int>> trajectory_generation(map < int,vector <vector<int> > > predictions,
  string possible_state);
  
  float collision_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory);
  
  float buffer_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory);
  
  float traffic_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory);
  
  void increment(int dt);

  void realize_state(map<int, vector < vector<int> > > predictions);

  void realize_constant_speed(map<int,vector<vector<int> > > predictions, int lane, int s);

  int _max_accel_for_lane(map<int,vector<vector<int> > > predictions, int lane, int s);

  void realize_keep_lane(map<int, vector< vector<int> > > predictions);

  void realize_lane_change(map<int,vector< vector<int> > > predictions, string direction);

  void realize_prep_lane_change(map<int,vector< vector<int> > > predictions, string direction);

};

// constants to balance
const int horizon = 5;
const float TRAFFIC = 100;
const float CHANGE = 15;
const float COLLISION = 1000;
const float UNSAFE = 300;
const float STRESS = 100;


/**
 * Initializes BPVehicle
 */
BPVehicle::BPVehicle(int lane, int s, int v, int a, int target_lane) {

    this->lane = lane;
    this->s = s;
    this->v = v;
    this->a = a;
	this->target_lane = target_lane;
    state = "CS";

}

BPVehicle::~BPVehicle() {}



void BPVehicle::update_state(map<int,vector < vector<int> > > predictions, int lane, int s, int v) {
	/*
    Updates the "state" of the BPVehicle by assigning one of the
    following values to 'self.state':

    "KL" - Keep Lane
     - The BPVehicle will attempt to drive its target speed, unless there is 
       traffic in front of it, in which case it will slow down.

    "LCL" or "LCR" - Lane Change Left / Right
     - The BPVehicle will IMMEDIATELY change lanes and then follow longitudinal
       behavior for the "KL" state in the new lane.

    "PLCL" or "PLCR" - Prepare for Lane Change Left / Right
     - The BPVehicle will find the nearest BPVehicle in the adjacent lane which is
       BEHIND itself and will adjust speed to try to get behind that BPVehicle.

    INPUTS
    - predictions 
    A dictionary. The keys are ids of other BPVehicles and the values are arrays
    where each entry corresponds to the BPVehicle's predicted location at the 
    corresponding timestep. The FIRST element in the array gives the BPVehicle's
    current position. Example (showing a car with id 3 moving at 2 m/s):

    {
      3 : [
        {"lane": 0, "s" : 4},
        {"lane": 0, "s" : 6},
        {"lane": 0, "s" : 8},
        {"lane": 0, "s" : 10},
      ]
    }
    */
    
	/*
	Find the most flowing lane
	If difference is not clear, prefer the target lane selected in the previous step.
	*/
	
	//update sensor measure
	this->lane = lane;
    this->s = s;
    this->v = v;
	
	// remember current lane
	int lane_memory = this->lane;
	
	int target_lane_temp;
	float min_cost = 9999999;
	
	for (int l =0; l < lanes_available; l++)
    {
		// generate a rough idea of what trajectory we would
        // follow IF we chose this lane.
		this->lane = l;
		vector<vector<int>> trajectory = trajectory_generation(predictions, "KL");
		
		// calculate the "cost" associated with that trajectory.
        float cost = 0;
		
		// traffic cost to measure the flowability of the lanes
		cost += TRAFFIC * traffic_cost(predictions, trajectory);
		
		// change cost to boost the previous lane
        if(l != this->target_lane){cost += CHANGE;}
		
		if (cost < min_cost)
        {
            min_cost = cost;
            target_lane_temp = l;
        }
	}
	
	this->target_lane = target_lane_temp;
	// restore lane from memory
	this->lane = lane_memory;
	
	
	/*
	Then according to target lane, update the possible states list. And choose the best one.
	*/
	this->state = "KL";
	min_cost = 500;
	
	if(this->lane != this->target_lane)
	{
		vector<string> possible_states;
		if(this->lane < this->target_lane){possible_states = {"PLCR", "LCR"};}
		if(this->lane > this->target_lane){possible_states = {"PLCL", "LCL"};}
		
		// generate a rough idea of what trajectory we would
        // follow IF we chose this state.
		for (int i =0; i < possible_states.size(); i++)
		{
			vector<vector<int>> trajectory = trajectory_generation(predictions, possible_states[i]);
			
			// calculate the "cost" associated with that trajectory.
			float cost = 0;
			
			// stress cost - prefer to change lane if possible
			cost += STRESS * abs(trajectory[horizon][0] - this->target_lane);
			
			// collision cost 
			cost += COLLISION * collision_cost(predictions, trajectory);
			
			// buffer cost - security distance with the other vehicle
			cost += UNSAFE * buffer_cost(predictions, trajectory);
			
			if (cost < min_cost)
			{
				min_cost = cost;
				this->state = possible_states[i];
			}
		}
	}
}


vector<vector<int>> BPVehicle::trajectory_generation(map<int,vector < vector<int> > > predictions, string possible_state){
	/*
    Generate the trajectory according to state 
    */
    map<int, vector<vector<int> > >::iterator it;  
    
    // remember current state
    vector<int> memories = {this->lane, this->s, this->v, this->a};
    string state_memory = this->state;
    
	// generate a rough idea of what trajectory we would
	// follow IF we chose this state.
	vector<vector<int>> trajectory;
	this->state = possible_state;
	map<int,vector < vector<int> > > predictions_copy = predictions;
	for (int h =0; h <= horizon; h++)
	{
		if (possible_state.compare("LCL") == 0 | possible_state.compare("LCR") == 0){if( h > 0){this->state = "KL";}}
		trajectory.push_back({this->lane, this->s});
		realize_state(predictions_copy);
		increment(1);
		// need to remove first prediction for each other vehicles for each step.
		for (it=predictions_copy.begin(); it!=predictions_copy.end(); ++it)
		{
			it->second.erase(it->second.begin());
		}
    }    
	
    
	// restore state from memories
	this->state = state_memory;
	this->lane = memories[0];
	this->s = memories[1];
	this->v = memories[2];
	this->a = memories[3];
	
	return trajectory;
}          

float BPVehicle::collision_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory) {
	map<int, vector<vector<int> > >::iterator it; 
	
	// collision cost
	float cost = 0;
	
	for (it=predictions.begin(); it!=predictions.end(); ++it)
	{
		bool behind = false;
		bool ahead = false;
		bool collision = false;
		for (int h =0; h <= horizon; h++)
		{
			// ego and other car are on the same position
			if(trajectory[h][0] == it->second[h][0]
			&&  trajectory[h][1] == it->second[h][1]){collision = true;}
			// other car is ahead, then ego pass over in the same lane 
			if(trajectory[h][0] == it->second[h][0] 
			&& trajectory[h][1] < it->second[h][1]){ahead = true;}
			if(ahead && trajectory[h][0] == it->second[h][0] 
			&& trajectory[h][0] > it->second[h][0]){collision = true;}
			// other car is behind , then ego pass behind in the same lane 
			if(trajectory[h][0] == it->second[h][0] 
			&& trajectory[h][1] > it->second[h][1]){behind = true;}
			if(behind && trajectory[h][0] == it->second[h][0]
			&& trajectory[h][1] < it->second[h][1]){collision = true;}
		}
		if (collision){cost += 1;}
	}
	return cost;
}

float BPVehicle::buffer_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory) {
	map<int, vector<vector<int> > >::iterator it; 
	
	// buffer cost - cost if an other vehicle is less than 20m far
	float cost = 0;
	
	int dist = preferred_buffer;
	for (it=predictions.begin(); it!=predictions.end(); ++it)
	{
		for (int h =0; h <= horizon; h++)
		{
			if(trajectory[h][0] == it->second[h][0]
			&&  abs(trajectory[h][1] - it->second[h][1]) < dist)
			{dist = abs(trajectory[h][1] - it->second[h][1]);}
		}  
	}
	if(dist < preferred_buffer){cost += (float)(preferred_buffer - dist) / preferred_buffer;}
	return cost;
}

float BPVehicle::traffic_cost(map<int,vector < vector<int> > > predictions, vector<vector<int>> trajectory) {
	map<int, vector<vector<int> > >::iterator it; 
	
	// to measure the flowability of the lanes
	float cost = 0;
	
	int dist = 70;
	for (it=predictions.begin(); it!=predictions.end(); ++it)
	{
		if(trajectory[0][0] == it->second[0][0]
			&&  (it->second[0][1] - trajectory[0][1]) >= 0 && (it->second[0][1] - trajectory[0][1]) < dist)
			{dist = it->second[0][1] - trajectory[0][1];} 
	}
	cost += (float)(70 - dist) / 70;
	return cost;
}

void BPVehicle::increment(int dt = 1) {

	this->s += this->v * dt;
    this->v += this->a * dt;
	
}


void BPVehicle::realize_state(map<int,vector < vector<int> > > predictions) {
   
	/*
    Given a state, realize it by adjusting acceleration and lane.
    Note - lane changes happen instantaneously.
    */
    string state = this->state;
    if(state.compare("CS") == 0)
    {
    	realize_constant_speed(predictions, this->lane, this->s);
    }
    else if(state.compare("KL") == 0)
    {
    	realize_keep_lane(predictions);
    }
    else if(state.compare("LCL") == 0)
    {
    	realize_lane_change(predictions, "L");
    }
    else if(state.compare("LCR") == 0)
    {
    	realize_lane_change(predictions, "R");
    }
    else if(state.compare("PLCL") == 0)
    {
    	realize_prep_lane_change(predictions, "L");
    }
    else if(state.compare("PLCR") == 0)
    {
    	realize_prep_lane_change(predictions, "R");
    }

}

void BPVehicle::realize_constant_speed(map<int,vector<vector<int> > > predictions, int lane, int s) {
	int max_acc = 0;
	
	map<int, vector<vector<int> > >::iterator it = predictions.begin();
    vector<vector<vector<int> > > in_front;
    while(it != predictions.end())
    {
       
    	int v_id = it->first;
    	
        vector<vector<int> > v = it->second;
        
        if((v[0][0] == lane) && (v[0][1] >= s))
        {
        	in_front.push_back(v);
        }
        it++;
    }
    
    if(in_front.size() > 0)
    {
    	int min_s = 10000;
    	vector<vector<int>> leading = {};
    	for(int i = 0; i < in_front.size(); i++)
    	{
    		if(in_front[i][0][1] < min_s)
    		{
    			min_s = in_front[i][0][1];
    			leading = in_front[i];
    		}
    	}
    	
    	int next_pos = leading[1][1];   
    	int my_next = s + this->v;
    	int separation_next = next_pos - my_next;
    	int available_room = separation_next - preferred_buffer;
    	max_acc = min(max_acc, available_room);
    }
    
    max_acc = max(max_acc, -this->max_acceleration);
	max_acc = max(max_acc, -this->v);
	
    this->a = max_acc;
}

int BPVehicle::_max_accel_for_lane(map<int,vector<vector<int> > > predictions, int lane, int s) {

	int delta_v_til_target = target_speed - v;
    int max_acc = min(max_acceleration, delta_v_til_target);

    map<int, vector<vector<int> > >::iterator it = predictions.begin();
    vector<vector<vector<int> > > in_front;
    while(it != predictions.end())
    {
       
    	int v_id = it->first;
    	
        vector<vector<int> > v = it->second;
        
        if((v[0][0] == lane) && (v[0][1] >= s))
        {
        	in_front.push_back(v);
        }
        it++;
    }
    
    if(in_front.size() > 0)
    {
    	int min_s = 10000;
    	vector<vector<int>> leading = {};
    	for(int i = 0; i < in_front.size(); i++)
    	{
    		if(in_front[i][0][1] < min_s)
    		{
    			min_s = in_front[i][0][1];
    			leading = in_front[i];
    		}
    	}
    	
    	int next_pos = leading[1][1];   
    	int my_next = s + this->v;
    	int separation_next = next_pos - my_next;
    	int available_room = separation_next - preferred_buffer;
    	max_acc = min(max_acc, available_room);
    }
    
    max_acc = max(max_acc, -this->max_acceleration);
	max_acc = max(max_acc, -this->v);
	
    return max_acc;

}

void BPVehicle::realize_keep_lane(map<int,vector< vector<int> > > predictions) {
	this->a = _max_accel_for_lane(predictions, this->lane, this->s);
}

void BPVehicle::realize_lane_change(map<int,vector< vector<int> > > predictions, string direction) {
	int delta = 1;
    if (direction.compare("L") == 0)
    {
    	delta = -1;
    }
    this->lane += delta;
    int lane = this->lane;
    int s = this->s;
    this->a = _max_accel_for_lane(predictions, lane, s);
}

void BPVehicle::realize_prep_lane_change(map<int,vector<vector<int> > > predictions, string direction) {
	int delta = 1;
    if (direction.compare("L") == 0)
    {
    	delta = -1;
    }
    int lane_togo = this->lane + delta;

    map<int, vector<vector<int> > >::iterator it = predictions.begin();
    vector<vector<vector<int> > > at_behind;
    while(it != predictions.end())
    {
    	int v_id = it->first;
        vector<vector<int> > v = it->second;

        if((v[0][0] == lane_togo) && (v[0][1] <= this->s))
        {
        	at_behind.push_back(v);

        }
        it++;
    }
    if(at_behind.size() > 0)
    {
    	int max_s = -1000;
    	vector<vector<int> > nearest_behind = {};
    	for(int i = 0; i < at_behind.size(); i++)
    	{
    		if(at_behind[i][0][1] > max_s)
    		{
    			max_s = at_behind[i][0][1];
    			nearest_behind = at_behind[i];
    		}
    	}
    	int target_vel = nearest_behind[1][1] - nearest_behind[0][1];
    	int delta_v = this->v - target_vel;
    	int delta_s = this->s - nearest_behind[0][1];
    	if(delta_v != 0)
    	{

    		int time = -2 * delta_s/delta_v;
    		int a;
    		if (time == 0)
    		{
    			a = this->a;
    		}
    		else
    		{
    			a = delta_v/time;
    		}
    		a = min(a, this->max_acceleration);
			a = max(a, -this->max_acceleration);
			a = max(a, -this->v); 
			
    		this->a = a;
    	}
    	else
    	{
    		int my_min_acc = max(-this->max_acceleration,-delta_s);
			my_min_acc = max(my_min_acc, -this->v);
    		this->a = my_min_acc;
    	}

    }
	else
	{
		this->a = _max_accel_for_lane(predictions, this->lane, this->s);
	}

}


#endif