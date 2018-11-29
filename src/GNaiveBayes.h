#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>



using namespace std;

// Trained data during Lesson 3
vector<double> p_y = {0.236, 0.52, 0.244};
vector<vector<double>> mean_y = {{1.18644, 2.03147, 9.90349, -1.15607},
		  {0.946154, 1.99702, 9.97501, 0.00544296},
		  {0.759563, 2.07172, 9.98015, 1.1326}};
vector<vector<double>> variance_y = {{0.434166, 0.999321, 0.93146, 0.324336},
			  {0.691972, 0.125624, 1.09383, 0.0248785},
			  {0.477709, 1.06481, 0.977555, 0.311976}};

vector<string> possible_labels = {"left","keep","right"};



string predict(vector<double> sample)
{
	/*
		Once trained, this method is called and expected to return 
		a predicted behavior for the given observation.

		INPUTS

		observation - a 4 tuple with s, d, s_dot, d_dot.
		  - Example: [3.5, 0.1, 8.5, -0.2]

		OUTPUT

		A label representing the best guess of the classifier. Can
		be one of "left", "keep" or "right".
		"""
		# TODO - complete this
	*/
	
    // data engeneering
    
    sample[0] = int(sample[1]/4);
	sample[1] = fmod(sample[1], 4);
    
    // return the argmax y
    int argmax = 1;
    double max = 0;
    for (int l=0; l<possible_labels.size(); l++)
    {
        double y = 1;
        for (int f=0; f<sample.size(); f++)
        {
            // implement the product of P(xi|y)
            y *= exp(-1 * (sample[f] - mean_y[l][f]) * (sample[f] - mean_y[l][f])
                           / (2 * variance_y[l][f]))
                       / sqrt(2 * M_PI * variance_y[l][f]);
        }
        y *= p_y[l];
        if (y > max)
        {
            max = y;
            argmax = l;
        }
    }
	
	return possible_labels[argmax];

}