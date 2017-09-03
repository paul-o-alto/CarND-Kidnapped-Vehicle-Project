/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <array>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles 
        // to first position (based on estimates of x, y, theta and 
        // their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.

    num_particles = 1000;

    std::default_random_engine gen;

    std::normal_distribution<double> N_x(x, std[0]);
    std::normal_distribution<double> N_y(y, std[1]);
    std::normal_distribution<double> N_theta(theta, std[2]);

    for (int i=0; i<num_particles; i++) {
        Particle particle;
        particle.id = i;
        particle.x = N_x(gen);
        particle.y = N_y(gen);
        particle.theta = N_theta(gen);
        particle.weight = 1;
        
        particles.push_back(particle);
        weights.push_back(1);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and 
        // std::default_random_engine useful.
	// http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	// http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    for (int i = 0; i < num_particles; i++) {
        double new_x;
        double new_y;
        double new_theta;

        if (yaw_rate == 0) {

            new_x = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            new_y = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            new_theta = particles[i].theta;
            
        } else {

            new_x = particles[i].x + velocity/yaw_rate*(
                sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta));
            new_y = particles[i].y + velocity/yaw_rate*(
                cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
            new_theta = particles[i].theta + yaw_rate*delta_t;

        }

        normal_distribution<double> N_x(new_x, std_pos[0]);
        normal_distribution<double> N_y(new_y, std_pos[1]);
        normal_distribution<double> N_theta(new_theta, std_pos[2]);
        
        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_theta(gen);
    }

   

}

void ParticleFilter::dataAssociation(std::vector<Particle> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed 
	//       measurement and assign the observed measurement to this particular 
	//       landmark.
	// NOTE: this method will NOT be called by the grading code. But you will 
	//       probably find it useful to implement this method and use it as a 
	//       helper during the updateWeights phase.

        double min_delta_x, min_delta_y;
        double dist, min_dist;
        double x_pred, y_pred, theta, cos_temp, sin_temp;
	double obs_map_x, obs_map_y;
	double sense_x, sense_y;
        double x_diff, y_diff;
        double gauss_norm, exponent;


	for (int i = 0; i < predicted.size(); i++) {
	    int association; 
	    min_dist = NAN;
	    x_pred = predicted[i].x;
	    y_pred = predicted[i].y;
	    theta  = predicted[i].theta;
	  
	    for(int j = 0; j < observations.size(); j++) {
		cos_temp = cos(theta);
		sin_temp = sin(theta);
                obs_map_x = observations[i].x + cos_temp*x_pred - sin_temp*y_pred;
		obs_map_y = observations[i].y + sin_temp*x_pred + cos_temp*y_pred;

	        x_diff = abs(x_pred - obs_map_x);
	        y_diff = abs(y_pred - obs_map_y);
	   
	        dist = sqrt((x_diff*x_diff + y_diff*y_diff));
	        if (dist < min_dist) {
	            min_dist = dist;
	            min_delta_x = x_diff;
		    min_delta_y = y_diff;
	        }
	    }

	    // calculate normalization term
	    gauss_norm = (1/(2 * 3.14 * 0.3 * 0.3)); 
	    // calculate exponent
	    exponent= ((pow(min_delta_x,2))/(2*pow(0.3,2)) 
		     + (pow(min_delta_y,2))/(2*pow(0.3,2)));
	    //calculate weight using normalization terms and exponent
            weights[i] = gauss_norm * exp(-exponent); //predicted[i].weight;
	    particles[i].weight = weights[i];
	}

	

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian 
	// distribution. You can read more about this distribution here: 
	//   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. 
	//        Your particles are located according to the MAP'S coordinate 
	//        system. You will need to transform between the two systems.
	//        Keep in mind that this transformation requires both rotation 
	//        AND translation (but no scaling).
	// The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	// and the following is a good resource for the actual equation to 
	// implement (look at equation 3.33):
	//   http://planning.cs.uiuc.edu/node99.html
	
	//Particle predicted;
	double x, y;
	double obs_map_x, obs_map_y;
        double theta;
	Map::single_landmark_s l_mark;							      

	dataAssociation(particles, observations);
	std::vector<Map::single_landmark_s> lm_list = map_landmarks.landmark_list;

	double best_weight = -NAN;
	int best_idx = 0;
        for (int a = 0; a < num_particles; a++) {
            if (particles[a].weight > best_weight) {
                best_weight = particles[a].weight;
		best_idx = a;
	    }
        }
	Particle best_p = particles[best_idx];
	double sin_val = sin(best_p.theta);
	double cos_val = cos(best_p.theta);
	std::vector<int> associations;
	std::vector<double> sense_x;
	std::vector<double> sense_y;

	for (int i = 0; i < observations.size(); i++) {
            x = observations[i].x;
	    y = observations[i].y;

	    for (int j = 0; j < lm_list.size(); j++) {
		l_mark = lm_list[j];
		obs_map_x = x + cos_val*best_p.x - sin_val*best_p.y;
	        obs_map_y = y + sin_val*best_p.x + cos_val*best_p.y;	

	        associations.push_back(l_mark.id_i);
	        sense_x.push_back(obs_map_x - l_mark.x_f);
	        sense_y.push_back(obs_map_y - l_mark.y_f); 
	              
	    }
	}

        for (int k=0; k < num_particles; k++) {
	    SetAssociations(particles[k], associations, 
	                    sense_x, sense_y);
	}
	    
	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
        std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<int> d(weights.begin(), weights.end()); 
	
        int particle_index;
    	Particle sampled_particle;    

	for(int n=0; n<num_particles; ++n) {
	    particle_index = d(gen); 
	    sampled_particle = particles[particle_index];
            particles[n].x = sampled_particle.x;
	    particles[n].y = sampled_particle.y;
	    particles[n].theta = sampled_particle.theta;
	    particles[n].weight = sampled_particle.weight;	
            weights[n] = sampled_particle.weight;
	    
        }
        
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
