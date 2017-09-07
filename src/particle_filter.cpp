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

    num_particles = 100;

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
        double delta_x, delta_y;
        double gauss_norm, exponent, product;

        
	 
	
	
    
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
	double x, y, theta;
	Map::single_landmark_s l_mark;							      
	std::vector<Map::single_landmark_s> lm_list = map_landmarks.landmark_list;

	double gauss_norm = (1/(2 * 3.14 * std_landmark[0] * std_landmark[1]));
	double x_div = 2*pow(std_landmark[0],2);
	double y_div = 2*pow(std_landmark[1],2);
	double mu_x, mu_y;
	double joint_prob, exponent;
	double dist, delta_x, delta_y;

	for (int i = 0; i < num_particles; i++) {
        	
            std::vector<int> associations;
            std::vector<double> sense_x;
            std::vector<double> sense_y;

	    x = particles[i].x;
	    y = particles[i].y;
	    theta = particles[i].theta;

	    std::vector<LandmarkObs> map_obs;
            LandmarkObs new_obs, old_obs;

            for(int j = 0; j < observations.size(); j++) {
                old_obs = observations[j];
                new_obs.x = particles[i].x + (old_obs.x*cos(theta)-old_obs.y*sin(theta));
                new_obs.y = particles[i].y + (old_obs.x*sin(theta)+old_obs.y*cos(theta));
                map_obs.push_back(new_obs); 
            }

            particles[i].weight = 1.0;

	    for(int k = 0; k < map_obs.size(); k++) {
                double closest = sensor_range;
                int association = 0;
 
	        for (int l = 0; l < lm_list.size(); l++) {
		    l_mark = lm_list[l];
      
  	            dist = sqrt(pow(l_mark.x_f-map_obs[k].x,2) + pow(l_mark.y_f-map_obs[k].y,2));
		    if (dist < closest) {
			closest = dist;
			association = l;
		    }		    
                }
		    mu_x = lm_list[association].x_f;
                    mu_y = lm_list[association].y_f;

                    delta_x = abs(mu_x - map_obs[k].x);
		    delta_y = abs(mu_y - map_obs[k].y);

	            exponent= ((pow(delta_x,2)/x_div)
	                     + (pow(delta_y,2)/y_div));
		    joint_prob = gauss_norm * pow(2.7, -exponent);
		 
		    if (joint_prob > 0 ) 
		    { 
		        particles[i].weight *= joint_prob;
		    }
                
                associations.push_back(association+1);
                sense_x.push_back(map_obs[k].x);
                sense_y.push_back(map_obs[k].y);
            }
	    particles[i] = SetAssociations(particles[i], 
		   		           associations,
				           sense_x, sense_y);
	    weights[i] = particles[i].weight;
	}
	    	
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
        default_random_engine gen;
	discrete_distribution<int> d(weights.begin(), weights.end()); 
	
        int particle_index;
    	vector<Particle> resampled_particles;    

	for(int n=0; n<num_particles; ++n) {
	    particle_index = d(gen); 
	    resampled_particles.push_back(particles[particle_index]);
        }
	particles = resampled_particles;
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
