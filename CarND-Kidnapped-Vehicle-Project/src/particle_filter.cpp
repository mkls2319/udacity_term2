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

#include "particle_filter.h"

using namespace std;



void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all w to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;

	// sensor noise x, y, and theta
	normal_distribution<double> N_x(0, std[0]);
	normal_distribution<double> N_y(0, std[1]);
	normal_distribution<double> N_theta(0, std[2]);

	for (int i = 0; i < num_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = x + N_x(gen);
		p.y = y + N_y(gen);
		p.theta = theta + N_theta(gen);
		p.weight = 1.0;

		particles.push_back(p);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	normal_distribution<double> N_x(0, std_pos[0]);
	normal_distribution<double> N_y(0, std_pos[1]);
	normal_distribution<double> N_theta(0, std_pos[2]);

	for (int i = 0; i < num_particles; i++)
	{
		if (fabs(yaw_rate) < 0.000001)
		{
      particles[i].x += velocity * delta_t * cos(particles[i].theta) + N_x(gen);
      particles[i].y += velocity * delta_t * sin(particles[i].theta) + N_y(gen);
    }
    else
		{
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta)) + N_x(gen);
      particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t)) + N_y(gen);
      particles[i].theta += yaw_rate * delta_t;
    }
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < observations.size(); i++)
	{

    double m_dist = numeric_limits<double>::max();

    int map_id = -1;

    for (unsigned int j = 0; j < predicted.size(); j++)
		{

      double c_dist;
			c_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      if (c_dist < m_dist)
			{
        m_dist = c_dist;
        map_id = predicted[j].id;
      }
    }

    observations[i].id = map_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the w of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	for (int i = 0; i < num_particles; i++) {

     vector<LandmarkObs> pred;

     for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

       if (fabs(map_landmarks.landmark_list[j].x_f -  particles[i].x) <= sensor_range && fabs(map_landmarks.landmark_list[j].y_f -  particles[i].y) <= sensor_range)
			 {
         pred.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f });
       }
     }

     vector<LandmarkObs> t_os;
     for (unsigned int j = 0; j < observations.size(); j++)
		 {
       double t_x = cos( particles[i].theta)*observations[j].x - sin( particles[i].theta)*observations[j].y +  particles[i].x;
       double t_y = sin( particles[i].theta)*observations[j].x + cos( particles[i].theta)*observations[j].y +  particles[i].y;
       t_os.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
     }


     dataAssociation(pred, t_os);

     particles[i].weight = 1.0;

     for (unsigned int j = 0; j < t_os.size(); j++)
		 {

       double ob_x, ob_y, pred_x, pred_y;
       ob_x = t_os[j].x;
       ob_y = t_os[j].y;

       int a_pred = t_os[j].id;

       for (unsigned int k = 0; k < pred.size(); k++)
			 {
         if (pred[k].id == a_pred)
				 {
           pred_x = pred[k].x;
           pred_y = pred[k].y;
         }
       }


       double std_x = std_landmark[0];
       double std_y = std_landmark[1];
       double ob_w = ( 1/(2*M_PI*std_x*std_y)) * exp( -( pow(pred_x-ob_x,2)/(2*pow(std_x, 2)) + (pow(pred_y-ob_y,2)/(2*pow(std_y, 2))) ) );

       particles[i].weight *= ob_w;
     }
   }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> np;

  vector<double> w;
  for (int i = 0; i < num_particles; i++)
	{
    w.push_back(particles[i].weight);
  }


  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto ind = uniintdist(gen);

  double m_weight = *max_element(w.begin(), w.end());

  uniform_real_distribution<double> unirealdist(0.0, m_weight);

  double beta = 0.0;

  for (int i = 0; i < num_particles; i++)
	{
    beta += unirealdist(gen) * 2.0;
    while (beta > w[ind]) {
      beta -= w[ind];
      ind = (ind + 1) % num_particles;
    }
    np.push_back(particles[ind]);
  }

  particles = np;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
