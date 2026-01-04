#include "Particle.hpp"

using namespace std;

Particle::Particle(){
}


Particle::Particle(shared_ptr<Particle> particle2use){
    this->reinitializeUsingParticle(particle2use);
}

void Particle::reinitializeUsingParticle(shared_ptr<Particle> particle2use){
    double* pos2use = particle2use->getPosition();
    double* vel2use = particle2use->getVelocity();
    for( int i = 0; i < 6; i++ ){
        pos[i] = pos2use[i];
        vel[i] = vel2use[i];
    }
    type = particle2use->getType();
    charge = particle2use->getCharge();
}

void Particle::reinitializeUsingParticle(Particle* particle2use){
    double* pos2use = particle2use->getPosition();
    double* vel2use = particle2use->getVelocity();
    for( int i = 0; i < 6; i++ ){
        pos[i] = pos2use[i];
        vel[i] = vel2use[i];
    }
    type = particle2use->getType();
    charge = particle2use->getCharge();
}


void Particle::serialize(double * objects, int shift){
    
    for( int i = 0; i < 6; i++ ){
        objects[shift+i]   = pos[i];
        objects[shift+i+6] = vel[i];
    }
    objects[shift+12] = type;
    objects[shift+13] = charge;
}


void Particle::deserialize(double * objects, int shift){
    for( int i = 0; i < 6; i++ ){
        pos[i] = objects[shift+i];
        vel[i] = objects[shift+6+i];
    }

    type = (int) objects[shift+12];
    charge = objects[shift+13];
}

double* Particle::getPosition(){
    return pos;
}

void Particle::setPosition(double input[6]){
    for( int i = 0; i < 6; i++ ){
        pos[i] = input[i];
    }
}

void Particle::setPosition(int idx, double input){
    pos[idx] = input;
}

double* Particle::getVelocity(){
    return vel;
}

void Particle::setVelocity(double input[6]){
    for( int i = 0; i < 6; i++ ){
        vel[i] = input[i];
    }
}

void Particle::setVelocity(int idx, double input){
    vel[idx] = input;
}


int Particle::getType(){
    return type;
}

void Particle::setType(int input){
    type = input;
}

double Particle::getCharge(){
    return charge;
}

void Particle::setCharge(double newCharge){
    charge = newCharge;
}
