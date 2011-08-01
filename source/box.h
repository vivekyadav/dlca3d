#ifndef _BOX_H_
#define _BOX_H_
#endif

#ifndef _SPACE_H_
#include "space.h"
#endif
class box
{
      int x,y,z,box_length;
	  
	  
	  public:
             float particles_contained;

      box();

      void setx(int);
	  void sety(int);
	  void setbox_length(int);
	  int getx();int gety();
	  int getbox_length();
	  float getparticles_contained();
	  friend void find_particles_contained(box*);
	  int count_particles(space* sp);
	  void find_particles_contained(space* sp);
	  
	  
};