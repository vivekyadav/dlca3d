#include "box.h"
#include<stdlib.h>

box::box()
{
	x=0;y=0;z=0;
	particles_contained=0;
	box_length=0;
}

void box::setx(int a)
{
	x=a;
}
void box::sety(int a)
{
	y=a;
}
void box::setbox_length(int l)
{
     box_length=l;
}
float box::getparticles_contained()
{
      return particles_contained;
}

int box::count_particles(space* sp)
{
	int i,j,sum=0,p,k;
	for( i=x;i<x+box_length;i++)
	{
         for( j=y;j<y+box_length;j++)
         {
			 for(k=z;k<z+box_length;k++)
			 {
				 p=sp->sites[sp->period(i)][sp->period(j)][sp->period(k)];
				 if(p!=-1 && sp->molecules[p].mol_type!="polymer")
					 sum++;
			 }
         }
    }
    return sum;
}

void box::find_particles_contained(space* sp)
{
	int times=(sp->l^2)/(box_length^2),i,particle_sum=0,m=0;
	for(i=0;i<sp->num_nano;i++)
	{
		m=int(   (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*sp->num_nano);
		x=sp->period(sp->molecules[m].x[0]-box_length/2);		
		y=sp->period(sp->molecules[m].x[1]-box_length/2);
		z=sp->period(sp->molecules[m].x[2]-box_length/2);
		
		particle_sum+=count_particles(sp);
	}
	particles_contained=(particle_sum/sp->num_nano);
}
