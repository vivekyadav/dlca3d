/* 
 * File:   main.cpp
 * Author: Vivek Yadav
 *
 * Created on October 31, 2010, 6:18 PM
 */
//Back-uped on 16/2/2011
//using namespace std;

#include<cstdlib>
#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include<math.h>
#include "space.h"
#include "box.h"





int main(int argc, char** argv) 
{
	int repeat=0,old_n=0,i,l,avg_num=1,run=0,movement=1;
	char* eptr;
	float con_poly=0,con_par=0,sticking_p = 1;
	double *rdf1,*rdf2,*rdf3,*rdf4;
	std::map<long,long> box_count_particles;
    
	l=200;
	l=strtod(argv[1],&eptr);
	std::cout<<"\nl = "<<l;
	con_par=strtod(argv[2],&eptr);
    
	if(argc>3)
	{
		con_poly=strtod(argv[3],&eptr);
		if(argc>4)
			{
                avg_num = strtol(argv[4],&eptr,10);
                if(argc>5)
                    sticking_p = strtod(argv[5],&eptr);
            }
	}
    
	std::cout<<"\nparticle concentration = "<<con_par;
	std::cout<<"\npolymer concentration = "<<con_poly;
	std::cout<<"\nNumber of runs = "<<avg_num;
    std::cout<<"\nSticking Probability = "<<sticking_p;
    
	space space1(con_poly,con_par,l,sticking_p);
    
	space1.r_max = l/3;
	space1.dr = 1;
    
	rdf1 = new double [int(space1.r_max/space1.dr)];
	rdf2 = new double [int(space1.r_max/space1.dr)];
	rdf3 = new double [int(space1.r_max/space1.dr)];
	rdf4 = new double [int(space1.r_max/space1.dr)];
    
	while(run++ < avg_num)
	{
		std::cout<<"\n\n\nSimulation "<<run;
		repeat = 0;movement=1;
	//std::cout<<"\nMin clusters = "<<int(space1.num_nano/1000);
	space1.initialize();
	space1.place_mols();
	printf("\nNumber of molecules is %d",space1.num_nano+space1.num_plmr);
	printf("\n\nStarting Simulation with %d clusters\n",space1.n_clusters);
	std::ofstream moment_ratios("moment_ratios.txt");
	
	while((old_n=space1.n_clusters)>0)
	{
		movement = space1.simulate();
		moment_ratios<<(space1.find_moment(2)/space1.find_moment(1))<<"\n";
		if(old_n==space1.n_clusters)
			repeat++;
		else
			repeat=0;
		//std::cout<<"\n"<<space1.n_clusters<<" clusters left....";
		//std::cout<<"\t movement = "<<movement; 		

		if(repeat >20 && movement == -1 && space1.n_clusters <= (space1.num_nano+space1.num_plmr)/20+20)
			break;
	}
	moment_ratios.close();
    space1.store();
    printf("\nSimulation %d Done",run);

	//fractal Dimension
	/*box b;
	for(i=5;i<space1.l/2;i+=5)
	{
		b.setbox_length(i);
		b.find_particles_contained(&space1);
		std::cout<<"\n particles in a box of length "<<i<<" "<<b.getparticles_contained();//<<" FDim= "<<log(float(b.getparticles_contained()))/log(float(i));
		if(box_count_particles.find(i) == box_count_particles.end())
		{
			box_count_particles[i] = 0;
		}
		box_count_particles[i] = box_count_particles[i] + b.getparticles_contained();
	}*/

	//radius of gyration
	space1.find_all_rG();
	space1.rdf(0,space1.num_nano,0,space1.num_nano,space1.conc_nano,rdf1);
	space1.rdf(0,space1.num_nano,space1.num_nano,space1.num_nano+space1.num_plmr,space1.conc_plmr,rdf2);
	space1.rdf(space1.num_nano,space1.num_nano+space1.num_plmr,0,space1.num_nano,space1.conc_nano,rdf3);
	space1.rdf(space1.num_nano,space1.num_nano+space1.num_plmr,space1.num_nano,space1.num_nano+space1.num_plmr,space1.conc_plmr,rdf4);
	


	//space1.store();
    space1.clear();
	}
	space1.store_rdf("rdf_n_n",rdf1);
	space1.store_rdf("rdf_n_p",rdf2);
	space1.store_rdf("rdf_p_n",rdf3);
	space1.store_rdf("rdf_p_p",rdf4);
	delete[] rdf1;
	delete[] rdf2;
	delete[] rdf3;
	delete[] rdf4;
	std::cout<<"\nAll simulations done\n";
	std::cout.flush();
	/*std::stringstream ss;
	std::string astring;
	ss<<"simulation."<<run<<".LvsN.txt";
	astring = ss.str();
	
	std::map<long,long>::iterator it = box_count_particles.begin();
	std::ofstream data("LvsN.txt");
	while(it != box_count_particles.end())
	{
		data<<it->first<<" "<<float(it->second)/run<<"\n";
		it++;
	}
	data.close();*/
	return 0;
}

