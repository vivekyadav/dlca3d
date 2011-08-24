#ifndef _SPACE_H_
#define _SPACE_H_


#include <string>
#include<map>
#include<vector>

class point
{
public:
    int x[3];           //the co-ordiantes of the point  DEPENDS on D = DIMENSION

};

class mol: public point
{
public:
    int con[6];         // DEPENDS on D = DIMENSION, con(0,1,2,3,4,5)=>(-y,-x,-z,y,x,z)
    int cID;
    std::string mol_type;

    mol();
};

class cluster
{
	public:
	std::vector<int> members;
	float cM[3];			//centre of mass
	float rG;				//Radius of Gyration
	int mass;

	cluster();
	~cluster();
};
class space
{

    int n_cluster_add,dx[2],sim_num;
    int * clusters_to_add;    
    
	std::map<int,cluster*> clusters;
    int D;
public:
	mol *molecules;
	int num_plmr,num_nano;
	float conc_plmr,conc_nano,sticking_p;
	int ***sites;
    int n_clusters,l;//count starts from 0;
    double r_max,dr;
    
    void initialize();
    space();    
    space(float ,float ,int,float );   
	void set_sim_num(int);

    void place_nano(int start,int num,std::string type);
    void place_plmr(int start,int num,std::string type);
    void place_mols();
    int period(int);
    int period_direction(int);
    void make_cluster_of_nano();
    void particleIdentity(int,int);
    void clear_clusters();
    void clear_molecules();
    void clear_sites();
    int simulate();
    //void join_mols();
    void add_two_clusters(int,int);
    void add_cluster(int);
    void add_all_clusters(int);
    int check(int,int,int);
    void store();
	void store_rdf(char * ,const double *);
	int check_for_duplicates();
	void clean_unused_clusters();
	void find_cM(int);
	float find_rG(int);
	void find_all_rG();
	void connections(int,int,int);
	void rdf(int, int, int, int,double,double *);
	void find_interior_particles(std::vector<int> &,double,const int&,const int&);
	double distance(int,int);
	double find_moment(int n);	//finds the 'n'th moment 
	int find_mass(int);	
	int can_cluster_move(int cID,int direction,int);
	void move_cluster(int cID,int direction,int);
	void clear();
};
#endif
