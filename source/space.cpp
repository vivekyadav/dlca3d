#include "space.h"
#include<string>
#include<fstream>
#include<ctime>
#include<cstdlib>
#include<math.h>
#include<iostream>
#include<sstream>
mol::mol()
{
	cID=-1;
	int i;
	for(i=0;i<6;i++)
	{
		con[i]=-9;
	}
}

cluster::cluster()
{
	mass=0;
}
cluster::~cluster()
{
	members.clear();
}
int space::find_mass(int c_id)
{
	int i;
	clusters[c_id]->mass=0;
	for(i=0;i<clusters[c_id]->members.size();i++)
	{
		if(molecules[clusters[c_id]->members[i]].mol_type == "nano")
		clusters[c_id]->mass++;
	}
	return clusters[c_id]->mass;
}
double space::find_moment(int n)	//finds the 'n'th moment 
{
	std::map<int,cluster*>::iterator it;
	double nth_moment=0;

	it = clusters.begin();
	while(it!=clusters.end())
	{
		nth_moment += pow(find_mass(it->first),n);
		it++;
	}
	return nth_moment;
}

void space::find_cM(int c_id)
{
	unsigned int i,j;
	for(j=0;j<D;j++)
	{
		clusters[c_id]->cM[j]=0;			
	}
	clusters[c_id]->mass=0;

	for(i=0;i<clusters[c_id]->members.size();i++)
	{
		if(molecules[clusters[c_id]->members[i]].mol_type == "nano")
		{
			for(j=0;j<D;j++)
			{
				clusters[c_id]->cM[j]+=molecules[clusters[c_id]->members[i]].x[j];			
			}
			clusters[c_id]->mass++;
		}
	}
	for(j=0;j<D;j++)
	{
		clusters[c_id]->cM[j]/=clusters[c_id]->mass;
	}
}

float space::find_rG(int c_id)
{
	find_cM(c_id);
	int i,j;
	clusters[c_id]->rG=0;
	for(i=0;i<clusters[c_id]->members.size();i++)
	{
		if(molecules[clusters[c_id]->members[i]].mol_type == "nano")
		for(j=0;j<D;j++)
                    clusters[c_id]->rG+=pow(float(molecules[clusters[c_id]->members[i]].x[j]-clusters[c_id]->cM[j]),2);
	}

	clusters[c_id]->rG/=clusters[c_id]->mass;
	clusters[c_id]->rG=pow(clusters[c_id]->rG,float(0.5));
	return clusters[c_id]->rG;
}

void space::find_all_rG()
{
	std::map<int,cluster*>::iterator it;
	std::fstream data("RvsN.txt",std::fstream::out | std::fstream::app);
	it=clusters.begin();
	float rg=0;
	while(it!=clusters.end())
	{
            rg=find_rG(it->first);
		/*std::stringstream ss;
		ss<<"./findD "<<(it->second->mass)<<" "<<rg<<"\n";
		std::string str(ss.str());
		const char* cstr = str.c_str();
		//std::cout<<"\n"<<cstr;
		std::system(cstr);*/

                data<<rg<<" "<<it->second->mass<<"\n";
		it++;
	}

	data.close();
}

space::space()
{
	conc_plmr=0.01;conc_nano=0.05;l=50;dx[0]=-1;dx[1]=1;
	initialize();
}

space::space(float concen_plmr,float concen_nano,int len)
{
	conc_plmr=concen_plmr;conc_nano=concen_nano;l=len;dx[0]=-1;dx[1]=1;
	initialize();
}

void space::store()
{
    int i,j,n;
    std::ofstream store_file1,store_file2;
    //char* filename="space.txt";

    n=num_plmr+num_nano;
    store_file1.open("nano.txt");
	store_file2.open("plmr.txt");
	
    for(i=0;i<n;i++)
    {
		if(molecules[i].mol_type=="nano")
		{
			for(j=0;j<D;j++)
				store_file1<<molecules[i].x[j]<<" ";
			store_file1<<"\n";            
		}
		else
		{
			for(j=0;j<D;j++)
				store_file2<<molecules[i].x[j]<<" ";
			store_file2<<"\n";
		}
    }

    store_file1.close();
    store_file2.close();
}

void space::store_rdf(char * filename,const double * rdf_arr)
{
	int n = r_max/dr,i;
	std::ofstream store_file;
	store_file.open(filename);
	
	for(i =0; i<n; i++)
		store_file<<dr*(i+1)<<" "<<rdf_arr[i]<<"\n";
	
	store_file.close();
}	
void space::initialize()
{
        D=3;
        
        num_plmr=int(pow(double(l),D)*conc_plmr);
        num_nano=int(pow(double(l),D)*conc_nano);
		sites=new int**[l];
        int j,i,total;
        
		total=num_plmr+num_nano;
        for(i=0;i<l;i++)
        {
            sites[i]=new int*[l];
			for(j=0;j<l;j++)
				sites[i][j]=new int[l];
        }
        molecules = new mol[total];
        for( i=0;i<total;i++)
		{
			clusters[i]=new cluster;
			//printf("\n%d, %d",i,clusters.find(i)->second->n_members);
		}
}
void space::clear_clusters()
{
	std::map<int,cluster*>::iterator it;
	for(it=clusters.begin(); it!=clusters.end(); it++)
			delete it->second;

	clusters.clear();
}
void space::clean_unused_clusters()
{	
	std::map<int,cluster*>::iterator it;
	it=clusters.find(n_clusters+1);
	clusters.erase(it,clusters.end());
}
void space::clear_molecules()
{
    delete[] molecules;
}
void space::clear_sites()
{
	int i,j;
	for(i=0;i<l;i++)
	{
		for(j=0;j<l;j++)
			delete this->sites[i][j];
			
		delete this->sites[i];
	}
	delete this->sites;
}
void space::clear()
{
	this->clear_clusters();
	this->clear_molecules();
	this->clear_sites();	
}
int space::period(int n)
{
    if( n>=l)
    n-=l;
    else if( n<0)
    n+=l;

    return n;
}
int space::period_direction(int n)
{
	int m=2*D;
    if( n>=m)
    n-=m;
    else if( n<0)
    n+=m;

    return n;
}

void space::place_mols()
{
    int i,j,k;
    std::string type;

    for(i=0;i<l;i++)
        for(j=0;j<l;j++)
			for(k=0;k<l;k++)
				sites[i][j][k]=-1;

    type="nano";
    place_nano(0,num_nano,type);
    make_cluster_of_nano();
    type="polymer";
    place_plmr(num_nano,num_plmr,type);

	//clearing extra cluster space
	clean_unused_clusters();
}
void space::place_nano(int start,int num,std::string type)
{
    int i,j,x_temp[3];
    srand(time(0));

    for(i=start;i<start+num;i++)
    {
        for(j=0;j<D;j++)
            x_temp[j]=(rand()%l);
        if(sites[x_temp[0]][x_temp[1]][x_temp[2]]!=-1)
        {
            i--;continue;
        }
        for(j=0;j<D;j++)
            molecules[i].x[j]=x_temp[j];

        molecules[i].mol_type=type;
        sites[x_temp[0]][x_temp[1]][x_temp[2]]=i;
    }
}
void space::connections(int i,int a,int p)			// con(0,1,2,3,4,5)=>(-y,-x,-z,y,x,z)
{
	int j,k;

	molecules[i].con[a]=p;
	molecules[i].con[period_direction(a+D)]=0;
	k=-1;
	
	for(j=1;j<D;j++)
	{
		molecules[i].con[period_direction(a+j)]=k;
		k--;
		molecules[i].con[period_direction(a+j+D)]=k;
		k--;
	}
	
}
void space::place_plmr(int start,int num,std::string type)
{
    int i,j,x_temp[3],px,py,pz,n,k,dx[2]={-1,1},is_cID_assigned;
    srand(time(0));
	std::map<int,cluster*>::iterator it;

    for(i=start;i<start+num;i++)
    {
        for(j=0;j<D;j++)
            x_temp[j]=(rand()%l);
		px=x_temp[0];py=x_temp[1];pz=x_temp[2];

        if(sites[px][py][pz]!=-1)
        {
            i--;continue;
        }
		
        for(j=0;j<D;j++)
            molecules[i].x[j]=x_temp[j];

        molecules[i].mol_type=type;
        sites[x_temp[0]][x_temp[1]][x_temp[2]]=i;

        //checking if any nano particle adjacent to plmr, joining it to that nano
        is_cID_assigned=0;
                   for( k=0;k<2;k++)
                   {
                        n=period(py+dx[k]);
                        if( sites[px][n][pz] != -1 && molecules[ sites[px][n][pz] ].mol_type == "nano")
                        {
                            if(k==0)
                            {
                                connections(i,0,sites[px][n][pz]);
                            }
                            else
                            {
                                connections(i,0+D,sites[px][n][pz]);
                            }
                            molecules[i].cID=molecules[ sites[px][n][pz] ].cID;
							//printf("\nPushing Mol. %d",i);
                            clusters[molecules[i].cID]->members.push_back(i);
							//it=clusters.find(molecules[i].cID);
							//it->second->members.push_back(i);
							//it->second->n_members++;
                            is_cID_assigned=1;
                            goto END;
                        }

                   }

                   for( k=0;k<2;k++)
                   {
                        n=period(px+dx[k]);
                        if( sites[n][py][pz] != -1 && molecules[ sites[n][py][pz] ].mol_type != "polymer")
                        {
                            if(k==0)
                            {
                                connections(i,1,sites[n][py][pz]);
                            }
                            else
                            {
                               connections(i,1+D,sites[n][py][pz]);
                            }
                            molecules[i].cID=molecules[ sites[n][py][pz] ].cID;
							//printf("\nPushing Mol. %d",i);
                            clusters[molecules[i].cID]->members.push_back(i);
							//it=clusters.find(molecules[i].cID);
							//it->second->members.push_back(i);
							//it->second->n_members++;
                            is_cID_assigned=1;
                            goto END;
                        }

                   }

				   for( k=0;k<2;k++)
                   {
                        n=period(pz+dx[k]);
                        if( sites[px][py][n] != -1 && molecules[ sites[px][py][n] ].mol_type != "polymer")
                        {
                            if(k==0)
                            {
                                connections(i,2,sites[px][py][n]);
                            }
                            else
                            {
                               connections(i,2+D,sites[px][py][n]);
                            }
                            molecules[i].cID=molecules[ sites[px][py][n] ].cID;
							//printf("\nPushing Mol. %d",i);
                            clusters[molecules[i].cID]->members.push_back(i);
							//it=clusters.find(molecules[i].cID);
							//it->second->members.push_back(i);
							//it->second->n_members++;
                            is_cID_assigned=1;
                            goto END;
                        }

                   }

        END:;
        if(is_cID_assigned==0)
        {
			space::n_clusters++;
            molecules[i].cID=space::n_clusters;
            clusters[space::n_clusters]->members.push_back(i);
            //space::n_clusters++;
        }
    }

}
void space::particleIdentity(int p,int c_id)
{
     molecules[p].cID=c_id;
     clusters[c_id]->members.push_back(p);

	 int i=0,px=molecules[p].x[0],py=molecules[p].x[1],pz=molecules[p].x[2],n,dx[2]={-1,1};
     for( i=0;i<2;i++)
     {
          n=period(py+dx[i]);
          if( sites[px][n][pz] != -1 && molecules[ sites[px][n][pz] ].cID != c_id)
              particleIdentity(sites[px][n][pz],c_id);

     }

     for( i=0;i<2;i++)
     {
          n=period(px+dx[i]);
          if( sites[n][py][pz] != -1 && molecules[ sites[n][py][pz] ].cID != c_id)
              particleIdentity(sites[n][py][pz],c_id);

     }
	 for( i=0;i<2;i++)
     {
          n=period(pz+dx[i]);
          if( sites[px][py][n] != -1 && molecules[ sites[px][py][n] ].cID != c_id)
              particleIdentity(sites[px][py][n],c_id);

     }
}
void space::make_cluster_of_nano()
{
    int i,cId=0;
    for( i=0;i<num_nano;i++)
     {
          if( molecules[i].cID == -1)
          {

              particleIdentity(i,cId);
              cId++;
          }
     }
    n_clusters=cId-1;
}

int space::can_cluster_move(int cID,int direction,int v)
{
	int i,x[3],j,p=0;
	for(i=0; i<clusters[cID]->members.size(); i++)
	{
		p = clusters[cID]->members[i];
		for(j=0;j<D;j++)
			x[j] = molecules[p].x[j];
		
		x[direction] = period(x[direction]+dx[v]);
		if(sites[x[0]][x[1]][x[2]]!=-1 && molecules[sites[x[0]][x[1]][x[2]]].cID != cID)
			{//std::cout<<"\nCluster "<<cID<<" cant move";
			return 0;}
	}
	
	return 1;
}

void space::move_cluster(int cID,int direction, int v)
{
	int i,p=0;
	for(i=0; i<clusters[cID]->members.size(); i++)
	{
		p = clusters[cID]->members[i];
		sites[molecules[p].x[0]][molecules[p].x[1]][molecules[p].x[2]]=-1;		
	}
	for(i=0; i<clusters[cID]->members.size(); i++)
	{
		p = clusters[cID]->members[i];
		molecules[p].x[direction] = period(molecules[p].x[direction]+dx[v]);
		sites[molecules[p].x[0]][molecules[p].x[1]][molecules[p].x[2]]=p;
	}
	/*std::cout<<"\nMoved CLuster "<<cID<<" in direction "<<direction<<", "<<v;
	std::cout.flush();*/
}

int space::simulate()
{
    std::map<int,cluster*>::iterator it;
    int d,direction,v,i,cID,p,px,py,pz,k,n,q,if_stop=0,n_molecules,skip,t,move_failed=0;
	n_molecules=num_nano+num_plmr;

    for( q=0 ;q<n_molecules ;q++)
    {
		if(move_failed >= int(n_molecules)-1)
			return -1;

        n_cluster_add=0;skip=0;

              //d=int(   (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )*n_clusters); // c_id is the cluster to be shifted
		      d=rand()%(n_clusters+1);	//generate numbers from 0 to n_clusters
              direction=rand()%(D);
              v=rand()%2;   //if shifting is to be done in positive or negative direction

	    it=clusters.begin();
	    for(i=0;i<d;i++)
		{
		  it++;
		}
	    cID=(*it).first;	//the cluster that we wish to move
			/*if(check_for_duplicates())
					{printf("Overlap Found before cluster movement");return 0;} */
			if(can_cluster_move(cID,direction,v) == 0)
				{move_failed++;continue;}
			else
				move_cluster(cID,direction,v);
			/*if(check_for_duplicates())
				{printf("Overlap Found after cluster movement");return 0;} */
			  //Combining clusters as neccesary ......
			  // clusters_to_add has the list of clusters that have come to contact and need to be merged
			  if(3*clusters[cID]->members.size() < n_clusters)
			  clusters_to_add = new int[3*clusters[cID]->members.size()];
			  else
				   clusters_to_add = new int[n_clusters];

              for( i=0;i<clusters[cID]->members.size();i++)
              {
                   p=clusters[cID]->members[i];
				   px=molecules[p].x[0];py=molecules[p].x[1];pz=molecules[p].x[2];
                   for( k=0;k<2;k++)
                   {
                        n=period(py+dx[k]);
                        if( sites[px][n][pz] != -1 && molecules[ sites[px][n][pz] ].cID != cID)
                        {
                            //add_cluster(molecules[ sites[px][n] ].cID);
                            if_stop=check(p,sites[px][n][pz],D*k);
							if(if_stop==0)
								goto after_add_all_clusters;
							else if(if_stop==-1)
							{
								move_cluster(cID,direction, abs(v-1));
								move_failed++;
								goto after_add_all_clusters;
							}
						}

                   }

                   for( k=0;k<2;k++)
                   {
                        n=period(px+dx[k]);
                        if( sites[n][py][pz] != -1 && molecules[ sites[n][py][pz] ].cID != cID)
						{
                            //add_cluster(molecules[ sites[n][py] ].cID);
                            if_stop=check(p,sites[n][py][pz],D*k+1);
							if(if_stop==0)
								goto after_add_all_clusters;
							else if(if_stop==-1)
							{
								move_cluster(cID,direction, abs(v-1));
								move_failed++;
								goto after_add_all_clusters;
							}
						}

                   }

				   for( k=0;k<2;k++)
                   {
                        n=period(pz+dx[k]);
                        if( sites[px][py][n] != -1 && molecules[ sites[px][py][n] ].cID != cID)
						{
                            //add_cluster(molecules[ sites[n][py] ].cID);
                            if_stop=check(p,sites[px][py][n],D*k+2);
							if(if_stop==0)
								goto after_add_all_clusters;
							else if(if_stop==-1)
							{
								move_cluster(cID,direction, abs(v-1));
								move_failed++;
								goto after_add_all_clusters;
							}
						}

                   }
              }
              
              if(n_cluster_add == 0)
				move_failed++;
			  add_all_clusters(cID);			   

			  after_add_all_clusters:;		//LABEL "after_add_all_clusters"

			  delete []clusters_to_add;
			  

			  end_of_for:;					//LABEL "end_of_for"

         }
	return move_failed;
}

int space::check(int mol1, int mol2,int direction)
{
    std::string type1=molecules[mol1].mol_type,type2=molecules[mol2].mol_type;
    std::map<int,cluster*>::iterator it;
    int cluster1=molecules[mol1].cID,cluster2=molecules[mol2].cID;
    if(type1=="nano")
    {
        if(type2=="nano")
        {
            add_cluster(molecules[mol2].cID);
        }
        else
        {
            direction+=D;
            if(direction>2*D-1)
                direction-=2*D;

            if(molecules[mol2].con[direction]==-9)
            {
                connections(mol2,direction,mol1);

                molecules[mol2].cID=cluster1;
                clusters[cluster1]->members.push_back(mol2);
                it = clusters.find(cluster2);
				delete it->second;
				clusters.erase(it);
				n_clusters--;
            }
			else
			{				
				return -1;
			}
        }
    }
    else
    {
        if(type2=="nano")
        {
            if(molecules[mol1].con[direction]==-9)
            {
                connections(mol1,direction,mol2);

                molecules[mol1].cID=cluster2;
                clusters[cluster2]->members.push_back(mol1);
                it = clusters.find(cluster1);
				delete it->second;
				clusters.erase(it);
				n_clusters--;
				return 0;
            }
			else
			{
				return -1;
			}
        }
		else		//when both neighbouring molecules are polymers
		{
			//printf("\nCluster taken back");
			return -1;
		}
    }
	return 1;
}

void space::add_two_clusters(int cluster1, int cluster2)
{

	int i=0,j=0;
     /*if( cluster1>cluster2)
     {
         i=cluster1;
         cluster1=cluster2;
         cluster2=i;
     } */
	 //printf("\nAdding Clusters  %d and %d ",cluster1,cluster2);
     int n_cluster1=clusters[cluster1]->members.size(),n_cluster2=clusters[cluster2]->members.size();

     j=0;
     for( i=n_cluster1,j=0;i<n_cluster1+n_cluster2;i++,j++)
     {
          clusters[cluster1]->members.push_back(clusters[cluster2]->members[j]);
          molecules[clusters[cluster1]->members[i]].cID=cluster1;
     }

}

//Iterate over the list of all the clusters to be added an dadd each one
void space::add_all_clusters(int c_id)
{
	int i=0;
	for( i=0;i<n_cluster_add;i++)
	{
		add_two_clusters(c_id,clusters_to_add[i]);
	}
	std::map<int,cluster*>::iterator it;
	for( i=0;i<n_cluster_add;i++)
	{
		it = clusters.find(clusters_to_add[i]);
		delete it->second;
		clusters.erase(it);
		n_clusters--;
	}
}

//This Function is to make sure that there are no repeatations in  the list of clusters to be added.
void space::add_cluster(int p)
{
	int i,flag = 1;
	for(i=0;i<n_cluster_add;i++)
	{
		if( clusters_to_add[i] == p)
		{
			flag = 0;
			break;
		}
	}

	if( flag)
	{
		clusters_to_add[n_cluster_add++] = p;
	}
}

int space::check_for_duplicates()
{
	int i,j,n_molecules=num_nano+num_plmr;
	for(i=0;i<n_molecules-1;i++)
	{
		for(j=i+1;j<n_molecules;j++)
			if(molecules[i].x[0]==molecules[j].x[0] && molecules[i].x[1]==molecules[j].x[1] && molecules[i].x[2]==molecules[j].x[2])
				{
					printf("\nOverlap of molecules %d and %d",i,j);
					printf("of cluster %d and %d",molecules[i].cID,molecules[j].cID);
					return 1;
				}
	}
	return 0;
}

double space::distance(int i,int j)
{
	return pow(pow(double(this->molecules[i].x[0]-this->molecules[j].x[0]),2)+pow(double(this->molecules[i].x[1]-this->molecules[j].x[1]),2)+pow(double(this->molecules[i].x[2]-this->molecules[j].x[2]),2),0.5);
}

void space::rdf(const int p1_start,const int p1_len,const int p2_start,const int p2_len,const double concentration,double * rdf_arr)
{
	int i,j,num_increments=r_max/dr;
	double d,pi=3.14159265,volume=0;

	std::vector<int> interior_particles;
	find_interior_particles(interior_particles, r_max,p1_start,p1_len);
	//printf("\n Num int particles = %d R_MAX = %f  num_increments = %d", interior_particles.size(),r_max,num_increments);
	
	double *particles_in_dr = new double[int(r_max/dr)];
	for(i=0; i<int(this->r_max/this->dr); i++)
		particles_in_dr[i] = 0;

	for(i=0; i<interior_particles.size(); i++)
	{
		for(j=p2_start; j<p2_start+p2_len; j++)
		{
			d=this->distance(i,j);
			//printf("\n D = %f",d);
			if(d<r_max && d!=0)
				particles_in_dr[int(d/dr)]++;
		}
	//	particles_in_dr[0]--;	//to exclude particle 'i' itself
	}

	for(i=0; i<num_increments; i++)
	{
		volume= 4*pi*(pow(dr*(i+1),3)-pow(dr*i,3))/3;
		//printf("\nparticle_in_dr = %f ,volume = %f",particles_in_dr[i],volume);		
		rdf_arr[i] = particles_in_dr[i]/(concentration*interior_particles.size()*volume);
		//printf(" particle_in_dr = %f",particles_in_dr[i]);
	}
	delete[] particles_in_dr;
}

//find all those particles on which we can fit a sphere of radius 'r_max' without going out of the box
void space::find_interior_particles(std::vector<int> &ip,double r_max,const int & p1_start,const int & p1_len)
{
	int i,x,y,z;
	for(i=p1_start; i<p1_start+p1_len; i++)
	{
		x=this->molecules[i].x[0];y=this->molecules[i].x[1];z=this->molecules[i].x[2];
		if(x>r_max && y>r_max && z>r_max && x<(this->l-r_max) && y<(this->l-r_max) && z<(this->l-r_max))
			{ip.push_back(i);/* printf("\n pushed "); */}
	}
}

