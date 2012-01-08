#include<iostream>
#include<cstdlib>
#include<string>
#include<sstream>

int main(int argc, char** argv)
{
	float c1,c2,inc,sticking_p;
	int l,i,runs;
	std::string str;
    float concentrations[] = {0.001,0.003,0.0055,0.01,0.03,0.072,0.1,0.139,0.2,0.24,0.3,0.4,0.469,0.5};

	inc= 0.002;
	c1 = 0.001;
    c2 = 0;
	l = 50;
	runs = 3;
    sticking_p = 1;
	for(i=0;i<14;i++)
	{
        c1 = concentrations[i];
		std::stringstream ss1,ss2,ss3,ss4;
		str.clear();

		ss3<<"./simulation "<<l<<" "<<c1<<" "<<c2<<" "<<runs<<" "<<sticking_p<<"\n";
		str.append(ss3.str());
		system(str.c_str());
		str.clear();

		ss1<<"mkdir result_"<<l<<"_"<<c1<<"_"<<c2<<"_"<<runs<<"\n";
		str.append(ss1.str());
		system(str.c_str());
		str.clear();
		
		ss2<<"cp -r System_States moment_ratios.txt RvsN.txt rdf_n_n rdf_n_p rdf_p_n rdf_p_p result_"<<l<<"_"<<c1<<"_"<<c2<<"_"<<runs<<"\n";
		str.append(ss2.str());
		system(str.c_str());
		str.clear();

		ss4<<"rm -r System_States moment_ratios.txt RvsN.txt rdf_n_n rdf_n_p rdf_p_n rdf_p_p \n";
		str.append(ss4.str());
		system(str.c_str());
		str.clear();
        
	}
	return 0;
}
