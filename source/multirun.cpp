#include<iostream>
#include<cstdlib>
#include<string>
#include<sstream>

int main(int argc, char** argv)
{
	float c1,c2,inc,sticking_p;
	int l,i,runs;
	std::string str;

	inc=0.0001;
	c1 = 0.0005;
	l = 100;
	runs = 50;
    sticking_p = 0.1;
	for(i=0;i<10;i++)
	{
		c2=inc*i;

		std::stringstream ss1,ss2,ss3,ss4;
		str.clear();

		ss3<<"./simulation "<<l<<" "<<c1<<" "<<c2<<" "<<runs<<" "<<sticking_p<<"\n";
		str.append(ss3.str());
		system(str.c_str());
		str.clear();

		ss1<<"mkdir result_"<<l<<"_"<<c1<<"_"<<c2<<"_"<<runs<<"\n";
		str.append(ss1.str());
		//std::cout<<str;
		system(str.c_str());
		str.clear();
		
		ss2<<"cp corrected_FD.txt nano.txt moment_ratios.txt plmr.txt LvsN.txt RvsN.txt rdf_n_n rdf_n_p rdf_p_n rdf_p_p result_"<<l<<"_"<<c1<<"_"<<c2<<"_"<<runs<<"\n";
		str.append(ss2.str());
		//std::cout<<str;
		system(str.c_str());
		str.clear();

		ss4<<"rm corrected_FD.txt nano.txt moment_ratios.txt plmr.txt LvsN.txt RvsN.txt rdf_n_n rdf_n_p rdf_p_n rdf_p_p \n";
		str.append(ss4.str());
		//std::cout<<str;
		system(str.c_str());
		str.clear();

		//std::cout<<cstr2;
		//system(cstr2);
	}
	return 0;
}
