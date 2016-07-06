#include<iostream>
#include<vector>
#include<fstream>
#include<vector>
#include<cstdlib>
#include<cmath>
#include<math.h>
#include<iomanip>

using namespace std;


struct point
{
        int go_ref;
        int go_alt;
        int ba_ref;
        int ba_alt;
        float likeAA;
        float likeAB;
        float likeBB;
};

struct marker
{
        string chro;
        int position;
        char ref;
        char alt;
};


float fplus=0.6;
float error=0.01;
void read_list(char *filelist,vector<string> &biglist)
{
        ifstream infile(filelist,std::ios::in);
        if(!infile.is_open())
        {
                cout<<"Can't open file "<<filelist<<"\n";
                exit(1);
        }
        while(!infile.eof())
        {
                string temp;
                getline(infile,temp);
                if(temp.length() < 1)
                        continue;
                biglist.push_back(temp);
        }
        infile.close();
}

void read_marker(string filelist,vector<struct marker> &bigmarker)
{
        ifstream infile(filelist.c_str(),std::ios::in);
        if(!infile.is_open())
        {
                cout<<"Can't open file "<<filelist<<"\n";
                exit(1);
        }
        while(!infile.eof())
        {
                string temp;
                getline(infile,temp);
                if(temp.length() < 1)
                        continue;
                vector<int> seppos;
                for(int i=0;i<temp.length();i++)
                {
                        if(temp[i] == '\t' || temp[i] == ' ')
                        {
                                seppos.push_back(i);
                                if(seppos.size() >= 3)
                                        break;
                        }
                }
                struct marker temp_marker;
                temp_marker.chro=temp.substr(0,seppos[0]);
                temp_marker.position=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0] -1).c_str());
                temp_marker.ref=temp[seppos[1]+1];
		temp_marker.alt=temp[seppos[2]+1];
                bigmarker.push_back(temp_marker);
        }
        infile.close();
}

void construct_hash(char *filepos,vector<vector<struct marker> > &chosen_marker,vector<struct marker> &marker_buck,int &keyslot)
{
	//vector<struct marker> marker_buck;
	ifstream infile(filepos,std::ios::in);
	if(!infile.is_open())
	{
		cout<<"Can't open file "<<filepos<<"\n";
		exit(1);
	}
	while(!infile.eof())
	{
		string temp;
		getline(infile,temp);
		if(temp.length() < 1)
			continue;
		//cout<<temp<<"\n";
		vector<int> seppos;
		for(int i=0;i<temp.length();i++)
		{
			if(temp[i] == '\t' || temp[i] == ' ')
			{
				seppos.push_back(i);
				if(seppos.size() >= 3)
					break;
			}
		}
		if(seppos.size() == 0)
		{
			cerr<<"No position in the chosen marker file \n";
			exit(1);
		}
		struct marker temp_ma;
		temp_ma.chro=temp.substr(0,seppos[0]);
		temp_ma.position=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0] -1).c_str());
		temp_ma.ref=temp[seppos[1]+1];
		temp_ma.alt=temp[seppos[2]+1];
		marker_buck.push_back(temp_ma);
	}
	infile.close();
	keyslot=int(marker_buck.size()/3);
	chosen_marker.resize(keyslot);
	//cout<<"Ok here\n";
	long total_add_count=0;
	float ten_plus=1.0;
	for(int i=0;i<marker_buck.size();i++)
	{
		long total_add_count=0;
		int ten_plus=1;
		for(int k=0;k<marker_buck[i].chro.length();k++)
		{
			total_add_count += int(marker_buck[i].chro[k])*ten_plus*(keyslot+1);
			ten_plus =ten_plus*2;
		}
		total_add_count += marker_buck[i].position*(keyslot+1);
		int key_value=abs(total_add_count)% keyslot;
		chosen_marker[key_value].push_back(marker_buck[i]);
	//	cout<<key_value<<"\t"<<marker_buck[i].chro<<"\t"<<marker_buck[i].position<<"\n";
	}
}


void read_big_data(vector<string> &filelist,vector<vector<struct marker> > &chosen_marker,vector<vector<struct point> > &super_point,int keyslot)
{
	for(int i=0;i<filelist.size();i++)
	{
		vector<struct point> big_point;
		ifstream infile(filelist[i].c_str(),std::ios::in);
		if(!infile.is_open())
		{
			cout<<"Can't open file "<<filelist[i]<<"\n";
			exit(1);
		}
		cout<<"Reading file "<<filelist[i]<<"\n";
		//continue;
		while(!infile.eof())
		{
			string temp;
			getline(infile,temp);
			if(temp.length() < 1)
				continue;
			vector<int> seppos;
			for(int k=0;k<temp.length();k++)
			{
				if(temp[k] == ' ' || temp[k] == '\t')
				{
					seppos.push_back(k);
					if(seppos.size() >= 4)
						break;
				}
			}
			//if(seppos.size() >= 4)
			////	break;
			if(seppos.size() != 4)
			{
				cerr<<"the pileup input format is wrong "<<temp<<"\n";
				exit(1);
			}
			struct point temp_point;
			temp_point.go_ref=0;
			temp_point.go_alt=0;
			temp_point.ba_ref=0;
			temp_point.ba_alt=0;

	
			string chro=temp.substr(0,seppos[0]);
			int position=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0] -1).c_str());
		//	cout<<chro<<"\t"<<position<<"\t"<<temp<<"\t"<<keyslot<<endl;
			//exit(1);
			 long total_add_count=0;
                int ten_plus=1;
                for(int k=0;k<chro.length();k++)
                {
                        total_add_count += int(chro[k])*ten_plus*(keyslot+1);
                        ten_plus =ten_plus*2;
                }
                total_add_count += position*(keyslot+1);
                int key_value=abs(total_add_count)% keyslot;

		//	cout<<key_value<<endl;
		//	exit(1);
			//struct point temp_point;
			for(int m=0;m<chosen_marker[key_value].size();m++)
			{
				if(chosen_marker[key_value][m].chro.compare(chro) == 0 && chosen_marker[key_value][m].position == position)
				{
				//	cout<<chosen_marker[key_value][m].chro<<"\t"<<chro<<"\t"<<chosen_marker[key_value][m].position<<"\t"<<position<<"\n";
					for(int k=seppos[2]+1,l=seppos[3]+1;k<seppos[3];k++,l++)
					{
						if(temp[k] == chosen_marker[key_value][m].ref)
						{
							if(temp[l] - 33 >= 20)
							{
								++temp_point.go_ref;
							}
							else
							{
								++temp_point.ba_ref;
							}
						}
						else if(temp[k] == chosen_marker[key_value][m].alt)
						{
							if(temp[l] - 33 >=20)
							{
								++temp_point.go_alt;
							}
							else
							{
								++temp_point.ba_alt;
							}
						}
						else if(temp[k] == '*' && temp[l] == '*')
						{
							temp_point.go_ref=0;
							temp_point.go_alt=0;
							temp_point.ba_ref=0;
							temp_point.ba_alt=0;
							break;
						}
					}
					big_point.push_back(temp_point);
					break;
				}
			}
		
		}
		super_point.push_back(big_point);
	}
/*	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			cout<<super_point[i][k].go_ref<<"\t"<<super_point[i][k].ba_ref<<"\t"<<super_point[i][k].go_alt<<"\t"<<super_point[i][k].ba_alt<<"\n";
		}
	}*/
}
	

void calculate_likelihood(vector<vector<struct point> > &super_point)
{
	for(int i=0;i<super_point[0].size();i++)
	{
		for(int k=0;k<super_point.size();k++)
		{
			float probAA=0.0,probAB=0.0,probBB=0.0;
			float readA=super_point[k][i].go_ref + fplus*super_point[k][i].ba_ref;
			float readB=super_point[k][i].go_alt + fplus*super_point[k][i].ba_alt;
			if(readA > 10)
			{
				readA=10;
			}
			if(readB > 10)
			{
				readB=20;
			}
			if(readA != 0.0 || readB != 0.0)
			{
				probAA=pow(0.6,fplus*readA)*pow(error,readB);
				probAB=pow((1-fplus),readA+readB);
				probBB=pow(error,readA)*pow(0.6,readB);
			}
			else if(readA == 0.0 && readB == 0.0)
			{
				probAA=0.3333;
				probAB=0.3333;
				probBB=0.3333;
			}
			else
			{
				cerr<<"abormal happens\n";
				exit(1);
			}
			super_point[k][i].likeAA=log10(probAA/(probAA + probAB + probBB));
			super_point[k][i].likeAB=log10(probAB/(probAA+probAB+probBB));
			super_point[k][i].likeBB=log10(probBB/(probAA+probAB+probBB));
		}
	}
	/*for(int i=0;i<super_point[0].size();i++)
	{
		for(int k=0;k<super_point.size();k++)
		{
			cout<<super_point[k][i].likeAA<<"\t"<<super_point[k][i].likeAB<<"\t"<<super_point[k][i].likeBB<<"\n";
		}
	}*/
}


void output_file(char *fileout,vector<struct marker> &bigmarker,vector<vector<struct point> > &super_point,vector<string> &biglist)
{
	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cout<<"Can't open file "<<fileout<<"\n";
		exit(1);
	}
	outfile<<"##fileforamt=VCFv4.0\n##source=Hetmap:calculate\n##Format<ID=GT,number,Type=String,Description=\"genotype\">\n##Format=<ID=GL,number=3,Type=float,Description=\"three log10-scaled likelihoods for RR RA AA genotypes\">\n";
	outfile<<"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	for(int i=0;i<biglist.size();i++)
	{
		int start=biglist[i].find_first_of(".");
		outfile<<"\t"<<biglist[i].substr(0,start);
	}
	outfile<<"\n";
	
	for(int i=0;i<super_point[0].size();i++)
	{
		outfile<<bigmarker[i].chro<<"\t"<<bigmarker[i].position<<"\t.\t"<<bigmarker[i].ref<<"\t"<<bigmarker[i].alt<<"\t100\tPASS\t.\tGT:GL";
		for(int k=0;k<super_point.size();k++)
		{
			outfile<<"\t./.:"<<super_point[k][i].likeAA<<","<<super_point[k][i].likeAB<<","<<super_point[k][i].likeBB;
		}
		outfile<<"\n";
	}
	outfile.close();
}



int main(int argc,char *ARGV[])
{
	if(argc != 4)
	{
		cerr<<"The parameter number is incorrect, it should be ./usage filelist file_position fileout\n";
		exit(1);
	}
	char *filelist=ARGV[1];
	char *filepos=ARGV[2];
	char *fileout=ARGV[3];
	vector<string> biglist;
	read_list(filelist,biglist);
	if(biglist.size() == 0)
	{
		cerr<<"there is no input file in the filelist\n";
		exit(1);
	}
	vector<struct marker> bigmarker;
	vector<vector<struct point> > super_point;
	vector<vector<struct marker> > chosen_marker;
	int keyslot=0;
	construct_hash(filepos,chosen_marker,bigmarker,keyslot);
	//read_marker(biglist[0],bigmarker);
	/*for(int i=0;i<chosen_marker.size();i++)
	{
		for(int k=0;k<chosen_marker[i].size();k++)
		{
			cout<<chosen_marker[i][k].chro<<"\t"<<chosen_marker[i][k].position<<"\t"<<chosen_marker[i][k].alt<<"\t"<<chosen_marker[i][k].ref<<"\n";
		}
	}*/
	//exit(1);

	read_big_data(biglist,chosen_marker,super_point,keyslot);
	//cout<<"ok here\n";
//	exit(1);
	//int serial=-1;
/*	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cout<<"Can't open file "<<fileout<<"\n";
		exit(1);
	}*/
	/*for(int i=0;i<super_point[0].size();i++)
	{
		for(int k=0;k<super_point.size();k++)
		{
			//++serial;
			outfile<<bigmarker[i].chro<<"\t"<<bigmarker[i].position<<"\t"<<bigmarker[i].alt<<"\t"<<bigmarker[i].ref<<"\t"<<super_point[k][i].go_ref<<"\t"<<super_point[k][i].ba_ref<<"\t"<<super_point[k][i].go_alt<<"\t"<<super_point[k][i].ba_alt<<"\n";
		}
	}*/
	calculate_likelihood(super_point);
	output_file(fileout,bigmarker,super_point,biglist);

}
