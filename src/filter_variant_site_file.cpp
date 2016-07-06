#include<iostream>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<map>


using namespace std;

struct base_point
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
				if(seppos.size() >= 2)
					break;
			}
		}
		struct marker temp_marker;
		temp_marker.chro=temp.substr(0,seppos[0]);
		temp_marker.position=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0] -1).c_str());
		temp_marker.ref=temp[seppos[1]+1];
		bigmarker.push_back(temp_marker);
	}
	infile.close();
}

void count_allele(vector<string> &biglist,vector<vector<int> > &big_count)
{
	map<char,int> conver;
	conver['A']=0;
	conver['T']=1;
	conver['G']=2;
	conver['C']=3;
	for(int i=0;i<biglist.size();i++)
	{
		cout<<"Reading file "<<biglist[i]<<"\n";
		ifstream infile(biglist[i].c_str(),std::ios::in);
		if(!infile.is_open())
		{
			cerr<<"Can't open file "<<biglist[i]<<"\n";
			exit(1);
		}
		int serial=-1;
		while(!infile.eof())
		{
			string temp;
			getline(infile,temp);
			if(temp.length() < 1)
				continue;
			++serial;
			vector<int> seppos;
			for(int k=0;k<temp.length();k++)
			{
				if(temp[k] == ' ' || temp[k] == '\t')
				{
					seppos.push_back(k);
					if(seppos.size() == 4)
						break;
				}
			}
			map<char,int> count;
			for(int k=seppos[2]+1,l=seppos[3]+1;k<seppos[3];k++,l++)
			{
				if(temp[k] == '*' && temp[l] == '*')
				{
					++big_count[serial][4];
					break;
				}
				if(temp[l] - 33 >= 20)
				{
					if(count.count(temp[k]) > 0)
					{
						continue;
					}
					else
					{
						if(conver.count(temp[k]) > 0)
						{
							++big_count[serial][conver[temp[k]]];
							count[temp[k]]=1;
						}
					}
				}
			}
		}
		infile.close();
	}
}

int main(int argc,char *ARGV[])
{
	if(argc != 5)
	{
		cerr<<"Incorrect number of parameters. it should be ./usage filelist maf_cretria miss_cretria fileout\n";
		exit(1);
	}
	char *filelist=ARGV[1];
	float flt_maf=atof(ARGV[2]);
	float flt_miss=atof(ARGV[3]);
	char *fileout=ARGV[4];
	vector<string> biglist;
	read_list(filelist,biglist);
	vector<struct marker> bigmarker;
	vector<vector<struct base_point> > super_point;
	if(biglist.size() == 0)
	{
		cerr<<"there is no input file in the file list "<<filelist<<"\n";
		exit(1);
	}
	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cout<<"can't open file "<<fileout<<"\n";
		exit(1);
	}
	read_marker(biglist[0],bigmarker);
	vector<vector<int> > big_count;
	//big_count.resize(bigmarker.size());
	for(int i=0;i<bigmarker.size();i++)
	{
		vector<int> test;
		test.resize(5,0);
		big_count.push_back(test);
	}
	count_allele(biglist,big_count);
	map<int,char> re_conver;
	re_conver[0]='A';
	re_conver[1]='T';
	re_conver[2]='G';
	re_conver[3]='C';

	map<char,int> conver;
        conver['A']=0;
        conver['T']=1;
        conver['G']=2;
        conver['C']=3;

	for(int i=0;i<big_count.size();i++)
	{
		int max=0;
		int id=0;
		for(int k=0;k<4;k++)
		{
			if(k == conver[bigmarker[i].ref])
			{
				continue;
			}
			if(big_count[i][k] > max)
			{
				max=big_count[i][k];
				id=k;
			}
		}
		if(max == 0)
		{
			continue;
		}
		float maf_rate=0.0;
		int ref_count=big_count[i][conver[bigmarker[i].ref]];
		if(ref_count > max)
		{
			maf_rate=1.0*max/(max + ref_count);
		}
		else
		{
			maf_rate=1.0*ref_count/(max + ref_count);
		}
		float miss_rate=1.0*big_count[i][4]/(biglist.size());
		if(miss_rate > flt_miss || maf_rate < flt_maf)
		{
			continue;
		}
	//	cout<<bigmarker[i].alt<<"\n";
		bigmarker[i].alt=re_conver[id];
	//	cout<<big_count[i][0]<<"\t"<<big_count[i][1]<<"\t"<<big_count[i][2]<<"\t"<<big_count[i][3]<<"\t"<<id<<"\t"<<max<<"\n";
	//	cout<<bigmarker[i].alt<<"\n";
		outfile<<bigmarker[i].chro<<"\t"<<bigmarker[i].position<<"\t"<<bigmarker[i].ref<<"\t"<<bigmarker[i].alt<<"\t"<<ref_count<<"\t"<<max<<"\t"<<maf_rate<<"\t"<<miss_rate<<"\n";
	}
}
