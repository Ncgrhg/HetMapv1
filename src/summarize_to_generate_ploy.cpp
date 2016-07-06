#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<map>
//#include<bam.h>
//#include<faidx.h>
#include<vector>
#include<string>
#include<sstream>
#include<fstream>
#define WHERE fprintf(stderr,"%d\n",__LINE__)

using namespace std;

struct point
{
	string chro;
	int pos;
	char ref;
	//vector<char>base;
	vector<int>number;
};
struct i_point
{
	int chro;
	int pos;
	char ref;
	vector<int> number;
};


void read_genome(char *fileref,map<string,int> &conver,int &total_len)
{
	ifstream infile(fileref,std::ios::in);
	if(!infile.is_open())
	{
		cout<<"Can't open file "<<fileref<<"\n";
		exit(1);
	}
	vector<string> chro_name;
	while(!infile.eof())
	{
		string temp;
		getline(infile,temp);
		if(temp.length() < 1)
			continue;
		if(temp[0] == '>')
		{
			vector<int> seppos;
			for(int i=1;i<temp.length();i++)
			{
				if(temp[i] == ' ' || temp[i] == '\t')
				{
					seppos.push_back(i);
					break;
				}
			}
			if(seppos.size() == 0)
			{
				seppos.push_back(temp.length()-1);
			}
			chro_name.push_back(temp.substr(1,seppos[0]));
		}
		else
		{
			total_len += temp.length();
		}
	}
	for(int i=0;i<chro_name.size();i++)
	{
		conver[chro_name[i]]=i;
	}
	infile.close();
}

void read_file_list(char *filebiglist,map<string,int> &conver,int keyslot,int total_len,vector<vector<struct point> > &super_point)
{
	vector<string> filelist;
	ifstream infile(filebiglist,std::ios::in);
	if(!infile.is_open())
	{
		cout<<"Can't open file "<<filebiglist<<"\n";
		exit(1);
	}
	while(!infile.eof())
	{
		string temp;
		getline(infile,temp);
		if(temp.length() < 1)
			continue;
		filelist.push_back(temp);
	}
	infile.close();
	map<char,int> base_change;
	base_change['A']=0;
	base_change['T']=1;
	base_change['G']=2;
	base_change['C']=3;
	int random_seed=int(total_len/23);
	for(int i=0;i<filelist.size();i++)
	{
		ifstream readin(filelist[i].c_str(),std::ios::in);
		if(!readin.is_open())
		{
			cout<<"Can't open file "<<filelist[i]<<"\n";
			exit(1);
		}
		cout<<"reading file "<<filelist[i]<<"\n";
		while(!readin.eof())
		{
			string temp;
			getline(readin,temp);
			if(temp.length() < 1)
				continue;
			vector<int> seppos;
			for(int k=0;k<temp.length();k++)
			{
				if(temp[k] == '\t' || temp[k] == ' ')
				{
					seppos.push_back(k);
				}
			}
			//cout<<temp<<"\t"<<seppos.size()<<endl;
			if(seppos.size() == 0)
			{
				cerr<<"wrong input file format\n";
			//	cout<<temp<<endl;
				exit(1);
			}
			if(seppos.size() != 4)
			{
				cout<<temp<<endl;
				continue;
			}
			string chro=temp.substr(0,seppos[0]);
			int position=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0] -1).c_str());
			char refbase=temp[seppos[1]+1];
			string seq_base=temp.substr(seppos[2]+1, seppos[3] -seppos[2]-1);
			string seq_qual=temp.substr(seppos[3]+1,temp.length() - seppos[3] -1);
			long total_add_count=0;
			int ten_plus=1;
			for(int k=0;k<chro.length();k++)
			{
				total_add_count += int(chro[k])*ten_plus*(random_seed+1);
				ten_plus=ten_plus*2;
			}
			total_add_count += position*(random_seed+1);
			int key_value=abs(total_add_count) % keyslot;
			int flag=0;
			for(int l=0;l<super_point[key_value].size();l++)
			{
				if(super_point[key_value][l].chro.compare(chro) == 0 && super_point[key_value][l].pos == position)
				{
					map<char,int> big_hash;
					for(int m=0;m<seq_base.length();m++)
					{
						if (int(seq_qual[m]) - 33 < 20)
						{
							continue;
						}
						if(base_change.count(seq_base[m]) > 0)
						{
							//++super_point[key_value][l].number[base_change[seq_base[m]]];
							if(big_hash.count(seq_base[m]) > 0)
							{
								continue;
							}
							else
							{
								++super_point[key_value][l].number[base_change[seq_base[m]]];
								big_hash[seq_base[m]]=1;
							}
						}
					}
					flag=1;
				}
			}
			if(flag == 0)
			{
				struct point temp_point;
				temp_point.chro=chro;
				temp_point.pos=position;
				temp_point.ref=refbase;
				temp_point.number.resize(4,0);
				map<char,int> big_hash;
				for(int m=0;m<seq_base.length();m++)
				{
					if(int(seq_qual[m]) - 33 < 20)
					{
						continue;
					}
					if(base_change.count(seq_base[m]) > 0)
					{
						if(big_hash.count(seq_base[m]) > 0)
						{
							continue;
						}
						else
						{
							++temp_point.number[base_change[seq_base[m]]];
							 big_hash[seq_base[m]]=1;
						}
					}
				}
				super_point[key_value].push_back(temp_point);
			}
			
			
		}
	}
//	cout<<"it's ok here\n";
//	exit(1);
}
	

void count_poly_site(vector<vector<struct point> > &super_point,int &total_poly,int number)
{
	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			int flag=0;
			for(int m=0;m<4;m++)
			{
				if(super_point[i][k].number[m] >= number)
				{
					flag=1;
					break;
				}
			}
			if(flag == 1)
			{
				++total_poly;
			}
		}
	}
}

void swap(int *real_num,int start,int end)
{
	/*struct point temp;
	temp=super_point[start];
	super_point[start]=super_point[end];
	super_point[end]=temp;*/
	//cout<<real_num[start]<<"\t"<<real_num[end]<<"\n";
	int temp;
	temp=real_num[start];
	real_num[start]=real_num[end];
	real_num[end]=temp;
}

void qsort(struct i_point *super_point,int *real_num,int start,int end)
{
	if(!(start < end))
		return;
	swap(real_num,start,(start + end)/2);
	int current=start;
	for(int i=start+1;i<=end;i++)
	{
		//cout<<super_point[real_num[i]].chro<<"\t"<<super_point[real_num[start]].chro<<endl;
		//exit(1);
		if(super_point[real_num[i]].chro <super_point[real_num[start]].chro || (super_point[real_num[i]].chro == super_point[real_num[start]].chro && super_point[real_num[i]].pos < super_point[real_num[start]].pos))
		{
			swap(real_num,++current,i);
		}
	}
	swap(real_num,start,current);
	qsort(super_point,real_num,start,current-1);
	qsort(super_point,real_num,current+1,end);
}

int main(int argc,char *ARGV[])
{
	if(argc != 5)
	{
		cerr<<"incorrect parameter number, it should be ./programme referece Pile_filelist fileout snp_number\n";
		exit(1);
	}
	char *fileref=ARGV[1];
	char *filelist=ARGV[2];
	char *fileout=ARGV[3];
	int number=atoi(ARGV[4]);
	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cerr<<"Can't open file "<<fileout<<"\n";
		exit(1);
	}
	map<string,int> conver;
	map<int,string> re_conver;
	map<string,int>::iterator it;
	//for(it=conver.begin();it!=conver.end();it++)
	//{
	//	re_conver[it->second]=it->first;
	//}
	int total_len=0;
	read_genome(fileref,conver,total_len);
	for(it=conver.begin();it!=conver.end();it++)
        {
                re_conver[it->second]=it->first;
        }
	//map<string,int>::iterator it;
	/*for(it=conver.begin();it!=conver.end();it++)
	{
		cout<<it->first<<"\t"<<it->second<<"\n";
	}
	//int keyslot=int(total_len/30);
	cout<<total_len<<"\n";*/
	int keyslot=int(total_len/30);
	vector<vector<struct point> > super_point;
	super_point.resize(keyslot);
	read_file_list(filelist,conver,keyslot,total_len,super_point);
/*	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			outfile<<super_point[i][k].chro<<"\t"<<super_point[i][k].ref<<"\t"<<super_point[i][k].pos<<"\t"<<super_point[i][k].number[0]<<"\t"<<super_point[i][k].number[1]<<"\t"<<super_point[i][k].number[2]<<"\t"<<super_point[i][k].number[3]<<"\n";
		}
	}*/
	//exit(1);
	int total_poly=0;
	count_poly_site(super_point,total_poly,number);
	struct i_point *big_point;
	big_point=new struct i_point[total_poly];
	int serial=0;
	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			//big_point[serial]=super_point[i][k];
			int flag=0;
			for(int m=0;m<4;m++)
			{
				if(super_point[i][k].number[m] >= number)
				{
					flag=1;
					break;
				}
			}
			if(flag == 0)
			{
				continue;
			}
			big_point[serial].chro=conver[super_point[i][k].chro];
			big_point[serial].pos=super_point[i][k].pos;
			big_point[serial].ref=super_point[i][k].ref;
			big_point[serial].number=super_point[i][k].number;
			++serial;
		}
	}
	int *real_num;
	real_num=new int[total_poly];
	for(int i=0;i<total_poly;i++)
	{
		real_num[i]=i;
		//outfile<<big_point[i].chro<<"\t"<<big_point[i].pos<<"\t"<<big_point[i].ref<<"\t"<<big_point[real_num[i]].number[0]<<"\t"<<big_point[real_num[i]].number[1]<<"\t"<<big_point[real_num[i]].number[2]<<"\t"<<big_point[real_num[i]].number[3]<<"\n";
	}
	//exit(1);
	qsort(big_point,real_num,0,total_poly-1);
//	cout<<"ok here"<<endl;
	for(int i=0;i<total_poly;i++)
	{
		outfile<<re_conver[big_point[real_num[i]].chro]<<"\t"<<big_point[real_num[i]].pos<<"\t"<<big_point[real_num[i]].ref<<"\t"<<big_point[real_num[i]].number[0]<<"\t"<<big_point[real_num[i]].number[1]<<"\t"<<big_point[real_num[i]].number[2]<<"\t"<<big_point[real_num[i]].number[3]<<"\n";
	}
		
}

