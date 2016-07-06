#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<bam.h>
#include<map>
#include<faidx.h>
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
        vector<char> base;
        vector<char> qual;
};

struct i_point
{
	int chro;
	int pos;
	char ref;
	vector<char> base;
	vector<char> qual;
};

int construct_hash(char *filepoly,vector<vector<struct point> > &super_point,int big_keyslot)
{
	int poly_number=0;
	ifstream infile(filepoly,std::ios::in);
	if(!infile.is_open())
	{
		cerr<<"Can't open file "<<filepoly<<"\n";
		exit(1);
	}
	
	while(!infile.eof())
	{
		string temp;
		getline(infile,temp);
		if(temp.length() < 1)
			continue;
		++poly_number;
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
		//cout<<temp<<endl;
		string chro=temp.substr(0,seppos[0]);
		int pos=atoi(temp.substr(seppos[0]+1,seppos[1] - seppos[0]-1).c_str());
		long total_add_count=0;
		int ten_plus=1;
		for(int k=0;k<chro.length();k++)
		{
			total_add_count += int(chro[k])*ten_plus*(big_keyslot+1);
			ten_plus=ten_plus*2;
		}
		total_add_count += pos*(big_keyslot+1);
		int key_value=abs(total_add_count) % big_keyslot;
		int flag=0;
		//for(int l=0;l<super_point[key_value].size();l++)
		//{
			struct point temp_point;
			temp_point.chro=chro;
			temp_point.pos=pos;
			temp_point.ref=temp[seppos[1]+1];
			super_point[key_value].push_back(temp_point);
		//}
	}
	infile.close();
	//cout<<"poly number is "<<poly_number<<"\n";
	return(poly_number);
}

void swap(int *real_num,int start,int end)
{
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
		if(super_point[real_num[i]].chro <super_point[real_num[start]].chro || (super_point[real_num[i]].chro == super_point[real_num[start]].chro && super_point[real_num[i]].pos < super_point[real_num[start]].pos))
		{
			swap(real_num,++current,i);
		}
	}
	 swap(real_num,start,current);
        qsort(super_point,real_num,start,current-1);
        qsort(super_point,real_num,current+1,end);
}



int main(int argc, char *ARGV[])
{
	if(argc != 5)
	{
		cerr<<"Incorrect number of parameters, the command should be ./programme bamfile reference_file poly_site_file output_file\n";
		exit(1);
	}
	bam1_t* bam_big=bam_init1();
	bamFile infile=bam_open(ARGV[1],"r");
	char *fileref=ARGV[2];
	char *filepoly=ARGV[3];
	char *fileout=ARGV[4];
	faidx_t *genomefaid;
	genomefaid=fai_load(fileref);
	if(genomefaid == 0)
	{
		cerr<<"Can't open the index file for the reference "<<fileref<<"\n";
		exit(1);
	}
	if(infile ==NULL)
	{
		cerr<<"Can't open the bam file for reading "<<ARGV[1];
		exit(1);
	}
/*	ifstream inpoly(filepoly,std::ios::in);
	if(!inpoly.is_open())
	{
		cout<<"Can't open polymorphism file "<<filepoly<<"\n";
		exit(1);
	}*/
	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cout<<"Can't open file for writting "<<fileout<<"\n";
		exit(1);
	}
	bam_header_t *header;
	header=bam_header_read(infile);
	int total_genome_len=0;
	for(int i=0;i<header->n_targets;i++)
	{
		total_genome_len += header->target_len[i];
	}


	map<string,int> conver;
	map<int,string> re_conver;
	for(int i=0;i<header->n_targets;i++)
	{
		
		//cout<<header->target_name[i]<<"\n";
		conver[string(header->target_name[i])]=i;
		re_conver[i]=string(header->target_name[i]);
	}
	

	int keyslot=int(total_genome_len/60);
	vector<vector<struct point> > super_point;
	super_point.resize(keyslot);
	int poly_number=0;
	//construct_hash(filepoly,super_point,keyslot);
	poly_number=construct_hash(filepoly,super_point,keyslot);




	
	/*for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			outfile<<super_point[i][k].chro<<"\t"<<super_point[i][k].pos<<"\t"<<super_point[i][k].ref<<"\n";
		}
	}*/
	//exit(1);

//	cout<<"ok here"<<endl;
	//exit(1);


	while(bam_read1(infile,bam_big) >= 0)
	{
		const bam1_core_t *bam_core=&bam_big->core;
		uint8_t *align_seq=bam1_seq(bam_big),*align_qual=bam1_qual(bam_big);
		char temp_seq[bam_core->l_qseq];
		char temp_qual[bam_core->l_qseq];
		char temp_check[bam_core->l_qseq];
	for(int i=0;i<bam_core->l_qseq-1;++i)
	{
		temp_seq[i]=bam_nt16_rev_table[bam1_seqi(align_seq,i)];
		temp_check[i]=char(align_qual[i]+33);
	}
	string str_seq=string(temp_seq);
	string str_qual=string(temp_check);
	
	char *chro=header->target_name[bam_core->tid];
	int map_quality=bam_core->qual;
	int left_pos=bam_core->pos +1;
	if(map_quality < 60)
	{
		continue;
	}
	int match_total=0;
	uint32_t *cigar=bam1_cigar(bam_big);
	for(int k=0;k<bam_core->n_cigar;k++)
	{
		int op=cigar[k] & BAM_CIGAR_MASK;
		int cl=cigar[k] >> BAM_CIGAR_SHIFT;
		if(op == 0)
		{
			match_total += cl;
		}
	}
	float ratio=1.0*match_total/bam_core->l_qseq;
	if(ratio < 0.92)
	{
		continue;
	}
	int end_point=left_pos + bam_core->l_qseq -1;
	stringstream ss;
	string str_chro;
	ss <<chro;
	ss >> str_chro;
	char reg_buffer[100];
	sprintf(reg_buffer,"%s:%d-%d",chro,left_pos,end_point + int(bam_core->l_qseq*0.5));
	const char *region=reg_buffer;
	int str_length=int(bam_core->l_qseq * 1.5);
	char *fetch=fai_fetch(genomefaid,region,&str_length);
	string str_fetch=string(fetch);
	int start_point=0;
	int indel_number=0;
	for(int i=0;i<bam_core->n_cigar;i++)
	{
		int op=cigar[i] & BAM_CIGAR_MASK;
                int cl=cigar[i] >> BAM_CIGAR_SHIFT;
		if(op == 4)
		{
			for(int k=0;k<cl;k++)
			{
				str_seq[k+start_point]='-';
			}
			start_point += cl;
		}
		else if(op == 0)
		{
			start_point += cl;
		}
		else if(op == 2)
		{
			str_seq.insert(start_point,cl,'N');
			str_qual.insert(start_point,cl,'#');
			start_point += cl;
		}
		else if(op == 1)
		{
			str_seq.erase(start_point,cl);
			str_qual.erase(start_point,cl);
			start_point -= cl;
		}
		else
		{
			cerr<<"abormal happens\n";
			exit(1);
		}
	}
	int add_ref=-1;
	int pos_cord=left_pos -1;
	for(int i=0;i<str_seq.length();i++)
	{
		if(str_seq[i] == '-')
		{
			continue;
		}
		++pos_cord;
		++add_ref;
		if(str_seq[i] == 'N')
		{
			continue;
		}
		long total_add_count=0;
		int ten_plus=1;
		for(int k=0;k<str_chro.length();k++)
		{
			total_add_count += int(str_chro[k])*ten_plus*(keyslot+1);
			ten_plus=ten_plus*2;
		}
		total_add_count += pos_cord*(keyslot+1);
		int key_value=abs(total_add_count)%keyslot;
		int flag=0;
		for(int l=0;l<super_point[key_value].size();l++)
		{
			if(super_point[key_value][l].chro.compare(str_chro) == 0 && super_point[key_value][l].pos == pos_cord)
			{
				super_point[key_value][l].base.push_back(str_seq[i]);
				super_point[key_value][l].qual.push_back(str_qual[i]);
				super_point[key_value][l].ref=str_fetch[add_ref];
				flag=1;
				break;
			}
		}
	}
	/*for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			outfile<<super_point[i][k].chro<<"\t"<<super_point[i][k].pos<<"\t"<<super_point[i][k]
	*/





	}//this is the end of the while loop
	
	struct i_point *big_point;
	big_point=new struct i_point[poly_number];
	//cout<<poly_number<<"\n";
	int *num_buck;
	num_buck=new int[poly_number];

	int serial_num=0;
	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			//big_point[serial_num]=super_point[i][k];
			big_point[serial_num].chro=conver[super_point[i][k].chro];
		//	cout<<big_point[serial_num].chro<<"\n";
			big_point[serial_num].pos=super_point[i][k].pos;
			big_point[serial_num].ref=super_point[i][k].ref;
			big_point[serial_num].base=super_point[i][k].base;
			big_point[serial_num].qual=super_point[i][k].qual;
			num_buck[serial_num]=serial_num;
			++serial_num;
		}
	}
	//exit(1);
	qsort(big_point,num_buck,0,poly_number-1);
	
	for(int i=0;i<poly_number;i++)
	{
		outfile<<re_conver[big_point[num_buck[i]].chro]<<"\t"<<big_point[num_buck[i]].pos<<"\t"<<big_point[num_buck[i]].ref<<"\t";
		if(big_point[num_buck[i]].base.size() == 0)
		{
			outfile<<"*";
		}
		else
		{
			for(int l=0;l<big_point[num_buck[i]].base.size();l++)
			{
				outfile<<big_point[num_buck[i]].base[l];
			}
		}
		outfile<<"\t";
		if(big_point[num_buck[i]].qual.size() == 0)
		{
			outfile<<"*";
		}
		else
		{
			for(int l=0;l<big_point[num_buck[i]].qual.size();l++)
			{
				outfile<<big_point[num_buck[i]].qual[l];
			}
		}
		outfile<<"\n";
	}
/*	for(int i=0;i<super_point.size();i++)
        {
                for(int k=0;k<super_point[i].size();k++)
                {
                        *if(super_point[i][k].ref || super_point[i][k].ref == ' ')
 *                         {
 *                                                         continue;
 *                                                                                 }

                        outfile<<super_point[i][k].chro<<"\t"<<super_point[i][k].pos<<"\t"<<super_point[i][k].ref<<"\t";
			if(super_point[i][k].base.size() == 0)
			{
				outfile<<"*";
			}
			else
			{
                        	for(int l=0;l<super_point[i][k].base.size();l++)
                        	{
                                	outfile<<super_point[i][k].base[l];
                        	}
			}
                        outfile<<"\t";
			if(super_point[i][k].qual.size() == 0)
			{
				outfile<<"*";
			}
			else
			{
                        	for(int l=0;l<super_point[i][k].qual.size();l++)
                        	{
                                	outfile<<super_point[i][k].qual[l];
                        	}
			}
                        outfile<<"\n";
                }
        }*/
        outfile.close();
        bam_header_destroy(header);
        bam_close(infile);
        bam_destroy1(bam_big);
}


