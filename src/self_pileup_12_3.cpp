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

int main(int argc,char *ARGV[])
{
	if(argc != 4)
	{
		cerr<<"incorrect number of parameters, the command should be ./usage bamfile reference_file pile_output_file\n";
		exit(1);
	}
	bam1_t* bam_big=bam_init1();//bam file alignment init
	bamFile infile=bam_open(ARGV[1],"r");// ARGV[1] is the bam filehandle
	char *fileref=ARGV[2];// this is the reference sequence and it must be indexed using faidx
	char *fileout=ARGV[3];// this is output pileup file
	faidx_t *genomefaid;
	genomefaid=fai_load(fileref);
	if(genomefaid == 0)
	{
		cerr<<"Can't open the index file for the reference file "<<fileref<<"\n";
		exit(1);
	}
	bam_header_t *header;// drefine bam header 
	if(infile == NULL)
	{
		cerr<<"Can't open the bam file for reading "<<ARGV[1]<<"\n";
		exit(1);
	}
	if( bam_big == NULL)
	{
		cerr<<"Can't open inint the bam file \n";
		exit(1);
	}
	ofstream outfile(fileout,std::ios::out);
	if(!outfile.is_open())
	{
		cout<<"Can't open file "<<fileout<<"\n";
		exit(1);
	}
	header=bam_header_read(infile);
	int total_genome_len=0;
	for(int i=0;i<header->n_targets;i++)
	{
		total_genome_len += header->target_len[i];
	}
	//cout<<"the total len is "<<total_genome_len<<"\n";
	int keyslot=int(total_genome_len/60);
	//ut<<keyslot<<"\n";
	//it(1);
	vector<vector<struct point> > super_point;
	super_point.resize(keyslot);
	while(bam_read1(infile,bam_big) >= 0)
	{
		const bam1_core_t *bam_core=&bam_big->core;
		uint8_t *align_seq=bam1_seq(bam_big),*align_qual=bam1_qual(bam_big);
		//cout<<bam1_qname(bam_big)<<"\n";
		char temp_seq[bam_core->l_qseq];
		int temp_qual[bam_core->l_qseq];
		char temp_check[bam_core->l_qseq];
		//cout<<bam_core->l_qseq<<"\n";
		for(int i=0;i< bam_core->l_qseq-1;++i)
		{
			//fputc(bam_nt16_rev_table[bam1_seqi(align_seq,i)],stdout);
			//fputc('\n',stdout);
			temp_seq[i]=bam_nt16_rev_table[bam1_seqi(align_seq,i)];
			//temp_qual[i]=int(align_qual[i]);
			temp_check[i]=char(align_qual[i]+33);
		}
		//fputc('\n',stdout);
		//cout<<temp_seq<<"\n"<<temp_qual<<"\n"<<temp_check<<"\n\n";
		string str_seq=string(temp_seq);
		string str_qual=string(temp_check);
		//string str_qual=string(temp_check);
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
			int cl=cigar[k]>>BAM_CIGAR_SHIFT;
			if(op == 0)
			{
				match_total += cl;
			}
		}
		//cout<<bam_core->l_qseq<<"\n";
		float ratio=1.0*match_total/bam_core->l_qseq;
		//cout<<match_total<<"\n";
		if(ratio < 0.92)
		{
			continue;
		}
		//cout<<str_seq<<"\n";
		//cout<<bam1_qname(bam_big)<<"\n"<<header->target_name[bam_core->tid]<<"\n"<<str_seq<<"\n"<<str_qual<<"\n"<<bam_core->pos<<"\n"<<bam_core->qual<<"\n\n";
		//int start_point=0;
		//int indel_number=0;
		int end_point=left_pos + bam_core->l_qseq -1;
		stringstream ss;
	//	cout<<left_pos<<"\t"<<end_point<<"\n";
		//string merge=string(chro) + ":" +string(left_pos) + "-" + string(end_point);
		//cout<<merge<<"\n";
		//chro >> ss;
		//string str_chro,str_start,str_end;
		string str_chro;
		ss << chro;
		ss >> str_chro;
		/*ostringstream oss;
		ss << left_pos;
		str_start=oss.str();
		//str_end=to_string(end_point);
		ss << end_point;
		str_end=oss.str();
		string region_str=str_chro + ":" + str_start + "-" + str_end;
		cout<<region_str<<"\n";*/
		char reg_buffer[100];
		sprintf(reg_buffer,"%s:%d-%d",chro,left_pos,end_point + int(bam_core->l_qseq*0.5));
		const char *region=reg_buffer;
		int str_length=int(bam_core->l_qseq * 1.5);
		//cout<<chro<<"\t"<<left_pos<<"\t"<<end_point<<endl;
		//sprintf(region,"%s:%d-%d",chro,left_pos,end_point);
		//int qseq_number=bam_core->l_qseq;
		//cout<<region<<"\n"<<qseq_number<<endl;
		//cout<<region<<"\n";
		char *fetch=fai_fetch(genomefaid,region,&str_length);
		//string refer_seq=string(fetch);
		//cout<<str_seq<<"\n"<<refer_seq<<"\n";
	//	cout<<str_seq<<"\n"<<fetch<<"\n";
/*		if(bam_core->n_cigar == 1)
		{
			continue;
		}*/
		//cout<<str_seq<<"\n"<<fetch<<"\n";
		string str_fetch=string(fetch);
		//cout<<bam_core->n_cigar<<"\n";
		int start_point=0;
		int indel_number=0;
		for(int i=0;i<bam_core->n_cigar;i++)
		{
			int op=cigar[i] & BAM_CIGAR_MASK;
			int cl=cigar[i]>>BAM_CIGAR_SHIFT;
			//cout<<"cigar "<<cl<<op<<"\n";
		}
		//cout<<"\n";
		//int pos_cord=bam_core->pos;
		for(int i=0;i<bam_core->n_cigar;i++)
		{
			int op=cigar[i] & BAM_CIGAR_MASK;
                        int cl=cigar[i]>>BAM_CIGAR_SHIFT;
		//	cout<<cl<<op<<"\n";
			if(op == 4)
			{
				for(int k=0;k<cl;k++)
				{
					str_seq[k + start_point]='-';
					//str_seq.erase(start_point,cl);
					//str_qual.erase(start_point,cl);
					
				}
				start_point +=cl;
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
				cerr<<"abormal happens"<<"\n";
				exit(1);
			}
		}
	//	cout<<"\n";
		/*for(int i=0;i<str_seq.length();i++)
		{
			if(str_seq[i] != str_fetch[i])
			{
				cout<<str_seq[i]<<"\t"<<str_fetch[i]<<"\t"<<str_qual[i]<<"*\n";
			}
			else
			{
				cout<<str_seq[i]<<"\t"<<str_fetch[i]<<"\t"<<str_qual[i]<<"\n";
			}
		}*/
		int add_ref=-1;
		int pos_cord=left_pos-1;
		for(int i=0;i<str_seq.length();i++)
		{
			//pos_cord;
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
			if(str_seq[i] != str_fetch[add_ref])
			{
			//	cout<<str_seq<<"\n"<<fetch<<"\n";
			//	cout<<region<<"\n";
			//	cout<<pos_cord<<"\t"<<left_pos<<"\t"<<str_seq[i]<<"\t"<<str_fetch[add_ref]<<"\n";
				long total_add_count=0;
				int ten_plus=1;
				//cout<<str_chro.length()<<"\t"<<str_chro<<endl;
				for(int k=0;k<str_chro.length();k++)
				{
					total_add_count += int(str_chro[k])*ten_plus*(keyslot+1);
					ten_plus=ten_plus*2;
				}
				total_add_count += pos_cord*(keyslot+1);
		//		cout<<total_add_count<<"\n";
				int key_value=abs(total_add_count) % keyslot;
				int flag=0;
				//cout<<str_chro<<"\t"<<pos_cord<<"\t"<<key_value<<"\n";
				//cout<<key_value<<endl;
				//cout<<super_point[key_value].size()<<"\n";
				for(int l=0;l<super_point[key_value].size();l++)
				{
				//	cout<<super_point[key_value][l].chro<<"\t"<<str_chro<<"\t"<<super_point[key_value][l].pos<<"\t"<<pos_cord<<"\n";
					if(super_point[key_value][l].chro.compare(str_chro) == 0 && super_point[key_value][l].pos == pos_cord)
					{
						//cout<<"ok\n";
						
						super_point[key_value][l].base.push_back(str_seq[i]);
						super_point[key_value][l].qual.push_back(str_qual[i]);
						flag=1;
						break;
					}
				}
				if(flag ==0)
				{
					struct point temp_point;
					temp_point.chro=str_chro;
					temp_point.pos=pos_cord;
				//	cout<<temp_point.pos<<" this is the pos_cord\n";
					temp_point.ref=str_fetch[add_ref];
					temp_point.base.push_back(str_seq[i]);
					temp_point.qual.push_back(str_qual[i]);
					super_point[key_value].push_back(temp_point);
				}
			}
			//++pos_cord;
		}
		//exit(1);
	}
	for(int i=0;i<super_point.size();i++)
	{
		for(int k=0;k<super_point[i].size();k++)
		{
			/*if(super_point[i][k].ref || super_point[i][k].ref == ' ')
			{
				continue;
			}*/
			
			outfile<<super_point[i][k].chro<<"\t"<<super_point[i][k].pos<<"\t"<<super_point[i][k].ref<<"\t";
			for(int l=0;l<super_point[i][k].base.size();l++)
			{
				outfile<<super_point[i][k].base[l];
			}
			outfile<<"\t";
			for(int l=0;l<super_point[i][k].qual.size();l++)
                        {
                                outfile<<super_point[i][k].qual[l];
                        }
			outfile<<"\n";
		}
	}
	outfile.close();
	bam_header_destroy(header);
	bam_close(infile);
	bam_destroy1(bam_big);
}

		
