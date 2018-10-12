/*
Copyright (C) 2017-2018 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/

#ifndef FASTQ_H_
#define FASTQ_H_

#include <iostream>
#include <getopt.h>
#include <cmath>

#include "data_container.h"
#include "algorithm/compression/compression_manager.h"
#include "algorithm/compression/rans_static.h"
#include "algorithm/compression/qual_codec.h"
#include "algorithm/open_hashtable.h"
#include "algorithm/timer.h"

const uint8_t BASE_MAP[256] =
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,3,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
 	 0,0,0,0,0,0,0,0};

const uint8_t BASE_MAP_COMP[5] = {1,0,3,2,4};

const static char BASE_INV_MAP[5] = {'A','T','G','C','N'};

#define BLOCK_SIZE 50000
#define KMER_SHIFT 2
#define KMER_SHIFT_MASK ((1 << KMER_SHIFT) - 1)
#define KMER_EXTRACT_POS(DATA,P,NBITS) ( (DATA >> (NBITS - ((P*KMER_SHIFT) + KMER_SHIFT)) ) & ((1 << KMER_SHIFT) - 1) )
#define KMER_EXTRACT_REV_POS(DATA,P,BASES) ( DATA >> (((BASES-1)-P)*2) & ((1 << KMER_SHIFT) - 1) )

/****************************
*  Core fastq block
****************************/
struct yon1_fq_block_t {
	yon1_fq_block_t() : n(0){
		bases_packed.header.data_header.controller.type = tachyon::YON_TYPE_8B;
		quals.header.data_header.controller.type = tachyon::YON_TYPE_CHAR;
	}
	typedef tachyon::yon1_dc_t dc;

	uint64_t n;
	dc bases_matches;
	dc bases_matches_contig;
	dc bases_packed;
	dc quals;
	dc streams[8];
};

/****************************
*  Core compression
****************************/
struct stats_test {
	stats_test() : bases(0), quals(0), names(0), bases_c(0), quals_c(0), names_c(0){}
	uint64_t bases, quals, names;
	uint64_t bases_c, quals_c, names_c;
};

static bool CompressFqz(RangeCoder& rc, fqz& fqz_ins, tachyon::yon1_dc_t& dc){
	dc.data.reset();
	dc.data.resize(dc.data_uncompressed.size()*4 + 65536);
	rc.output(dc.data.data());
	rc.StartEncode();
	fqz_ins.encode_qual(&rc,
	                     dc.data_uncompressed.data(),
	                     dc.data_uncompressed.size());

	rc.FinishEncode();
	size_t sz3 = rc.size_out();
	dc.data.n_chars_ = sz3;
	//std::cerr << "fqz->" << sz3 << std::endl;
	//std::cerr << "quals="<< dc.data_uncompressed.size() << "->" << dc.data.size() << "\t" << (float)dc.data_uncompressed.size()/dc.data.size() << std::endl;

	//total_filesize += container.bases.data.size() + container.bases.strides.size();
	//stats.quals_c += dc.data.size() + dc.strides.size();
	//total_filesize += dc.data.size() + dc.strides.size();
	//dc.reset();

	return true;
}

template <class uint_t, class int_t>
bool DeltaEncodeInternal(tachyon::yon1_dc_t& dc){
	if(dc.data_uncompressed.size() == 0)
		return true;

	dc.data.resize(dc.data_uncompressed.size() + 65536);

	const uint_t* data = reinterpret_cast<const uint_t*>(dc.data_uncompressed.data());
	const uint32_t n_e = dc.data_uncompressed.size() / sizeof(int_t);
	uint_t* out = reinterpret_cast<uint_t*>(dc.data.data());

	out[0] = data[0];
	for(int i = 1; i < n_e; ++i){
		int_t delta = (int_t)data[i] - (int_t)data[i-1];
		uint_t zig = (delta << 1) ^ (delta >> 31);
		//std::cerr << "," << zig;
		out[i] = zig;
	}
	//std::cerr << std::endl;

	dc.data.n_chars_ = dc.data_uncompressed.size();
	dc.data_uncompressed = std::move(dc.data);
	dc.data.reset();

	return true;
}

static bool DeltaEncode(tachyon::yon1_dc_t& dc){
	if(dc.data_uncompressed.size() == 0){
		std::cerr << "is empty return" << std::endl;
		return true;
	}

	tachyon::TACHYON_CORE_TYPE type = dc.GetDataPrimitiveType();
	switch(type){
	case(tachyon::YON_TYPE_8B):  return(DeltaEncodeInternal<uint8_t,  int8_t>(dc));  break;
	case(tachyon::YON_TYPE_16B): return(DeltaEncodeInternal<uint16_t, int16_t>(dc)); break;
	case(tachyon::YON_TYPE_32B): return(DeltaEncodeInternal<uint32_t, int32_t>(dc)); break;
	case(tachyon::YON_TYPE_64B): return(DeltaEncodeInternal<uint64_t, int64_t>(dc)); break;
	default: std::cerr << "unknown type: " << (int)type << std::endl; break;
	}

	return false;
}

bool CompressBlock(yon1_fq_block_t& container,
		fqz& fqz_ins,
		RangeCoder& rc,
		tachyon::algorithm::CompressionManager& cm,
		stats_test& stats,
		const int32_t rid = -1)
{
	container.n = 0;
	if(rid != -1){
		std::cerr << tachyon::utility::timestamp("DEBUG") << "contig-" << rid << std::endl;
	} else {
		std::cerr << tachyon::utility::timestamp("DEBUG") << "default Hits=" << std::endl;
	}
	uint32_t total_new = 0;
	uint32_t total_new_comp = 0;

	if(container.quals.data_uncompressed.size()){
		container.quals.UpdateContainer(true, true);
		CompressFqz(rc, fqz_ins, container.quals);
		std::cerr << "quals="<< container.quals.data_uncompressed.size() << "->" << container.quals.data.size() << "\t" << (float)container.quals.data_uncompressed.size()/container.quals.data.size() << std::endl;
		stats.quals_c += container.quals.data.size() + container.quals.strides.size();
		container.quals.reset();
		container.quals.header.data_header.controller.type = tachyon::YON_TYPE_CHAR;
	} else
		return true; // no data if no quality available most likely

	// If packed bases
	if(container.bases_packed.data_uncompressed.size()){
		container.bases_packed.UpdateContainer(false, true);
		cm.zstd_codec.Compress(container.bases_packed);
		std::cerr << "bases_packed="<< container.bases_packed.data_uncompressed.size() << "->" << container.bases_packed.data.size() << "\t" << (float)container.bases_packed.data_uncompressed.size()/container.bases_packed.data.size() << std::endl;
		stats.bases_c += container.bases_packed.data.size() + container.bases_packed.strides.size();
		container.bases_packed.reset();
		container.bases_packed.header.data_header.controller.type = tachyon::YON_TYPE_8B;
	}

	//
	if(container.bases_matches.data_uncompressed.size()){
		container.bases_matches.UpdateContainer(true, true);
		if(DeltaEncode(container.bases_matches) == false){
			std::cerr << "failed to compress matches" << std::endl;
		}
		cm.zstd_codec.Compress(container.bases_matches);

		//cm.zstd_codec.Compress(container.bases_matches);
		std::cerr << "bases_matches="<< container.bases_matches.data_uncompressed.size() << "->" << container.bases_matches.data.size() << "\t" << (float)container.bases_matches.data_uncompressed.size()/container.bases_matches.data.size() << std::endl;
		stats.bases_c += container.bases_matches.data.size() + container.bases_matches.strides.size();
		container.bases_matches.reset();
	}

	if(container.bases_matches_contig.data_uncompressed.size()){
		container.bases_matches_contig.UpdateContainer(true, true);
		cm.zstd_codec.Compress(container.bases_matches_contig);

		//cm.zstd_codec.Compress(container.bases_matches);
		std::cerr << "bases_matches_contig="<< container.bases_matches_contig.data_uncompressed.size() << "->" << container.bases_matches_contig.data.size() << "\t" << (float)container.bases_matches_contig.data_uncompressed.size()/container.bases_matches_contig.data.size() << std::endl;
		stats.bases_c += container.bases_matches_contig.data.size() + container.bases_matches_contig.strides.size();
		container.bases_matches_contig.reset();
	}

	//total_filesize += container.bases_packed.data.size() + container.bases_packed.strides.size();

	for(int p = 0; p < 8; ++p){
		if(container.streams[p].data_uncompressed.size() == false){
			container.streams[p].reset();
			continue;
		}

		container.streams[p].UpdateContainer(true, true);

		if(p == 1 || p == 5 || p == 6){
			//std::cerr << "before delta encode" << std::endl;
			if(DeltaEncode(container.streams[p]) == false){
				std::cerr << "failed" << std::endl;
				exit(1);
			}
		}

		total_new += container.streams[p].data_uncompressed.size() + container.streams[p].strides_uncompressed.size();
		cm.zstd_codec.Compress(container.streams[p]);

		std::cerr << "stream-" << p << "\t" <<
				container.streams[p].data_uncompressed.size() << "\t" << container.streams[p].strides_uncompressed.size() << "\t" <<
				container.streams[p].data.size() << "\t" << container.streams[p].strides.size() << std::endl;

		total_new_comp += container.streams[p].data.size() + container.streams[p].strides.size();
		container.streams[p].reset();
	}
	stats.names_c += total_new_comp;
	//std::cerr << "total for name=" << total_new << " vs " << total_name << "\t" << total_new_comp << "\t" << (float)total_name/total_new_comp << "-fold" << std::endl;

	std::cout << tachyon::utility::timestamp("LOG") << tachyon::utility::ToPrettyDiskString(stats.bases_c) << "\t" << (float)stats.bases/stats.bases_c
			<< "\t" << tachyon::utility::ToPrettyDiskString(stats.quals_c) << "\t" << (float)stats.quals/stats.quals_c << "\t"
			<< "\t" << tachyon::utility::ToPrettyDiskString(stats.names_c) << "\t" << (float)stats.names/stats.names_c << std::endl;

	return true;
}


/****************************
*  Core fastq
****************************/
struct fastq_t {
    fastq_t() : revcmp(false), rid(-1), pos_s(0), pos_e(0){}
    ~fastq_t() = default;

    fastq_t(const fastq_t& other) : revcmp(other.revcmp), rid(other.rid), pos_s(other.pos_s), pos_e(other.pos_e){
        for(int i = 0; i < 4; ++i) fields[i] = other.fields[i];
    }

    fastq_t(fastq_t&& other) : revcmp(other.revcmp), rid(other.rid), pos_s(other.pos_s), pos_e(other.pos_e){
		for(int i = 0; i < 4; ++i) fields[i] = std::move(other.fields[i]);
	}

    fastq_t& operator=(const fastq_t& other){
		revcmp = other.revcmp; rid = other.rid; pos_s = other.pos_s; pos_e = other.pos_e;
		for(int i = 0; i < 4; ++i) fields[i] = other.fields[i];
	}

    fastq_t& operator=(fastq_t&& other) noexcept{
    	revcmp = other.revcmp; rid = other.rid; pos_s = other.pos_s; pos_e = other.pos_e;
        for(int i = 0; i < 4; ++i) fields[i] = std::move(other.fields[i]);
    }

    void reset(void){
        revcmp = false; rid = -1; pos_s = 0; pos_e = 0;
        //for(int i = 0; i < 4; ++i) fields[i].clear();
    }

    bool operator<(const fastq_t& other) const{
    	if(fields[1].size() < other.fields[1].size()) return true;
    	if(other.fields[1].size() < fields[1].size()) return false;

    	if(rid >= 0 && other.rid >= 0){
    		if(pos_s < other.pos_s) return true;
    		if(other.pos_s < pos_s) return false;
    		if(pos_e < other.pos_e) return true;
    		if(other.pos_e < pos_e) return false;
    	}

    	return false;
    }

    bool revcmp;
    int32_t rid;
    int64_t pos_s, pos_e;
	std::string fields[4];
};

bool FastqPrefixSort(const fastq_t& f1, const fastq_t& f2){
	if(f1.fields[1].size() < f2.fields[1].size()) return true;
	if(f2.fields[1].size() < f1.fields[1].size()) return false;

	for(int i = 0; i < 32; ++i){
		if(BASE_MAP[f1.fields[1][i]] < BASE_MAP[f2.fields[1][i]]) return true;
		if(BASE_MAP[f2.fields[1][i]] < BASE_MAP[f1.fields[1][i]]) return false;
	}

	return false;

}

bool FastqPrefixSortQual(const fastq_t& f1, const fastq_t& f2){
	if(f1.fields[1].size() < f2.fields[1].size()) return true;
	if(f2.fields[1].size() < f1.fields[1].size()) return false;

	for(int i = 0; i < 32; ++i){
		if(BASE_MAP[f1.fields[3][i]] < BASE_MAP[f2.fields[3][i]]) return true;
		if(BASE_MAP[f2.fields[3][i]] < BASE_MAP[f1.fields[3][i]]) return false;
	}

	return false;

}

class FastqContainer {
public:
    FastqContainer() : n(0), m(0), rcd(nullptr){}
    ~FastqContainer(){ delete [] rcd; }

    inline void operator+=(const fastq_t& rec){
        if(n == m) this->capacity(std::max(m*2,(size_t)10000));
        rcd[n++] = rec;
    }

    void capacity(const uint32_t new_size){
        if(new_size < m) return;
        fastq_t* temp = rcd;
        rcd = new fastq_t[new_size];
        for(int i = 0; i < n; ++i) rcd[i] = std::move(temp[i]);
        delete [] temp;
        m = new_size;
    }

    void reset(void){
        for(int i = 0; i < n; ++i)
            rcd[i].reset();
        n = 0;
    }

    bool Process(yon1_fq_block_t& block, stats_test& stats, const FastqContainer& fq, const int32_t rid){
    	if(rid >= 0){
			//std::sort(&fq.rcd[0], &fq.rcd[fq.n]);
		} else {

			//std::cerr << "prefix sort" << std::endl;
			//std::sort(&fq.rcd[0], &fq.rcd[fq.n], FastqPrefixSortQual);

			/*
			std::cerr << fq.rcd[0].fields[1] << " " << 0 << "/" << 0 << std::endl;
			for(int z = 1; z < fq.n; ++z){
				// Prefix matching
				int i = 0; int matches = 0;
				for(; i < fq.rcd[z].fields[1].size(); ++i){
					if(fq.rcd[z].fields[1].at(i) != fq.rcd[z-1].fields[1].at(i)) std::cerr << fq.rcd[z].fields[1].at(i);
					else { std::cerr << "-"; ++matches; }
				}
				std::cerr << " " << matches << "/" << (float)matches/fq.rcd[z].fields[1].size() << std::endl;
			}
			*/

		}
		//std::cerr << "done sorted" << std::endl;

    	uint32_t total_name = 0;
    	for(int i = 0; i < n; ++i){
			if(rcd[i].rid == -1){
				//if(i == 0){ // first one
					//std::cerr << "rid-1 first" << std::endl;
					//++n_no_hits;
					//std::cerr << "packing bases" << std::endl;
					uint8_t packed_bases = 0;
					int p = 0;
					for(; p + 4 < rcd[i].fields[1].size(); p += 4){
						packed_bases = ((rcd[i].fields[1].at(p) & 3) << 6) | ((rcd[i].fields[1].at(p+1) & 3) << 4) | ((rcd[i].fields[1].at(p+2) & 3) << 2) | ((rcd[i].fields[1].at(p+3) & 3) << 0);
						//std::cerr << "," << (int)packed_bases;
						block.bases_packed.AddLiteral((uint8_t)packed_bases);
						packed_bases = 0;
					}

						// Add residual
						uint32_t diff = rcd[i].fields[1].size() - p;
						assert(diff < 4);
						for(int k = 0; k < diff; ++k){
							packed_bases |= ((rcd[i].fields[1].at(p+k) & 3) << k*2);
						}
						block.bases_packed.AddLiteral((uint8_t)packed_bases);
						++block.bases_packed;

						// Add stride
						block.bases_packed.AddStride(ceil((float)rcd[i].fields[1].size()/4));
					/*
					} else { // not first
						uint8_t pval = 0; uint8_t poff = 0;
						for(int j = 0; j < fq.rcd[i].fields[1].size(); ++j){
							if(fq.rcd[i].fields[1].at(j) != fq.rcd[i-1].fields[1].at(j)){
								//std::cerr << fq.rcd[i].fields[1].at(j);
								pval = (pval << 4) | (BASE_MAP[fq.rcd[i].fields[1].at(j)]);
								if(++poff == 2){
									block.bases_packed.AddLiteral((uint8_t)pval);
									++block.bases_packed;
									poff = 0; pval = 0;
								}
							}
							else {
								//std::cerr << "-";
								pval = (pval << 4) | 6;
								if(++poff == 2){
									block.bases_packed.AddLiteral((uint8_t)pval);
									++block.bases_packed;
									poff = 0; pval = 0;
								}
							}
						}
						if(poff != 0){
							block.bases_packed.AddLiteral((uint8_t)pval);
							++block.bases_packed;
							poff = 0; pval = 0;
						}
					}
					*/

				} else {
					block.bases_matches_contig.Add(rcd[i].rid);
					block.bases_matches.Add(rcd[i].pos_s);
					++block.bases_matches;
					++block.bases_matches_contig;
				}

				// Quals
				//yon1_fq_block_t::dc* target = &block.quals;
				for(int p = 0; p < rcd[i].fields[3].size(); ++p){
					//uint8_t zig = rcd[i].fields[3].at(p) - 33;
					block.quals.AddLiteral(rcd[i].fields[3].at(p));
				}
				block.quals.header.n_additions += rcd[i].fields[3].size();
				++(block.quals);
				block.quals.AddStride(rcd[i].fields[3].size());

			std::vector<std::string> tokens = tachyon::utility::split(rcd[i].fields[0],' ');
			//@ERR174310.14839
			std::vector<std::string> p1 = tachyon::utility::split(tokens[0],'.');
			//std::cerr.write(&p1[0][4], p1[0].size() - 4);

			// All names copied. Todo: check uniformiy of strings
			block.streams[0].AddString(p1[0]);
			++block.streams[0];
			block.streams[0].AddStride(p1[0].size());

			// Todo: guranteed to be ordered. check arithmetic progression
			block.streams[1].Add(atoi(p1[1].c_str()));
			++block.streams[1];

			// HSQ1008_141:5:1101:11962:4395/1
			std::vector<std::string> p2 = tachyon::utility::split(tokens[1], ':');
			std::string end = p2[4].substr(0, p2[4].size() - 2);

			// Todo: guranteed to be ordered. check arithmetic progression
			block.streams[2].AddString(p2[0]);
			++block.streams[2];
			block.streams[2].AddStride(p2[0].size());

			block.streams[3].Add(atoi(p2[1].c_str()));
			++block.streams[3];
			block.streams[4].Add(atoi(p2[2].c_str()));
			++block.streams[4];
			block.streams[5].Add(atoi(p2[3].c_str()));
			++block.streams[5];
			block.streams[6].Add(atoi(end.c_str()));
			++block.streams[6];
			block.streams[7].Add((int8_t)p2[4].back());
			++block.streams[7];

			total_name += rcd[i].fields[0].size();

			stats.bases += rcd[i].fields[1].size();
			stats.names += rcd[i].fields[0].size();
			stats.quals += rcd[i].fields[3].size();
    	}

    	return true;
    }

    bool ProcessQuals(yon1_fq_block_t& block);
    bool ProcessBases(yon1_fq_block_t& block);
    bool ProcessNames(yon1_fq_block_t& block);

public:
    size_t n, m;
    fastq_t* rcd;
};


/****************************
*  Core reference sequence
****************************/
struct fasta_b {
	static const uint32_t n_bases_unit = 64 / 2;
	static const uint32_t b_off_start  = n_bases_unit - 1;
	static const uint64_t lookup_mask  = 3L << 2*(n_bases_unit-1);

	fasta_b(void) : b_off(b_off_start), n_d(0), m_d(0), n_b(0), d(nullptr){}
	fasta_b(const uint32_t n) : b_off(b_off_start), n_d(0), m_d(n), n_b(0), d(new uint64_t[n]){ memset(d, 0, m_d*sizeof(uint64_t)); }
	~fasta_b(){ delete [] d; }

	fasta_b& Clone(const fasta_b& other){
		delete [] d;
		n_b = other.n_b;
		m_d = other.n_d;
		n_d = other.n_d;
		d = new uint64_t[m_d];
		memcpy(d, other.d, other.n_d*sizeof(uint64_t));
		return(*this);
	}

	void resize(const uint32_t n){ n_d = 0; m_d = n; d = new uint64_t[n]; memset(d, 0, m_d*sizeof(uint64_t)); }

	void operator+=(const char base){
		if(b_off == -1){ ++n_d; b_off = b_off_start; }
		d[n_d] |= ((uint64_t)BASE_MAP[base] & 3) << (2*b_off--);
		++n_b;
	}

	inline char operator[](const uint32_t p) const{
		return(BASE_INV_MAP[(d[p/n_bases_unit] >> ((n_bases_unit - (p%n_bases_unit) - 1)*2)) & 3]);
	}

	void reset(){ b_off = b_off_start; n_d = 0; memset(d, 0, m_d*sizeof(uint64_t)); }

	inline const uint32_t& size() const{
		//return(this->n_d*n_bases_unit + (n_bases_unit - (b_off+1)));
		return n_b;
	}

	uint64_t b_off;
	uint32_t  n_d, m_d, n_b;
	uint64_t* d;
};

// Idea: compute (k,w)-minimizers for every read. Store tuple (read offset, match)
//       for every unique minimizer. Whenever matching, search to every minimizer
//       position and try to extend left and right until no longer perfect matching
//       occurs or some N in-order mismatches have been found.
template <uint32_t n_bases, class int_t = uint64_t>
struct kmer_base {
	const static int_t word_size    = sizeof(int_t);
	const static int_t bit_size     = word_size * 8;
	const static int_t bases_word   = bit_size / KMER_SHIFT;
	const static int_t bits_used    = n_bases % bases_word == 0 ? bit_size : (n_bases % bases_word) * KMER_SHIFT;
	const static int_t offset_start = bit_size - KMER_SHIFT;
	const static int_t trim_mask    = (bit_size == bits_used ? ~(int_t)0 : ((int_t)1 << (bits_used)) - 1);
    const static int_t n_words      = ceil((float)n_bases / bases_word);
    const static int_t upper4_mask  = (int_t)15 << (bits_used - 4);

    // Function pointer definition.
    typedef void (kmer_base::*add_func)(const char& value);

    kmer_base() :
        d_fwd(new int_t[n_words]), d_rev(new int_t[n_words]),
        f_fwd(n_words == 1 ? &kmer_base::AddFwdSingle : &kmer_base::AddFwdMultiple),
        f_rev(n_words == 1 ? &kmer_base::AddRevSingle : &kmer_base::AddRevSingle)
    {
    	memset(d_fwd, 0, word_size*n_words);
    	memset(d_rev, 0, word_size*n_words);
    }
	virtual ~kmer_base(){ delete [] d_fwd; delete [] d_rev; }

	inline void operator+=(const char v){ this->Add(v); }
	inline char AtFwd(const uint32_t p) const{ return(BASE_INV_MAP[KMER_EXTRACT_POS(d_fwd[0],p,bits_used)]); }
	inline char AtFwd(const uint32_t p, const uint32_t word, const uint32_t bits = bit_size) const{ return(BASE_INV_MAP[KMER_EXTRACT_POS(d_fwd[word],p,bits)]); }
	inline char AtRev(const uint32_t p) const{ return(BASE_INV_MAP[KMER_EXTRACT_REV_POS(d_rev[0],p,n_bases)]); }
	inline char AtRev(const uint32_t p, const uint32_t word) const{ return(BASE_INV_MAP[KMER_EXTRACT_REV_POS(d_rev[word],p,n_bases%bases_word)]); }

	inline char AtRevFwd(const uint32_t p) const{ return(BASE_INV_MAP[KMER_EXTRACT_POS(d_rev[0],p,bits_used)]); }

    inline void Add(const char v){ (this->*f_fwd)(v); (this->*f_rev)(v); }

	bool operator<(const kmer_base& other) const{
		// smaller if:
		//    1) forward is smaller than other forward AND smaller than other rev
		//    2) rev is smaller than other rev AND smaller than other fwd
		if(d_fwd[0] < other.d_fwd[0] && d_fwd[0] < other.d_rev[0] && ((d_fwd[0] & upper4_mask) >> (bits_used - 4)) != 0) return true;
		if(d_rev[0] < other.d_rev[0] && d_rev[0] < other.d_fwd[0] && ((d_rev[0] & upper4_mask) >> (bits_used - 4)) != 0) return true;
		return false;
	}

	inline void operator=(const kmer_base& other){ d_fwd[0] = other.d_fwd[0]; d_rev[0] = other.d_rev[0]; }

    std::ostream& Print(std::ostream& stream){
    	//std::cerr << "extra=" << (n_bases%bases_word) << std::endl;
    	const uint32_t n_extra = (n_bases%bases_word) == 0 ? bases_word : (n_bases%bases_word);
    	for(int i = 0; i < n_extra; ++i)
			stream << this->AtFwd(i, n_words-1, bits_used);

    	for(int j = (int)n_words - 2; j >= 0; --j){
			for(int i = 0; i < bases_word; ++i)
				stream << this->AtFwd(i, j);
    	}

    	return(stream);
    }
	//virtual void PrintRevComp(std::ostream& stream);

	inline void reset(){ memset(d_fwd, 0, word_size*n_words); memset(d_rev, 0, word_size*n_words); }

    inline void AddFwdSingle(const char& value){ d_fwd[0] = ((d_fwd[0] << KMER_SHIFT) | (BASE_MAP[value] & KMER_SHIFT_MASK)) & trim_mask; }
	inline void AddRevSingle(const char& value){ d_rev[0] = (d_rev[0] >> KMER_SHIFT) | (((int_t)BASE_MAP_COMP[BASE_MAP[value] & KMER_SHIFT_MASK]) << (bit_size - KMER_SHIFT)); }

	void AddFwdMultiple(const char& value){
		/*
		 *            2     1     0
		 *  start  | A B | A B | A B |
		 *  shift  | B   | A B | A B |
		 *  move   | B A | B A | B   |
		 *  add    | B A | B A | B C |
		 *
		 */

		// Shift value off first
		d_fwd[n_words-1] <<= KMER_SHIFT;

		for(int i = (int)n_words - 2; i >= 0; --i){
			d_fwd[i+1]  |= (d_fwd[i] & ((int_t)3 << (bit_size - 2))) >> (bit_size - 2);
			d_fwd[i]   <<= KMER_SHIFT;
		}

		// Shift in lowest
		d_fwd[0] |= BASE_MAP[value];
		d_fwd[n_words-1] &= trim_mask;
	}

	int_t* d_fwd, *d_rev;
    add_func f_fwd, f_rev;
};

struct hash_support_t {
	uint64_t contig: 6, pos: 58;
};

static uint64_t identity_bins[11];

void BuildMinimizers(tachyon::hash::HashTable<uint64_t, uint64_t>& htable, fasta_b* fasta){
	std::ifstream f; f.open("/home/mk21/Downloads/fa_ref/Homo_sapiens.GRCh38.dna.toplevel.b.fa");
	std::string line;

	kmer_base<64,uint64_t> x; // running kmer
	kmer_base<64,uint64_t> x2; // best kmer

	uint64_t position = 0;
	uint32_t rel_pos  = 0;
	uint64_t abs_pos  = 0;
	uint32_t window   = 0;
	const uint32_t window_size = 100;

	uint64_t n_minimizers_c = 0;
	uint32_t n_contigs = 0;

	uint64_t wins_a = 0, wins_b = 0;
	uint64_t b_data = 0;
	uint64_t bit2_fasta_size = 0;

	fasta_b fasta_buffer(300000000);

	getline(f, line);
	if(f.good() == false){
		std::cerr << "could not get first line" << std::endl;
		return;
	}
	std::cerr << "newline: " << line << std::endl;
	b_data += line.size();

	while(getline(f, line)){
		b_data += line.size();
		if(line[0] == '>'){
			std::cerr << "minimizers=" << n_minimizers_c << " " << wins_a << "," << wins_b << std::endl;
			std::cerr << "occupied=" << htable.occupied() << std::endl;
			std::cerr << "fasta=" << fasta_buffer.size() << std::endl;
			bit2_fasta_size += (fasta_buffer.n_d+1)*sizeof(uint64_t);
			fasta[n_contigs].Clone(fasta_buffer);
			std::cerr << "Cloning into: " << n_contigs << std::endl;
			fasta_buffer.reset();
			position = 0; abs_pos = 0;
			++n_contigs;
			if(n_contigs == 32) break;
			std::cerr << "newline: " << line << std::endl;
			continue;
		}

		std::cerr << line << std::endl;

		//std::cerr << line << std::endl;
		//if(n_contigs == 22)
		//	std::cerr << n_contigs << "\t" << position << std::endl;

		uint32_t start_pos = position;
		for(int i = 0; i < line.size(); ++i, ++position, ++window){
			//std::cerr << "adding: " << line.at(i) << "->" << (int)BASE_MAP[line.at(i)] << std::endl;
			x += line.at(i);
			fasta_buffer += line.at(i);


			if(position >= 63){
				x.Print(std::cerr);
				std::cerr << std::endl;
			}


			if(window == 63){
				rel_pos = window;
				abs_pos = position;
				x2 = x;
			}

			if(x < x2){
				x2 = x;
				rel_pos = window;
				abs_pos = position;
			}

			if(window == window_size){
				++n_minimizers_c;

				//wins_a += x2.data[0] <= x2.data_rev_cmp[0];
				//wins_b += x2.data_rev_cmp[0] < x2.data[0];

				// Hash smallest
				hash_support_t store; store.contig = n_contigs; store.pos = abs_pos - 63;
				uint64_t* store_int = reinterpret_cast<uint64_t*>(&store);

				if(x2.d_fwd[0] <= x2.d_rev[0]){
					htable.SetItem(x2.d_fwd, *store_int, sizeof(uint64_t));
					/*
					std::cerr << rel_pos << " " << abs_pos-32 << " hit= ";
					x2.Print(std::cerr) << " ";
					for(int p = store.pos + 1; p < store.pos + 32; ++p){
						std::cerr << fasta_buffer[p];
					}
					std::cerr << std::endl;
					*/
					//if(x2.data[0] != 0) exit(1);
				}
				else htable.SetItem(x2.d_rev, *store_int, sizeof(uint64_t));

				window = 0; rel_pos = 0;
				x2.reset();
			}
		}

		// Debug
		/*
		std::cerr << line << " ";
		for(int p = start_pos; p < position; ++p){
			std::cerr << fasta_buffer[p];
		}
		std::cerr << std::endl;
		*/

		// Assert
		/*
		if(line[0] != 'N'){
			if(line[0] != fasta_buffer[position-line.size()]){
				std::cerr << "fail=" << line[0] << "!=" << fasta_buffer[position-line.size()] << std::endl;
			}
		}
		*/
	}
	fasta[n_contigs].Clone(fasta_buffer);
	std::cerr << "minimizers=" << n_minimizers_c << std::endl;
	std::cerr << "data=" << tachyon::utility::ToPrettyDiskString(b_data) << std::endl;
	std::cerr << "fasta2bit=" << bit2_fasta_size << "b -> " << tachyon::utility::ToPrettyDiskString(bit2_fasta_size) << std::endl;

	uint64_t actual_bytes = 0;
	for(int i = 0; i < 32; ++i){
		actual_bytes += fasta[i].n_d*sizeof(uint64_t);
	}
	std::cerr << "actual=" << actual_bytes << std::endl;

	return;
}

bool ExtendRevComp(fastq_t& fq,
                   uint32_t p,
                   const hash_support_t& hit,
                   kmer_base<32,uint64_t>& kmer,
                   const fasta_b* fasta)
{
	uint64_t right_end  = hit.pos + 32;
	uint64_t left_start = hit.pos;

	// Check matching
	int32_t b_matches = 0;
	for(int k = left_start, g = 0; k < right_end; ++k, ++g){
		b_matches += (kmer.AtRev(g) == fasta[hit.contig][k]);
	}
	assert(b_matches == 32);

	/*
	std::cerr << "rev-comp= ";
	kmer.Print(std::cerr);
	std::cerr << " ";

	for(int k = left_start; k < right_end; ++k){
		std::cerr << fasta[hit.contig][k];
	}
	std::cerr << std::endl;
	*/

	if(hit.pos <= fq.fields[1].size() || hit.pos + fq.fields[1].size() >= fasta[hit.contig].n_b){
		//std::cerr << "oobunds return" << std::endl;
		fq.rid = -1;
		return false;
	}

	uint32_t end_left = 0;
	uint8_t  skipped_left = 0;
	for(int k = p + 1; k < fq.fields[1].size(); ++k){
		if(BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][left_start - 1 - end_left]]]] == fq.fields[1].at(k)){
			//std::cerr << "match left: " << k << "->" << BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][left_start - 1 - end_left]]]] << "==" << bases.at(k) << std::endl;

			++end_left;
		} else {
			//std::cerr << "left: no match: " << k << " this=" << bases.at(k) << " ref=" << BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][left_start - 1 - end_left]]]] << std::endl;
			if(skipped_left++ == 0) ++end_left;
			else break;
		}
	}
	++end_left; // make left end non-inclusive [from, to)

	//
	uint32_t end_right = 0;
	uint8_t  skipped_right = 0;
	bool right_anchored = true;
	if(p-31 != 0){
		right_anchored = false;
		for(int k = p-32; k != 0; --k){
			if(BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][right_end + end_right]]]] == fq.fields[1].at(k)){
				//std::cerr << "match right: " << k << "->" << BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][right_end + end_right]]]] << "==" << bases.at(k) << std::endl;

				++end_right;
			} else {
				//std::cerr << "right: no match: " << (p-31) - k << " this=" << bases.at(k) << " ref=" << BASE_INV_MAP[BASE_MAP_COMP[BASE_MAP[fasta[hit.contig][right_end + end_right]]]] << std::endl;
				if(skipped_right++ == 0) ++end_right;
				else break;
			}
		}
	}

	//std::cerr << "range=" << (p-32)-end_right << "->" << p + end_left << std::endl;

	// Add one to end because inclusive
	const uint32_t range = ((p + end_left)) - ((p-(32-right_anchored)) - end_right);
	const float identity = (float)range/fq.fields[1].size();
	//std::cerr << "revcmp identity=" << identity << " skip-left=" << (int)skipped_left << " skipped-right=" << (int)skipped_right << std::endl;

	++identity_bins[(uint32_t)(identity*10)];

	if(identity < 0.9){
		//std::cerr  << bases << std::endl;
		fq.rid = -1;
		return false;
	}

	fq.revcmp = true;
	// Relative offsets
	fq.pos_e  = p + end_left;
	fq.pos_s  = (p-(32-right_anchored)) - end_right;

	/*
	std::cerr << "raw= l=" << end_left << " r=" << end_right << " p=" << p << " -> " << p+end_left << "," << (((int32_t)p-(32-right_anchored)) - end_right) << std::endl;
	std::cerr << "rev=id" << identity << "@" << fq.fields[1].size() << " rid=" << fq.rid << " rel_pos=" << fq.pos_s << "-" << fq.pos_e << " " << hit.pos << "@" << p << "->" << hit.pos+fq.pos_s << "-" << hit.pos+fq.pos_e << " skips=" << (int)skipped_left << "," << (int)skipped_right << " anchor=" << right_anchored << std::endl;
	std::cerr << "left=" << (p - fq.pos_s) << " right=" << fq.pos_e - p << std::endl;
	std::cerr << "l=" << hit.pos - ((p-right_anchored) - fq.pos_s) << " r=" << hit.pos + (fq.pos_e - p) << std::endl;
	*/

	// Absolute FASTA offsets
	fq.pos_s = hit.pos - ((p-right_anchored) - fq.pos_s);
	fq.pos_e = hit.pos + (fq.pos_e - p);

	assert(identity <= 1 && identity >= 0);

	return true;
}

bool Extend(fastq_t& fq,
	        uint32_t p,
			const hash_support_t& hit,
			kmer_base<32,uint64_t>& kmer,
			const fasta_b* fasta)
{
	uint64_t right_end  = hit.pos + 32;
	uint64_t left_start = hit.pos;

	// Check matching
	int32_t b_matches = 0;
	for(int k = left_start, g = 0; k < right_end; ++k, ++g){
		b_matches += (kmer.AtFwd(g) == fasta[hit.contig][k]);
	}

	if(b_matches != 32){
		return ExtendRevComp(fq,p,hit,kmer,fasta);
		//return false;
	}

	if(hit.pos <= fq.fields[1].size() || hit.pos + fq.fields[1].size() >= fasta[hit.contig].n_b){
		//std::cerr << "oobunds return" << std::endl;
		fq.rid = -1;
		return false;
	}

	/*
	for(int k = left_start; k < right_end; ++k){
		std::cerr << fasta[hit->contig][k];
	}
	std::cerr << std::endl;
	*/

	uint32_t end_right = 0;
	uint8_t  skipped_right = 0;
	//std::cerr << "expand right: " << std::endl;
	for(int k = p + 1; k < fq.fields[1].size(); ++k){
		if(fasta[hit.contig][right_end + end_right] == fq.fields[1].at(k)){
			//std::cerr << "match right: " << k << "->" << fasta[hit->contig][right_end + end_right]<< "==" << fq.fields[1].at(k) << std::endl;
			++end_right;
		} else {
			if(skipped_right++ == 0) ++end_right;
			//std::cerr << "right: no match: " << k - p << " this=" << fq.fields[1].at(k) << " ref=" << fasta[hit->contig][right_end + end_right] << std::endl;
			else break;
		}
	}
	++end_right; // make right end non-inclusive [from, to)

	uint32_t end_left = 0;
	uint8_t skipped_left = 0;
	bool left_anchored = true;
	if(p-31 != 0){
		left_anchored = false;
		for(int k = p-32; k != 0; --k){
			if(fasta[hit.contig][left_start - 1 - end_left] == fq.fields[1].at(k)){
				//std::cerr << "match left: " << k << "->" << fasta[hit->contig][left_start - 1 - end_left]<< "==" << fq.fields[1].at(k) << std::endl;

				++end_left;
			} else {
				//std::cerr << "left: no match: " << (p-31) - k << " this=" << fq.fields[1].at(k) << " ref=" << fasta[hit->contig][left_start - 1 - end_left] << std::endl;
				if(skipped_left++ == 0) ++end_left;
				else break;
			}
		}
	}

	//std::cerr << "range=" << (p-32)-end_left << "->" << p + end_right << std::endl;

	// Add one to end because inclusive
	const uint32_t range = (p + end_right) - ((p-(32-left_anchored))-end_left);
	const float identity = (float)range/fq.fields[1].size();
	//std::cerr << "fwd identity=" << identity << " skip-left=" << (int)skipped_left << " skipped-right=" << (int)skipped_right << std::endl;

	++identity_bins[(uint32_t)(identity*10)];

	if(identity < 0.9){
		//std::cerr  << bases << std::endl;
		fq.rid = -1;
		return false;
	}

	fq.revcmp = false;
	// Relative offsets
	fq.pos_s  = (p-(32-left_anchored))-end_left;
	fq.pos_e  = p + end_right;

	//std::cerr << "fwd=id" << identity << "@" << fq.fields[1].size() << " " << fq.rid << " " << fq.pos_s << "-" << fq.pos_e << " " << hit.pos << "@" << p << "->" << hit.pos+fq.pos_s << "-" << hit.pos+fq.pos_e << " skips=" << (int)skipped_left << "," << (int)skipped_right << " anchor=" << left_anchored << std::endl;
	//std::cerr << "left=" << (p - fq.pos_s) << " right=" << fq.pos_e - p << std::endl;
	//std::cerr << "l=" << hit.pos - ((p-left_anchored) - fq.pos_s) << " r=" << hit.pos + (fq.pos_e - p) << std::endl;

	// Absolute FASTA offsets
	fq.pos_s = hit.pos - ((p-left_anchored) - fq.pos_s);
	fq.pos_e = hit.pos + (fq.pos_e - p);

	assert(identity <= 1 && identity >= 0);
	return true;
}

bool GetBlock(std::istream& stream, const uint32_t n_bytes, tachyon::yon_buffer_t& buffer){
	assert(buffer.capacity() >= n_bytes);

	uint64_t prev_pos = stream.gcount();
	stream.read(buffer.data(), n_bytes); // read first n_bytes
	// search from tail to next new line
	uint64_t cur_pos = (uint64_t)stream.gcount() - prev_pos;
	buffer.n_chars_ = cur_pos;
	std::cerr << "done reading: " << buffer.size() << std::endl;

	uint32_t pos[5]; uint32_t n_pos = 0;
	for(int i = buffer.size() - 1; i != 0; --i){
		if(buffer[i] == '\n'){
			pos[n_pos] = i;
			if(++n_pos == 5) break;
		}
	}
	std::cerr << "done get pos" << std::endl;

	for(int i = 4; i != 0; --i){
		std::cerr << "pos " << pos[i] << "->" << pos[i-1] << ": " << pos[i-1] - pos[i] << std::endl;
		std::cerr << std::string(&buffer[pos[i] + 1], pos[i-1] - (pos[i] + 1)) << std::endl;
	}

	return true;
}

int fastq(int argc, char** argv){
	tachyon::hash::HashTable<uint64_t, uint64_t> htable(100000000);
	fasta_b* fasta = new fasta_b[32];
	BuildMinimizers(htable, fasta);

	std::string line;
	fastq_t fq;
	FastqContainer fqc[33]; // 33 is unaligned
	for(int i = 0; i < 33; ++i)
		fqc[i].capacity(BLOCK_SIZE + 128);

	uint32_t i = 0;
	yon1_fq_block_t container;
	//yon1_fq_block_t blocks[32];

	uint32_t j = 0;
	uint32_t total_name = 0;

	tachyon::algorithm::CompressionManager cm;
	cm.zstd_codec.SetCompressionLevel(3);

	stats_test stats;

	fqz fqz_ins;
	RangeCoder rc;

	kmer_base<32,uint64_t> kmer;

	uint64_t n_hits = 0, n_no_hits = 0;
	uint64_t hits_fwd = 0, hits_rev = 0;

	tachyon::algorithm::Timer t; t.Start();

	//tachyon::yon_buffer_t buff(12000000);
	//GetBlock(std::cin,10000000,buff);

	//delete [] fasta;
	//return 0;
	memset(identity_bins, 0, sizeof(uint64_t)*11);
	//bool qual_seen[128];

	while(getline(std::cin, line)){
		fq.fields[i%4] = line;
		if(i%4 == 0 && i != 0){
			//memset(qual_seen, 0, 128);

			uint64_t* ret_target = nullptr;

			// Iterate over seuqence bases.
			//std::cerr << fq.fields[1] << std::endl;

			hash_support_t* kmer_hit = nullptr;
			for(uint32_t p = 0; p < fq.fields[1].size(); ++p){
				kmer += fq.fields[1].at(p);

				if(p >= 31){
					if(kmer.d_fwd[0] <= kmer.d_rev[0]){
						htable.GetItem(kmer.d_fwd, ret_target, sizeof(uint64_t));
					}
					else htable.GetItem(kmer.d_rev, ret_target, sizeof(uint64_t));

					if(ret_target != nullptr){
						++n_hits;
						kmer_hit = reinterpret_cast<hash_support_t*>(ret_target);
						assert(kmer_hit->contig < 24);
						if(kmer_hit->pos > fasta[kmer_hit->contig].n_b){
							std::cerr << "out of bounds: " << kmer_hit->contig << ":" << kmer_hit->pos << "/" << fasta[kmer_hit->contig].n_b << std::endl;
							exit(1);
						}

						fq.rid = kmer_hit->contig;
						Extend(fq, p, *kmer_hit, kmer, fasta);

						// if extension did not work out well we keep trying
						// otherwise we break
						if(fq.rid >= 0) break;
						else {
							//std::cerr << "continue; not worked out well: " << p << std::endl;
							ret_target = nullptr;
						}
					}
				}
			}
			kmer.reset(); //kmer_best.clear();

			// Set target fq container
			FastqContainer* target_fq = &fqc[0];
			//if(fq.rid >= 0) target_fq = &fqc[fq.rid];
			//else ++identity_bins[0];

			if(target_fq->n == BLOCK_SIZE){
				target_fq->Process(container, stats, *target_fq, 0);
				CompressBlock(container,fqz_ins,rc,cm,stats,0);

				/*
				if(target_fq->rcd[0].rid >= 0){
					std::cerr << target_fq->rcd[0].rid << ":" << target_fq->rcd[0].pos_s << "-" << target_fq->rcd[0].pos_e << " " << target_fq->rcd[0].revcmp << ":" << target_fq->rcd[0].fields[1] << std::endl;
					for(int z = 1; z < target_fq->n; ++z){
						std::cerr << target_fq->rcd[z].rid << ":" << target_fq->rcd[z].pos_s << "-" << target_fq->rcd[z].pos_e << " " << target_fq->rcd[z].revcmp << ":" << target_fq->rcd[z].fields[1] << std::endl;
					}
				}
				*/
				target_fq->reset();
			}
			*target_fq += fq;
			fq.reset();
		}


		++i; ++j;
	}

	if(fqc[0].n){
		fqc[0].Process(container, stats, fqc[0], 0);
		CompressBlock(container,fqz_ins,rc,cm,stats,0);
	}

	/*
	for(int i = 0; i < 33; ++i){
		if(fqc[i].n == 0)
			continue;

		fqc[i].Process(container, stats, fqc[i], i);
		CompressBlock(container,fqz_ins,rc,cm,stats,i);
	}
	*/

	std::cout << "Final cost=" << tachyon::utility::ToPrettyDiskString(stats.names_c+stats.bases_c+stats.quals_c) << std::endl;
	std::cout << "names\t" << stats.names << "\t" << stats.names_c << "\t" << (float)stats.names_c/stats.names << std::endl;
	std::cout << "bases\t" << stats.bases << "\t" << stats.bases_c << "\t" << (float)stats.bases_c/stats.bases << std::endl;
	std::cout << "qual\t" << stats.quals << "\t" << stats.quals_c << "\t" << (float)stats.quals_c/stats.quals << std::endl;

	delete [] fasta;

	for(int i = 0; i < 11; ++i)
		std::cerr << i << "\t" << identity_bins[i] << std::endl;

	std::cerr << tachyon::utility::timestamp("DONE") << t.ElapsedString() << std::endl;
	return 0;
}

#endif /* STATS_H_ */
