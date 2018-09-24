/*
Copyright (C) 2017-current Genome Research Ltd.
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
#ifndef TACHYON_VARIANT_WRITER_H_
#define TACHYON_VARIANT_WRITER_H_

#include <cassert>
#include <fstream>
#include <thread>

#include "index.h"
#include "tachyon.h"
#include "data_container.h"
#include "header_footer.h"
#include "variant_container.h"

namespace tachyon {

struct VariantWriteSettings {
public:
	VariantWriteSettings();
	~VariantWriteSettings() = default;

	inline void SetThreads(const int32_t n_threads){ this->n_threads = n_threads; }
	inline void SetPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void SetCompressionLevel(const int32_t compression_level){ this->compression_level = compression_level; }
	inline void SetCheckpointBases(const int32_t bases){ this->checkpoint_bases = bases; }
	inline void SetCheckpointVariants(const int32_t variants){ this->checkpoint_n_snps = variants; }
	inline void SetPermuteGenotypes(const bool yes = true){ this->permute_genotypes = yes; }

public:
	bool permute_genotypes; // permute GT flag
	bool encrypt_data; // encryption flag
	int32_t checkpoint_n_snps; // number of variants until checkpointing
	int32_t checkpoint_bases; // number of bases until checkpointing
	int32_t n_threads; // number of parallel importer threads
	int32_t compression_level; // compression level sent to ZSTD

	std::string output_prefix; // output file prefix
};

class VariantWriterInterface {
public:
	typedef VariantWriterInterface  self_type;
	typedef VariantWriteSettings    settings_type;
	typedef yon1_dc_t               container_type;

public:
	VariantWriterInterface();
	virtual ~VariantWriterInterface();

	void WriteIndex(void);
	virtual bool open(const std::string output) =0;

	/**<
	 * Write the mandatory Tachyon archive header information. The user needs
	 * to update the header with provenance tracking inforamtion priro to writing
	 * it out.
	 * @param header Src global header to be written.
	 * @return       Return TRUE upon success or FALSE otherwise.
	 */
	bool WriteFileHeader(yon_vnt_hdr_t& header);

	yon_vnt_hdr_t& UpdateHeaderImport(yon_vnt_hdr_t& header, const std::string command_string);
	yon_vnt_hdr_t& UpdateHeaderView(yon_vnt_hdr_t& header, const std::string command_string);

	bool Write(yon1_vb_t& container, yon1_idx_rec& index_entry);
	bool UpdateIndex(yon1_idx_rec& index_entry);

	/**<
	 * Wrapper function for writing a valid VariantBlock to a
	 * Tachyon archive. Writes the block itself then the block footer
	 * and end-of-block marker.
	 * @param block       Source VariantBlock for writing.
	 * @param index_entry Source VariantIndexEntry to update the variant index with.
	 * @return            Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteBlock(yon1_vb_t& block, yon1_idx_rec& index_entry);

	void operator+=(const yon1_vnt_t& rec);
	void operator+=(const yon1_vc_t& container);

	bool close(void);

private:
	/**<
	 * Supportive function to inject the essential Tacyon variant block footer into
	 * a byte stream.
	 * @param footer Src data container housing the footer/header to be written.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteBlockFooter(const container_type& footer);

	/**<
	 * Supportive function to inject the Tachyon end-of-block marker following the
	 * emission of a block.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteEndOfBlock(void);

	/**<
	 * Finalize the writing of a Tachyon archive. Writes the
	 * footer, the index, and the checksums.
	 * @param checksums Source checksum container for archive integrity checks.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool WriteFinal();

public:
	uint64_t n_blocks_written;
	uint64_t n_variants_written;
	std::ostream* stream;
	yon_index_t index;

	// Local variant container to add to.
	uint32_t n_s;
	settings_type settings;
	yon1_vc_t     variant_container;

	class VariantWriterInterfaceImpl;
	std::unique_ptr<VariantWriterInterfaceImpl> mImpl;
};

class VariantWriterFile : public VariantWriterInterface{
public:
	typedef VariantWriterFile self_type;

public:
	VariantWriterFile();
	~VariantWriterFile();
	bool open(const std::string output);

private:
	void CheckOutputNames(const std::string& input);

public:
	// Stream information
	std::string filename;
	std::string basePath;
	std::string baseName;
};

class VariantWriterStream : public VariantWriterInterface{
public:
	typedef VariantWriterStream self_type;

public:
	VariantWriterStream();
	~VariantWriterStream();
	bool open(const std::string output){ return true; }

};

}

#endif /* IO_VARIANT_IMPORTWRITER_H_ */
