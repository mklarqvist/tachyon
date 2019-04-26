#include <strings.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <openssl/evp.h>

#include "variant_writer.h"
#include "utility.h"
#include "tachyon.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "algorithm/compression/compression_manager.h"

namespace tachyon {

class VariantWriterInterface::VariantWriterInterfaceImpl {
public:
	VariantWriterInterfaceImpl() = default;
	~VariantWriterInterfaceImpl() = default;

public:
	algorithm::VariantDigestManager digest;
};

VariantWriterInterface::VariantWriterInterface() :
	n_blocks_written(0),
	n_variants_written(0),
	stream(nullptr),
	n_s(0),
	mImpl(new VariantWriterInterfaceImpl)
{}

VariantWriterInterface::~VariantWriterInterface()
{

}

bool VariantWriterInterface::WriteFileHeader(yon_vnt_hdr_t& header) {
	if (stream == nullptr)
		return false;

	// Write basic header prefix.
	stream->write(&TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH);

	// Pack header into a byte-stream, compress it, and write
	// it out.
	yon_buffer_t temp(500000);
	yon_buffer_t temp_cmp(temp);
	temp << header;
	algorithm::CompressionManager manager;
	manager.zstd_codec.Compress(temp, temp_cmp, 20);
	uint32_t l_data   = temp.size();
	uint32_t l_c_data = temp_cmp.size();
	utility::SerializePrimitive(l_data,   *stream);
	utility::SerializePrimitive(l_c_data, *stream);
	stream->write(temp_cmp.data(), l_c_data);
	return(stream->good());
}

yon_vnt_hdr_t& VariantWriterInterface::UpdateHeaderImport(yon_vnt_hdr_t& header, const std::string command_string) {
	VcfExtra e;
	e.key = "tachyon_importVersion";
	e.value = tachyon::TACHYON_PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  tachyon::TACHYON_PROGRAM_NAME + '-' + tachyon::TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importCommand";
	e.value = tachyon::LITERAL_COMMAND_LINE;
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_importSettings";
	e.value = command_string;
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);

	return(header);
}

yon_vnt_hdr_t& VariantWriterInterface::UpdateHeaderView(yon_vnt_hdr_t& header, const std::string command_string) {
	VcfExtra e;
	e.key = "tachyon_viewVersion";
	e.value = tachyon::TACHYON_PROGRAM_NAME + "-" + VERSION + ";";
	e.value += "libraries=" +  tachyon::TACHYON_PROGRAM_NAME + '-' + tachyon::TACHYON_LIB_VERSION + ","
			+   SSLeay_version(SSLEAY_VERSION) + ","
			+  "ZSTD-" + ZSTD_versionString()
			+  "; timestamp=" + tachyon::utility::datetime();
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);
	e.key = "tachyon_viewCommand";
	e.value = tachyon::LITERAL_COMMAND_LINE;
	header.literals_ += "##" + e.key + "=" + e.value + '\n';
	header.extra_fields_.push_back(e);

	return(header);
}


bool VariantWriterInterface::Write(yon1_vb_t& container, yon1_idx_rec& index_entry)
{
	this->WriteBlock(container, index_entry); // write block
	this->UpdateIndex(index_entry); // Update index.
	return(this->stream->good());
}

bool VariantWriterInterface::UpdateIndex(yon1_idx_rec& index_entry) {
	assert(this->stream != nullptr);
	index_entry.block_id        = this->n_blocks_written;
	index_entry.byte_offset_end = this->stream->tellp();
	this->index                += index_entry;

	index_entry.Print(std::cerr);
	std::cerr << std::endl;

	++this->n_blocks_written;
	this->n_variants_written += index_entry.n_variants;
	index_entry.reset();

	return true;
}

bool VariantWriterInterface::WriteBlock(yon1_vb_t& block, yon1_idx_rec& index_entry) {
	index_entry.byte_offset = this->stream->tellp();
	block.write(*this->stream);
	this->WriteBlockFooter(block.footer_support);
	this->WriteEndOfBlock(); // End-of-block marker.

	return(this->stream->good());
}

bool VariantWriterInterface::close() {
	if (this->stream == nullptr) return false;
	if (this->stream->good() == false) return false;

	this->WriteFinal();
	this->stream->flush();

	return true;
}

bool VariantWriterInterface::WriteFinal() {
	// Write last if any.
	if (variant_container.size()) {
		this->variant_container.PrepareWritableBlock(this->n_s, this->settings.compression_level);
		this->index.IndexContainer(this->variant_container, this->n_blocks_written);
		this->Write(variant_container.block_, this->index.GetCurrent());
		this->variant_container.clear();
	}

	// Done importing
	this->stream->flush();

	// Write global footer.
	yon_ftr_t footer;
	footer.offset_end_of_data = this->stream->tellp();
	footer.n_blocks           = this->n_blocks_written;
	footer.n_variants         = this->n_variants_written;
	assert(footer.n_blocks == this->index.GetLinearSize());

	uint64_t last_pos = this->stream->tellp();
	this->WriteIndex(); // Write index.
	std::cerr << utility::timestamp("PROGRESS") << "Index size: " << utility::ToPrettyDiskString((uint64_t)this->stream->tellp() - last_pos) << "..." << std::endl;
	last_pos = this->stream->tellp();
	this->mImpl->digest.finalize();       // Finalize SHA-512 digests.
	*this->stream << this->mImpl->digest;
	std::cerr << utility::timestamp("PROGRESS") << "Checksum size: " << utility::ToPrettyDiskString((uint64_t)this->stream->tellp() - last_pos) << "..." << std::endl;
	last_pos = this->stream->tellp();
	*this->stream << footer; // Write global footer and EOF marker.
	std::cerr << utility::timestamp("PROGRESS") << "Footer size: " << utility::ToPrettyDiskString((uint64_t)this->stream->tellp() - last_pos) << "..." << std::endl;

	this->stream->flush();
	return(this->stream->good());
}

void VariantWriterInterface::WriteIndex(void) {
	*this->stream << this->index;
	this->stream->flush();
}

bool VariantWriterInterface::WriteBlockFooter(const container_type& footer) {
	if (this->stream == nullptr) return false;
	const uint64_t start_footer_pos = this->stream->tellp();
	utility::SerializePrimitive(footer.header.data_header.uLength, *this->stream);
	utility::SerializePrimitive(footer.header.data_header.cLength, *this->stream);
	this->stream->write(reinterpret_cast<const char*>(&footer.header.data_header.crc[0]), MD5_DIGEST_LENGTH);
	*this->stream << footer.data;
	return(this->stream->good());
}

bool VariantWriterInterface::WriteEndOfBlock(void) {
	if (this->stream == nullptr) return false;
	utility::SerializePrimitive(TACHYON_BLOCK_EOF, *this->stream);
	return(this->stream->good());
}

void VariantWriterInterface::operator+=(const yon1_vnt_t& rec) {
	// Make sure we have enough space.
	if (variant_container.capacity() < settings.checkpoint_n_snps + 10)
		variant_container.reserve(settings.checkpoint_n_snps + 10);

	if (variant_container.size() == settings.checkpoint_n_snps) {
		// write
		this->variant_container.PrepareWritableBlock(this->n_s, this->settings.compression_level);
		this->index.IndexContainer(this->variant_container, this->n_blocks_written);
		this->Write(variant_container.block_, this->index.GetCurrent());
		this->variant_container.clear();
	}

	if (variant_container.size()) {
		if (rec.rid != variant_container.front().rid) {
			//std::cerr << "different rid: write" << std::endl;
			// write
			this->variant_container.PrepareWritableBlock(this->n_s, this->settings.compression_level);
			this->index.IndexContainer(this->variant_container, this->n_blocks_written);
			this->Write(this->variant_container.block_, this->index.GetCurrent());
			this->variant_container.clear();
		} else if (rec.pos < this->variant_container.front().pos) {
			std::cerr << "unsorted file: " << rec.pos << " < " << variant_container.front().pos << std::endl;
			exit(1);
		} else if (variant_container.back().pos - variant_container.front().pos > settings.checkpoint_bases) {
			//std::cerr << "length cutoff=" << (variant_container.back().pos - variant_container.front().pos) << "/" << settings.checkpoint_bases << std::endl;
			// write
			this->variant_container.PrepareWritableBlock(this->n_s, this->settings.compression_level);
			this->index.IndexContainer(this->variant_container, this->n_blocks_written);
			this->Write(variant_container.block_, this->index.GetCurrent());
			this->variant_container.clear();
		}
	}

	variant_container += rec;
}

bool VariantWriterFile::open(const std::string output) {
	if (output.size() == 0) {
		std::cerr << utility::timestamp("ERROR", "WRITER") << "No output file/file prefix provided!" << std::endl;
		return false;
	}
	std::ofstream* ostream = reinterpret_cast<std::ofstream*>(this->stream);

	this->filename = output;
	this->CheckOutputNames(output);
	ostream->open(this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX, std::ios::out | std::ios::binary);

	// Check streams
	if (!this->stream->good()) {
		std::cerr << utility::timestamp("ERROR", "WRITER") << "Could not open: " << this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX << "!" << std::endl;
		return false;
	}

	if (!SILENT) {
		std::cerr << utility::timestamp("LOG", "WRITER") << "Opening: " << this->basePath + this->baseName + '.' + TACHYON_OUTPUT_SUFFIX << "..." << std::endl;
	}

	return true;
}

void VariantWriterFile::CheckOutputNames(const std::string& input) {
	std::vector<std::string> paths = utility::FilePathBaseExtension(input);
	this->basePath = paths[0];
	if (this->basePath.size() > 0)
		this->basePath += '/';

	if (paths[3].size() == TACHYON_OUTPUT_SUFFIX.size() && strncasecmp(&paths[3][0], &TACHYON_OUTPUT_SUFFIX[0], TACHYON_OUTPUT_SUFFIX.size()) == 0)
		this->baseName = paths[2];
	else this->baseName = paths[1];
}

VariantWriterFile::VariantWriterFile() { this->stream = new std::ofstream; }

VariantWriterFile::~VariantWriterFile() {
	this->stream->flush();
	delete this->stream;
}

VariantWriterStream::VariantWriterStream() { this->stream = &std::cout; }

VariantWriterStream::~VariantWriterStream() {
	this->stream->flush();;
}

} /* namespace Tachyon */
