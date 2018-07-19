#ifndef CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_
#define CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_

namespace tachyon{
namespace containers{

class FormatContainerInterface{
private:
    typedef FormatContainerInterface self_type;
    typedef std::size_t              size_type;

public:
    FormatContainerInterface() : primitive_type( YON_TYPE_32B), n_entries(0){}
    FormatContainerInterface(const size_t n_entries) : primitive_type(YON_TYPE_32B), n_entries(n_entries){}
    virtual ~FormatContainerInterface(){}

    // Capacity
	inline bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

    virtual std::ostream& to_vcf_string(std::ostream& stream, const U32 position, const U64 sample_number) const =0;
    virtual io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const U32 position, const U64 sample) const =0;
    virtual io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const U32 position, const U64 sample) const =0;
    virtual bool emptyPosition(const U32& position) const =0;
    virtual bool emptyPosition(const U32& position, const U64& sample) const =0;

protected:
    TACHYON_CORE_TYPE primitive_type;
	size_t  n_entries;
};

}
}



#endif /* CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_ */
