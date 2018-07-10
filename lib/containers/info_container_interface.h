#ifndef CONTAINERS_INFO_CONTAINER_INTERFACE_H_
#define CONTAINERS_INFO_CONTAINER_INTERFACE_H_

namespace tachyon{
namespace containers{

class InfoContainerInterface{
private:
    typedef InfoContainerInterface self_type;
    typedef std::size_t            size_type;

public:
    InfoContainerInterface() : primitive_type(YON_TYPE_32B), n_entries(0), n_capacity(0){}
    InfoContainerInterface(const size_t n_entries) : primitive_type(YON_TYPE_32B), n_entries(n_entries), n_capacity(n_entries){}
    virtual ~InfoContainerInterface(){}

    // Capacity
	inline const bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }
	inline const size_type& capacity(void) const{ return(this->n_capacity); }

    virtual std::ostream& to_vcf_string(std::ostream& stream, const U32 position) const =0;
    virtual io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer, const U32 position) const =0;
    virtual io::BasicBuffer& to_json_string(io::BasicBuffer& buffer, const U32 position) const =0;
    virtual const bool emptyPosition(const U32& position) const =0;

protected:
    TACHYON_CORE_TYPE primitive_type;
	size_t  n_entries;
	size_t  n_capacity;
};

}
}



#endif /* CONTAINERS_INFO_CONTAINER_INTERFACE_H_ */
