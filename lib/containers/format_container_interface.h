#ifndef CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_
#define CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_

namespace tachyon{
namespace containers{

class FormatContainerInterface{
private:
    typedef FormatContainerInterface self_type;
    typedef std::size_t              size_type;

public:
    FormatContainerInterface() : primitive_type( YON_TYPE_32B), n_entries(0), n_capacity(0){}
    FormatContainerInterface(const size_t n_entries) : primitive_type(YON_TYPE_32B), n_entries(n_entries), n_capacity(0){}
    virtual ~FormatContainerInterface(){}

	/**<
	 * Convert a FormatContainer into a DataContainer. This is primarily
	 * done for writing out a Tachyon archive.
	 * @return Returns a DataContainer with the contextual representation of the data into this container.
	 */
	virtual DataContainer ToDataContainer(void) =0;

	/**<
	 * Add data from a FormatContainer into an already existing DataContainer.
	 * @param container Destination DataContainer.
	 * @return          Returns a reference to the input DataContainer with the contextual representation of the data in this container added to it.
	 */
	virtual DataContainer& UpdateDataContainer(DataContainer& container) =0;

    // Capacity
	inline bool empty(void) const{ return(this->n_entries == 0); }
	inline const size_type& size(void) const{ return(this->n_entries); }

    virtual std::ostream& ToVcfString(std::ostream& stream, const uint32_t position, const uint64_t sample_number) const =0;
    virtual io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint32_t position, const uint64_t sample) const =0;
    virtual io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const uint32_t position, const uint64_t sample) const =0;
    virtual bool emptyPosition(const uint32_t& position) const =0;
    virtual bool emptyPosition(const uint32_t& position, const uint64_t& sample) const =0;

    virtual bcf1_t* UpdateHtslibVcfRecord(const uint32_t position,
                                          bcf1_t* rec,
	                                      bcf_hdr_t* hdr,
	                                      const std::string& tag) const =0;

protected:
    TACHYON_CORE_TYPE primitive_type;
	size_t  n_entries;
	size_t  n_capacity;
};

}
}



#endif /* CONTAINERS_FORMAT_CONTAINER_INTERFACE_H_ */
