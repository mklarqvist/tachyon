#ifndef IO_FORWARDINGSTREAM_H_
#define IO_FORWARDINGSTREAM_H_

namespace Tachyon{
namespace IO{

#include <iostream>
#include <streambuf>

template<typename CharType, typename Traits = std::char_traits<CharType> >
class ForwardingStreamBuf : public std::basic_streambuf<CharType, Traits>
{
public:
    typedef Traits traits_type;
    typedef typename traits_type::int_type int_type;
    typedef typename traits_type::pos_type pos_type;
    typedef typename traits_type::off_type off_type;

    ForwardingStreamBuf(std::basic_streambuf<CharType, Traits> *baseStreamBuf)
        : _baseStreamBuf(baseStreamBuf)
    {
    }

protected:
    virtual int_type overflow(int_type c = traits_type::eof()){
        if(this->_baseStreamBuf == NULL)
            return traits_type::eof();

        if(traits_type::eq_int_type(c, traits_type::eof()))
        	return traits_type::not_eof(c);
        else {
        	// DO SOMETHING HERE
            CharType ch = traits_type::to_char_type(c);
            if( ch >= 'A' && ch <= 'z' )
                ch++; // Do some meaningless transformation
            return _baseStreamBuf->sputc(ch);
        }
    }

    virtual int sync(){
        if( _baseStreamBuf == NULL )
            return -1;
        else
            return _baseStreamBuf->pubsync();
    }
private:
    std::basic_streambuf<CharType, Traits>* _baseStreamBuf;
};

template<typename CharType, typename Traits = std::char_traits<CharType> >
class ForwardingStream : public std::basic_ostream<CharType, Traits>{
public:
    ForwardingStream(std::basic_ostream<CharType, Traits> &stream)
        : std::basic_ostream<CharType, Traits>(NULL), _buffer(stream.rdbuf())
    {
        this->init(&_buffer);
    }

    ForwardingStreamBuf<CharType, Traits>* rdbuf() const{ return &_buffer; }
private:
    ForwardingStreamBuf<CharType, Traits> _buffer;
};

}

}



#endif /* IO_FORWARDINGSTREAM_H_ */
