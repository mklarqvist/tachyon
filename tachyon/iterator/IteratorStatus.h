#ifndef ITERATORSTATUS_H_
#define ITERATORSTATUS_H_

namespace tachyon{
namespace iterator{

enum TACHYON_ITERATOR_STATUS {
	// All OK statuses are positive
	YON_IT_START = 1,
	YON_IT_GOOD = 2,
	YON_IT_END = 3,
	// All errors are negative
	YON_IT_ERROR_GENERAL = -1,
	YON_IT_ERROR_UNINITIALIZED = -2,
	YON_IT_ERROR_MALFORMED = -3,
	YON_IT_ERROR_UNKNOWN_TYPE = -4
};

}
}

#endif /* ITERATORSTATUS_H_ */
