#include "numbers.h"

Numbers::Numbers() {}

Numbers::~Numbers()
{
    if (this->_num_array != nullptr) {
        delete this->_num_array;
        this->_num_array = nullptr;
    }
}

void Numbers::_Create()
{
    if (this->_num == 0 && this->_num_array == nullptr) {
        this->_num_array = new int[_N_];
        for (int i = 0; i < _N_; i++)
            this->_num_array[i] = (i + 1);
    }
}
