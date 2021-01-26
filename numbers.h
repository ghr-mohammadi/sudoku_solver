#ifndef NUMBERS_H
#define NUMBERS_H

#define _N_ 9

class Numbers
{
public:
    Numbers();
    ~Numbers();

    int _num = 0;
    int * _num_array = nullptr;

    void _Create();
};

#endif //NUMBERS_H