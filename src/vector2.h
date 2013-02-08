#ifndef pta_vector2_h
#define pta_vector2_h
#include <vector>

template<typename T>
class vector2 {
    private:
        int rows;
        int columns;
        std::vector<T> data;

    public:
        void assign(int rows_, int columns_, T value) {
            rows = rows_;
            columns = columns_;
            data.assign(rows * columns, value);
        }

        T operator()(int i, int j) const {
            return data[i * columns + j];
        }

        T& operator()(int i, int j) {
            return data[i * columns + j];
        }
};

#endif
