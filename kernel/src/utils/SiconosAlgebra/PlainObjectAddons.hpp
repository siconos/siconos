inline Index size(unsigned int index) const
{
    if (index == 0)
    {
        return this->cols();
    }
    else if (index == 1)
    {
        return this->rows();
    }
    else
    {
        return 0;
    }   
}

inline Index size() const {
    assert(this->cols() == 1);
    return this->rows();
}

void display() const {
    
}


inline void zero() {
    this->setZero();
}

inline Scalar norm2() {
    return this->norm();
}

// inline Scalar normInf() {
//     return this->cwiseAbs().rowwise().sum().maxCoeff();
// }

size_t nnz(double tol = 1e-14) {
    size_t nnz = 0;
    double* arr = this->data();
    for (size_t i = 0; i < this->size(0) * this->size(1); ++i) {
      if (fabs(arr[i]) > tol) {
        nnz++;
      }
    }
    return nnz;
}

inline Scalar vector_sum() {
    assert(this->cols() == 1);
    return this->sum();
}

/** copy the vector into an array
   *
   *  \param data the memory where to copy the data
   *  \return the number of element written (size of the vector)
   */
inline unsigned int copyData(double* data) {
    assert(this->cols() == 1);
    unsigned int size = this->size();
    if constexpr(std::is_same_v<Scalar, double>) {
        std::memcpy(data, this->data(), size * sizeof(Scalar)); // WARNING : Remember that Eigen::Scalar has to be double
    } else {
        []<bool flag = false>() {
            static_assert(flag, "Eigen::Scalar has to be double");
        } ();
    }
    return size;
}