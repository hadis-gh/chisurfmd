
template<typename T, int M, int N>
T& Mat<T, M, N>::operator()(int i, int j) {
    return elements[i][j];
}
template<typename T, int M, int N>
const T& Mat<T, M, N>::operator()(int i, int j) const {
    return elements[i][j];
}

template<typename T, int M, int N>
Vec<T, N>& Mat<T, M, N>::operator()(int i){
    return elements[i];
}

template<typename T, int M, int N>
const Vec<T, N>& Mat<T, M, N>::operator()(int i) const{
    return elements[i];
}




template<typename T, int M, int N>
Mat<T, M, N> Mat<T, M, N>::transpose() const {
    Mat<T, M, N> T_mat;
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            T_mat(j, i) = elements[i][j]; 
        }
    }
    return T_mat;
}

template<typename T, int M, int N>
template <int MM, int NN>
Mat<T, M, N> Mat<T, M, N>::operator* (const Mat<T, MM, NN> &other) const {
    Mat<T, M, N> result;
    static_assert(N == MM, "matrix multiplication is not allowded!");

    for (int i=0; i<M; i++){
        for (int k=0; k<NN; k++){
            for (int j=0; j<N; j++){
                result(i, j) += elements[i][k] * other(k, j);
            }
        }
    }
    return result;
}

template<typename T, int M, int N>
Mat<T, M, N> & Mat<T, M, N>::operator* (const T &scaler) {
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            elements[i][j] *= scaler;
        }
    }
    return *this;
}

template<typename T, int M, int N>
template <int MM, int NN>
Mat<T, M, N> Mat<T, M, N>::inner_product (const Mat<T, MM, NN> &other){
    return ((*this).transpose()) * other;
}
template<typename T, int M, int N>
template <int MM, int NN>
Mat<T, M, N> Mat<T, M, N>::outer_product (Mat<T, MM, NN> &other) const{
    return (*this) * (other.transpose());
}
template<typename T, int M, int N>
T Mat<T, M, N>::get_determinant() const{
    T determinant = elements[0][0] *elements[1][1] - elements[0][1] *elements[1][0];
    // static_assert(determinant != 0, "the matrix does not have inverse! (zero determinant)"); //How can do this?
    return determinant;
}
template<typename T, int M, int N>
Mat<T, M, N> Mat<T, M, N>::inverse () const{
    static_assert(N == 2 && M == 2, "only implemented for 2x2!");
    Mat<T, M, N> inverse_mat;

    inverse_mat (0, 0) = elements[1][1];
    inverse_mat (0, 1) = - elements[0][1];
    inverse_mat (1, 0) = - elements[1][0];
    inverse_mat (1, 1) = elements[0][0];

    return inverse_mat*(1/inverse_mat.get_determinant());        
}
template<typename T, int M, int N>
template <int MM, int NN>
bool Mat<T, M, N>::operator==(const Mat<T, MM, NN> &other) const{
    if (M != MM || N != NN){
        return false;
    }
    for (int i=0; i<M; i++){
        for (int j=0; j<N; j++){
            if (elements[i][j] != other(i, j)){
                return false;
            }
        }
    }
    return true;
}
template<typename T, int M, int N>
template <int MM, int NN>
bool Mat<T, M, N>::operator!=(const Mat<T, MM, NN> &other) const{
    return !(*this == other);
}

template<typename T, int M, int N>
std::ostream& operator<<(std::ostream& COUT, const Mat<T, M, N>& ref){
COUT << "{" << std::endl;
for (int i=0; i<M; i++){
    for (int j=0; j<N; j++){
        COUT << ref(i, j) << ", ";
    }
    COUT << std::endl;
}
COUT << "}";
return COUT;    
};