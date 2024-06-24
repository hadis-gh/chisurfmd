
template<typename T , int D> 
T& Vec<T, D>::operator[](int index) {
    return elements[index];
}

template<typename T , int D> 
const T& Vec<T, D>::operator[] (int index) const{
    return elements[index];
}

template<typename T , int D> 
Vec<T, D> Vec<T, D>::operator+(const Vec<T, D>& other) const {
    Vec <T, D> result(*this);
    result += other;
    return result;
}

template<typename T , int D>
template<typename TT>
Vec<T, D>& Vec<T, D>::operator+=(const Vec<TT, D>& other) {
    for (int i=0; i < D; i++)
        elements[i] = elements[i] + other.elements[i];
    return *this;
}

template<typename T , int D> 
Vec<T, D> Vec<T, D>::operator-(const Vec<T, D>& other) const {
    Vec <T, D> result(*this);
    result -= other;
    return result;
}

template<typename T , int D>
template<typename TT> 
Vec<T, D>& Vec<T, D>::operator-=(const Vec<TT, D>& other){
    for (int i=0; i < D; i++)
        elements[i] = elements[i] - other[i];
    return *this;
}

template<typename T , int D> 
Vec<T, D> Vec<T, D>::operator-() const {
    Vec<T, D> result;
    for (int i=0; i<D; i++)
        result.elements[i] = -elements[i];
    return result;
}

template<typename T , int D> 
T Vec<T, D>::operator* (const Vec<T, D>& other) const {
    T result=0;
    for (int i=0; i<D; i++)
        result += elements[i]* other.elements[i];
    return result;
}

template<typename T , int D> 
Vec<T, D> Vec<T, D>::operator*(const T& number) const {
    Vec<T, D> result(*this);
    result *= number;
    return result;
}

template<typename T , int D> 
Vec<T, D> Vec<T, D>::operator/(const T& number) const {
    Vec<T, D> result(*this);
    result /= number;
    return result;
}

template<typename T , int D>
Vec<T, D>& Vec<T, D>::operator *= (const T& number){
    for (int i=0; i< D; i++)
        elements[i] = number* elements[i];
    return *this;
}

template<typename T , int D>
Vec<T, D>& Vec<T, D>::operator /= (const T& number){
    for (int i=0; i< D; i++)
        elements[i] = elements[i]/ number;
    return *this;
}


template<typename T , int D>
template<typename TT>    
bool Vec<T, D>::operator==(const Vec<TT, D>& other)const{
    for (int i=0; i < D; i++){
        if (elements[i] != other[i])
            return false;
    }
    return true;
}
template<typename T , int D>
bool Vec<T, D>::operator>=(const T& number) const{ 
    for (int i=0; i<D; i++){
        if (elements[i]< number){
            return false;
        }
    }
    return true;
}

template<typename T , int D>
bool Vec<T, D>::operator<=(const T& number)const{ 
        return !(*this > number);
    }

template<typename T , int D>
bool Vec<T, D>::operator>(const T& number)const{ 
    for (int i=0; i<D; i++){
        if (elements[i]<= number){
            return false;
        }
    }
    return true;
}

template<typename T , int D>
bool Vec<T, D>::operator<(const T& number)const{ 
    for (int i=0; i<D; i++){
        if (elements[i]>= number){
            return false;
        }
    }
    return true;
    }

template<typename T , int D>
Vec<T, D> Vec<T, D>::operator- (const T& number)const{
    for (int i=0; i< D; i++)
        elements[i] -= number;
    return *this;
}

template<typename T , int D>
Vec<T, D> Vec<T, D>::operator+ (const T& number)const{
    for (int i=0; i< D; i++)
        elements[i] += number;
    return *this;
}
template<typename T , int D>
T Vec<T, D>::sum()const{
    T result = elements[0];
    for (int i=1; i<D; i++)
        result += elements[i];
    return result;    
}

template<typename T , int D>
T Vec<T, D>::dot(const Vec<T, D>& other)const{
    Vec<T, D> result;
    for (int i=0; i<D; i++)
        result.elements[i] = elements[i]* other.elements[i];
    return result.sum();
}

template<typename T , int D>
T Vec<T, D>::abs2()const{
    return dot(*this);
}

template<typename T , int D> 
Vec<T, D> operator*(const T& number, const Vec<T, D>& vec){
    Vec<T, D> result(vec);
    result *= number;
    return result;
}

template<typename T, int D>
std::ostream& operator<<(std::ostream& COUT, const Vec<T, D>& ref) {
    COUT << "(" ;
    COUT << ref[0];
    for (int i=1; i<D; i++)
        COUT << " , " << ref[i];
    COUT << ")";
    return COUT;
}