#include <iostream>
#include <array>

template<typename T, int D = 2>
class Vector
{
	std::array<T, D> m_elements;

	public:

	T& operator[] (int i)
	{
		return m_elements[i];
	}

	const T& operator[] (int i) const
	{
		return m_elements[i];
	}

	auto operator-() const
	{
		Vector<T, D> ret;
		for(int i = 0; i < D; ++i)
			ret[i] = -m_elements[i];
		return ret;
	}

	auto operator+(const Vector<T, D> b) const
	{	
		Vector<T, D> ret;
		for(int i = 0; i < D; ++i)
			ret[i] = m_elements[i] + b[i];
		return ret;
	}

};

int main()
{
	Vector<int> v1;
	Vector<int> v2;
	v1 [0] = {1};
	v1 [1] = 3;
	v2 [0] = {2};
	v2 [1] = 3;

	const auto v_neg = -v1;
	const auto v_sum = v1 + v2;

	std::cout << "-v1:" << v_neg[0] << std::endl;
	std::cout << "v1 + v2:" << v_sum[0] << std::endl;


};
