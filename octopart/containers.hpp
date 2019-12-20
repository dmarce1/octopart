#pragma once

#include <array>
#include <cassert>
#include <vector>

template<class T>
class vector: public std::vector<T> {
	using base_type = std::vector<T>;
	using base_type::base_type;

public:
	T operator[](int i) const {
		assert(i >= 0);
		assert(i < base_type::size());
		return base_type::operator [](i);
	}
	T& operator[](int i) {
		assert(i >= 0);
		assert(i < base_type::size());
		return base_type::operator [](i);
	}

};

template<class T, auto N>
class array: public std::array<T, N> {
	using base_type = std::array<T,N>;
	using base_type::base_type;

public:
	T operator[](int i) const {
		assert(i >= 0);
		assert(i < base_type::size());
		return base_type::operator [](i);
	}
	T& operator[](int i) {
		assert(i >= 0);
		assert(i < base_type::size());
		return base_type::operator [](i);
	}

};
