// Copyright (c) 2025 Iviesever

#pragma once

template<typename T>
struct computation_type
{
	using type = T;
};

template<> struct computation_type<short> { using type = double; };
template<> struct computation_type<int> { using type = double; };
template<> struct computation_type<long> { using type = double; };
template<> struct computation_type<long long> { using type = double; };
template<> struct computation_type<float> { using type = double; };

template<> struct computation_type<double> { using type = double; };
template<> struct computation_type<long double> { using type = long double; };

template<typename T>
using computation_type_t = typename computation_type<T>::type;
