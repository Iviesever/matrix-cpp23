// Copyright (c) 2025 Iviesever

#pragma once

#include <concepts>
#include <complex>

template<typename T>
concept is_matrix_element =
std::default_initializable<T> &&
std::copy_constructible<T> &&
	requires(T a, T b)
{
	{ a + b } -> std::same_as<T>;
	{ a - b } -> std::same_as<T>;
	{ a * b } -> std::same_as<T>;
	{ a / b } -> std::same_as<T>;
} &&
	requires(T a)
{
	{ -a } -> std::same_as<T>;
	a <=> 0;
};

template<typename T>
concept is_matrix = requires(T m)
{
	m.rows();
	m.cols();
};

template<typename T>
struct is_complex_impl: std::false_type {}; 

template<typename U>
struct is_complex_impl<std::complex<U>>: std::true_type {};

template<typename T>
concept is_complex = is_complex_impl<T>::value;
