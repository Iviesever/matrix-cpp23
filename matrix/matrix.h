// Copyright (c) 2025 Iviesever

#pragma once

#include "matrix_concepts.h"
#include "matrix_traits.h"
#include "matrix_math.h"
#include <format>
#include <vector>
#include <variant>
#include <array>
#include <mdspan>
#include <print>
#include <memory>
#include <ranges>
#include <algorithm>
#include <utility>
#include <span>
#include <functional>
#include <numeric>

struct _Op_swap_rows
{
	size_t row1, row2;
};

template<typename U>
struct _Op_add_row_multiple
{
	size_t row1, row2;
	U scalar;
};

template<typename U>
struct _Op_scale_row
{
	size_t row1;
	U scalar;
};

template<typename U>
using operation_type = std::variant<_Op_swap_rows, _Op_add_row_multiple<U>, _Op_scale_row<U>>;

template<is_matrix_element T, size_t ...Dims>
	requires (sizeof...(Dims) == 2 || sizeof...(Dims) == 0)
class matrix
{
public:

	static constexpr bool is_static = bool(sizeof...(Dims) > 0);

	using data_type = decltype([]
		{
			if constexpr(is_static)
			{
				return std::array<T, (Dims * ...)>{};
			}
			else
			{
				return std::vector<T>{};
			}
		}());

	using extents_type = decltype([]
		{
			if constexpr(is_static)
			{
				return std::extents<size_t, Dims...>{};
			}
			else
			{
				return std::dextents<size_t, 2>{};
			}
		}());

	using view_type = std::mdspan<T, extents_type>;

	using computation_t = computation_type_t<T>;

	using computation_matrix = matrix<computation_t, Dims...>;

	using complex_t = std::complex<double>;

	using operation_t = std::variant<
		_Op_swap_rows, _Op_add_row_multiple<computation_t>, _Op_scale_row<computation_t>>;

private:

	data_type data;
	view_type view;

	mutable std::vector<operation_t> _History;

	mutable std::unique_ptr<computation_matrix>
		_Row_echelon_matrix,
		_Reduced_Row_Echelon_matrix,
		_Inverse_matrix,
		_Adjugate_matrix;

	mutable std::unique_ptr<size_t> _Rand;

	mutable size_t swap_times{};

	template<is_matrix Matrix>
	void transform(const Matrix & other)
	{
		std::ranges::transform
		(
			other.data,
			data.begin(),
			[](const Matrix::value_type & value)
			{
				return static_cast<value_type>(value);
			}
		);
	}

	template<is_matrix Matrix>
	void _QR_decomposition(Matrix & Q, Matrix & R)
	{
		const auto n = Q.rows();
		if(n != Q.cols())
		{
			throw std::logic_error("QR decomposition requires a square matrix.");
		}

		std::vector<double> 
			x_data(n),
			v_data(n),
			vT_R_sub_data(n),
			Q_sub_v_data(n);

		for(auto j : std::views::iota(0ull, n))
		{
			const auto sub_dim = n - j;
			if(sub_dim <= 1) continue;

			auto x = std::span(x_data.data(), sub_dim);
			for(auto i : std::views::iota(0ull, sub_dim))
			{
				x[i] = R[j + i, j];
			}

			auto x_norm_sq = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
			auto x_norm = std::sqrt(x_norm_sq);

			if(is_zero(x_norm)) continue;

			auto sign = (x[0] >= 0) ? 1.0 : -1.0;
			x[0] += sign * x_norm;

			auto v_norm_sq = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
			auto v_norm = std::sqrt(v_norm_sq);

			if(is_zero(v_norm)) continue;

			auto v = std::span(v_data.data(), sub_dim);

			std::ranges::for_each(std::views::iota(0ull, sub_dim), [&](auto i)
				{
					v[i] = x[i] / v_norm;
				});

			const auto vT_R_sub_size = R.cols() - j;
			auto vT_R_sub = std::span(vT_R_sub_data.data(), vT_R_sub_size);
			for(auto col : std::views::iota(0ull, vT_R_sub.size()))
			{
				double sum{};
				for(auto row : std::views::iota(0ull, sub_dim))
				{
					sum += v[row] * R[j + row, j + col];
				}
				vT_R_sub[col] = sum;
			}

			for(auto row : std::views::iota(0ull, sub_dim))
			{
				for(auto col : std::views::iota(0ull, vT_R_sub.size()))
				{
					R[j + row, j + col] -= 2.0 * v[row] * vT_R_sub[col];
				}
			}

			auto Q_sub_v = std::span(Q_sub_v_data.data(), n);
			for(auto row : std::views::iota(0ull, n))
			{
				double sum{};
				for(auto col : std::views::iota(0ull, sub_dim))
				{
					sum += Q[row, j + col] * v[col];
				}
				Q_sub_v[row] = sum;
			}

			for(auto row : std::views::iota(0ull, n))
			{
				for(auto col : std::views::iota(0ull, sub_dim))
				{
					Q[row, j + col] -= 2.0 * Q_sub_v[row] * v[col];
				}
			}
		}
	}

	template<is_matrix_element, size_t...D>
		requires (sizeof...(D) == 2 || sizeof...(D) == 0) friend class matrix;


public:

	auto solve_homogeneous()
	{
		auto rank = r();

		auto A = rref();

		std::vector<std::vector<computation_t>> basis;

		if(rank == cols())
		{
			std::vector<computation_t> trivial_solution(cols(), 0);
			basis.push_back(trivial_solution);
		}
		else
		{
			auto s = cols() - rank;

			// 找到主元列
			std::vector<size_t> pivot_columns;
			for(size_t i = 0; i < rank; ++i)
			{
				for(size_t j = 0; j < cols(); ++j)
				{
					if(A[i, j] == 1)
					{
						pivot_columns.push_back(j);
						break;
					}
				}
			}

			// 找到自由列
			std::vector<size_t> free_columns;
			for(size_t j = 0; j < cols(); ++j)
			{
				if(std::find(pivot_columns.begin(), pivot_columns.end(), j) == pivot_columns.end())
				{
					free_columns.push_back(j);
				}
			}

			// 为每个自由变量创建一个解向量
			for(auto free_col_index : free_columns)
			{
				std::vector<computation_t> solution_vector(cols(), 0);

				// 将当前自由变量设为1
				solution_vector[free_col_index] = 1;

				// 根据RREF矩阵求解主元变量
				for(size_t i = 0; i < rank; ++i)
				{
					auto pivot_col = pivot_columns[i];
					solution_vector[pivot_col] = - A[i, free_col_index];
				}
				basis.push_back(solution_vector);
			}

		}
			return basis;

	}

	template<is_matrix Matrix>
	auto solve_non_homogeneous(const Matrix & other)
	{
		if(other.cols() != 1 || other.rows() != rows())
			throw std::invalid_argument("can't solve equations");

		auto r1 = r();


		auto A = rref();
		auto b = computation_matrix(other);
		apply_history_to(b);

		for(size_t i = r1; i < rows(); ++i)
		{
			if(!is_zero(b[i, 0]))
			{
				return std::vector<std::vector<computation_t>>{};
			}
		}

		std::vector<computation_t> particular_solution(cols(), 0);
		std::vector<size_t> pivot_columns;
		pivot_columns.reserve(r1);
		for(size_t i = 0; i < r1; ++i)
		{
			for(size_t j = 0; j < cols(); ++j)
			{
				if(A[i, j] == 1)
				{
					pivot_columns.push_back(j);
					// 主元变量的值由 b 决定
					particular_solution[j] = b[i, 0];
					break;
				}
			}
		}

		if(r1 == cols()) // 情况1: 唯一解
		{
			// 没有自由变量，齐次解只有零向量，通解就是该特解
			std::vector<std::vector<computation_t>> result;
			result.push_back(particular_solution);
			return result;
		}
		else // r1 < cols(), 情况2: 无限解
		{
			// 通解 = 特解 + 齐次解
			auto basis = solve_homogeneous(); // 获取齐次解的基

			basis.push_back(particular_solution); // 将特解添加到返回集合中

			return basis;
		}

	}

	auto lambda(size_t maxIter = 1000)
	{
		const size_t n = rows();
		if(n == 0) return std::vector<complex_t> {};
		if(n != cols())
		{
			throw std::logic_error("rows != cols");
		}

		auto A = matrix<double>(*this);

		matrix<double> Q(n, n), R(n, n);

		for(auto iter : std::views::iota(0ull, maxIter))
		{
			double max_off_diag{};
			for(auto i : std::views::iota(1ull, n))
			{
				for(auto j : std::views::iota(0ull, i))
				{
					max_off_diag = std::max(max_off_diag, std::abs((double)A[i, j]));
				}
			}
			if(is_zero(max_off_diag)) break;

			R = A;
			Q = Q.generate_e();

			_QR_decomposition(Q, R);

			A = R * Q;
		}

		std::vector<complex_t> eigenvalues;
		size_t i{};
		while(i < n)
		{
			if(i == n - 1 || is_zero(A[i + 1, i]))
			{
				eigenvalues.emplace_back(A[i, i], 0.0);
				++i;
			}
			else
			{
				double a = A[i, i], b = A[i, i + 1];
				double c = A[i + 1, i], d = A[i + 1, i + 1];

				double trace = a + d;
				double det = a * d - b * c;
				double discriminant = trace * trace - 4 * det;

				if(discriminant >= 0)
				{
					double sqrt_disc = std::sqrt(discriminant);
					eigenvalues.emplace_back((trace + sqrt_disc) / 2.0, 0.0);
					eigenvalues.emplace_back((trace - sqrt_disc) / 2.0, 0.0);
				}
				else
				{
					double sqrt_abs_disc = std::sqrt(-discriminant);
					eigenvalues.emplace_back(trace / 2.0, sqrt_abs_disc / 2.0);
					eigenvalues.emplace_back(trace / 2.0, -sqrt_abs_disc / 2.0);
				}
				i += 2;
			}
		}

		std::ranges::transform(eigenvalues, eigenvalues.begin(), [&](const auto & c)
			{
				return complex_t{ try_cast_int(c.real()), try_cast_int(c.imag()) };
			});

		std::ranges::sort(eigenvalues, [](const auto & c1, const auto & c2)
			{ 
				if(c1.real() != c2.real())
					return c1.real() < c2.real();
				return c1.imag() < c2.imag();
			});

		return eigenvalues;
	}

	using value_type = T;

	matrix() requires(is_static) : data{}, view(data.data()) {}

	matrix(size_t rows, size_t cols) requires(!is_static) : data(rows * cols), view(data.data(), rows, cols) {}

	matrix(const matrix & other): data(other.data), view(data.data(), other.rows(), other.cols()) {}

	template<is_matrix Matrix>
		requires(is_static)
	explicit matrix(const Matrix & other): view(data.data(), other.rows(), other.cols())
	{
		transform(other);
	}

	template<is_matrix Matrix>
		requires(!is_static)
	explicit matrix(const Matrix & other)
	{
		data.resize(other.rows() * other.cols());
		view = view_type(data.data(), other.rows(), other.cols());

		transform(other);
	}

	void fill(const value_type & value)
	{
		std::ranges::fill(data, value);
	}

	constexpr matrix(std::initializer_list<std::initializer_list<value_type>> list)
		requires (!is_static)
	{
		size_t rows = list.size();
		size_t cols = (rows > 0) ? list.begin()->size() : 0;

		bool is_valid = std::ranges::all_of(list, [&cols](const auto & row)
			{
				return row.size() == cols;
			});

		if(!is_valid)
		{
			throw std::invalid_argument("Invalid cols");
		}

		data.assign(rows * cols, T{});
		view = view_type(data.data(), rows, cols);

		for(size_t r{}; auto & row_list : list)
		{
			std::ranges::copy(row_list, &view[r, 0]);
			++r;
		}
	}

	constexpr matrix(std::initializer_list<std::initializer_list<value_type>> list)
		requires (is_static) : data{}, view(data.data())
	{
		constexpr size_t static_rows = view.static_extent(0);
		constexpr size_t static_cols = view.static_extent(1);

		if(list.size() != static_rows)
		{
			throw std::invalid_argument("Invalid rows");
		}

		bool is_valid = std::ranges::all_of(list, [&static_rows](const auto & row)
			{
				return row.size() == static_cols;
			});

		if(!is_valid)
		{
			throw std::invalid_argument("Invalid cols");
		}

		for(size_t r{}; const auto & row_list : list)
		{
			std::ranges::copy(row_list, &view[r, 0]);
			++r;
		}

	}

	auto operator=(const matrix & other)
	{
		if(this == &other)
		{
			return *this;
		}

		data = other.data;
		view = view_type(data.data(), other.rows(), other.cols());

		return *this;
	}

	decltype(auto) operator[](this auto && self, size_t r, size_t c)
	{
		return std::forward_like<decltype(self)>(self.view[r, c]);
	}

	constexpr auto rows() const
	{
		if constexpr(is_static)
		{
			return view.static_extent(0);
		}
		return view.extent(0);
	}

	constexpr auto cols() const
	{
		if constexpr(is_static)
		{
			return view.static_extent(1);
		}
		return view.extent(1);
	}

	void swap_rows(size_t row_1, size_t row_2)
	{
		for(size_t j{}; j < cols(); ++j)
			std::swap(view[row_1, j], view[row_2, j]);
	}

	void swap_cols(size_t col_1, size_t col_2)
	{
		for(size_t i{}; i < rows(); ++i)
			std::swap(view[i, col_1], view[i, col_2]);
	}

	void scale_row(size_t row, computation_t k)
	{
		for(size_t j{}; j < cols(); ++j)
			if(auto & num{ view[row, j] }; num != 0)
				num *= k;
	}

	void scale_col(size_t col, computation_t k)
	{
		for(size_t i{}; i < rows(); ++i)
			if(auto & num{ view[i, col] }; num != 0)
				num *= k;
	}

	void add_row_multiple(size_t row_1, size_t row_2, computation_t k)
	{
		for(size_t j{}; j < cols(); ++j)
			view[row_1, j] += k * view[row_2, j];
	}

	void add_col_multiple(size_t col_1, size_t col_2, computation_t k)
	{
		for(size_t i{}; i < rows(); ++i)
			view[i, col_1] += k * view[i, col_2];
	}



	const auto & ref() const
	{
		if(_Row_echelon_matrix)
			return *_Row_echelon_matrix;

		_History.clear();

		auto temp{ computation_matrix(*this) };

		size_t i{}, j{};

		while(i < rows() && j < cols())
		{
			auto higher_i{ i };

			for(auto next_i : std::views::iota(i, rows()) | std::views::drop(1))
			{
				if(abs(temp[next_i, j]) > abs(temp[higher_i, j]))
					higher_i = next_i;
			}

			if(is_zero(temp[higher_i, j]))
			{
				temp[higher_i, j] = 0;
				++j;
				continue;
			}

			if(i != higher_i)
			{
				temp.swap_rows(i, higher_i);

				_History.emplace_back(_Op_swap_rows{ i, higher_i });

				++swap_times;
			}

			for(auto next_i : std::views::iota(i, rows()) | std::views::drop(1))
			{
				if(is_zero(temp[next_i, j]))
				{
					temp[next_i, j] = 0;
					continue;
				}

				computation_t multiple_num = -1 * temp[next_i, j] / temp[i, j];

				temp.add_row_multiple(next_i, i, multiple_num);

				_History.emplace_back(_Op_add_row_multiple<computation_t>{ next_i, i, multiple_num });
			}

			++i, ++j;
		}

		_Row_echelon_matrix = std::make_unique<computation_matrix>(std::move(temp));
		return *_Row_echelon_matrix;
	}

	const auto & rref() const
	{
		if(_Reduced_Row_Echelon_matrix)
			return *_Reduced_Row_Echelon_matrix;

		auto temp = [&]
			{
				if(!_Row_echelon_matrix)
					return ref();
				else
					return *_Row_echelon_matrix;
			}();

		for(size_t j{}; auto i : std::views::iota(0ull, rows()) | std::views::reverse)
		{
			while(j < cols() && is_zero(temp[i, j]))
			{
				temp[i, j] = 0;
				++j;
			}

			if(j == cols())
			{
				j = 0;
				continue;
			}

			computation_t scale_num = 1 / temp[i, j];

			_History.emplace_back(_Op_scale_row<computation_t>{ i, scale_num });

			temp.scale_row(i, scale_num);


			for(auto k : std::views::iota(0ull, i) | std::views::reverse)
			{
				if(is_zero(temp[k, j]))
				{
					temp[k, j] = 0;
					continue;
				}

				auto multiple_num = -temp[k, j];

				_History.emplace_back(_Op_add_row_multiple<computation_t>{ k, i, multiple_num });

				temp.add_row_multiple(k, i, multiple_num);
			}

			j = 0;
		}

		_Reduced_Row_Echelon_matrix = std::make_unique<decltype(temp)>(std::move(temp));
		return *_Reduced_Row_Echelon_matrix;
	}

	auto det() const
	{
		if(rows() != cols())
			throw std::logic_error("There is no determinant.");

		const auto & temp = [&]
			{
				if(!_Row_echelon_matrix)
					return ref();
				else
					return *_Row_echelon_matrix;
			}();

		auto num = std::ranges::fold_left(std::views::iota(0ull, rows()) | std::views::transform([&](auto i)
			{
				return temp[i, i];
			}), computation_t{ 1 }, std::multiplies<>{});

		if(swap_times % 2 == 1)
			num *= -1;

		if(num == -0) num = 0;

		return num;
	}

	const auto & inv() const
	{
		if(det() == 0)
			throw std::logic_error("det = 0");

		if(!_Reduced_Row_Echelon_matrix)
			rref();

		computation_matrix E(generate_e());

		apply_history_to(E);

		_Inverse_matrix = std::make_unique<computation_matrix>(std::move(E));

		return *_Inverse_matrix;
	}

	template<is_matrix Matrix>
	void apply_history_to(Matrix & result) const
	{
		if(result.rows() != rows())
			throw std::invalid_argument("matrix's size isn't same");

		if(!_Reduced_Row_Echelon_matrix)
			rref();

		for(const auto & op_variant : _History)
		{
			std::visit([&](auto && op)
				{
					using OpType = std::decay_t<decltype(op)>;
					if constexpr(std::is_same_v<OpType, _Op_swap_rows>)
					{
						//std::print("swap_rows {} {}\n", op.row1, op.row2);
						result.swap_rows(op.row1, op.row2);
					}
					else if constexpr(std::is_same_v<OpType, _Op_add_row_multiple<computation_t>>)
					{
						//std::print("add_row_multiple {} {} {}\n", op.row1, op.row2, op.scalar);
						result.add_row_multiple(op.row1, op.row2, op.scalar);
					}
					else if constexpr(std::is_same_v<OpType, _Op_scale_row<computation_t>>)
					{
						//std::print("scale_row {} {}\n", op.row1, op.scalar);
						result.scale_row(op.row1, op.scalar);
					}
				}, op_variant);
		}
	}

	const auto & adj() const
	{
		if(_Adjugate_matrix)
			return *_Adjugate_matrix;

		auto temp = [&]
			{
				if(!_Inverse_matrix)
					return inv();
				else
					return *_Inverse_matrix;
			}();

		auto d = det();

		for(auto i : std::views::iota(0ull, rows()))
			temp.scale_row(i, d);

		_Adjugate_matrix = std::make_unique<computation_matrix>(std::move(temp));
		return *_Adjugate_matrix;
	}

	const auto & r() const
	{
		if(_Rand)
			return *_Rand;

		_Rand = std::make_unique<size_t>(0ull);

		const auto & temp = [&]
			{
				if(!_Reduced_Row_Echelon_matrix)
					return rref();
				else
					return *_Reduced_Row_Echelon_matrix;
			}();

		for(size_t i{}, j{}; i < rows() && j < cols(); ++i, ++j)
		{
			while(j < cols() && temp[i, j] == 0)
				++j;

			if(j == cols())
				break;

			++(*_Rand);
		}

		return *_Rand;
	}

	auto tr()const
	{
		return std::ranges::fold_left(std::views::iota(0ull, rows()) | std::views::transform([&](auto i)
			{
				return view[i, i];
			}), computation_t{}, std::plus<>{});
	}

	auto generate_o()const
	{
		if constexpr(is_static)
			return matrix{};
		else
			return matrix{ rows(), cols() };
	}

	auto generate_e()const
	{
		auto e = [&]
			{
				if constexpr(is_static)
					return matrix{};
				else
					return matrix{ rows(), cols() };
			}();

		std::ranges::for_each(std::views::iota(0ull, rows()), [&](auto i)
			{
				e[i, i] = 1;
			});

		return e;
	}

};

template<is_matrix Matrix>
auto operator+(const Matrix & m1, const Matrix & m2)
{
	auto temp = m1.generate_o();

	for(size_t i{}; i < temp.rows(); ++i)
		for(size_t j{}; j < temp.cols(); ++j)
		{
			temp[i, j] = m1[i, j] + m2[i, j];
		}

	return temp;
}

template<is_matrix Matrix>
auto operator-(const Matrix & m1, const Matrix & m2)
{
	auto temp = m1.generate_o();

	for(size_t i{}; i < temp.rows(); ++i)
		for(size_t j{}; j < temp.cols(); ++j)
		{
			temp[i, j] = m1[i, j] - m2[i, j];
		}

	return temp;
}

template<is_matrix Matrix_1, is_matrix Matrix_2>
auto operator*(const Matrix_1 & m1, const Matrix_2 & m2)
{
	if(m1.cols() != m2.rows())
		throw std::invalid_argument("Matrix_1 cols != Matrix_2 rows");

	using _Value_type = std::common_type_t<typename Matrix_1::value_type, typename Matrix_2::value_type>;

	auto temp = matrix <_Value_type>{ m1.rows(), m2.cols() };

	for(size_t i{}; i < temp.rows(); ++i)
		for(size_t j{}; j < temp.cols(); ++j)
			for(size_t k{}; k < temp.cols(); ++k)
			{
				temp[i, k] += m1[i, j] * m2[j, k];
			}

	return temp;
}

template<is_matrix Matrix>
struct std::formatter<Matrix>
{
	constexpr auto parse(std::format_parse_context & ctx)
	{
		return ctx.begin();
	}

	auto format(const Matrix & m, std::format_context & ctx)const
	{
		auto out = ctx.out();

		if(m.rows() == 0)
			return std::format_to(out, "[]");

		for(size_t i{}; i < m.rows(); ++i)
		{
			out = std::format_to(out, "[ ");

			for(size_t j{}; j < m.cols(); ++j)
			{
				out = std::format_to(out, "{}", m[i, j]);

				if(j < m.cols() - 1)
				{
					out = std::format_to(out, "\t");
				}
			}

			if(i < m.rows() - 1)
				out = std::format_to(out, " ]\n");
			else
				out = std::format_to(out, " ]");
		}

		return out;
	}

};

template<is_complex Complex>
struct std::formatter<Complex>
{
	constexpr auto parse(std::format_parse_context & ctx)
	{
		return ctx.begin();
	}

	auto format(const Complex & c, std::format_context & ctx)const
	{

		if(c.real() == 0 && c.imag() != 0)
			return std::format_to(ctx.out(), "{}i", c.imag());

		if(c.imag() > 0)
			return std::format_to(ctx.out(), "{} + {}i", c.real(), c.imag());
		else if(c.imag() < 0)
		{
			auto t = -c.imag();
			return std::format_to(ctx.out(), "{} - {}i", c.real(), t);
		}
		else
			return std::format_to(ctx.out(), "{}", c.real());
	}

};

template<is_complex Complex>
struct std::formatter<std::vector<Complex>>
{
	constexpr auto parse(std::format_parse_context & ctx)
	{
		auto it = ctx.begin();
		if(it != ctx.end() && *it != '}')
		{
			sep_ = *it; // 例如格式字符串为"{:;}"，则分隔符为';'
			++it;
		}
		return it;
	}

	auto format(const std::vector<Complex> & vec, std::format_context & ctx) const
	{
		auto out = ctx.out();
		*out++ = '[';

		std::formatter<Complex> elem_formatter;

		bool first = true;
		for(const auto & elem : vec)
		{
			if(!first)
			{
				*out++ = sep_;
				if(sep_ != ' ') *out++ = ' ';
			}

			out = elem_formatter.format(elem, ctx);
			first = false;
		}

		*out++ = ']'; 
		return out;
	}

private:
	char sep_ = ',';
};
