// Copyright (c) 2025 Iviesever

#pragma once
#include <format>
#include <stdexcept>
#include <numeric>

struct rational
{
public:

	long long num{}, den{ 1 };

	rational() = default;

	rational(long long numerator, long long denominator = 1): num{ numerator }, den{ denominator }
	{
		if(den == 0)
		{
			throw std::invalid_argument("denominator can't be 0");
		}

		reduce();
	}

	explicit operator double() const
	{
		return static_cast<double>(num) / den;
	}

private:

	void reduce()
	{
		auto common = std::gcd(num, den);

		num /= common;
		den /= common;

		if(den < 0)
		{
			num = -num;
			den = -den;
		}
	}

};

inline auto operator+(const rational & a, const rational & b)
{
	return rational(a.num * b.den + a.den * b.num, a.den * b.den);
}

inline auto operator+=(rational & a, const rational & b)
{
	return a = rational(a.num * b.den + a.den * b.num, a.den * b.den);
}

inline auto operator-(const rational & a)
{
	return rational(a.num * -1, a.den);
}

inline auto operator-(const rational & a, const rational & b)
{
	return rational(a.num * b.den - a.den * b.num, a.den * b.den);
}

inline auto operator*(const rational & a, const rational & b)
{
	return rational(a.num * b.num, a.den * b.den);
}

inline auto operator*=(rational & a, const rational & b)
{
	return a = rational(a.num * b.num, a.den * b.den);
}

inline auto operator/(const rational & a, const rational & b)
{
	return rational(a.num * b.den, a.den * b.num);
}

inline auto operator<=>(const rational & a, const rational & b)
{
	return (a.num * b.den) <=> (a.den * b.num);
}

inline auto operator==(const rational & a, const rational & b)
{
	return (a.num * b.den) == (a.den * b.num);
}

template<>
struct std::formatter<rational>
{
	constexpr auto parse(std::format_parse_context & context)
	{
		return context.begin();
	}

	auto format(const rational & r, std::format_context & context)const
	{
		if(r.den == 1)
			return std::format_to(context.out(), "{}", r.num);
		else
			return std::format_to(context.out(), "{}/{}", r.num, r.den);
	}
};

template <std::ranges::range Range>
	requires (!std::is_same_v<Range, std::string> &&
!std::is_same_v<Range, std::string_view>)
struct std::formatter<Range>
{

	constexpr auto parse(std::format_parse_context & ctx)
	{
		return ctx.begin();
	}

	template <typename FormatContext>
	auto format(const Range & r, FormatContext & ctx) const
	{
		auto out = ctx.out();

		*out++ = '[';

		bool first = true;
		for(const auto & val : r)
		{
			if(!first)
			{
				out = std::format_to(out, ", ");
			}
			first = false;

			out = std::format_to(out, "{}", val);
		}

		// 4. 打印闭括号
		*out++ = ']';

		return out;
	}
};
