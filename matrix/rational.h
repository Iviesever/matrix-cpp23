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
// 添加一个 requires 子句以避免与 std::string 和 std::string_view 的默认格式化器冲突
	requires (!std::is_same_v<Range, std::string> &&
!std::is_same_v<Range, std::string_view>)
struct std::formatter<Range>
{

	// parse 函数可以保持简单，因为它不处理特殊的格式说明符
	constexpr auto parse(std::format_parse_context & ctx)
	{
		return ctx.begin();
	}

	// format 函数执行实际的格式化工作
	template <typename FormatContext>
	auto format(const Range & r, FormatContext & ctx) const
	{
		auto out = ctx.out();

		// 1. 打印开括号
		*out++ = '[';

		// 2. 迭代范围内的元素
		bool first = true;
		for(const auto & val : r)
		{
			if(!first)
			{
				// 在元素之间打印分隔符
				out = std::format_to(out, ", ");
			}
			first = false;

			// 3. 递归地调用 std::format_to 来格式化当前元素
			//    这会自动为 val 找到正确的格式化器（可能是另一个范围，或是 rational，或是 int 等）
			out = std::format_to(out, "{}", val);
		}

		// 4. 打印闭括号
		*out++ = ']';

		return out;
	}
};