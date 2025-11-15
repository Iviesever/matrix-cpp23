#pragma once
#include <concepts>
#include <type_traits>
#include <cmath>
#include <limits>

auto try_cast_int(std::floating_point auto value)
{
	const double nearest_int = std::round(value);

	const double diff = std::abs(value - nearest_int);

	constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();
	constexpr double relative_factor = 1000.0;

	const double absolute_tolerance = std::sqrt(std::numeric_limits<double>::epsilon());
	const double tolerance = std::max(absolute_tolerance, std::abs(value) * relative_factor * machine_epsilon);


	if(diff < tolerance)
	{
		return nearest_int;
	}

	return value;
};

auto is_zero(auto value)
{
	using U = computation_type<decltype(value)>::type;

	if constexpr(std::is_floating_point_v<U>)
	{
		U absolute_tolerance = std::sqrt(std::numeric_limits<U>::epsilon());
		return std::abs(value) < absolute_tolerance;
	}
	else
	{
		return value == U{};
	}
}

auto abs(auto value)
{
	if(value < 0)
		value = -value;
	return value;
}