#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <array>
#include <stdexcept>
#include <string>
#include <utility>

namespace quadrature
{

namespace detail
{

template <int order>
void validate(int n)
{
    static_assert(order > 0, "quadrature order must be a positive integer");
    if (n <= 0) {
        throw std::invalid_argument(
            "\nthe number of quadrature subintervals must be a positive integer "
            "(n > 0) but n = " +
            std::to_string(n));
    }
}

}  // namespace detail

template <int order>
struct gauss_legendre_rule;

template <>
struct gauss_legendre_rule<2> {
    static constexpr std::array<double, 2> abscissa = {
        -0.5773502691896258,
        0.5773502691896258};
    static constexpr std::array<double, 2> weights = {1.0, 1.0};
};

template <>
struct gauss_legendre_rule<3> {
    static constexpr std::array<double, 3> abscissa = {
        -0.7745966692414834,
        0.0,
        0.7745966692414834};
    static constexpr std::array<double, 3> weights = {
        5.0 / 9.0,
        8.0 / 9.0,
        5.0 / 9.0};
};

template <>
struct gauss_legendre_rule<4> {
    static constexpr std::array<double, 4> abscissa = {
        -0.8611363115940526,
        -0.3399810435848563,
        0.3399810435848563,
        0.8611363115940526};
    static constexpr std::array<double, 4> weights = {
        0.34785484513745385,
        0.6521451548625461,
        0.6521451548625461,
        0.34785484513745385};
};

template <int order, typename T = double, typename integrand>
T gauss_legendre(integrand&& f, double a, double b, int n)
{
    detail::validate<order>(n);
    using rule = gauss_legendre_rule<order>;
    double const dx = (b - a) / n;
    T result{};
    for (int i = 0; i < n; ++i) {
        double const mid = a + (i + 0.5) * dx;
        for (int j = 0; j < order; ++j) {
            result += f(mid + 0.5 * dx * rule::abscissa[j]) * rule::weights[j];
        }
    }
    return result * (0.5 * dx);
}

}  // namespace quadrature

#endif
