//////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   fc.cpp
/// @brief  计算阶乘数和组合数
///
/// 计算小整数范围的阶乘和组合数
///
/// @version 1.1
/// @author  曹业
/// @date    2016年8月1日
///
///
/// 修订说明：第一个版本
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cstddef>

///      阶乘数
///
///      使用直接连乘，计算小整数的阶乘数
///
///      @param m                阶乘截止整数
///      @return                 阶乘数
///      @see
///      @note
size_t Factorial(size_t m)
{
    size_t f = 1;
    
    for (size_t i = 1; i <= m;       i++) { f *= i; }
    
    return f;
}

///      组合数
///
///      使用组合数公式，c(n,m)=n!/(m!*(n-m)!)，计算小整数范围的组合数
///
///      @param m                阶乘截止整数
///      @return                 阶乘数
///      @see
///      @note
size_t Combination(size_t n, size_t m)
{

	size_t bound = std::max(n - m, m);
    size_t c = 1;
    
    for (size_t i = n; i >  bound;   i--) { c *= i; }

    for (size_t i = 1; i <= n-bound; i++) { c /= i; }
    
    return c;
}
