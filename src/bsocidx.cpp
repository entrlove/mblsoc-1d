//////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   bsocidx.cpp
/// @brief  生成基底的索引形式
///
/// 生成给定基底的索引形式
///
/// @version 1.1
/// @author  曹业
/// @date    2016年8月1日
///
///
/// 修订说明：第一个版本
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include "dec.h"

///      基底的索引形式
///
///      根据约化和保留位置，计算基底的索引形式
///
///      @param site_num            晶格格点数目
///      @param basis_without_index 基底
///      @param position_reduce     约化部分的位置
///      @param position_remain     保留部分的位置
///      @return                    基底的索引形式
///      @see
///      @note

basis_soc_with_index Basis_soc_index(size_t site_num,
                                     std::string &basis_without_index,
                                     std::vector<size_t> &position_reduce,
                                     std::vector<size_t> &position_remain)
{
    /// 定义局部变量
    basis_soc_with_index b_soc_with_index;
    std::string s_reduce(2*position_reduce.size(), '0');
    std::string s_remain(2*position_remain.size(), '0');
    size_t particle_reduce = 0;
    size_t particle_remain = 0;
    
    /// 设定基底索引类对象的基底值
    b_soc_with_index.basis = basis_without_index;
    
    /// 提取约化部分的图案和粒子数
    for (size_t i = 0; i < position_reduce.size(); i++)
    {
        if ((s_reduce[i] = basis_without_index[position_reduce[i]]) == '1') {particle_reduce++;};
        if ((s_reduce[i+position_reduce.size()] = basis_without_index[position_reduce[i]+site_num]) == '1')
        {particle_reduce++;};
    }
    /// 提取保留部分的图案和粒子数
    for (size_t i = 0; i < position_remain.size(); i++)
    {
        if ((s_remain[i] = basis_without_index[position_remain[i]]) == '1') {particle_remain++;};
        if ((s_remain[i+position_remain.size()] = basis_without_index[position_remain[i]+site_num]) == '1')
        {particle_remain++;};
    }
    /// 记录基底索引值
    b_soc_with_index.index_of_basis = index_val(s_reduce, s_remain);
    /// 记录约化和保留部分的粒子数
    b_soc_with_index.particle_of_subspace = particle_num_subspace(particle_reduce, particle_remain);
    
    /// 返回基底的索引形式
    return b_soc_with_index;
} 
