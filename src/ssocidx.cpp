//////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   ssocidx.cpp
/// @brief  生成系统基底空间的索引形式
///
/// 生成自旋1/2的自旋轨道耦合费米子基底空间的索引形式
///
/// @version 1.1
/// @author  曹业
/// @date    2016年8月1日
///
///
/// 修订说明：第一个版本
//////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "dec.h"

///      基底空间的索引
///
///      计算自旋轨道耦合的基底空间索引形式及约化与保留子空间维度
///
///      @param site_num         晶格格点数目
///      @param particle_num     粒子数目
///      @param position_reduce  约化位置
///      @param position_remain  保留位置
///      @param dim_of_space     约化和保留的子空间维度（输出）
///      @param space_with_index 基底空间的索引形式（输出）
///      @return                 不返回有效值
///      @see
///      @note

void Space_soc_index(size_t site_num,
                     size_t particle_num,
                     std::vector<size_t> &position_reduce,
                     std::vector<size_t> &position_remain,
                     std::vector<dim_soc_subspace> &dim_of_subspace,
                     std::vector<basis_soc_with_index> &space_with_index)
{
    /// 清空存储
    space_with_index.clear();
    dim_of_subspace.clear();
    
    /// 定义并初始化约化与保留位置粒子的数目
    size_t particle_num_remain = std::min(particle_num, position_remain.size()*2);
    size_t particle_num_reduce = particle_num - particle_num_remain;
    /// 定义约化与保留的基底空间，具有确定粒子数分配的基底空间索引形式
    std::vector<std::string> space_remain;
    std::vector<std::string> space_reduce;
    std::vector<basis_soc_with_index> space;
    
    /// 遍历粒子数分配
    do
    {
        /// 生成保留部分的基底空间
        Space_soc(position_remain.size(), particle_num_remain, space_remain);
        /// 生成约化部分的基底空间
        Space_soc(position_reduce.size(), particle_num_reduce, space_reduce);
        /// 存储保留部分和约化部分的基底空间的维度
        dim_of_subspace.push_back( dim_soc_subspace(space_reduce.size(), space_remain.size(),
                                                    particle_num_reduce, particle_num_remain)
                                  );
        /// 直积保留部分和约化部分的基底空间
        Space_soc_combin_subspace(
            site_num,
            particle_num_reduce,
            particle_num_remain,
            position_reduce,
            position_remain,
            space_reduce,
            space_remain,
            space
        );
        /// 对具有确定粒子数分配的基底空间按照索引排序
        std::sort(space.begin(), space.end());
        /// 将具有确定粒子数分配的基底空间插入到基底空间的最后
        space_with_index.insert(space_with_index.end(), space.begin(), space.end());
        /// 计算下一个粒子数分配方案
        particle_num_reduce++;
        particle_num_remain--;
        
    }
    while ( particle_num_reduce <= std::min(particle_num, position_reduce.size()*2) );
    
    return;
}

///      直积约化与保留的基底空间
///
///      对于确定的粒子数分配，合并保留和约化的基底空间
///
///      @param site_num             晶格格点数目
///      @param particle_num_reduce  约化位置粒子数
///      @param particle_num_remain  保留位置粒子数
///      @param position_reduce      约化位置
///      @param position_remain      保留位置
///      @param space_reduce         约化部分的基底空间
///      @param space_remain         保留部分的基底空间
///      @param space                基底空间的索引形式（输出）
///      @return                     不返回有效值
///      @see
///      @note
void Space_soc_combin_subspace(size_t site_num,
                               size_t particle_num_reduce,
                               size_t particle_num_remain,
                               std::vector<size_t> &position_reduce,
                               std::vector<size_t> &position_remain,
                               std::vector<std::string> &space_reduce,
                               std::vector<std::string> &space_remain,
                               std::vector<basis_soc_with_index> &space)
{
    /// 定义基底，位置变量
    std::string basis(2*site_num, '0');
    size_t pos = 0;
    
    /// 清空基底空间存储，重新分配内存
    space.clear();
    space.resize(space_reduce.size()*space_remain.size());
    
    /// 遍历约化部分的基底空间
    for(std::string s_reduce: space_reduce)
    {
        /// 设定约化位置的值
        for (size_t i = 0; i < position_reduce.size(); i++)
        {
            basis[position_reduce[i]] = s_reduce[i];
            basis[position_reduce[i]+site_num] = s_reduce[i+position_reduce.size()];
        }
        /// 遍历保留部分的基底空间
        for (std::string s_remain : space_remain)
        {
            /// 设定保留位置的值
            for (size_t i = 0; i < position_remain.size(); i++)
            {
                basis[position_remain[i]] = s_remain[i];
                basis[position_remain[i]+site_num] = s_remain[i+position_remain.size()];
            }
            /// 记录基底
            space[pos].basis = basis;
            /// 记录基底索引值
			space[pos].index_of_basis = index_val(s_reduce, s_remain);
            /// 记录基底子空间粒子数
            space[pos].particle_of_subspace = particle_num_subspace( particle_num_reduce,
                                                                     particle_num_remain );
            /// 自增基底位置
            pos++;
        }
    }

    return;
}



