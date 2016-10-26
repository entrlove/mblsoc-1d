//////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   sp.cpp
/// @brief  生成系统的基底空间
///
/// 生成自旋1/2的费米子基底空间
/// 自旋守恒：计算具有固定粒子数且固定自旋数的单自旋组份，双自旋组份基底空间
/// 自旋轨道耦合：计算具有固定粒子数的基底空间
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

///      单组份自旋基底空间
///
///      计算自旋守恒的单组份费米子基底空间
///
///      @param site_num         晶格格点数目
///      @param particle_num     粒子数目
///      @param space            基底（输出）
///      @return                 不返回有效值
///      @see                    Factorial, Combination, next_permutation
///      @note
void Space_one_component(size_t site_num, size_t particle_num, std::vector<std::string> &space)
{
    /// 判断是否粒子数大于晶格格点数
    if (particle_num > site_num)
    {
        throw std::runtime_error("Particle number must be less than site number.");
    }
 
    /// 定义并初始化基底（所有自旋从右侧填充）
    std::string basis = std::string(site_num-particle_num, '0') + std::string(particle_num, '1');
    size_t pos = 0;
    
    /// 清空基底空间存储，重新分配内存
    space.clear();
    space.resize( Combination(site_num, particle_num) );
    
    /// 生成所有基底
    do
    {
        space[pos++] = basis;
    }
    while ( std::next_permutation(basis.begin(), basis.end()) );
    
    return;
}

///      双组份自旋基底空间
///
///      计算自旋守恒的双组份费米子基底空间
///
///      @param site_num         晶格格点数目
///      @param particle_num_up  自旋向上粒子数目
///      @param particle_num_dn  自旋向下粒子数目
///      @param space            基底（输出）
///      @return                 不返回有效值
///      @see                    Space_one_component
///      @note
void Space_two_component(size_t site_num, size_t particle_num_up, size_t particle_num_dn,
                         std::vector<std::string> &space)
{
    /// 判断是否任一自旋组份的粒子数大于晶格格点数
    if ( (particle_num_up > site_num) || (particle_num_dn > site_num) )
    {
        throw std::runtime_error("Particle number of each component must be less than \
                                  site number.");
    }
    
    /// 定义两个自旋组份的基底空间
    std::vector<std::string> space_up;
    std::vector<std::string> space_dn;
	size_t pos = 0;

    /// 清空基底空间存储，重新分配内存
    space.clear();
    space.resize(Combination(site_num, particle_num_up)*Combination(site_num, particle_num_dn));
    
    /// 生成两个自旋组份的基底空间
    Space_one_component(site_num, particle_num_up, space_up);
    Space_one_component(site_num, particle_num_dn, space_dn);
    
    /// 直积两个自旋组份的基底空间生成基底空间
    for (std::string basis_up : space_up)
    {
        for (std::string basis_dn : space_dn)
        {
			/// 拼接两个自旋组分的图案，自旋向上在低位，自旋向下在高位
            space[pos++] = basis_up + basis_dn;
        }
    }
    
    return;
}

///      自旋轨道耦合的基底空间
///
///      计算自旋轨道耦合的费米子基底空间
///
///      @param site_num         晶格格点数目
///      @param particle_num     粒子数目
///      @param space            基底（输出）
///      @return                 不返回有效值
///      @see                    Space_two_component
///      @note
void Space_soc(size_t site_num, size_t particle_num, std::vector<std::string> &space)
{
    /// 判断是否粒子数大于晶格格点数目的两倍
    if (particle_num > site_num * 2)
    {
        throw std::runtime_error("Particle number of SOC Fermions must be less than \
                                  twice of the site number.");
    }
    
    /// 定义并初始化两个自旋组份的粒子数，定义自旋守恒的基底空间
    size_t particle_num_up = std::min(site_num, particle_num);
    size_t particle_num_dn = particle_num - particle_num_up;
    std::vector<std::string> space_spin_conservation;
    
    /// 清空基底空间存储
    space.clear();
    
    /// 对所有可能的自旋组份分配计算自旋守恒的基底空间，对这些空间求和得到基底空间
    do
    {
        Space_two_component(site_num, particle_num_up--, particle_num_dn++, space_spin_conservation);
        space.insert(space.end(), space_spin_conservation.begin(), space_spin_conservation.end());
    }
    while(particle_num_dn <= std::min(site_num, particle_num));
    
    return;
}
