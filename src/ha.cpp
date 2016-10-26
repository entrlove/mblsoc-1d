/////////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   ha.cpp
/// @brief  系统的哈密顿量
///
/// 动能，自旋轨道耦合（Rashba，Dresselhaus），对角（hz，相互作用，无序，边界）哈密顿量
/// 系统哈密顿量
///
/// @version 1.6
/// @author  曹业
/// @date    2016年10月4日
///
///
/// 修订说明：第六个版本
/////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>
/*** 使用PETSc库 *********************************/
#include <petsc.h>
/************************************************/
/*** 使用std::complex<double>代替MKL_Complex16 ***/
#include <ccomplex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
/***********************************************/
#include "dec.h"

///      动能及Dresselhaus自旋轨道耦合哈密顿量
///
///      动能及Dresselhaus自旋轨道耦合哈密顿量
///
///      @param space            基底
///      @param position_reduce  约化位置
///      @param position_remain  保留位置
///      @param A                哈密顿量（输出）
///      @param row_beg          当前进程的起始行
///      @param row_end          当前进程的结束行（+1）
///      @param site_num         晶格格点数目
///      @param particle_num     粒子数
///      @param hopping_diatance 跳跃距离
///      @param t                动能强度
///      @param lambda_y         Dresselhaus自旋轨道耦合强度
///      @param boundary         布尔值，0表示周期性边界条件，1表示开边界条件
///      @return                 不返回有效值
///      @see
///      @note
void Hamiltonian_hd(std::vector<basis_soc_with_index> &space,
                    std::vector<size_t> &position_reduce,
                    std::vector<size_t> &position_remain,
                    Mat A,
                    int row_beg,
                    int row_end,
                    size_t site_num,
                    size_t particle_num,
                    size_t hopping_distance,
                    double t,
                    double lambda_y,
                    bool boundary
                    )
{
    /// 哈密顿量的行，列位置
    size_t row, col;
    /// 晶格格点当前位置，跳跃目标位置（向左，向右）
    size_t pos_cur, pos_lft, pos_rgt, pos_up_cur, pos_up_lft, pos_up_rgt, pos_dn_cur, pos_dn_lft, pos_dn_rgt;
    /// 当前基底（列），目标基底（行）
    std::string basis, basis_next;
    /// 目标基底迭代器
    std::vector<basis_soc_with_index>::iterator iter_basis_next;
    /// 定义lambda表达式用于插入矩阵元
    auto Insert_element =
        [&](size_t pos_curr, size_t pos_targ, std::complex<double> val)
        {
            if (basis[pos_curr] == '1' && basis[pos_targ] == '0')
            {
                basis_next = basis;
                std::swap(basis_next[pos_curr], basis_next[pos_targ]);
                /// 查找目标基底在基底空间中的位置
                iter_basis_next = std::lower_bound(space.begin(), space.end(),
                    Basis_soc_index(site_num, basis_next, position_reduce, position_remain));
                /// 得到目标基底列位置
                col = iter_basis_next - space.begin();
                /// 插入矩阵元
                MatSetValue(A, static_cast<int>(row), static_cast<int>(col), val, INSERT_VALUES);
            }
        };
    /// 遍历哈密顿量行，<row|H_hop|col>
    for (row = row_beg; row < row_end; row++)
    {
        /// 取出当前行所对应基底
        basis = space[row].basis;
        /// 遍历晶格格点，计算向右跳跃哈密顿量
        for (pos_cur = 0; pos_cur < site_num; pos_cur++)
        {
            /// 计算向右向左跳跃的晶格位置
            pos_rgt =  (pos_cur+hopping_distance) % site_num;
            pos_lft = ((pos_cur-hopping_distance) + site_num) % site_num;
            /// 计算两个自旋组份当前位置和跳跃的位置
            pos_up_cur = pos_cur;
            pos_up_lft = pos_lft;
            pos_up_rgt = pos_rgt;
            pos_dn_cur = pos_cur + site_num;
            pos_dn_lft = pos_lft + site_num;
            pos_dn_rgt = pos_rgt + site_num;
            /// 右向跳跃
            if (pos_rgt > pos_cur || !(boundary))
            {
                Insert_element(pos_up_cur, pos_up_rgt, std::complex<double>(-t, -lambda_y));
                Insert_element(pos_dn_cur, pos_dn_rgt, std::complex<double>(-t, +lambda_y));
            }
            /// 左向跳跃
            if (pos_lft < pos_cur || !(boundary))
            {
                Insert_element(pos_up_cur, pos_up_lft, std::complex<double>(-t, +lambda_y));
                Insert_element(pos_dn_cur, pos_dn_lft, std::complex<double>(-t, -lambda_y));
            }
        }
    }
    return;
}

///      Rashba自旋轨道耦合哈密顿量
///
///      Rashba自旋轨道耦合哈密顿量
///
///      @param space            基底
///      @param position_reduce  约化位置
///      @param position_remain  保留位置
///      @param A                哈密顿量（输出）
///      @param row_beg          当前进程的起始行
///      @param row_end          当前进程的结束行（+1）
///      @param site_num         晶格格点数目
///      @param particle_num     粒子数
///      @param hopping_diatance 跳跃距离
///      @param lambda_z         Rashba自旋轨道耦合强度
///      @param boundary         布尔值，0表示周期性边界条件，1表示开边界条件
///      @return                 不返回有效值
///      @see
///      @note
void Hamiltonian_sr(std::vector<basis_soc_with_index> &space,
                    std::vector<size_t> &position_reduce,
                    std::vector<size_t> &position_remain,
                    Mat A,
                    int row_beg,
                    int row_end,
                    size_t site_num,
                    size_t particle_num,
                    size_t hopping_distance,
                    double lambda_z,
                    bool boundary
                    )
{
    /// 哈密顿量的行，列位置
    size_t row, col;
    /// 晶格格点当前位置，跳跃目标位置（向左，向右）
    size_t pos_cur, pos_lft, pos_rgt, pos_up_cur, pos_dn_cur, pos_up_lft, pos_up_rgt, pos_dn_lft, pos_dn_rgt;
    /// 当前基底（列），目标基底（行）
    std::string basis, basis_next;
    /// 目标基底的迭代器
    std::vector<basis_soc_with_index>::iterator iter_basis_next;
    /// 定义lambda表达式用于插入矩阵元
    auto Insert_element =
        [&](size_t pos_curr, size_t pos_targ, std::complex<double> val)
        {
            if (basis[pos_curr] == '1' && basis[pos_targ] == '0')
            {
                basis_next = basis;
                std::swap(basis_next[pos_curr], basis_next[pos_targ]);
                /// 查找目标基底在基底空间中的位置
                iter_basis_next = std::lower_bound(space.begin(), space.end(),
                Basis_soc_index(site_num, basis_next, position_reduce, position_remain));
                /// 得到目标基底列位置
                col = iter_basis_next - space.begin();
                /// 插入矩阵元
                MatSetValue(A, static_cast<int>(row), static_cast<int>(col), val, INSERT_VALUES);
            }
        };
    /// 遍历哈密顿量的行，<row|H_hop|col>
    for (row = row_beg; row < row_end; row++)
    {
        /// 取出当前行所对应基底
        basis = space[row].basis;
        /// 遍历晶格格点，计算向右跳跃哈密顿量
        for (pos_cur = 0; pos_cur < site_num; pos_cur++)
        {
            /// 计算向右向左跳跃的晶格位置
            pos_rgt =  (pos_cur+hopping_distance) % site_num;
            pos_lft = ((pos_cur-hopping_distance) + site_num) % site_num;
            /// 计算两个自旋组份当前位置和跳跃的位置
            pos_up_cur = pos_cur;
            pos_up_lft = pos_lft;
            pos_up_rgt = pos_rgt;
            pos_dn_cur = pos_cur + site_num;
            pos_dn_lft = pos_lft + site_num;
            pos_dn_rgt = pos_rgt + site_num;
            /// 右向跳跃
            if (pos_rgt > pos_cur || !(boundary))
            {
                Insert_element(pos_up_cur, pos_dn_rgt, std::complex<double>(+lambda_z, 0));
                Insert_element(pos_dn_cur, pos_up_rgt, std::complex<double>(-lambda_z, 0));
            }
            /// 左向跳跃
            if (pos_lft < pos_cur || !(boundary))
            {
                Insert_element(pos_up_cur, pos_dn_lft, std::complex<double>(-lambda_z, 0));
                Insert_element(pos_dn_cur, pos_up_lft, std::complex<double>(+lambda_z, 0));
            }
        }
    }
    return;
}

///      对角哈密顿量
///
///      生成无序，相互作用，磁场及边界哈密顿量
///
///      @param space            基底
///      @param A                对角哈密顿量（输出）
///      @param row_beg          当前进程的起始行
///      @param row_end          当前进程的结束行（+1）
///      @param site_num         晶格格点数目
///      @param dis              格点无序强度
///      @param hz               磁场强度
///      @param U                相互作用强度
///      @param hb               边界哈密顿量强度（左侧）
///      @param ub               边界哈密顿量强度（右侧）
///      @return                 不返回有效值
///      @see
///      @note
void Hamiltonian_dg(std::vector<basis_soc_with_index> &space,
                    Mat A,
                    int row_beg,
                    int row_end,
                    size_t site_num,
                    std::vector<double> &dis,
                    double hz,
                    double U,
                    double hb,
                    double ub
                    )
{
    /// 哈密顿量的行位置
    size_t row;
    /// 晶格格点当前位置
    size_t pos_cur;
    /// 当前基底
    std::string basis;
    /// 哈密顿量矩阵元值
    double val;
   
    /// 遍历哈密顿量的行
    for (row = row_beg; row < row_end; row++)
    {
        /// 取出当前列所对应基底
        basis = space[row].basis;
        /// 初始化矩阵元值
        val = 0;
        /// 累加边界哈密顿量矩阵元
        if (basis[0] == '1') val += hb;
        if (basis[site_num] == '1') val -= hb;
        if (basis[site_num - 1] == '1') val += ub;
        if (basis[site_num - 1 + site_num] == '1') val += ub;
        /// 遍历晶格格点
        for (pos_cur = 0; pos_cur < site_num; pos_cur++)
        {
            /// 累加磁场矩阵元
            if (basis[pos_cur] == '1') { val += hz; }
            if (basis[pos_cur+site_num] == '1') { val += (-hz); }
            /// 累加无序矩阵元
            if (basis[pos_cur] == '1') { val += dis[pos_cur]; }
            if (basis[pos_cur+site_num] == '1') { val += dis[pos_cur]; }
            /// 累加相互作用矩阵元
            if (basis[pos_cur] == '1' && basis[pos_cur+site_num] == '1') { val += (-U); }
        }
        /// 插入对角矩阵元的COO形式
        MatSetValue(A, static_cast<int>(row), static_cast<int>(row), val, INSERT_VALUES);
    }
    return;
}

