#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <stdexcept>
#include <vector>
/*** 使用std::complex<double>代替MKL_Complex16 ***/
#include <ccomplex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
/***********************************************/
#include "dec.h"

///      约化密度矩阵
///
///      对一个纯态，求保留位置的约化密度矩阵
///
///      @param dim_of_subspace  每一个约化、保留子空间对的维度
///      @param vector           系统的一个纯态
///      @param density_matrix   约化密度矩阵，按照子空间对顺序排列（dm1, dm2, ...）（输出）
///      @return                 不返回有效值
///      @see
///      @note
void Density_matrix_reduce(std::vector<dim_soc_subspace> &dim_of_subspace, MKL_Complex16 *vector,
                           MKL_Complex16 *density_matrix_remain)
{
    /// 密度矩阵起始位置
    MKL_Complex16 *dm_current = density_matrix_remain;
    /// 纯态起始位置
    MKL_Complex16 *vc_current = vector;
    /// 遍历所有子空间对
    for (dim_soc_subspace d : dim_of_subspace)
    {
        /// 初始化密度矩阵的子块
        std::fill(dm_current, dm_current+d.dim_subspace_remain*d.dim_subspace_remain, std::complex<double>(0, 0));
        /// 对密度矩阵子块赋值
        for (size_t i = 0; i < d.dim_subspace_reduce; i++)
        {
            cblas_zher (CblasColMajor, CblasUpper, static_cast<int>(d.dim_subspace_remain), 1.0,
                        vc_current, 1, dm_current, static_cast<int>(d.dim_subspace_remain));
            vc_current += d.dim_subspace_remain;
        }
        /// 移动指针指向下一个密度矩阵子块
        dm_current += d.dim_subspace_remain*d.dim_subspace_remain;
    }
    return;
}
///      纠缠熵
///
///      对一个纯态，求保留位置的纠缠熵
///
///      @param dim_of_subspace  每一个约化、保留子空间对的维度
///      @param vector           系统的一个纯态
///      @param ee               纠缠熵（输出）
///      @return                 不返回有效值
///      @see
///      @note
void Entanglement_entropy(std::vector<dim_soc_subspace> &dim_of_subspace, MKL_Complex16 *vector,
                          double &ee)
{
    /// 初始化纠缠熵
    ee = 0;
    /// 计算约化密度矩阵子块的长度和
    MKL_INT dm_length = 0;
    for (dim_soc_subspace d : dim_of_subspace) {dm_length += d.dim_subspace_remain*d.dim_subspace_remain;}
    /// 给约化密度矩阵申请内存
    MKL_Complex16 *density_matrix_remain = (MKL_Complex16 *)mkl_malloc(dm_length*sizeof(MKL_Complex16), 64);
    /// 计算密度矩阵
    Density_matrix_reduce(dim_of_subspace, vector, density_matrix_remain);
    /// 密度矩阵起始位置
    MKL_Complex16 *dm_current = density_matrix_remain;
    /// 分块对角化密度矩阵
    for (dim_soc_subspace d : dim_of_subspace)
    {
        char jobz        = 'N';
        char range       = 'A';
        char uplo        = 'U';
        MKL_INT n        = static_cast<MKL_INT>(d.dim_subspace_remain);
        MKL_INT lda      = n;
        MKL_INT ldz      = n;
        MKL_INT il       = 0;
        MKL_INT iu       = 0;
        MKL_INT m;
        MKL_INT info;
        double abstol    = -1;
        double vl        = 0;
        double vu        = 0;
        MKL_INT *isuppz  = (MKL_INT       *)mkl_malloc(n*2 *sizeof(MKL_INT),       64);
        double *w        = (double        *)mkl_malloc(n   *sizeof(double),        64);
        MKL_Complex16 *z = (MKL_Complex16 *)mkl_malloc(n*n *sizeof(MKL_Complex16), 64);
        MKL_Complex16 *a = dm_current;
        
        info = LAPACKE_zheevr(LAPACK_COL_MAJOR, jobz, range, uplo, n, a, lda,
                              vl, vu, il, iu, abstol, &m, w, z, ldz, isuppz);
        ee += std::accumulate(w, w+m, 0.0, [](double s, double eg){return s - std::abs(eg)*log(std::abs(eg));});
        /// 移动指针指向下一个密度矩阵子块
        dm_current += d.dim_subspace_remain*d.dim_subspace_remain;
        /// 释放内存
        mkl_free(isuppz);
        mkl_free(w);
        mkl_free(z);
    }
    /// 释放内存
    mkl_free(density_matrix_remain);
    return;
}



