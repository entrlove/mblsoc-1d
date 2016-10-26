//////////////////////////////////////////////////////////////////////////
///     COPYRIGHT NOTICE
///     Copyright (c) 2016, Swinburne University of Technology
///     All rights reserved.
///
/// @file   dec.h
/// @brief  自定义的类和函数声明
///
/// 由于要对基底空间进行索引排序以简化约化密度矩阵的计算复杂度，需计算基底对应的索引值。
/// 对于每一个基底，根据需要约化和需要保留的晶格格点位置，可以视为两个部分，即约化部分
/// 和保留部分之和，这两个部分构成两个子晶格。在晶格中填充一定粒子数的图案对应一个唯一
/// 的索引值。因此，对于一个基底，我们把两个子晶格的索引值作为这个基底的索引。在排序过
/// 程中，对粒子数分配方案进行首要排序，索引值进行次要排序，为避免在每一次排序过程中都
/// 计算两个子晶格的粒子数，将这两个量做为单独变量置于类中
///
/// @version 1.0
/// @author  曹业
/// @date    2016年8月1日
///
///
/// 修订说明：第一个版本
//////////////////////////////////////////////////////////////////////////

#ifndef dec_h
#define dec_h

#include <cstddef>
#include <string>
#include <vector>
/*** 使用PETSc库 *********************************/
#include <petsc.h>
/************************************************/
/*** 使用std::complex<double>代替MKL_Complex16 ***/
#include <ccomplex>
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
/***********************************************/
/*** 使用mpir库 *********************************/
#include <mpir.h>
/***********************************************/

///     本类的功能：存储两个索引值
///
///     约化部分和保留部分的索引值
class index_val
{
public:
    mpz_t index_val_reduce; /// 约化部分的索引值
    mpz_t index_val_remain; /// 保留部分的索引值
    /// 构造函数
	index_val(std::string &i_reduce, std::string &i_remain)
	{
		/// 初始化
		mpz_init(index_val_reduce);
		mpz_init(index_val_remain);
        /// 设定值
		mpz_set_str(index_val_reduce, i_reduce.c_str(), 2);
		mpz_set_str(index_val_remain, i_remain.c_str(), 2);
	};
    /// 构造函数（默认）
	index_val()
	{
		/// 初始化
		mpz_init(index_val_reduce);
		mpz_init(index_val_remain);
	};
	/// 拷贝构造
	index_val(const index_val &orig)
	{
		/// 初始化
		mpz_init(index_val_reduce);
		mpz_init(index_val_remain);
		/// 设定值
		mpz_set(index_val_reduce, orig.index_val_reduce);
		mpz_set(index_val_remain, orig.index_val_remain);
	}
	/// 重载赋值运算符（Assignment operator）
	index_val &operator=(const index_val &rhs)
	{
		mpz_set(index_val_reduce, rhs.index_val_reduce);
		mpz_set(index_val_remain, rhs.index_val_remain);
		return *this;
	}
	/// 析构函数
	~index_val()
	{
		mpz_clear(index_val_reduce);
		mpz_clear(index_val_remain);
	}
};

///     本类的功能：存储两个粒子数
///
///     约化部分和保留部分的粒子数
class particle_num_subspace
{
public:
    size_t particle_num_subspace_reduce; /// 约化部分的粒子数
    size_t particle_num_subspace_remain; /// 保留部分的粒子数
    /// 构造函数
    particle_num_subspace(size_t p_reduce, size_t p_remain)
        : particle_num_subspace_reduce(p_reduce),particle_num_subspace_remain(p_remain) {};
    /// 构造函数（默认）
    particle_num_subspace() : particle_num_subspace_reduce(), particle_num_subspace_remain() {};
};

///     本类的功能：存储基底，索引以及粒子数
///
///     存储基底和其对应的唯一的索引，显式存储了两个子晶格的粒子数
class basis_soc_with_index
{
public:
    index_val index_of_basis; /// 基底索引
    particle_num_subspace particle_of_subspace; /// 两部分的粒子数
    std::string basis; /// 基底
    /// 构造函数
    basis_soc_with_index(std::string &i_reduce, std::string &i_remain, 
		                 size_t p_reduce, size_t p_remain, 
						 std::string b) : index_of_basis(i_reduce, i_remain),
                                          particle_of_subspace(p_reduce, p_remain), 
										  basis(b)
    {};
    /// 构造函数（默认）
    basis_soc_with_index() : index_of_basis(), particle_of_subspace(), basis() {};
};

///     本类的功能：存储两个子空间维度和粒子数
///
///     存储具有确定粒子数的两个子晶格维度及其的粒子数
class dim_soc_subspace
{
public:
    size_t dim_subspace_reduce; /// 约化部分的维度
    size_t dim_subspace_remain; /// 保留部分的维度
    particle_num_subspace particle_of_subspace; /// 两部分的粒子数
    /// 构造函数
    dim_soc_subspace(size_t d_reduce, size_t d_remain, size_t p_num_reduce, size_t p_num_remain):
        dim_subspace_reduce(d_reduce),
        dim_subspace_remain(d_remain),
        particle_of_subspace(p_num_reduce, p_num_remain)
    {};
    /// 构造函数（默认）
    dim_soc_subspace() : dim_subspace_reduce(), dim_subspace_remain(), particle_of_subspace() {};
};

///     本类的功能：存储COO形式的矩阵元
///
///     存储矩阵元的值及其对应的行和列值
class matrix_element_with_index
{
public:
    double val; /// 矩阵元
    size_t row; /// 行值
    size_t col; /// 列值
    /// 构造函数
    matrix_element_with_index(double v, size_t r, size_t c) : val(v), row(r), col(c)
    {};
    /// 构造函数（默认）
    matrix_element_with_index() : val(), row(), col() {};
};

///      index_val类重载运算符<
///
///      内联函数，用于重载index_val类的<运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    mpz_cmp
///      @note
inline bool
operator<(const index_val &lhs, const index_val &rhs)
{
	int flag_reduce = mpz_cmp(lhs.index_val_reduce, rhs.index_val_reduce);
	int flag_remain = mpz_cmp(lhs.index_val_remain, rhs.index_val_remain);
	/// 优先比较index_val_reduce，若相等二者则进一步比较index_val_remain
	return ( flag_reduce < 0 || (flag_reduce == 0 && flag_remain < 0) );
}

///      index_val类重载运算符==
///
///      内联函数，用于重载index_val类==运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    mpz_cmp
///      @note
inline bool
operator==(const index_val &lhs, const index_val &rhs)
{
	int flag_reduce = mpz_cmp(lhs.index_val_reduce, rhs.index_val_reduce);
	int flag_remain = mpz_cmp(lhs.index_val_remain, rhs.index_val_remain);
	/// 当且仅当index_val_reduce和index_val_remain相等，返回真值
    return (flag_reduce == 0 && flag_remain == 0);
}

///      index_val类重载运算符！=
///
///      内联函数，用于重载index_val类！=运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    mpz_cmp
///      @note
inline bool
operator!=(const index_val &lhs, const index_val &rhs)
{
	int flag_reduce = mpz_cmp(lhs.index_val_reduce, rhs.index_val_reduce);
	int flag_remain = mpz_cmp(lhs.index_val_remain, rhs.index_val_remain);
	/// 若index_val_reduce或index_val_remain不等，返回真值
    return (flag_reduce != 0 || flag_remain != 0);
}

///      particle_num_subspace类重载运算符<
///
///      内联函数，用于重载particle_num_subspace类<运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    
///      @note
inline bool
operator<(const particle_num_subspace &lhs, const particle_num_subspace &rhs)
{
    return (lhs.particle_num_subspace_reduce <  rhs.particle_num_subspace_reduce);
}

///      particle_num_subspace类重载运算符==
///
///      内联函数，用于重载particle_num_subspace类==运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    
///      @note
inline bool
operator==(const particle_num_subspace &lhs, const particle_num_subspace &rhs)
{
    return (lhs.particle_num_subspace_reduce == rhs.particle_num_subspace_reduce);
}

///      particle_num_subspace类重载运算符!=
///
///      内联函数，用于重载particle_num_subspace类!=运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    
///      @note
inline bool
operator!=(const particle_num_subspace &lhs, const particle_num_subspace &rhs)
{
    return (lhs.particle_num_subspace_reduce != rhs.particle_num_subspace_reduce);
}

///     basis_soc_with_index类重载运算符<
///
///      内联函数，用于重载basis_soc_with_index类<运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    
///      @note
inline bool
operator<(const basis_soc_with_index &lhs, const basis_soc_with_index &rhs)
{
    return ( lhs.particle_of_subspace <  rhs.particle_of_subspace) ||
           ((lhs.particle_of_subspace == rhs.particle_of_subspace) &&
            (lhs.index_of_basis       <  rhs.index_of_basis));
}

///      matrix_element_with_index类重载运算符<
///
///      内联函数，用于重载matrix_element_with_index类<运算符
///
///      @param lhs              运算符左侧对象
///      @param rhs              运算符右侧对象
///      @return                 返回布尔值
///      @see                    
///      @note
inline bool
operator<(const matrix_element_with_index &lhs, const matrix_element_with_index &rhs)
{
    return ( lhs.row < rhs.row) || ((lhs.row == rhs.row) && (lhs.col < rhs.col));
}

///      模板函数，插入矩阵元
///
///      插入矩阵元及其对称位置矩阵元的COO形式（适用于实对称矩阵）
///
///      @param matrix           矩阵
///      @param row              矩阵的行
///      @param col              矩阵的列
///      @param val              插入元素值
///      @param row_num          插入元素所在行
///      @param col_num          插入元素所在列
///      @return                 不返回有效值
///      @see
///      @note
template<class T>
void Hamiltonian_val(std::vector<T> &matrix, std::vector<MKL_INT> &row, std::vector<MKL_INT> &col,
	T val, size_t row_num, size_t col_num)
{
	MKL_INT r = static_cast<MKL_INT>(row_num);
	MKL_INT c = static_cast<MKL_INT>(col_num);

	/// 插入矩阵元
	matrix.push_back(val);
	row.push_back(r);
	col.push_back(c);

	if (row_num == col_num) return;

	/// 如果非对角元，插入厄米共轭矩阵元
	matrix.push_back(val);
	row.push_back(c);
	col.push_back(r);

	return;
}

/// Function declaration
/// Generate basis
size_t Factorial(size_t m);
size_t Combination(size_t n, size_t m);
void Space_one_component(size_t site_num, size_t particle_num, std::vector<std::string> &space);
void Space_two_component(size_t site_num, size_t particle_num_up, size_t particle_num_dn,
                         std::vector<std::string> &space);
void Space_soc(size_t site_num, size_t particle_num, std::vector<std::string> &Space);
void Space_soc_combin_subspace
(
    size_t site_num,
    size_t particle_num_reduce,
    size_t particle_num_remain,
    std::vector<size_t> &position_reduce,
    std::vector<size_t> &position_remain,
    std::vector<std::string> &space_reduce,
    std::vector<std::string> &space_remain,
    std::vector<basis_soc_with_index> &space
 );
void Space_soc_index
(
    size_t site_num,
    size_t particle_num,
    std::vector<size_t> &position_reduce,
    std::vector<size_t> &position_remain,
    std::vector<dim_soc_subspace> &dim_of_subspace,
    std::vector<basis_soc_with_index> &space_with_index
 );

/// Density matrix
void Matrix_sum(MKL_Complex16 *matrix_lhs, MKL_Complex16 *matrix_rhs, size_t row_num, size_t col_num,
                size_t ldm_lhs, size_t ldm_rhs);
void Density_matrix_from_vector(size_t dim_vector, MKL_Complex16 *vector, MKL_Complex16 *density_matrix);
void Density_matrix_reduce(size_t dim_whole, size_t dim_remain,
                           std::vector<dim_soc_subspace> &dim_of_subspace,
                           MKL_Complex16 *density_matrix_whole, MKL_Complex16 *density_matrix_remain);

/// Generate basis with index
basis_soc_with_index Basis_soc_index(size_t site_num,
                                     std::string &basis_without_index,
                                     std::vector<size_t> &position_reduce,
                                     std::vector<size_t> &position_remain
                                     );
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
                    );
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
                    );
void Hamiltonian_dg(std::vector<basis_soc_with_index> &space,
                    Mat A,
                    int row_beg,
                    int row_end,
                    size_t site_num,
                    std::vector<double> &dis,
                    double hz,
                    double U,
                    double hb,
                    double mb
                    );
void Density_matrix_reduce(std::vector<dim_soc_subspace> &dim_of_subspace, MKL_Complex16 *vector,
                           MKL_Complex16 *density_matrix_remain);
void Entanglement_entropy(std::vector<dim_soc_subspace> &dim_of_subspace, MKL_Complex16 *vector,
                          double &ee);
void ravg(std::vector<double> &eigval, double &r);

#endif /* dec_h */
