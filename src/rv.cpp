#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
///      相邻能级差的平均比值
///
///      对输入的本征值序列进行排序，计算相邻能级差的比值，并对其进行平均
///
///      @param eigval           本征值序列，未排序
///      @param r                相邻能级差的平均比值
///      @return                 不返回有效值
///      @see
///      @note
void ravg(std::vector<double> &eigval, double &r)
{
    /// 获取本征值个数
    int num = static_cast<int>(eigval.size());
    /// 如果本征值个数小于3，设定r为NaN，返回
    if (num < 3) {r = std::numeric_limits<double>::quiet_NaN(); return; }
    /// 排序本征值序列
    std::sort(eigval.begin(),eigval.end());
    /// 申请空间
    std::vector<double> dn(num-1);
    std::vector<double> rn(num-2);
    /// 计算相邻能级差
    for (int i = 0; i < dn.size(); i++){
        dn[i] = eigval[i+1] - eigval[i];
    }
    /// 计算相邻能级差的比值
    for (int i = 0; i < rn.size(); i++){
        rn[i] = std::min(dn[i+1],dn[i])/std::max(dn[i+1],dn[i]);
    }
    /// 计算相邻能级差比值的平均值
    r  = std::accumulate(rn.begin(), rn.end(), 0.0);
    r /= rn.size();
    return;
}
