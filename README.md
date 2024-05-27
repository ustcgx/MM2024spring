# 烟尘扩散模型
二阶偏微分方程有限元方法模拟仿真

## 问题描述
一座高为h的烟囱向外以速率Q喷出烟尘，空间内有恒定风速v，建模描述空间各处烟尘浓度分布（虽然仅要求地面处，但实际模拟需要计算全空间）

## 扩散基本知识
扩散方程是二阶PDE，自由扩散基本方程为： 
$\frac{\partial c}{\partial t} = D\Delta c$

其中c是浓度，是空间和时间，即x,y,z,t的函数，是要求解的目标.

在本例中，我们有额外条件
1. 边界条件1：无穷远处浓度为0
2. 边界条件2：灰尘不能穿越地面，此处有一个一阶条件，即建模假设的3.
3. 初始条件：初始t=0时刻，c处处为0
4. 附加条件：仅在烟囱口(0,0,h)处，浓度以恒定速率增加（Q $kg/s$，考虑到离散化后单位体积是1 $m^3$，也可以直接认为是Q $ks/(m^3·s$)）
   
最终完整的问题描述见 report.pdf

## 文件说明
1. diffusion.m：扩散模拟函数，接收参数模拟扩散过程（也可以继承上次的扩散的时间和相应结果以继续模拟
2. diffusion_gravity.m：改自diffusion.m，只修改了两处地方，考虑重力影响
3. diffusion_reflection.m：<font color="red">TODO</font>
4. visualization.m：绘图函数，对扩散过程可视化，主要是接受当前空间内浓度函数C，选择某些特点的位置进行浓度可视化
5. simulation_instance：对一组给定参数的问题实例进行模拟可视化
6. simulation_h（殷）：探究 **烟囱高度** 对 **地面浓度最大值** 的影响，模拟多个实例，对最终各实例的结果绘制上述两者的关系图(此时由于我们不关心模拟的具体过程，可以选择不可视化，设置diffusion函数show_fig=false)
7. simulation_v（殷）：探究 **风速** 和 **地面浓度分布情况** 的影响，直接保存模拟过程中的图像，根据图像来比较说明风速会导致地面浓度中心偏移
8. simulation_r：<font color="red">TODO</font>
9. simulation_w（李）：增加 **重力项**，直接模拟绘图，观察 **地面浓度分布情况** ，直接保存模拟过程图像，根据图像来比较说明重力影响可能导致地面浓度增加更迅速.
10. ......