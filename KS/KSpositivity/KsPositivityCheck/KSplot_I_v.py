import matplotlib.pyplot as plt
import numpy as np

x = np.array([17*17,33*33,65*65,129*129])
y = np.array([6.57229,2.13234,0.570987,0.147639])
z = np.array([0.571766,0.156583,0.0401509,0.0101183])

plt.figure(figsize=(8,6))            #图片尺寸
plt.loglog(x,y,'r-o',label=r'$|| v - v_h ||_{\infty}$',ms=5)            #双对数坐标 红色标记圆圈  标记尺寸10
plt.loglog(x,z,'b-o',label=r'$|| v - v_h ||_{0}$',ms=5)


ln_x = np.log(x)                       #坐标取对数
ln_y = np.log(y)                       #坐标取对数
ln_z = np.log(z)
A = np.vstack([ln_x,np.ones(len(ln_x))]).T      # 合并数组-转置 构造系数矩阵
B = np.vstack([ln_x,np.ones(len(ln_x))]).T


m,c = np.linalg.lstsq(A,ln_y,rcond=None)[0]        # 最小二乘法拟合


m1,c1 = np.linalg.lstsq(B,ln_z,rcond=None)[0]

print('拟合斜率：%.5f'%m+'\n'+'拟合b值：%.5f'%c)
print('拟合斜率：%.5f'%m1+'\n'+'拟合b值：%.5f'%c1)
plt.plot(np.e**ln_x,np.e**(m*ln_x+c-0.6), 'r-',label=r'$CN^{-0.9402}$')  #round(-m,2) 四舍五入 -两位小数
plt.plot(np.e**ln_x,np.e**(m1*ln_x+c1-0.6), 'b-',label=r'$CN^{-0.9962}$')  #round(-m,2) 四舍五入 -两位小数


plt.legend()
plt.show()
