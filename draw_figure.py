# import matplotlib.pyplot as plt
# import numpy as np
#
# x, y = np.mgrid[-2 : 2 : 20j, -2 : 2 : 20j]
# z = 50 * np.sin(x + y)                     # 测试数据
# ax = plt.subplot(111, projection = '3d')   # 三维图形
# ax.plot_surface(x, y, z, rstride = 2, cstride = 1, cmap = plt.cm.Blues_r)
# ax.set_xlabel('x')                         # 设置坐标轴标签
# ax.set_xlabel('y')
# ax.set_xlabel('z')
# plt.show()

# import GPy
# # GPy.tests()

from GPy.models.gplvm import GPLVM  # observed data: simulate data
# from GPy.models.ss_mrd import IBPPrior_SSMRD
import numpy as np
Y=np.array([[1,2,3],
            [2,3,4]])
output=GPLVM(Y,1,init='PCA')
print(output)
