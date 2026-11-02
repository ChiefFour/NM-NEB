import os
import re
import numpy as np
with open('1.txt','r') as one:
    one1=one.read()
with open('2.txt','r') as two:
    two2=two.read()
# 打开一个文件以追加模式写入结果
with open('999', 'a') as f:
    with open('input1', 'r') as f1:
        lines = f1.readlines()
        for B1 in np.arange(1.0, 1.05, 0.02):
            for B2 in np.arange(1.0, 1.05, 0.02):
                for B3 in np.arange(1.0, 1.05, 0.02):
                    for A1 in np.arange(108, 110, 1):
                        fpath = "input"
                        f.write(f'109 -120 {A1}\n')
                        with open(fpath, 'w') as f2:
                            for eachline in lines:
                                content2 = re.sub('B1', '{:.4f}'.format(B1), eachline)
                                content3 = re.sub('B2', '{:.4f}'.format(B2), content2)
                                content4 = re.sub('B3', '{:.4f}'.format(B3), content3)
                                content5 = re.sub('A1', '{:.4f}'.format(A1), content4)
                                f2.writelines(content5)
                        os.system('./zmat_to_cart')
                        with open('output.txt','r') as f3:
                            contents=f3.read()
                        fpath11="{:.3f}-{:.3f}-{:.3f}-{:.0f}.inp".format(B1,B2,B3,A1)
                        with open(fpath11,'w') as f4:
                            f4.write(f'{one1}{contents}{two2}')
                        #os.system('rm input output.txt')
