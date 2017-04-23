# plot files sent by mike

col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
col7 = []
with open('/Users/danielpe/Dropbox/dark/leicester/data/test_civ_n23/l_s0_c4_nh10_fort.16') as f:
    for line in f:
        col1.append(float((line.split())[0]))
        col2.append(float((line.split())[1]))
        col3.append(float((line.split())[2]))
        col4.append(float((line.split())[3]))
        col5.append(float((line.split())[4]))
        col6.append(float((line.split())[5]))
        col7.append(float((line.split())[6]))
import matplotlib.pyplot as plt
plt.plot(col2,col6)
plt.axhline(y=42.2,color='red')
plt.show()
